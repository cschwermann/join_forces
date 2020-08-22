!! Author:  Christian Schwermann
!! E-mail:  c.schwermann@wwu.de
!! Date:    06/11/2019
!! Project: Joinforces
!! File:    joinforces.f90 
!! Copyright: © 2019 Christian Schwermann, ALL RIGHTS RESERVED
!!
!!***********************************************************************
!!
!! Analysis tool, takes the average forces from lammps TMD simulations
!! for some windows and joins them approximately in the middle,
!! minimizing the change in free energies induced by changing the joining point
!!
!!***********************************************************************
program Joinforces
   implicit none
   !! Precision
   integer, parameter :: DP = Selected_real_kind(15, 307)
   !! Input: filenames
   character(len=64)  :: filnam
   !! Output: text files
   character(len=64)  :: outnam, outnam2, outnam3
   !! Internal: #arguments, loop index, read & write error, index of least difference, #gridpoints
   integer            :: numargs, i, j, ierror, iout, jbest, npoints
   !! Internal: force +-dx, free energy +- dx
   real(kind=DP) :: force, forcemin, forcemax, free, freemin, freemax
   !! Internal: x value, x step, x min and max, joining x, difference of free energies, minimum diff
   real(kind=DP) :: x, dx, xmin, xmax, xj, diffj, diffmin
   !! Internal: factors, damping for tanh, shift of free energy
   real(kind=DP) :: fac1, fac2, damp, shift
   !! Internal: forces for both input files
   real(kind=DP), allocatable :: force1(:,:), force2(:,:)
   !! Internal: total force, free energy, old free energy
   real(kind=DP), allocatable :: force_full(:), free_ener(:), old_free(:)
   !!

   ! should be input, but lazy
   outnam = "differences.dat"
   outnam2 = "average_force.dat"
   outnam3 = "free_energy.dat"

   npoints = 200
   damp = 10.0_DP

   numargs = Iargc()

   if( numargs < 2 ) then
      write(*,*) "Error: need at least two input files!"
      error stop
   end if

   ! Open all input files in parallel
   call Getarg( 1, filnam )
   open( unit = 11, file = Trim( filnam ), status = "old", action = "read", iostat = ierror )
   if( ierror .ne. 0 ) then
      write(*,*) "Error reading file number ",1 ," ",trim(filnam),"!"
      error stop
   end if

   call Getarg( 2, filnam )
   open( unit = 12, file = Trim( filnam ), status = "old", action = "read", iostat = ierror )
   if( ierror .ne. 0 ) then
      write(*,*) "Error reading file number ",2 ," ",trim(filnam),"!"
      error stop
   end if

   ! count length of input files
   i = 0
   do
      read(11,*,iostat=ierror) x, force
      if( ierror /= 0 ) exit
      i = i + 1
   enddo
   allocate(force1(1:2,1:i))
   write(*,*) "Length of file 1: ", i
   rewind(11)

   i = 0
   do
      read(12,*,iostat=ierror) x, force
      if( ierror /= 0 ) exit
      i = i + 1
   enddo
   allocate(force2(1:2,1:i))
   write(*,*) "Length of file 2: ", i
   rewind(12)

   ! read input data
   do i = 1, size( force1(:,:), 2 )
      read(11,*,iostat=ierror) force1(1,i), force1(2,i)
   enddo
   do i = 1, size( force2(:,:), 2 )
      read(12,*,iostat=ierror) force2(1,i), force2(2,i)
   enddo

   close(11)
   close(12)

   ! sort arrays!
   call sort( force1(:,:) )
   call sort( force2(:,:) )

   ! find x range
   xmin= min( minval(force1(1,:)), minval(force2(1,:)) )
   xmax= max( maxval(force1(1,:)), maxval(force2(1,:)) )
   dx = ( xmax - xmin ) / real( npoints - 1, kind=DP )

   ! interpolate smoothly
   call interpolate( force1, npoints / 2 , xmin, xmax, .false. )
   call interpolate( force2, npoints / 2 , xmin, xmax, .false. )

   ! calculate derivatives
   do i=2, npoints / 2
      force1(3,i) = ( force1(2,i) - force1(2,i-1) ) / ( force1(1,i) - force1(1,i-1) )
      force2(3,i) = ( force2(2,i) - force2(2,i-1) ) / ( force2(1,i) - force2(1,i-1) )
   end do
   force1(3,1) = 0.0_DP
   force2(3,1) = 0.0_DP

   ! write linear interpolation
   open( unit = 42, file = "linear.dat", status = "replace", action = "write" )
   do i = 1, size( force1(:,:), 2 )
      write(42,*) force1(1,i), force1(2,i), force2(2,i) 
   end do
   close(42)

   call interpolate( force1, npoints, xmin, xmax, .true. )
   call interpolate( force2, npoints, xmin, xmax, .true. )

   ! write cubic interpolation
   open( unit = 43, file = "cubic.dat", status = "replace", action = "write" )
   do i = 1, size( force1(:,:), 2 )
      write(43,*) force1(1,i), force1(2,i), force2(2,i) 
   end do
   close(43)

   ! iterate over joining points
   open( unit = 13, file = Trim( outnam ), status = "replace", action = "write", iostat = iout )
   if( iout /= 0) then
      write(*,*) "Error opening file ", Trim( outnam ),"!"
      error stop
   end if

   allocate( force_full(1:npoints) )
   allocate( free_ener(1:npoints) )
   allocate( old_free(1:npoints) )

   old_free(1:npoints) = 0.0_DP
   diffmin = 1.0E20_DP
   jbest = -1
   do j = Int(npoints*0.2), Int(npoints*0.8)
      xj = ( j - 1 ) * dx + xmin

      ! calculate free energy for joining point xj
      free = 0.0_DP
      force_full(1:npoints) = 0.0_DP
      free_ener(1:npoints) = 0.0_DP
      do i=1, npoints
         x = ( i - 1 ) * dx + xmin
   
         fac1 = 0.5_DP + 0.5_DP * tanh( damp * ( xj - x ) )
         fac2 = 0.5_DP + 0.5_DP * tanh( damp * ( x - xj ) )
         force_full(i) = fac1 * force1(2,i) + fac2 * force2(2,i)
         free = free + force_full(i) * dx
         free_ener(i) = free
   
      end do

      ! shift minimum to zero
      !free_ener(1:npoints) = free_ener(1:npoints) - minval( free_ener(1:npoints) )

      ! calculate difference to previous free energy
      diffj = sum( ( free_ener(1:npoints) - old_free(1:npoints) ) ** 2 )

      ! alternative: try to put free energies as close as possible to each other
      ! brute force here
      do i = 1, npoints
         shift = minval( free_ener(1:npoints) ) * 2.0_DP * real( i - npoints / 2, kind=DP ) / real( npoints, kind=DP )
         diffj = min( diffj, sum( ( free_ener(1:npoints) - old_free(1:npoints) - shift ) ** 2) )
      end do

      ! update minimum difference
      if( diffj < diffmin ) then
         diffmin = diffj
         jbest = j
      end if
      ! update
      old_free(1:npoints) = free_ener(1:npoints)

      write(13,*) xj, diffj
   end do

   close(13)

   open( unit = 14, file = Trim( outnam2 ), status = "replace", action = "write", iostat = iout )
   if( iout /= 0) then
      write(*,*) "Error opening file ", Trim( outnam2 ),"!"
      error stop
   end if

   open( unit = 15, file = Trim( outnam3 ), status = "replace", action = "write", iostat = iout )
   if( iout /= 0) then
      write(*,*) "Error opening file ", Trim( outnam3 ),"!"
      error stop
   end if


   ! calculate optimal free energy, force and +- 1dx  for reference / as error bar
   free = 0.0_DP
   freemin = 0.0_DP
   freemax = 0.0_DP
   do i=1, npoints
      x = ( i - 1 ) * dx + xmin
   
      xj = ( jbest - 1 ) * dx + xmin
      fac1 = 0.5_DP + 0.5_DP * tanh( damp * ( xj - x ) )
      fac2 = 0.5_DP + 0.5_DP * tanh( damp * ( x - xj ) )
      force = fac1 * force1(2,i) + fac2 * force2(2,i)
      free = free + force * dx

      xj = ( jbest - 2 ) * dx + xmin
      fac1 = 0.5_DP + 0.5_DP * tanh( damp * ( xj - x ) )
      fac2 = 0.5_DP + 0.5_DP * tanh( damp * ( x - xj ) )
      forcemin = fac1 * force1(2,i) + fac2 * force2(2,i)
      freemin = freemin + forcemin * dx

      xj = ( jbest ) * dx + xmin
      fac1 = 0.5_DP + 0.5_DP * tanh( damp * ( xj - x ) )
      fac2 = 0.5_DP + 0.5_DP * tanh( damp * ( x - xj ) )
      forcemax = fac1 * force1(2,i) + fac2 * force2(2,i)
      freemax = freemax + forcemax * dx

      write(14,*) x, force, forcemin, forcemax
      write(15,*) x, free, freemin, freemax

   end do

   close(14)
   close(15)

   ! Short output
   write(*,*) "Number of points: ", npoints
   write(*,*) "Minimal difference ", diffmin, " at x = ", ( jbest - 1 ) * dx + xmin


contains

   subroutine interpolate( array, npoints, xmin, xmax , tmode )
      implicit none
      !! Input: 2D array with x and y values, also output after interpolation
      real(kind=DP), intent(inout), allocatable :: array(:,:)
      !! Input: #points of grid
      integer, intent(in) :: npoints
      !! Input: start and end x value
      real(kind=DP), intent(in) :: xmin, xmax
      !! Input: mode, either linear (0) or cubic (1)
      logical, intent(in) :: tmode
      !! Internal: temporary result array
      real(kind=DP) :: temp(3,npoints)
      !! Internal: loop index, indices of 2 closest points, size of array
      integer :: i, i1, i2, arsize
      !! Internal: x and y values of 2 closest points, derivative at 2 closest points
      real(kind=DP) :: x1, x2, y1, y2, k1, k2
      !! Internal: x value, x step, parameters (linear: slope and offset, cubic: more complicated)
      real(kind=DP) :: x, dx, t, a, b
      !! Internal: step size of array, min x of array
      real(kind=DP) :: dx2, xarmin, xarmax

      arsize = size( array(:,:), 2)
      dx = ( xmax - xmin ) / real( npoints - 1, kind=DP )

      do i = 1, npoints
         x = ( i - 1 ) * dx + xmin

         i1 = -1
         i2 = -1
         do j = 1, arsize
            xj = array(1,j)
            if ( x > array(1,j) .and. x < array(1,j+1) ) then
               i1 = j
               i2 = j + 1
               exit
            else if ( abs(x - xj) <= 1.E-8_DP ) then
               i1 = j 
               i2 = j
               exit
            else if ( x > array(1,arsize) ) then
               i1 = arsize-1 
               i2 = arsize
               exit
            end if
         end do

         ! only works for evenly spaced data
         i1 = Max( i1, 1 )
         i2 = Max( i2, 2 )
         i1 = Min( i1, arsize - 1 )
         i2 = Min( i2, arsize )

         x1 = array(1,i1)
         x2 = array(1,i2)
         y1 = array(2,i1)
         y2 = array(2,i2)

         if( i1 == i2 ) then
            temp(1,i) = x
            temp(2,i) = y1

            cycle
         end if

         if( tmode == .false. ) then
            ! linear
            t = ( y2 - y1 ) / ( x2 - x1 )
            b = y1 - t * x1

            temp(2,i) = t * x + b
            ! derivative, doesn't really work
            temp(3,i) = t
         else
            ! cubic
            if( size(array(:,:), 1 ) < 3 ) then
               write(*,*) "Can not execute cubic spline before linear interpolation!"
               error stop
            end if
            k1 = array(3,i1)
            k2 = array(3,i2)

            t = ( x - x1 ) / ( x2 - x1 )
            a = k1 * ( x2 - x1 ) - ( y2 - y1 )
            b = -1.0_DP * k2 * ( x2 - x1 ) + ( y2 - y1 )
            temp(2,i) = ( 1.0_DP - t ) * y1 + t * y2 + t * ( 1.0_DP - t ) * ( a * ( 1.0_DP - t )  + b * t )

            temp(3,i) = ( y2 - y1 ) / ( x2 - x1 ) & 
                    & + ( 1.0_DP - 2.0_DP * t ) * ( ( a * ( 1.0_DP -t ) ) + b * t ) / ( x2 - x1 ) &
                    & + t * ( 1.0_DP - t) * ( b - a ) / ( x2 - x1 )
         end if

         temp(1,i) = x
      end do

      deallocate( array )
      allocate( array(1:3,1:npoints) )

      array(1:3,1:npoints) = temp(1:3,1:npoints)

   end subroutine interpolate

   subroutine sort( array )
      implicit none
      !! Input: 2D array with x and y values, also output after interpolation
      real(kind=DP), intent(inout) :: array(:,:)
      !! Input: start and end x value
      real(kind=DP) :: xmin
      !! Internal: temporary result array
      real(kind=DP) :: temp(2)
      !! Internal: loop index, indices of 2 closest points, size of array
      integer :: i, arsize, imin, imin1(1:1)
      !! Internal: step size of array, min x of array
      real(kind=DP) :: xarmin

      arsize = size( array(:,:), 2)
      xarmin = minval( array(1,:))

      ! Selection sort, because I'm lazy
      do i = 1, arsize-1
         temp(1:2) = array(1:2,i)
         imin1 = minloc(array(1,i:)) + i - 1
         imin = imin1(1)
         array(1:2,i) = array(1:2,imin)
         array(1:2,imin) = temp(1:2)
      end do

   end subroutine sort

end program joinforces
