!! Author:  Christian Schwermann
!! E-mail:  c.schwermann@wwu.de
!! Date:    06/11/2019
!! Project: Joinforces
!! File:    aver_lambda.f90 
!! Copyright: © 2019 Christian Schwermann, ALL RIGHTS RESERVED
!!
!!***********************************************************************
!!
!! Analysis tool, takes the constraint forces from lammps TMD simulations
!! and outputs the average.
!!
!!***********************************************************************
program Aver_lambda
   implicit none
   !! Precision
   integer, parameter :: DP = Selected_real_kind(15, 307)
   !! #arguments, argument number, step number
   integer            :: numargs, ia, istep
   !! errors
   integer            :: ierror, iout
   !! reaction coordinate and average, constraint force and average and sigma^2
   real(kind=DP)      :: rstep, raver, lambda, laver, lsigma2
   !! input and output, dummy and chars for reading
   character(len=64)  :: filnam, outnam, dummy, intchar, rchar
   !! max number of steps, initialized to -1 in case no argument is given
   integer            :: nstep = -1
   !! location of last point in data, to reverse x axis in bw direction
   real(kind=DP)      :: endloc
   !! reverse direction?
   logical            :: treverse = .false.
   !! 

   ! should be input, but lazy
   outnam = "average.dat"

   numargs = Iargc()

   call Getarg( 1, filnam )
   open( unit = 11, file = Trim( filnam ), status = "old", action = "read", iostat = ierror )
   if( ierror /= 0 ) then
      write(*,*) "Error reading file number ", 11, " ", Trim( filnam ), "!"
      error stop
   end if

   if( numargs == 2 ) then
      call Getarg( 2, intchar )
      read(intchar, '(I64)') nstep
   end if
   if( numargs == 3) then
      call Getarg( 3, rchar )
      read(rchar,'(F15.8)') endloc
      treverse = .true.
   end if

   istep = 0
   ! headline
   read(11,*,iostat=ierror)
   ! initialize
   raver = 0.0_DP
   laver = 0.0_DP
   lsigma2 = 0.0_DP
   ! read in TMD output
   do
      istep = istep + 1
      read(11,*,iostat=ierror) dummy, rstep, dummy, dummy, dummy, lambda, dummy, dummy
      raver = raver + rstep
      laver = laver + lambda
      lsigma2 = lsigma2 + lambda ** 2
      if( ierror /= 0 ) exit
      if( istep == nstep ) exit
   end do
      
   raver = raver / Real( istep, kind=DP )
   laver = laver / Real( istep, kind=DP )
   lsigma2 = lsigma2 / Real( istep, kind=DP )

   if( treverse ) then
      ! in reverse direction, forces change sign (derivatives!)
      write(*,*) endloc - raver, - laver, Sqrt( lsigma2 - laver ** 2 )
   else
      write(*,*) raver, laver, Sqrt( lsigma2 - laver ** 2 )
   end if

end program
      
