module timers

   implicit none
   public :: second

contains

   function second() result(tic)
      real(4), dimension(2) :: time
      real(4) :: etime
      real(8) :: tic
      tic = dble(etime(time))
   end function second

end module timers
