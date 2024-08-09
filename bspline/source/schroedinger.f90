subroutine schroedinger(npnts,xint,mass,potential,energy,wavefunction)
!=======================================================================
!
!  Name:        schroedinger_1d
!  Description: Solves one-dimensional Schroedinger equation
!               using Fourier Grid Hamiltonian method
!
!  Input parameters:
!
!     npnts - number of grid points (should be power of two)
!     xint  - integration interval in Angstroems
!     mass  - mass of the particle in electron mass units
!     potential(1:npnts) - potential energy on the grid (kcal/mol)
!
!  Output parameters:
!
!     energy(1:npnts) - energy levels (kcal/mol)
!     wavefunction(i=1:npnts,j=1:npnts) - wavefunctions;
!                                         i - grid point;
!                                         j - quantum number.
!
!  Accuracy: integer*4 and real(kind=8)
!
!  Initial version: 11/26/2003
!
!-----------------------------------------------------------------------
!
!  $Author: souda $
!  $Date: 2009/04/12 00:14:37 $
!  $Revision: 1.2 $
!  $Log: schroedinger.f90,v $
!  Revision 1.2  2009/04/12 00:14:37  souda
!  initial version setup
!
!
!=======================================================================
   implicit none

   !-- input/output variables
   integer, intent(in)                          :: npnts
   real(kind=8),  intent(in)                          :: xint
   real(kind=8),  intent(in)                          :: mass
   real(kind=8),  intent(in),  dimension(npnts)       :: potential
   real(kind=8),  intent(out), dimension(npnts)       :: energy
   real(kind=8),  intent(out), dimension(npnts,npnts) :: wavefunction

   !-- parameters
   real(kind=8), parameter :: au2kcal = 627.5095d0
   real(kind=8), parameter :: kcal2au = 1.d0/au2kcal
   real(kind=8), parameter :: bohr2a  = 0.529177249d0
   real(kind=8), parameter :: a2bohr  = 1.d0/bohr2a
   real(kind=8), parameter :: pi      = 3.14159265358979d0
   real(kind=8), parameter :: hbar    = 1.d0

   !-- local variables
   integer :: i, j, ierr
   real(kind=8) :: xint_au
   real(kind=8), allocatable, dimension(:,:) :: hamiltonian    ! Hamiltonian matrix
   real(kind=8), allocatable, dimension(:,:) :: ekinetic       ! Kinetic energy matrix
   real(kind=8), allocatable, dimension(:)   :: work1, work2   ! work arrays

   !-- check whether npnts is an integer power of two

   if ( npnts.le.0 .or. .not.power2(npnts) ) then
      write(*,*) 'number of grid points (',npnts,') is invalid, should be power of 2'
      stop
   endif

   !-- Allocate arrays for the Hamiltonian matrix and kinetic energy matrix

   allocate (hamiltonian(npnts,npnts))
   hamiltonian = 0.d0
   allocate (ekinetic(npnts,npnts))
   ekinetic = 0.d0

   !-- Allocate work arrays
   allocate (work1(npnts), work2(npnts))
   work1 = 0.d0
   work2 = 0.d0

   !-- Calculate kinetic energy matrix (using FFT)
   xint_au = npnts * xint*a2bohr / (npnts-1)
   call gridke(npnts, xint_au, ekinetic)
   ekinetic = ekinetic/mass

   !-- add potential and kinetic energy matrices to form
   !   the full hamiltonian matrix

   do i=1,npnts
      hamiltonian(i,i) = potential(i)*kcal2au
   enddo

   hamiltonian = hamiltonian + ekinetic

   !-- diagonalize Hamiltonian matrix
   call rs(npnts,npnts,hamiltonian,energy,npnts,wavefunction,work1,work2,ierr)

   energy = energy*au2kcal

   if (ierr.ne.0) then
      write(*,*) 'Diagonalization error in schroedinger: exit...'
      call clean_memory
      stop
   endif

   call clean_memory
   return

contains

   !-- Checks whether the given integer (N) is a power of two
   logical function power2(n)

      implicit none
      integer, intent(in) :: n

      integer, parameter :: itwo = 2
      integer :: nw

      power2 = .false.

      nw = n
      if (nw.lt.2) then
         power2 = .false.
         return
      else
         power2 = .true.
      endif

      do while (nw.gt.1)
         if (mod(nw,itwo).ne.0) then
            power2 = .false.
            return
         else
            nw = nw/itwo
         endif
      enddo

      return
   end function power2

   subroutine clean_memory
      implicit none
      if (allocated(hamiltonian)) deallocate(hamiltonian)
      if (allocated(ekinetic   )) deallocate(ekinetic   )
      if (allocated(work1      )) deallocate(work1      )
      if (allocated(work2      )) deallocate(work2      )
   end subroutine clean_memory

   subroutine gridke(npnts_,xint_,hke_)
   !---------------------------------------------------------
   ! This subroutine calculates a kinetic energy matrix
   ! (npnts SHOULD BE AN INTEGER POWER OF TWO !!!)
   !---------------------------------------------------------
      implicit none
      integer, intent(in) :: npnts_
      real(kind=8),  intent(in) :: xint_
      real(kind=8),  intent(out), dimension(npnts_,npnts_) :: hke_

      real(kind=8), parameter :: zero=0.0d+00, one=1.0d+00, two=2.0d+00
      real(kind=8), parameter :: pt5=0.5d+00

      integer :: nn, i, j, k, kk
      real(kind=8) :: delk
      real(kind=8), allocatable :: pmom(:),delf(:)

      allocate (pmom(npnts))
      allocate (delf(2*npnts))

      delk = two*pi/xint_
      nn = npnts_/2
      do i=1,nn
         pmom(i) = (i - 1)*delk*hbar
         pmom(i+nn) = (i - 1 - nn)*delk*hbar
      enddo

      do i=1,npnts_

         do j=1,npnts_
            if (j.eq.i) then
               delf(2*j-1) = one
               delf(2*j)   = zero
            else
               delf(2*j-1) = zero
               delf(2*j)   = zero
            endif
         enddo

         !---- do fourier transform from the position ----
         !     representation to the momentum representation

         call four1(delf,npnts_,1)

         do k=1,2*npnts_
            kk = k/2 + mod(k,2)
            delf(k)=delf(k)*pmom(kk)**2/two
         enddo

         !---- now must do inverse fourier transform from the momentum ----
         !     representation to the postion representation

         call four1(delf,npnts_,-1)

         do j=1,npnts_
            hke_(j,i) = delf(2*j-1)/real(npnts_)
         enddo

      enddo

      deallocate (pmom,delf)
      return

   end subroutine gridke

   SUBROUTINE FOUR1(DATA,NN,ISIGN)
   !---1D FFT routine (Numerical Recipes)

      implicit none
      integer, intent(in) :: nn, isign
      real(kind=8),  intent(inout), dimension(2*nn) :: data

      integer :: n, j, i, m, mmax, istep
      real(kind=8) :: WR, WI, WPR, WPI, WTEMP, THETA, TEMPR, TEMPI

      N = 2*NN
      J = 1

      DO I=1,N,2
         IF (J.GT.I) THEN
            TEMPR = DATA(J)
            TEMPI = DATA(J+1)
            DATA(J) = DATA(I)
            DATA(J+1) = DATA(I+1)
            DATA(I) = TEMPR
            DATA(I+1) = TEMPI
         ENDIF
         M = N/2
         do while ((M.GE.2).AND.(J.GT.M))
            J = J - M
            M = M/2
         enddo
         J = J + M
      enddo

      MMAX = 2
  
      do while (N.GT.MMAX)
         ISTEP = 2*MMAX
         THETA = 6.28318530717959D0/(ISIGN*MMAX)
         WPR = -2.D0*DSIN(0.5D0*THETA)**2
         WPI = DSIN(THETA)
         WR = 1.D0
         WI = 0.D0
         DO M=1,MMAX,2
            DO I=M,N,ISTEP
               J = I + MMAX
               TEMPR = WR*DATA(J) - WI*DATA(J+1)
               TEMPI = WR*DATA(J+1) + WI*DATA(J)
               DATA(J) = DATA(I) - TEMPR
               DATA(J+1) = DATA(I+1) - TEMPI
               DATA(I) = DATA(I) + TEMPR
               DATA(I+1) = DATA(I+1) + TEMPI
            enddo
            WTEMP = WR
            WR = WR*WPR - WI*WPI + WR
            WI = WI*WPR + WTEMP*WPI + WI
         enddo
         MMAX = ISTEP
      enddo

      RETURN
 
   end subroutine four1

end subroutine schroedinger
