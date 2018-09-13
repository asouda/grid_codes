program fgh_bspline_q

!!------------------------------------------------------------------
!!
!! calculates vibrational energies and wavefunctions (1D)
!! using the b-spline for the given potentials
!! - using common shift for both potentials as given in the input
!! - calculates quadratic overlap attenuation parameters as well
!!
!!------------------------------------------------------------------
!!
!! $Author: souda $
!! $Date: 2015/04/30 15:45:00 $
!! $Revision: 1.2 $
!! $Log: fgh_bspline_abs.f90,v $
!!
!!------------------------------------------------------------------

   use bspline
   implicit none

   interface
      subroutine schroedinger(npnts,xint,mass,potential,energy,wavefunction)
      implicit none
      integer, intent(in)                          :: npnts
      real(8),  intent(in)                          :: xint
      real(8),  intent(in)                          :: mass
      real(8),  intent(in),  dimension(npnts)       :: potential
      real(8),  intent(out), dimension(npnts)       :: energy
      real(8),  intent(out), dimension(npnts,npnts) :: wavefunction
      end subroutine schroedinger
   end interface

   real(8), parameter :: ev2cm  = 8065.54477d0
   real(8), parameter :: au2ev  =   27.2113834d0
   real(8), parameter :: h2kcal =  627.5095d0              ! 1 Hartree = 627.5095 kcal/mol

   !-- b-spline order parameter
   integer, parameter :: korder=3

   integer :: n, npntspow2, i, j, k, npoints, ncoef, idot, nstates
   integer :: iargc, narg
   real(8) :: x_left, x_right, dx, xh, m, phase
   real(8) :: pot_min, xt, bt
   real(8) :: tstart, tfinish, xp, e0, e1, sij, sijup, sijdown
   real(8), allocatable, dimension(:)   :: x_data, reactant_pot_data, product_pot_data
   real(8), allocatable, dimension(:)   :: reactant_pot, product_pot, x
   real(8), allocatable, dimension(:)   :: reactant_en, product_en
   real(8), allocatable, dimension(:,:) :: reactant_wavef, product_wavef, overlap, alpha, qalpha

   ! b-spline arrays
   real(8), dimension(:), allocatable :: xknot, bscoef

   character(len=240) :: argument, data_file

   narg = iargc()
   if (narg.eq.0) then
      !write(*,*)
      !write(*,*) "=================================================================================================="
      !write(*,*) " Usage:"
      !write(*,*) " fgh_bspline_q.bin <name_of_the_file_with_the_potential> <n_grid_points> <particle_mass> <nstates>"
      !write(*,*) "=================================================================================================="
      stop
   endif

   !-- read the table with the potential data

   call getarg(1,argument)
   data_file = trim(argument)
   open(1,file=data_file,status='old')

   npoints = 0
   do
      read(1,*,end=11) xp, e0, e1
      npoints = npoints + 1
      if (npoints.eq.1) x_left = xp
   enddo
11 continue
   x_right = xp
   close(1)

   if (x_right.le.x_left) then
      write(*,*)
      write(*,*) "x_right should be greater than x_left... Sorry..."
      stop
   endif

   xh = x_right - x_left
   write(*,'( 1x,"Left  integration limit (Angstroems): ",f12.6)') x_left
   write(*,'( 1x,"Right integration limit (Angstroems): ",f12.6)') x_right
   write(*,'( 1x,"Length of the integration interval (Angstroems): ",f12.6)') xh

   !-- n - number of grid points (should be power of two)
   call getarg(2,argument)
   read(argument,*) n

   if (.not.power2(n)) then
      write(*,'(/1x,"*** The number of grid points (",I3,") is not a power of two.")') n
      npntspow2 = 1
      do while (n.gt.npntspow2.and.npntspow2.lt.4096)
         npntspow2 = npntspow2 * 2
      enddo
      n = 2 * npntspow2
      write(*,'(" Number of grid points is reset to ",i3/)') n
   endif
   write(*,'( 1x,"Number of grid points: ",i10)') n

   !-- Mass of the quantum particle (electronic mass units)

   call getarg(3,argument)
   read(argument,*) m
   write(*,'( 1x,"Mass of the quantum particle: ",f15.6," electronic mass units")') m

   !-- nstates - number of vibrational states to calculate (default = 10)
   call getarg(4,argument)
   read(argument,*) nstates
   if (nstates.le.0.or.nstates.gt.n) nstates = 10
   write(*,'( 1x,"Number of states to output: ",i10)') nstates

   !-- calculate potentials on the grid using the spline
   !   of the tabulated potentials

   !-- read the data points

   !-- allocate arrays
   allocate (x_data(npoints))
   allocate (reactant_pot_data(npoints))
   allocate (product_pot_data(npoints))

   open(1,file=data_file,status='old')
   do i=1,npoints
      read(1,*) x_data(i), reactant_pot_data(i), product_pot_data(i)
   enddo
   close(1)

   !-- allocate b-spline arrays
   allocate (bscoef(npoints), xknot(npoints+korder))

   !-- generate knot sequence
   call dbsnak(npoints, x_data, korder, xknot)

   allocate (x(n))
   dx = xh/(n-1)
   x(1) = x_data(1)
   do i=2,n-1
      x(i) = x(i-1) + dx
   enddo
   x(n) = x_data(npoints)

   !-- generate spline coefficients for reactant potential

   bscoef = 0.d0
   call dbsint (npoints, x_data, reactant_pot_data, korder, xknot, bscoef)

   !-- evaluate b-spline at the grid points

   allocate (reactant_pot(n), reactant_en(n), reactant_wavef(n,n))
   reactant_pot(1) = reactant_pot_data(1)
   pot_min = 999999.d0
   do i=2,n-1
      reactant_pot(i) = dbsval(x(i),korder,xknot,npoints,bscoef)
      if (reactant_pot(i).le.pot_min) pot_min = reactant_pot(i)
   enddo
   reactant_pot(n) = reactant_pot_data(npoints)

   !-- generate spline coefficients for product potential

   bscoef = 0.d0
   call dbsint (npoints, x_data, product_pot_data, korder, xknot, bscoef)

   !-- evaluate b-spline at the grid points

   allocate (product_pot(n), product_en(n), product_wavef(n,n))
   product_pot(1) = product_pot_data(1)
   do i=2,n-1
      product_pot(i) = dbsval(x(i),korder,xknot,npoints,bscoef)
      if (product_pot(i).le.pot_min) pot_min = product_pot(i)
   enddo
   product_pot(n) = product_pot_data(npoints)

   !-- adjust the lowest potential point to zero
   reactant_pot = (reactant_pot-pot_min)*h2kcal
   product_pot = (product_pot-pot_min)*h2kcal

   !-- write out the splined potentials

   idot = index(data_file,".")
   if (idot.gt.0) data_file = data_file(1:idot-1)

   open(unit=1,file=trim(data_file)//'_bspline.dat')
   do i=1,n
      write(1,'(3g20.10)') x(i), reactant_pot(i), product_pot(i)
   enddo
   close(1)

   open(unit=1,file=trim(data_file)//'_data.dat')
   do i=1,npoints
      write(1,'(3g20.10)') x_data(i), reactant_pot_data(i)*h2kcal, product_pot_data(i)*h2kcal
   enddo
   close(1)

   !-- calculate energy levels and wavefunctions

   call schroedinger(n,xh,m,reactant_pot,reactant_en,reactant_wavef)
   call schroedinger(n,xh,m,product_pot,product_en,product_wavef)
   write(*,*)
   write(*,*) 'Reactant eigenvalues (kcal/mol):'
   write(*,*) (reactant_en(k),k=1,nstates)
   write(*,*)
   write(*,*) 'Product eigenvalues (kcal/mol):'
   write(*,*) (product_en(k),k=1,nstates)

   !-- adjust phases of the wavefunctions

   phase = sign(1.d0,reactant_wavef(n/2,1))
   reactant_wavef = reactant_wavef*phase
   phase = sign(1.d0,product_wavef(n/2,1))
   product_wavef = product_wavef*phase

   !-- output the energies

   open(unit=1,file=trim(data_file)//"_reactant_en.dat")
   write(1,'(10g20.10)') (reactant_en(k),k=1,nstates)
   close(1)

   open(unit=1,file=trim(data_file)//"_product_en.dat")
   write(1,'(10g20.10)') (product_en(k),k=1,nstates)
   close(1)

   !-- output the wavefunctions

   open(unit=1,file=trim(data_file)//"_reactant_wf.dat")
   do i=1,n
      write(1,'(21g20.10)') x(i),(reactant_wavef(i,k),k=1,nstates)
   enddo
   close(1)

   open(unit=1,file=trim(data_file)//"_product_wf.dat")
   do i=1,n
      write(1,'(21g20.10)') x(i),(product_wavef(i,k),k=1,nstates)
   enddo
   close(1)

   !-- Overlap integrals and alpha parameters

   allocate (overlap(nstates,nstates))
   allocate (alpha(nstates,nstates))
   allocate (qalpha(nstates,nstates))
   overlap = 0.d0
   alpha = 0.d0
   qalpha = 0.d0

   do i=1,nstates
      do j=1,nstates

         sij = 0.d0
         sijup = 0.d0
         sijdown = 0.d0

         do k=1,n
            sij = sij + reactant_wavef(k,i)*product_wavef(k,j)
         enddo

         do k=2,n
            sijup = sijup + reactant_wavef(k,i)*product_wavef(k-1,j)
         enddo

         do k=1,n-1
            sijdown = sijdown + reactant_wavef(k,i)*product_wavef(k+1,j)
         enddo

         overlap(i,j) = sij
         alpha(i,j) = 0.5d0*(sijdown - sijup)/dx/sij

         qalpha(i,j) = -(sijup - 2.d0*sij + sijdown)/dx/dx/sij + &
                     &  (sijup-sijdown)*(sijup-sijdown)/4.d0/dx/dx/sij/sij

      enddo
   enddo

   open(unit=1,file=trim(data_file)//"_overlaps.dat")
   do i=1,nstates
      write(1,'(40g20.10)') (overlap(i,k),k=1,nstates)
   enddo
   close(1)

   open(unit=1,file=trim(data_file)//"_alphas.dat")
   do i=1,nstates
      write(1,'(40g20.10)') (alpha(i,k),k=1,nstates)
   enddo
   close(1)

   open(unit=1,file=trim(data_file)//"_qalphas.dat")
   do i=1,nstates
      write(1,'(40g20.10)') (qalpha(i,k),k=1,nstates)
   enddo
   close(1)

   !-- deallocate arrays
   deallocate (x_data, x)
   deallocate (bscoef, xknot)
   deallocate (reactant_pot_data, reactant_pot, reactant_en, reactant_wavef)
   deallocate (product_pot_data,  product_pot,  product_en,  product_wavef)
   deallocate (overlap, alpha, qalpha)

CONTAINS

   logical function power2(n)
   !-- Checks whether the given integer (N) is a power of two

      implicit none
      integer, intent(in) :: n

      integer, parameter :: itwo = 2
      integer :: nw

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

end program fgh_bspline_q
