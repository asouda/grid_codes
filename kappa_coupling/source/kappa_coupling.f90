program kappa_coupling

!!
!! Calculates the tunneling splittings, prefactors kappa in the Stuchebrukhov's expression,
!! vibronic couplings and overlap integrals for each pair of reactant and product
!! diabatic vibronic states calculated for given reactant and product proton potentials.
!!
!! [Eq. (1.7) in JCP 113, 10438 (2000)]
!!
!! diabatic profiles and electronic couplings are read from disk
!!

   use timers
   use bspline
   implicit none

   interface
      subroutine schroedinger(npnts,xint,mass,potential,energy,wavefunction)
      implicit none
      integer, intent(in)                          :: npnts
      real(kind=8),  intent(in)                          :: xint
      real(kind=8),  intent(in)                          :: mass
      real(kind=8),  intent(in),  dimension(npnts)       :: potential
      real(kind=8),  intent(out), dimension(npnts)       :: energy
      real(kind=8),  intent(out), dimension(npnts,npnts) :: wavefunction
      end subroutine schroedinger
   end interface

   real(kind=8), parameter :: ev2cm   = 8065.54477d0
   real(kind=8), parameter :: au2ev   =   27.2113834d0
   real(kind=8), parameter :: au2kcal =  627.5095d0
   real(kind=8), parameter :: kcal2au =  1.d0/627.5095d0
   real(kind=8), parameter :: bohr2a  = 0.529177249d0
   real(kind=8), parameter :: a2bohr  = 1.d0/bohr2a
   real(kind=8), parameter :: au2ps   = 2.418884326505d-5
   real(kind=8), parameter :: ps2au   = 1.d0/au2ps
   real(kind=8), parameter :: au2fs   = 2.418884326505d-2
   real(kind=8), parameter :: fs2au   = 1.d0/au2fs
   real(kind=8), parameter :: dalton  = 1822.8900409014022D0
   real(kind=8), parameter :: pi = 3.1415926535897932384626434D0

   !-- b-spline order parameter
   integer, parameter :: korder=3

   integer :: n, i, j, k, im, nstates, npoints, npntspow2, ixc, istate, jstate
   integer :: iargc, narg, filename_offset, split_level_plus, split_level_minus
   real(kind=8) :: x_left, x_right, dx, xh, mass, fcij, phase, dpot_left, dpot_right
   real(kind=8) :: xp, h11, h22, h12, v12x
   real(kind=8) :: sh, dh, shift_l, shift_r
   real(kind=8) :: viblevel1, viblevel2
   real(kind=8) :: tstart, tfinish
   real(kind=8), allocatable, dimension(:)   :: x_data, e0_data, e1_data, h11_data, h22_data, h12_data, x
   real(kind=8), allocatable, dimension(:)   :: pot_l, pot_r, xpot_l, xpot_r, pot_ground, pot_excited
   real(kind=8), allocatable, dimension(:)   :: el_coupling, en_l, en_r, en_ground, pot_diff
   real(kind=8), allocatable, dimension(:)   :: wavef_plus, wavef_minus, overlaps_plus, overlaps_minus
   real(kind=8), allocatable, dimension(:,:) :: wavef_l, wavef_r, wavef_ground, vib_overlap, vib_coupling_diab
   real(kind=8), allocatable, dimension(:,:) :: p_adiabaticity, tau_p, tau_e, kappa
   real(kind=8), allocatable, dimension(:,:) :: vib_coupling_kappa, vib_coupling_product, tunn_splitting

   real(kind=8) :: dgamma
   real(kind=8) :: xc, tunn_energy, df, vt, p

   character(len=80) :: argument, data_file
   character(len=40) :: particle, n_char
   character(len=5)  :: ovl_symb
   character(len=10) :: vib_symb

   !-- b-spline arrays
   real(kind=8), dimension(:), allocatable :: xknot, bscoef

   narg = iargc()
   if (narg.eq.0) then
      write(*,*)
      write(*,*) "=================================================================================================="
      write(*,*) " Usage:"
      write(*,*)
      write(*,*) "   xxx.bin <name_of_the_file_with_PT profiles> <n_grid_points> <particle_mass> <number_of_states>"
      write(*,*)
      write(*,*) "   Columns in file_with_PT_profiles (diabatic energies and couplings in kcal/mol):"
      write(*,*) "   x(A)   H11   H22   V12"
      write(*,*)
      write(*,*) "   Number of grid points should be power of two (default is 1024)"
      write(*,*) "   particle_mass should be in electron mass units (default is mass of a proton: 1836.15 emu)"
      write(*,*) "   default number of states in each diabatic potential is 1"
      write(*,*) "=================================================================================================="
      stop
   endif

   call getarg(1,argument)
   data_file = trim(argument)
   open(1,file=data_file,status='old')

   filename_offset = len(data_file)
   do i=len(data_file),1,-1
      if (data_file(i:i).eq.".") then
         filename_offset = i - 1
         exit
      endif
   enddo

   !-- scan input file and determine left and right integration limits

   npoints = 0
   x_left = 999.d0
   x_right = -999.d0
   do
      read(1,*,end=11) xp, h11, h22, h12
      npoints = npoints + 1
      if (xp.le.x_left) x_left = xp
      if (xp.ge.x_right) x_right = xp
   enddo
11 continue
   rewind 1

   !-- n - number of data points
   xh = x_right - x_left
   write(*,'(/1x,"Left  integration limit:            ",f12.6," A")') x_left
   write(*,'( 1x,"Right integration limit:            ",f12.6," A")') x_right
   write(*,'( 1x,"Length of the integration interval: ",f12.6," A")') xh
   write(*,'( 1x,"Number of data points:              ",i12)') npoints

   !-- n - number of grid points (should be power of two!)

   if (narg.gt.1) then
      call getarg(2,argument)
      read(argument,*) n
   else
      n = 1024
   endif

   if (.not.power2(n)) then
      write(*,'(/1x,"*** The number of grid points (",I3,") is not a power of two.")') n
      npntspow2 = 1
      do while (n.gt.npntspow2.and.npntspow2.lt.4096)
         npntspow2 = npntspow2 * 2
      enddo
      n = 2 * npntspow2
      write(*,'(" Number of grid points is reset to ",i3/)') n
   endif

   write(*,'(/1x,"Number of grid points for splined potentials: ",i10)') n
   write(n_char,'(i4.4,a2)') n, "pt"

   !-- Mass of the quantum particle (electronic mass units)

   if (narg.gt.2) then
      call getarg(3,argument)
      read(argument,*) mass
   else
      mass = 1836.152701d0
   endif

   !select case(im)
   !   case(0)
   !      mass = 1.d0
   !      particle = "electron"//"_"//trim(n_char)
   !   case(1)
   !      mass = 1836.152701d0
   !      particle = "proton"//"_"//trim(n_char)
   !   case(2)
   !      mass = 3670.482955d0
   !      particle = "deuteron"//"_"//trim(n_char)
   !   case(3)
   !      mass = 5497.9210703d0
   !      particle = "triton"//"_"//trim(n_char)
   !   case default
   !      write(*,'(/1x,"You must specify the particle kind! ")')
   !      stop "Come again..."
   !end select

   write(*,'(/1x,"Mass of the quantum particle: ",f15.6," electronic mass units")') mass
   write(*,'( 1x,"                              ",f15.6," atomic mass units")') mass/dalton

   !-- nstates - number of states of interest in the left and right potentials

   if (narg.gt.3) then
      call getarg(4,argument)
      read(argument,*) nstates
   else
      nstates = 1
   endif

   if (nstates.le.0) nstates = 1

   write(*,'(/1x,"Number of vibrational states in each diabatic potential: ",i5)') nstates

   !-- allocate arrays
   allocate (x_data(npoints))
   allocate (h11_data(npoints))
   allocate (h22_data(npoints))
   allocate (h12_data(npoints))

   allocate (x(n))
   allocate (pot_l(n), en_l(n), wavef_l(n,n))
   allocate (pot_r(n), en_r(n), wavef_r(n,n))
   allocate (el_coupling(n))

   !-- read data points
   do i=1,npoints
      read(1,*) x_data(i), h11_data(i), h22_data(i), h12_data(i)
   enddo
   close(1)

   !-- allocate b-spline arrays
   allocate (bscoef(npoints), xknot(npoints+korder))

   !-- generate knot sequence
   call dbsnak(npoints, x_data, korder, xknot)

   dx = xh/(n-1)
   x(1) = x_data(1)
   do i=2,n-1
      x(i) = x(i-1) + dx
   enddo
   x(n) = x_data(npoints)

   !-- generate spline coefficients for reactant potential
   bscoef = 0.d0
   call dbsint (npoints, x_data, h11_data, korder, xknot, bscoef)

   !-- evaluate b-spline at the grid points
   pot_l(1) = h11_data(1)
   do i=2,n-1
      pot_l(i) = dbsval(x(i),korder,xknot,npoints,bscoef)
   enddo
   pot_l(n) = h11_data(npoints)

   !-- generate spline coefficients for product potential
   bscoef = 0.d0
   call dbsint (npoints, x_data, h22_data, korder, xknot, bscoef)

   !-- evaluate b-spline at the grid points
   pot_r(1) = h22_data(1)
   do i=2,n-1
      pot_r(i) = dbsval(x(i),korder,xknot,npoints,bscoef)
   enddo
   pot_r(n) = h22_data(npoints)

   !-- generate spline coefficients for electronic coupling
   bscoef = 0.d0
   call dbsint (npoints, x_data, h12_data, korder, xknot, bscoef)

   !-- evaluate b-spline at the grid points
   el_coupling(1) = h12_data(1)
   do i=2,n-1
      el_coupling(i) = dbsval(x(i),korder,xknot,npoints,bscoef)
   enddo
   el_coupling(n) = h12_data(npoints)

   !-- calculate energy levels and wavefunctions
   !   in the diabatic electronic potentials

   !-- left diabatic potential (reactant)

   tstart = second()
   call schroedinger(n,xh,mass,pot_l,en_l,wavef_l)
   tfinish = second()
   write(*,*)
   write(*,*) ' Left potential: time elapsed --> ',tfinish-tstart,' seconds'
   write(*,*)
   write(*,*) 'Eigenvalues (kcal/mol):'
   write(*,*) (en_l(k),k=1,nstates)

   !-- adjust phase
   phase = sign(1.d0,wavef_l(n/2,1))
   wavef_l = wavef_l*phase


   !-- right diabatic potential (product)

   tstart = second()
   call schroedinger(n,xh,mass,pot_r,en_r,wavef_r)
   tfinish = second()
   write(*,*)
   write(*,*) 'Right potential: time elapsed --> ',tfinish-tstart,' seconds'
   write(*,*)
   write(*,*) 'Eigenvalues (kcal/mol):'
   write(*,*) (en_r(k),k=1,nstates)

   !-- adjust phase
   phase = sign(1.d0,wavef_r(n/2,1))
   wavef_r = wavef_r*phase

   !-- shift the diabatic potentials and
   !   the vibrational levels so that
   !   the ground vibrational states in the left
   !   and right potentials are degenerate
   !   and the ground state vibrational
   !   levels are zero

   shift_l = en_l(1)
   shift_r = en_r(1)

   do i=1,n
      pot_l(i) = pot_l(i) - shift_l
      en_l(i)  = en_l(i)  - shift_l
   enddo

   do i=1,n
      pot_r(i) = pot_r(i) - shift_r
      en_r(i)  = en_r(i)  - shift_r
   enddo

   !-- output the vibrational energy levels

   open(unit=1,file=data_file(1:filename_offset)//"_reactant_en.dat")
   write(1,'(30g20.10)') (en_l(k),k=1,nstates)
   close(1)

   open(unit=1,file=data_file(1:filename_offset)//"_product_en.dat")
   write(1,'(30g20.10)') (en_r(k),k=1,nstates)
   close(1)


   !-- output the diabatic potentials and vibrational wavefunctions

   open(unit=1,file=data_file(1:filename_offset)//"_left.dat")
   do i=1,n
      write(1,'(60g20.10)') x(i),pot_l(i),(en_l(k),en_l(k)+wavef_l(i,k)*30.d0,k=1,nstates)
   enddo
   close(1)

   open(unit=1,file=data_file(1:filename_offset)//"_right.dat")
   do i=1,n
      write(1,'(60g20.10)') x(i), pot_r(i), (en_r(k),en_r(k) + wavef_r(i,k)*30.d0,k=1,nstates)
   enddo
   close(1)

   !-- calculate Franck-Condon factors (overlaps)

   allocate (vib_overlap(nstates,nstates))

   do i=1,nstates
      do j=1,nstates
         fcij = 0.d0
         do k=1,n
            fcij = fcij + wavef_l(k,i)*wavef_r(k,j)
         enddo
         vib_overlap(i,j) = fcij
      enddo
   enddo

   open(unit=2,file=data_file(1:filename_offset)//"_overlaps.dat")
   do i=1,nstates
      write(2,'(40g20.10)') (vib_overlap(i,k),k=1,nstates)
   enddo
   close(2)

   !-- calculate vibronic couplings between locally diabatic vibronic states
   !   taking into account x-dependent electronic coupling: <mu|Vel(x)|nu>

   allocate (vib_coupling_diab(nstates,nstates))

   do i=1,nstates
      do j=1,nstates
         fcij = 0.d0
         do k=1,n
            fcij = fcij + wavef_l(k,i)*el_coupling(k)*wavef_r(k,j)
         enddo
         vib_coupling_diab(i,j) = fcij
      enddo
   enddo

   open(unit=2,file=data_file(1:filename_offset)//"_vibcouplings_diab.dat")
   do i=1,nstates
      write(2,'(40g20.10)') (vib_coupling_diab(i,k),k=1,nstates)
   enddo
   close(2)


   !-- calculate tunneling splittings and semiclassical vibronic couplings
   !   for each pair of vibronic states

   allocate (tunn_splitting(nstates,nstates))
   allocate (p_adiabaticity(nstates,nstates))
   allocate (tau_p(nstates,nstates))
   allocate (tau_e(nstates,nstates))
   allocate (kappa(nstates,nstates))
   allocate (vib_coupling_kappa(nstates,nstates))
   allocate (vib_coupling_product(nstates,nstates))

   allocate (xpot_l(n), xpot_r(n), pot_diff(n))
   allocate (pot_ground(n), pot_excited(n))
   allocate (en_ground(n), wavef_ground(n,n))
   allocate (wavef_plus(n), wavef_minus(n), overlaps_plus(2*nstates), overlaps_minus(2*nstates))

   write(*,*)
   write(*,'("#",160("="))')
   write(*,*) "# SEMICLASSICAL PARAMETERS"
   write(*,*) "# i - reactant state index (starting from 0)"
   write(*,*) "# j - product  state index (starting from 0)"
   write(*,*) "# i_ad - tunneling state index (starting from 0)"
   write(*,*) "# j_ad - tunneling  state index (starting from 0)"
   write(*,*) "# Sij - vibrational overlap integral"
   write(*,*) "# V^el - electronic coupling at the crossing point (kcal/mol)"
   write(*,*) "# tau_e - electronic transition time (fs)"
   write(*,*) "# tau_p - proton tunneling time (fs)"
   write(*,*) "# p = tau_e/tau_p - adiabaticity parameter"
   write(*,*) "# kappa - Stuchebrukhov prefactor"
   write(*,*) "# Delta - tunneling splitting in ground state adiabatic potential (kcal/mol)"
   write(*,*) "# V^el*<i|j> - vibronic coupling in electronically nonadiabatic limit (kcal/mol)"
   write(*,*) "# V^vib = kappa*Delta - semiclassical vibronic coupling (kcal/mol)"
   write(*,*) "# Delta/2 - vibronic coupling in electronically adiabatic limit"
   write(*,'("#",168("="))')
   write(*,'("#",t4,"i",t8,"j",t12,"i_a",t16,"j_a",t22,"<i|j>",t38,"V^el",t51,"tau_e",t67,"tau_p",t84,"p",t96,&
   &"kappa",t112,"Delta",t125,"V^el*|<i|j>|",t142,"V^(sc)",t156,"Delta/2")')
   write(*,'("#",168("-"))')

   do istate=1,nstates
      do jstate=1,nstates

      !-- istate - reactant state in the left diabatic potential
      !-- jstate - product state in the right diabatic potential

      !-- shift diabatic potentials to align the vibrational levels
      !   for a given pair of vibrational states

      do i=1,n
         xpot_l(i) = pot_l(i) - en_l(istate)
         xpot_r(i) = pot_r(i) - en_r(jstate)
      enddo

      !-- tunneling energy is zero by construction
      !   (vibrational levels are aligned at zero)
      tunn_energy = 0.d0


      !-- find the crossing point where pot_l(xc) = pot_r(xc)
      pot_diff = xpot_l - xpot_r
      do i=2,n
         if (pot_diff(i)*pot_diff(i-1).le.0.d0) then
            ixc = i
            exit
         endif
      enddo

      !-- electronic coupling at the crossing point
      v12x  = 0.5d0*(el_coupling(ixc)  + el_coupling(ixc-1))*kcal2au

      !-- the product of the electronic coupling at the crossing point and the overlap
      !   is the vibronic coupling in the electronically nonadiabatic limit

      vib_coupling_product(istate,jstate) = v12x*au2kcal*abs(vib_overlap(istate,jstate))

      !-- check if the tunneling energy is below the barrier

      if (xpot_l(ixc).le.0.d0) then
         write(*,*) "--------------- tunneling energy is above the barrier ------------"
         p_adiabaticity(istate,jstate) = 0.d0
         tau_p(istate,jstate) = 0.d0
         tau_e(istate,jstate) = 0.d0
         kappa(istate,jstate) = 0.d0
         vib_coupling_kappa(istate,jstate) = 0.d0
         tunn_splitting(istate,jstate) = 0.d0
         cycle
      endif

      !-- calculate symmetric (plus) and antisymmetric (minus) combinations
      !   of the reactant and product diabatic states

      do k=1,n
         wavef_plus(k)  = (wavef_l(k,istate) + wavef_r(k,jstate))/sqrt(2.d0)
         wavef_minus(k) = (wavef_l(k,istate) - wavef_r(k,jstate))/sqrt(2.d0)
      enddo

      !-- diagonalize 2x2 electronic Hamiltonian matrix
      !   to obtain ground and excited state adiabatic potentials

      pot_ground = 0.d0
      pot_excited = 0.d0
      do i=1,n
         h11 = xpot_l(i)
         h22 = xpot_r(i)
         sh = h11 + h22
         dh = h22 - h11
         h12 = el_coupling(i)
         pot_ground(i)  = 0.5d0*sh - 0.5d0*sqrt(dh*dh + 4.d0*h12*h12)
         pot_excited(i) = 0.5d0*sh + 0.5d0*sqrt(dh*dh + 4.d0*h12*h12)
      enddo

      !== calculate the p-, tau_p, tau_e, and kappa from Stuchebrukhov's expression

      !-- calculate the difference of the slopes of the
      !   diabatic profiles at the crossing point (pot_l = pot_r)
      !   (DeltaF in the paper)

      dpot_left  = (xpot_l(ixc) - xpot_l(ixc-1))/dx
      dpot_right = (xpot_r(ixc) - xpot_r(ixc-1))/dx

      df = abs(dpot_left - dpot_right)    ! in kcal/mol/A
      df = df*kcal2au/a2bohr              ! in a.u./Bohr

      !-- calculate tunneling velocity (Eq. 1.9)
      vt = sqrt(2.d0*(xpot_l(ixc)*kcal2au - tunn_energy)/mass)

      !-- adiabaticity parameter p (Eq. 1.8)
      p_adiabaticity(istate,jstate) = v12x*v12x/df/vt

      !-- tunneling times for electron and proton (Eqs. 4.12-4.13)
      tau_p(istate,jstate) = v12x/df/vt
      tau_e(istate,jstate) = 1.d0/v12x

      !-- kappa (Eq. 1.7)
      p = p_adiabaticity(istate,jstate)
      kappa(istate,jstate) = sqrt(2.d0*pi*p)*exp(p*log(p)-p)/dgamma(p+1.d0)

      !-- calculate tunneling splitting in the ground state adiabatic potential

      !-- calculate vibrational wavefunctions and energy levels in the ground state adiabatic potential
      en_ground = 0.d0
      wavef_ground = 0.d0
      call schroedinger(n,xh,mass,pot_ground,en_ground,wavef_ground)

      !-- find two states with the largest overlaps with symmetric and antisymmetric
      !   combinations of the reactant and product diabatic states;
      !   splitting between these two levels is the tunneling splitting
      !   for a given tunneling energy

      !-- absolute values of overlaps
      do i=1,2*nstates
         overlaps_plus(i)  = 0.d0
         overlaps_minus(i) = 0.d0
         do k=1,n
            overlaps_plus(i)  = overlaps_plus(i)  + wavef_plus(k) *wavef_ground(k,i)
            overlaps_minus(i) = overlaps_minus(i) + wavef_minus(k)*wavef_ground(k,i)
         enddo
      enddo
      overlaps_plus  = abs(overlaps_plus)
      overlaps_minus = abs(overlaps_minus)

      !-- find the indices of the largest overlaps
      split_level_plus  = 1
      split_level_minus = 1
      do i=1,2*nstates
         if (overlaps_plus(i) .gt.overlaps_plus(split_level_plus))  split_level_plus  = i
         if (overlaps_minus(i).gt.overlaps_minus(split_level_minus)) split_level_minus = i
      enddo

      !-- output the adiabatic potentials for a pair of diabatic potentials
      !   shifted to align ground state vibrational levels

      if (istate.eq.1.and.jstate.eq.1) then

         open(unit=1,file=data_file(1:filename_offset)//"_adiabatic.dat")
         do k=1,n
            write(1,'(60g20.10)') x(k), pot_ground(k), pot_excited(k)
         enddo
         close(1)

         open(unit=1,file=data_file(1:filename_offset)//"_plus-minus.dat")
         do k=1,n
            write(1,'(60g20.10)') x(k), wavef_plus(k), wavef_minus(k), &
            & wavef_ground(k,split_level_plus), wavef_ground(k,split_level_minus)
         enddo
         close(1)

      endif

      tunn_splitting(istate,jstate) = abs(en_ground(split_level_minus) - en_ground(split_level_plus))

      !-- half of this tunneling splitting times kappa
      !   is the semiclassical vibronic coupling

      vib_coupling_kappa(istate,jstate) = 0.5d0*kappa(istate,jstate)*tunn_splitting(istate,jstate)

      !-- output the tables

      write(*,'(4i4,10g15.6)')&
      & istate-1, jstate-1, split_level_plus-1, split_level_minus-1, vib_overlap(istate,jstate), v12x*au2kcal,&
      & tau_e(istate,jstate)*au2fs, tau_p(istate,jstate)*au2fs,&
      & p_adiabaticity(istate,jstate), kappa(istate,jstate), tunn_splitting(istate,jstate),&
      & vib_coupling_product(istate,jstate), vib_coupling_kappa(istate,jstate), 0.5d0*tunn_splitting(istate,jstate)

      enddo
   enddo

   write(*,'("#",168("="))')

   open(unit=2,file=data_file(1:filename_offset)//"_vibcouplings_product.dat")
   do i=1,nstates
      write(2,'(40g20.10)') (vib_coupling_product(i,k),k=1,nstates)
   enddo
   close(2)

   open(unit=2,file=data_file(1:filename_offset)//"_vibcouplings_kappa.dat")
   do i=1,nstates
      write(2,'(40g20.10)') (vib_coupling_kappa(i,k),k=1,nstates)
   enddo
   close(2)

   open(unit=2,file=data_file(1:filename_offset)//"_tunneling_splittings.dat")
   do i=1,nstates
      write(2,'(40g20.10)') (tunn_splitting(i,k),k=1,nstates)
   enddo
   close(2)

   open(unit=2,file=data_file(1:filename_offset)//"_tau_e.dat")
   do i=1,nstates
      write(2,'(40g20.10)') (tau_e(i,k)*au2fs,k=1,nstates)
   enddo
   close(2)

   open(unit=2,file=data_file(1:filename_offset)//"_tau_p.dat")
   do i=1,nstates
      write(2,'(40g20.10)') (tau_p(i,k)*au2fs,k=1,nstates)
   enddo
   close(2)

   open(unit=2,file=data_file(1:filename_offset)//"_p_adiabaticity.dat")
   do i=1,nstates
      write(2,'(40g20.10)') (p_adiabaticity(i,k),k=1,nstates)
   enddo
   close(2)

   open(unit=2,file=data_file(1:filename_offset)//"_kappa.dat")
   do i=1,nstates
      write(2,'(40g20.10)') (kappa(i,k),k=1,nstates)
   enddo
   close(2)

   !-- deallocate arrays

   deallocate (xpot_l, xpot_r, pot_diff)
   deallocate (pot_ground, pot_excited)

   deallocate (tunn_splitting)
   deallocate (p_adiabaticity)
   deallocate (tau_p, tau_e)
   deallocate (kappa)
   deallocate (vib_coupling_kappa)
   deallocate (vib_coupling_product)
   deallocate (en_ground, wavef_ground)
   deallocate (wavef_plus, wavef_minus, overlaps_plus, overlaps_minus)

   deallocate (x, pot_l, en_l, wavef_l, pot_r, en_r, wavef_r)
   deallocate (x_data)
   deallocate (h11_data)
   deallocate (h22_data)
   deallocate (h12_data)
   deallocate (el_coupling)

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


end program kappa_coupling

