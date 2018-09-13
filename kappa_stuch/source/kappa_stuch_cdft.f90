program kappa_stuch_cdft

!!
!! calculates the factor kappa in the Stuchebrukhov's expression
!! for the PCET coupling (Eq. 1.7 in JCP 113 (2000) 10438)
!!
!! Adiabatic and diabatic profiles are read from disk (CDFT calculations)
!!

   use timers
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
   real(kind=8), parameter :: pi = 3.1415926535897932384626434D0

   ! b-spline order parameter
   integer, parameter :: korder=3

   integer :: n, i, j, k, im, npoints, ixc, ixc0
   integer :: iargc, narg, filename_offset
   real(8) :: x_left, x_right, dx, xh, mass, fcij, fc00, phase, dpot_left, dpot_right
   real(8) :: reactant_pot_min, product_pot_min, xp, e0, e1, h11, h22, h12, s12, s12x0, s12x, v12x0, v12x
   real(8) :: sh, dh, e_ground, e_excited, eshift
   real(8) :: tstart, tfinish
   real(8), allocatable, dimension(:)   :: x_data, e0_data, e1_data, h11_data, h22_data, h12_data, s12_data
   real(8), allocatable, dimension(:)   :: pot_l, pot_r, x, en_l, en_r, pot_diff
   real(8), allocatable, dimension(:)   :: el_coupling, el_overlap
   real(8), allocatable, dimension(:,:) :: wavef_l, wavef_r, vib_coupling, overlap

   real(8) :: dgamma
   real(8) :: xc, tunn_energy, df, vt, p, tau_p, tau_e, kappa

   character(len=80) :: argument, data_file
   character(len=40) :: particle, n_char
   character(len=5)  :: ovl_symb
   character(len=10) :: vib_symb

   ! b-spline arrays
   real(8), dimension(:), allocatable :: xknot, bscoef

   narg = iargc()
   if (narg.eq.0) then
      write(*,*)
      write(*,*) "=============================================================================="
      write(*,*) " Usage:"
      write(*,*) "   xxx.bin <name_of_the_file_with_PT profiles> <n_grid_points> <particle_kind>"
      write(*,*) "   Columns in file_with_PT_profiles (energies in kcal/mol):"
      write(*,*) "   x(A)   E_adiab(1)  E_adiab(2)   H11   H12   H22"
      write(*,*) "   ..."
      write(*,*) "   Number of grid points should be power of two (32,64,128,etc)"
      write(*,*) "   Particle_kind: 0 - electron, 1 - proton, 2 - deuteron, 3 - triton"
      write(*,*) "=============================================================================="
      stop
   endif

   call getarg(1,argument)
   data_file = trim(argument)
   open(1,file=data_file,status='old')

   !-- scan input file and determine left and right integration limits

   npoints = 0
   x_left = 999.d0
   x_right = -999.d0
   do
      read(1,*,end=11) xp, e0, e1, h11, h12, h22, s12
      npoints = npoints + 1
      if (xp.le.x_left) x_left = xp
      if (xp.ge.x_right) x_right = xp
   enddo
11 continue
   rewind 1

   xh = x_right - x_left
   write(*,'(/1x,"Left  integration limit:            ",f12.6," A")') x_left
   write(*,'( 1x,"Right integration limit:            ",f12.6," A")') x_right
   write(*,'( 1x,"Length of the integration interval: ",f12.6," A")') xh
   write(*,'( 1x,"Number of data points:              ",i12)') npoints

   !-- n - number of grid points (should be power of two)
   call getarg(2,argument)
   read(argument,*) n
   write(*,'(/1x,"Number of grid points for splined potentials: ",i10)') n
   write(n_char,'(i4.4,a2)') n, "pt"

   !-- Mass of the quantum particle (electronic mass units)

   call getarg(3,argument)
   read(argument,*) im

   !write(*,'(/1x,"Quantum particle: ")')
   !write(*,'( 1x,"   (0) electron")')
   !write(*,'( 1x,"   (1) proton")')
   !write(*,'( 1x,"   (2) deuteron")')
   !write(*,'( 1x,"   (3) triton")')

   select case(im)
      case(0)
         mass = 1.d0
         particle = "electron"//"_"//trim(n_char)
      case(1)
         mass = 1836.152701d0
         particle = "proton"//"_"//trim(n_char)
      case(2)
         mass = 3670.482955d0
         particle = "deuteron"//"_"//trim(n_char)
      case(3)
         mass = 5497.9210703d0
         particle = "triton"//"_"//trim(n_char)
      case default
         write(*,'(/1x,"You must specify the particle kind! ")')
         stop "Come again..."
   end select

   write(*,'(/1x,"Mass of the quantum particle: ",f15.6," electronic mass units")') mass

   ! allocate arrays
   allocate (x_data(npoints))
   allocate (e0_data(npoints))
   allocate (e1_data(npoints))
   allocate (h11_data(npoints))
   allocate (h22_data(npoints))
   allocate (h12_data(npoints))
   allocate (s12_data(npoints))

   allocate (x(n))
   allocate (pot_l(n), en_l(n), wavef_l(n,n))
   allocate (pot_r(n), en_r(n), wavef_r(n,n))
   allocate (el_coupling(n))
   allocate (el_overlap(n))

   !-- read data points
   do i=1,npoints
      read(1,*) x_data(i), e0_data(i), e1_data(i), h11_data(i), h12_data(i), h22_data(i), s12_data(i)
   enddo
   close(1)

   ! allocate b-spline arrays
   allocate (bscoef(npoints), xknot(npoints+korder))

   ! generate knot sequence
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

   !-- generate spline coefficients for overlap integral
   bscoef = 0.d0
   call dbsint (npoints, x_data, s12_data, korder, xknot, bscoef)

   !-- evaluate b-spline at the grid points
   el_overlap(1) = s12_data(1)
   do i=2,n-1
      el_overlap(i) = dbsval(x(i),korder,xknot,npoints,bscoef)
   enddo
   el_overlap(n) = s12_data(npoints)


   !-- calculate energy levels and wavefunctions
   !   in the left potential well

   !-- left diabatic potential (reactant)

   tstart = second()
   call schroedinger(n,xh,mass,pot_l,en_l,wavef_l)
   tfinish = second()
   write(*,*)
   write(*,*) ' Left potential: time elapsed --> ',tfinish-tstart,' seconds'
   write(*,*)
   write(*,*) 'First 10  left eigenvalues (kcal/mol):'
   write(*,*) (en_l(k),k=1,10)

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
   write(*,*) 'First 10 right eigenvalues (kcal/mol):'
   write(*,*) (en_r(k),k=1,10)

   !-- adjust phase
   phase = sign(1.d0,wavef_r(n/2,1))
   wavef_r = wavef_r*phase

   !-- shift the right diabatic potential and
   !   its vibrational levels so that
   !   the ground vibrational states in the left
   !   and right potentials are degenerate

   eshift = en_r(1) - en_l(1)
   do i=1,n
      pot_r(i) = pot_r(i) - eshift
      en_r(i)  = en_r(i) - eshift
   enddo

   !-- tunneling energy is the energy of the ground state
   !   in the left well (E in the paper)
   tunn_energy = en_l(1)*kcal2au

   !-- output the wavefunctions

   filename_offset = len(data_file)
   do i=len(data_file),1,-1
      if (data_file(i:i).eq.".") then
         filename_offset = i - 1
         exit
      endif
   enddo

   open(unit=1,file=data_file(1:filename_offset)//"_"//trim(particle)//"_left.dat")
   do i=1,n
      write(1,'(40g20.10)') x(i),pot_l(i),(en_l(k),en_l(k)+wavef_l(i,k)*30.d0,k=1,10)
   enddo
   close(1)

   open(unit=1,file=data_file(1:filename_offset)//"_"//trim(particle)//"_right.dat")
   do i=1,n
      write(1,'(40g20.10)') x(i), pot_r(i), (en_r(k),en_r(k) + wavef_r(i,k)*30.d0, k=1,10)
   enddo
   close(1)


   !-- find the crossing point where pot_l(xc) = pot_r(xc) (unshifted)
   allocate (pot_diff(n))
   pot_diff = pot_l - (pot_r + eshift)
   do i=2,n
      if (pot_diff(i)*pot_diff(i-1).le.0.d0) then
         ixc0 = i
         exit
      endif
   enddo
   deallocate(pot_diff)


   !-- find the crossing point where pot_l(xc) = pot_r(xc)
   allocate (pot_diff(n))
   pot_diff = pot_l - pot_r
   do i=2,n
      if (pot_diff(i)*pot_diff(i-1).le.0.d0) then
         ixc = i
         exit
      endif
   enddo
   deallocate(pot_diff)

   !-- absolute value of the electronic overlap at the crossing points
   s12x0 = 0.5d0*(abs(el_overlap(ixc0)) + abs(el_overlap(ixc0-1)))
   s12x  = 0.5d0*(abs(el_overlap(ixc))  + abs(el_overlap(ixc-1)))

   !-- electronic coupling at the crossing points
   v12x0 = 0.5d0*(el_coupling(ixc0) + el_coupling(ixc0-1))*kcal2au
   v12x  = 0.5d0*(el_coupling(ixc)  + el_coupling(ixc-1))*kcal2au

   write(*,*)
   write(*,*) "---------------Unshifted diabatic potentials------------------------------"
   write(*,*) "Crossing point number: ", ixc
   write(*,*) "Crossing point coordinate value: ", x_left + (ixc-1)*dx
   write(*,*) "Electronic overlap  at the crossing point: ", s12x
   write(*,*) "Electronic coupling at the crossing point (eV): ", v12x*au2ev
   write(*,*) "                                       (cm^-1): ", v12x*au2ev*ev2cm
   write(*,*) "                                    (kcal/mol): ", v12x*au2kcal

   write(*,*)
   write(*,*) "---------------Shifted diabatic potentials--------------------------------"
   write(*,*) "Crossing point number: ", ixc0
   write(*,*) "Crossing point coordinate value: ", x_left + (ixc0-1)*dx
   write(*,*) "Electronic overlap  at the crossing point: ", s12x0
   write(*,*) "Electronic coupling at the crossing point (eV): ", v12x0*au2ev
   write(*,*) "                                       (cm^-1): ", v12x0*au2ev*ev2cm
   write(*,*) "                                    (kcal/mol): ", v12x0*au2kcal
   write(*,*) " Tunneling energy (kcal/mol):                   ", tunn_energy*au2kcal

   open(unit=1,file=data_file(1:filename_offset)//"_"//trim(particle)//"_profiles_vconst.dat")
   do i=1,n
      h11 = pot_l(i)
      h22 = pot_r(i)
      sh = h11 + h22
      dh = h22 - h11
      h12 = v12x*au2kcal
      e_ground  = 0.5d0*sh - 0.5d0*sqrt(dh*dh + 4.d0*h12*h12)
      e_excited = 0.5d0*sh + 0.5d0*sqrt(dh*dh + 4.d0*h12*h12)
      write(1,'(10g20.10)') x(i), h11, h22, h12, e_ground, e_excited, tunn_energy*au2kcal, s12x
   enddo
   close(1)

   open(unit=1,file=data_file(1:filename_offset)//"_"//trim(particle)//"_profiles.dat")
   do i=1,n
      h11 = pot_l(i)
      h22 = pot_r(i)
      sh = h11 + h22
      dh = h22 - h11
      h12 = el_coupling(i)
      s12 = el_overlap(i)
      e_ground  = 0.5d0*sh - 0.5d0*sqrt(dh*dh + 4.d0*h12*h12)
      e_excited = 0.5d0*sh + 0.5d0*sqrt(dh*dh + 4.d0*h12*h12)
      write(1,'(10g20.10)') x(i), h11, h22, h12, e_ground, e_excited, tunn_energy*au2kcal, s12
   enddo
   close(1)

   open(unit=1,file=data_file(1:filename_offset)//"_"//trim(particle)//"_profiles_vconst_unshifted.dat")
   do i=1,n
      h11 = pot_l(i)
      h22 = pot_r(i) + eshift
      sh = h11 + h22
      dh = h22 - h11
      h12 = v12x0*au2kcal
      s12 = s12x0
      e_ground  = 0.5d0*sh - 0.5d0*sqrt(dh*dh + 4.d0*h12*h12)
      e_excited = 0.5d0*sh + 0.5d0*sqrt(dh*dh + 4.d0*h12*h12)
      write(1,'(10g20.10)') x(i), h11, h22, h12, e_ground, e_excited, en_l(1), en_r(1) + eshift, s12
   enddo
   close(1)

   open(unit=1,file=data_file(1:filename_offset)//"_"//trim(particle)//"_profiles_unshifted.dat")
   do i=1,n
      h11 = pot_l(i)
      h22 = pot_r(i) + eshift
      sh = h11 + h22
      dh = h22 - h11
      h12 = el_coupling(i)
      s12 = el_overlap(i)
      e_ground  = 0.5d0*sh - 0.5d0*sqrt(dh*dh + 4.d0*h12*h12)
      e_excited = 0.5d0*sh + 0.5d0*sqrt(dh*dh + 4.d0*h12*h12)
      write(1,'(10g20.10)') x(i), h11, h22, h12, e_ground, e_excited, en_l(1), en_r(1) + eshift, s12
   enddo
   close(1)

   !== calculate the kappa from Stuchebrukhov's expression


   !-- calculate the difference of the slopes of the
   !   diabatic profiles at the crossing point (pot_l = pot_r)
   !   (DeltaF in the paper)

   dpot_left  = (pot_l(ixc) - pot_l(ixc-1))/dx
   dpot_right = (pot_r(ixc) - pot_r(ixc-1))/dx

   df = abs(dpot_left - dpot_right)    ! in kcal/mol/A
   df = df*kcal2au/a2bohr              ! in a.u./Bohr

   !-- calculate tunneling velocity (Eq. 1.9)
   vt = sqrt(2.d0*(pot_l(ixc)*kcal2au - tunn_energy)/mass)

   !-- adiabaticity parameter p (Eq. 1.8)
   p = v12x*v12x/df/vt

   !-- tunneling times for electron and proton (Eqs. 4.12-4.13)
   tau_p = v12x/df/vt
   tau_e = 1.d0/v12x

   write(*,*)
   write(*,*) " Proton   tunneling time (fs): ", tau_p*au2fs
   write(*,*) " Electron tunneling time (fs): ", tau_e*au2fs
   write(*,*) " Adiabaticity parameter: ", p

   open(unit=1,file=data_file(1:filename_offset)//"_"//trim(particle)//'_results.dat')

   write(1,*) "======================================================================"
   write(1,*) " UNSHIFTED (ORIGINAL) DIABATIC POTENTIALS"
   write(1,*) " Crossing point number: ", ixc0
   write(1,*) " Crossing point coordinate value: ", x_left + (ixc0-1)*dx
   write(1,*) " Electronic overlap  at the crossing point     : ", s12x0
   write(1,*) " Electronic coupling at the crossing point (eV): ", v12x0*au2ev
   write(1,*) "                                        (cm^-1): ", v12x0*au2ev*ev2cm
   write(1,*) "                                     (kcal/mol): ", v12x0*au2kcal
   write(1,*) "======================================================================"
   write(1,*) " SHIFTED DIABATIC POTENTIALS"
   write(1,*) " Crossing point number: ", ixc
   write(1,*) " Crossing point coordinate value: ", x_left + (ixc-1)*dx
   write(1,*) " Electronic overlap  at the crossing point     : ", s12x
   write(1,*) " Electronic coupling at the crossing point (eV): ", v12x*au2ev
   write(1,*) "                                        (cm^-1): ", v12x*au2ev*ev2cm
   write(1,*) "                                     (kcal/mol): ", v12x*au2kcal
   write(1,*) "----------------------------------------------------------------------"
   write(1,*) " Tunneling energy (kcal/mol):                    ", tunn_energy*au2kcal
   write(1,*) "======================================================================"
   write(1,*) " QUASICLASSICAL EXPRESSION PARAMETERS"
   write(1,*) " Proton   tunneling time (fs): ", tau_p*au2fs
   write(1,*) " Electron tunneling time (fs): ", tau_e*au2fs
   write(1,*) " Adiabaticity parameter: ", p

   !-- kappa (Eq. 1.7)

   kappa = sqrt(2.d0*pi*p)*exp(p*log(p)-p)/dgamma(p+1.d0)
   write(*,*)
   write(*,*) " kappa: ", kappa
   write(1,*) " kappa: ", kappa
   write(1,*) "====================================================================="


   !-- calculate Franck-Condon factors (overlaps) for first ten states

   write(*,*)
   write(*,*) "====================================================================="
   write(*,*) "FRANCK-CONDON OVERLAPS"
   write(*,*) "interval: ",x_left," - ",x_right
   write(*,*) "number of grid points: ",n
   write(*,*) "---------------------------------------------------------------------"

   write(1,*) "======================================================================"
   write(1,*) " FRANCK-CONDON OVERLAPS"
   write(1,*) " interval: ", x_left, " - ", x_right
   write(1,*) " number of grid points: ", n
   write(1,*) "----------------------------------------------------------------------"

   open(unit=2,file=data_file(1:filename_offset)//"_"//trim(particle)//"_overlaps.dat")
   do i=1,10
      do j=1,10
         fcij = 0.d0
         do k=1,n
            fcij = fcij + wavef_l(k,i)*wavef_r(k,j)
         enddo
         if (i.eq.1.and.j.eq.1) fc00 = fcij
         write(ovl_symb,'("<",i1,"|",i1,">")') i-1,j-1
         write(*,*) ovl_symb," = ", fcij
         write(1,*) ovl_symb," = ", fcij
         write(2,*) i-1, j-1, fcij
      enddo
      write(2,*)
   enddo
   close(2)
   write(*,*) "---------------------------------------------------------------------"
   write(1,*) "---------------------------------------------------------------------"
   write(*,*) " Nonadiabatic coupling V_el*<0|0> = ", v12x*fc00*au2ev*ev2cm," cm^-1"
   write(1,*) " Nonadiabatic coupling V_el*<0|0> = ", v12x*fc00*au2ev*ev2cm," cm^-1"

   !-- calculate vibronic couplings between locally diabatic vibronic states
   !   for first five states

   write(*,*)
   write(*,*) "====================================================================="
   write(*,*) "VIBRONIC COUPLINGS"
   write(*,*) "interval: ", x_left, " - ", x_right
   write(*,*) "number of grid points: ", n
   write(*,*) "---------------------------------------------------------------------"

   write(1,*) "#====================================================================="
   write(1,*) "FRANCK-CONDON OVERLAPS"
   write(1,*) "interval: ",x_left," - ",x_right
   write(1,*) "number of grid points: ",n
   write(1,*) "#---------------------------------------------------------------------"

   open(unit=2,file=data_file(1:filename_offset)//"_"//trim(particle)//"_vibcouplings.dat")
   do i=1,10
      do j=1,10
         fcij = 0.d0
         do k=1,n
            fcij = fcij + wavef_l(k,i)*el_coupling(k)*wavef_r(k,j)
         enddo
         if (i.eq.1.and.j.eq.1) fc00 = fcij
         write(vib_symb,'("<",i1,"|V_el|",i1,">")') i-1,j-1
         write(*,*) vib_symb," = ",fcij
         write(1,*) vib_symb," = ",fcij
         write(2,*) i-1, j-1, fcij
      enddo
      write(2,*)
   enddo
   close(2)
   write(*,*) "---------------------------------------------------------------------"
   write(1,*) "---------------------------------------------------------------------"
   write(*,*) " Nonadiabatic coupling <0|V_el|0> = ", v12x*fc00*au2ev*ev2cm," cm^-1"
   write(1,*) " Nonadiabatic coupling <0|V_el|0> = ", v12x*fc00*au2ev*ev2cm," cm^-1"
   write(*,*) "====================================================================="
   write(1,*) "====================================================================="
   close(1)

   ! deallocate arrays
   deallocate (x, pot_l, en_l, wavef_l, pot_r, en_r, wavef_r)
   deallocate (x_data)
   deallocate (e0_data)
   deallocate (e1_data)
   deallocate (h11_data)
   deallocate (h22_data)
   deallocate (h12_data)
   deallocate (s12_data)
   deallocate (el_coupling)
   deallocate (el_overlap)

end program kappa_stuch_cdft

