!==============================================================
!     test_sdme_plots.f90: Generate data files for SDME plots
!     Compile: gfortran -O2 -ffree-line-length-none -std=legacy \
!              -o test_sdme_plots test_sdme_plots.f90 diffrad_vpk.f90
!     (link against subroutines only, not main program)
!==============================================================
      program test_sdme_plots
      implicit none
      real*4  urand
      real*8  sdme_W, sdme_Wmax, calc_epsilon
      external sdme_W, sdme_Wmax, calc_epsilon
      real*8  r(15)
      integer iy, nevt, nbin
      parameter(nevt=2000000, nbin=80)

      real*8  h1(nbin), h2(nbin), h3(nbin)
      real*8  costh, phid, aPhi, eps, wmax
      real*8  aPI, binval, dcth, dpd, dpl
      integer ievt, i, ibin
      real*8  wint

      aPI = acos(-1d0)
      iy  = 54321

!     ══════════════════════════════════════════════════════
!     DATASET 1: SCHC+NPE (photoproduction-like), eps=0.8
!     ══════════════════════════════════════════════════════
      write(*,*) 'Generating SCHC+NPE dataset...'
      call init_sdme(r, 1)
      eps = 0.8d0
      wmax = sdme_Wmax(eps, r)

      do i=1,nbin
        h1(i)=0d0; h2(i)=0d0; h3(i)=0d0
      enddo

      do ievt = 1, nevt
        call sample_sdme_angles(costh, phid, aPhi, eps, r, wmax, iy)
        ibin = min(nbin, max(1, int((costh+1d0)/2d0*nbin) + 1))
        h1(ibin) = h1(ibin) + 1d0
        ibin = min(nbin, max(1, int(phid/(2d0*aPI)*nbin) + 1))
        h2(ibin) = h2(ibin) + 1d0
        ibin = min(nbin, max(1, int(aPhi/(2d0*aPI)*nbin) + 1))
        h3(ibin) = h3(ibin) + 1d0
      enddo

      open(20, file='sdme_schc_costh.dat')
      open(21, file='sdme_schc_phi.dat')
      open(22, file='sdme_schc_Phi.dat')
      do i = 1, nbin
        write(20,'(F10.4,1X,F12.1)') &
          -1d0+(dble(i)-0.5d0)*2d0/dble(nbin), h1(i)
        write(21,'(F10.4,1X,F12.1)') &
          (dble(i)-0.5d0)*360d0/dble(nbin), h2(i)
        write(22,'(F10.4,1X,F12.1)') &
          (dble(i)-0.5d0)*360d0/dble(nbin), h3(i)
      enddo
      close(20); close(21); close(22)

!     ══════════════════════════════════════════════════════
!     DATASET 2: r04_00 = 0.5 (longitudinal mix), eps=0.8
!     ══════════════════════════════════════════════════════
      write(*,*) 'Generating r04_00=0.5 dataset...'
      call init_sdme(r, 0)
      r(1) = 0.5d0
      eps = 0.8d0
      wmax = sdme_Wmax(eps, r)

      do i=1,nbin
        h1(i)=0d0; h2(i)=0d0; h3(i)=0d0
      enddo

      do ievt = 1, nevt
        call sample_sdme_angles(costh, phid, aPhi, eps, r, wmax, iy)
        ibin = min(nbin, max(1, int((costh+1d0)/2d0*nbin) + 1))
        h1(ibin) = h1(ibin) + 1d0
        ibin = min(nbin, max(1, int(phid/(2d0*aPI)*nbin) + 1))
        h2(ibin) = h2(ibin) + 1d0
        ibin = min(nbin, max(1, int(aPhi/(2d0*aPI)*nbin) + 1))
        h3(ibin) = h3(ibin) + 1d0
      enddo

      open(20, file='sdme_r04_costh.dat')
      open(21, file='sdme_r04_phi.dat')
      open(22, file='sdme_r04_Phi.dat')
      do i = 1, nbin
        write(20,'(F10.4,1X,F12.1)') &
          -1d0+(dble(i)-0.5d0)*2d0/dble(nbin), h1(i)
        write(21,'(F10.4,1X,F12.1)') &
          (dble(i)-0.5d0)*360d0/dble(nbin), h2(i)
        write(22,'(F10.4,1X,F12.1)') &
          (dble(i)-0.5d0)*360d0/dble(nbin), h3(i)
      enddo
      close(20); close(21); close(22)

!     ══════════════════════════════════════════════════════
!     DATASET 3: Isotropic (all SDMEs = 0), eps=0.8
!     ══════════════════════════════════════════════════════
      write(*,*) 'Generating isotropic dataset...'
      call init_sdme(r, 0)
      eps = 0.8d0
      wmax = sdme_Wmax(eps, r)

      do i=1,nbin
        h1(i)=0d0; h2(i)=0d0; h3(i)=0d0
      enddo

      do ievt = 1, nevt
        call sample_sdme_angles(costh, phid, aPhi, eps, r, wmax, iy)
        ibin = min(nbin, max(1, int((costh+1d0)/2d0*nbin) + 1))
        h1(ibin) = h1(ibin) + 1d0
        ibin = min(nbin, max(1, int(phid/(2d0*aPI)*nbin) + 1))
        h2(ibin) = h2(ibin) + 1d0
        ibin = min(nbin, max(1, int(aPhi/(2d0*aPI)*nbin) + 1))
        h3(ibin) = h3(ibin) + 1d0
      enddo

      open(20, file='sdme_iso_costh.dat')
      open(21, file='sdme_iso_phi.dat')
      open(22, file='sdme_iso_Phi.dat')
      do i = 1, nbin
        write(20,'(F10.4,1X,F12.1)') &
          -1d0+(dble(i)-0.5d0)*2d0/dble(nbin), h1(i)
        write(21,'(F10.4,1X,F12.1)') &
          (dble(i)-0.5d0)*360d0/dble(nbin), h2(i)
        write(22,'(F10.4,1X,F12.1)') &
          (dble(i)-0.5d0)*360d0/dble(nbin), h3(i)
      enddo
      close(20); close(21); close(22)

!     ══════════════════════════════════════════════════════
!     DATASET 4: r04_00=0.3, Re r04_10=-0.05, r04_1-1=-0.02
!                r1_11=0.4, r1_00=0.05, r5_00=0.15
!                "HERMES-like" realistic electroproduction
!     ══════════════════════════════════════════════════════
      write(*,*) 'Generating HERMES-like dataset...'
      call init_sdme(r, 0)
      r(1)  = 0.3d0      ! r04_00
      r(2)  = -0.05d0     ! Re r04_10
      r(3)  = -0.02d0     ! r04_1-1
      r(4)  = 0.4d0       ! r1_11
      r(5)  = 0.05d0      ! r1_00
      r(6)  = -0.03d0     ! Re r1_10
      r(7)  = 0.3d0       ! r1_1-1
      r(8)  = 0.01d0      ! Im r2_10
      r(9)  = -0.3d0      ! Im r2_1-1
      r(10) = 0.1d0       ! r5_11
      r(11) = 0.15d0      ! r5_00
      r(12) = -0.04d0     ! Re r5_10
      r(13) = 0.05d0      ! r5_1-1
      r(14) = 0.02d0      ! Im r6_10
      r(15) = -0.05d0     ! Im r6_1-1
      eps = 0.85d0
      wmax = sdme_Wmax(eps, r)

      do i=1,nbin
        h1(i)=0d0; h2(i)=0d0; h3(i)=0d0
      enddo

      do ievt = 1, nevt
        call sample_sdme_angles(costh, phid, aPhi, eps, r, wmax, iy)
        ibin = min(nbin, max(1, int((costh+1d0)/2d0*nbin) + 1))
        h1(ibin) = h1(ibin) + 1d0
        ibin = min(nbin, max(1, int(phid/(2d0*aPI)*nbin) + 1))
        h2(ibin) = h2(ibin) + 1d0
        ibin = min(nbin, max(1, int(aPhi/(2d0*aPI)*nbin) + 1))
        h3(ibin) = h3(ibin) + 1d0
      enddo

      open(20, file='sdme_hermes_costh.dat')
      open(21, file='sdme_hermes_phi.dat')
      open(22, file='sdme_hermes_Phi.dat')
      do i = 1, nbin
        write(20,'(F10.4,1X,F12.1)') &
          -1d0+(dble(i)-0.5d0)*2d0/dble(nbin), h1(i)
        write(21,'(F10.4,1X,F12.1)') &
          (dble(i)-0.5d0)*360d0/dble(nbin), h2(i)
        write(22,'(F10.4,1X,F12.1)') &
          (dble(i)-0.5d0)*360d0/dble(nbin), h3(i)
      enddo
      close(20); close(21); close(22)

      write(*,*) 'All data files written.'
      stop
      end
