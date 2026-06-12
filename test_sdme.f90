!==============================================================
!     test_sdme.f90: Test program for SDME angular distribution
!     Compile: gfortran -O2 -ffree-line-length-none -std=legacy \
!              -o test_sdme test_sdme.f90 diffrad_vpk.f90
!==============================================================
      program test_sdme
      implicit none
      real*4  urand
      real*8  sdme_W, sdme_Wmax, calc_epsilon
      external sdme_W, sdme_Wmax, calc_epsilon
      real*8  r(15)
      integer iy, nevt, nbin
      parameter(nevt=500000, nbin=50)

      real*8  hist_ct(nbin), hist_pd(nbin), hist_PL(nbin)
      real*8  costh, phid, aPhi, eps, wmax, wval
      real*8  aPI, avgbin, binval
      real*8  ebeam_t, Q2_t, xB_t, amp_t, anu_t, y_t, eps_t, wint
      integer ievt, i, ibin, nstars
      character(len=60) :: stars

      aPI = acos(-1d0)
      iy  = 12345
      stars = '************************************************************'

!     ──────────────────────────────────────────────────────
!     TEST 1: Isotropic (all SDMEs = 0)
!     ──────────────────────────────────────────────────────
      write(*,*) '========================================='
      write(*,*) ' TEST 1: Isotropic (all SDMEs = 0)'
      write(*,*) '========================================='
      call init_sdme(r, 0)
      eps = 0.8d0

      wmax = sdme_Wmax(eps, r)
      write(*,'(A,F12.6)') '  Wmax = ', wmax

      write(*,'(A,F12.6)') '  W(0, 0, 0)     = ', &
        sdme_W(0d0, 0d0, 0d0, eps, r)
      write(*,'(A,F12.6)') '  W(1, 0, 0)     = ', &
        sdme_W(1d0, 0d0, 0d0, eps, r)
      write(*,'(A,F12.6)') '  W(0, pi, pi)   = ', &
        sdme_W(0d0, aPI, aPI, eps, r)
      write(*,'(A,F12.6)') '  Expected 3/(16*pi^2) = ', &
        3d0/(16d0*aPI**2)
      write(*,*)

!     ──────────────────────────────────────────────────────
!     TEST 2: SCHC + NPE (photoproduction-like)
!     ──────────────────────────────────────────────────────
      write(*,*) '========================================='
      write(*,*) ' TEST 2: SCHC+NPE (photoproduction)'
      write(*,*) '  r04_00=0, r1_11=0.5, r1_1-1=0.5'
      write(*,*) '  Im r2_1-1=-0.5, eps=0.8'
      write(*,*) '========================================='
      call init_sdme(r, 1)
      eps = 0.8d0

      wmax = sdme_Wmax(eps, r)
      write(*,'(A,F12.6)') '  Wmax = ', wmax

      do i=1,nbin
        hist_ct(i) = 0d0
        hist_pd(i) = 0d0
        hist_PL(i) = 0d0
      enddo

      do ievt = 1, nevt
        call sample_sdme_angles(costh, phid, aPhi, eps, r, wmax, iy)

        ibin = min(nbin, max(1, int((costh+1d0)/2d0*nbin) + 1))
        hist_ct(ibin) = hist_ct(ibin) + 1d0

        ibin = min(nbin, max(1, int(phid/(2d0*aPI)*nbin) + 1))
        hist_pd(ibin) = hist_pd(ibin) + 1d0

        ibin = min(nbin, max(1, int(aPhi/(2d0*aPI)*nbin) + 1))
        hist_PL(ibin) = hist_PL(ibin) + 1d0
      enddo

      avgbin = dble(nevt)/dble(nbin)

      write(*,*)
      write(*,*) '  cos(Theta) distribution (~sin^2 shape):'
      write(*,*) '  ─────────────────────────────────────'
      do i = 1, nbin
        binval = -1d0 + (dble(i)-0.5d0)*2d0/dble(nbin)
        nstars = nint(hist_ct(i)/avgbin * 25d0)
        nstars = min(60, max(0, nstars))
        write(*,'(F7.3,A,I6,1X,A)') binval, ' |', &
          nint(hist_ct(i)), stars(1:nstars)
      enddo

      write(*,*)
      write(*,*) '  phi distribution (decay plane angle):'
      write(*,*) '  ─────────────────────────────────────'
      do i = 1, nbin
        binval = (dble(i)-0.5d0)*360d0/dble(nbin)
        nstars = nint(hist_pd(i)/avgbin * 25d0)
        nstars = min(60, max(0, nstars))
        write(*,'(F7.1,A,I6,1X,A)') binval, ' |', &
          nint(hist_pd(i)), stars(1:nstars)
      enddo

      write(*,*)
      write(*,*) '  Phi distribution (production vs lepton):'
      write(*,*) '  ─────────────────────────────────────'
      do i = 1, nbin
        binval = (dble(i)-0.5d0)*360d0/dble(nbin)
        nstars = nint(hist_PL(i)/avgbin * 25d0)
        nstars = min(60, max(0, nstars))
        write(*,'(F7.1,A,I6,1X,A)') binval, ' |', &
          nint(hist_PL(i)), stars(1:nstars)
      enddo

!     ──────────────────────────────────────────────────────
!     TEST 3: Large r04_00 = 0.5 (longitudinal component)
!     ──────────────────────────────────────────────────────
      write(*,*)
      write(*,*) '========================================='
      write(*,*) ' TEST 3: r04_00 = 0.5 (longitudinal mix)'
      write(*,*) '  All other SDMEs = 0, eps = 0.8'
      write(*,*) '========================================='
      call init_sdme(r, 0)
      r(1) = 0.5d0
      eps = 0.8d0

      wmax = sdme_Wmax(eps, r)
      write(*,'(A,F12.6)') '  Wmax = ', wmax

      do i=1,nbin
        hist_ct(i) = 0d0
      enddo

      do ievt = 1, nevt
        call sample_sdme_angles(costh, phid, aPhi, eps, r, wmax, iy)
        ibin = min(nbin, max(1, int((costh+1d0)/2d0*nbin) + 1))
        hist_ct(ibin) = hist_ct(ibin) + 1d0
      enddo

      write(*,*)
      write(*,*) '  cos(Theta) (sin^2+cos^2 mixture):'
      write(*,*) '  ─────────────────────────────────────'
      do i = 1, nbin
        binval = -1d0 + (dble(i)-0.5d0)*2d0/dble(nbin)
        nstars = nint(hist_ct(i)/avgbin * 25d0)
        nstars = min(60, max(0, nstars))
        write(*,'(F7.3,A,I6,1X,A)') binval, ' |', &
          nint(hist_ct(i)), stars(1:nstars)
      enddo

!     ──────────────────────────────────────────────────────
!     TEST 4: epsilon calculation
!     ──────────────────────────────────────────────────────
      write(*,*)
      write(*,*) '========================================='
      write(*,*) ' TEST 4: epsilon calculation'
      write(*,*) '========================================='
      amp_t = 0.93827d0
      ebeam_t = 10.6d0

      Q2_t = 2.0d0; xB_t = 0.3d0
      anu_t = Q2_t/(2d0*amp_t*xB_t)
      y_t = anu_t/ebeam_t
      eps_t = calc_epsilon(ebeam_t, y_t, Q2_t)
      write(*,'(A,F6.1,A,F5.2,A,F6.3,A,F7.4)') &
        '  E=',ebeam_t,' Q2=',Q2_t,' y=',y_t,' eps=',eps_t

      Q2_t = 1.0d0; xB_t = 0.15d0
      anu_t = Q2_t/(2d0*amp_t*xB_t)
      y_t = anu_t/ebeam_t
      eps_t = calc_epsilon(ebeam_t, y_t, Q2_t)
      write(*,'(A,F6.1,A,F5.2,A,F6.3,A,F7.4)') &
        '  E=',ebeam_t,' Q2=',Q2_t,' y=',y_t,' eps=',eps_t

      Q2_t = 5.0d0; xB_t = 0.4d0
      anu_t = Q2_t/(2d0*amp_t*xB_t)
      y_t = anu_t/ebeam_t
      eps_t = calc_epsilon(ebeam_t, y_t, Q2_t)
      write(*,'(A,F6.1,A,F5.2,A,F6.3,A,F7.4)') &
        '  E=',ebeam_t,' Q2=',Q2_t,' y=',y_t,' eps=',eps_t

      Q2_t = 0.001d0; xB_t = 0.001d0
      anu_t = Q2_t/(2d0*amp_t*xB_t)
      y_t = anu_t/ebeam_t
      eps_t = calc_epsilon(ebeam_t, y_t, Q2_t)
      write(*,'(A,F6.1,A,F7.4,A,F7.4,A,F7.4)') &
        '  E=',ebeam_t,' Q2=',Q2_t,' y=',y_t,' eps=',eps_t

!     ──────────────────────────────────────────────────────
!     TEST 5: Normalization check
!     ──────────────────────────────────────────────────────
      write(*,*)
      write(*,*) '========================================='
      write(*,*) ' TEST 5: Normalization (integral of W)'
      write(*,*) '========================================='

      call init_sdme(r, 0)
      eps = 0.8d0
      call integrate_W(r, eps, wint)
      write(*,'(A,F10.6,A)') '  Isotropic:  ', wint, '  (expect 1.0)'

      call init_sdme(r, 1)
      call integrate_W(r, eps, wint)
      write(*,'(A,F10.6,A)') '  SCHC+NPE:   ', wint, '  (expect 1.0)'

      call init_sdme(r, 0)
      r(1) = 0.5d0
      call integrate_W(r, eps, wint)
      write(*,'(A,F10.6,A)') '  r04_00=0.5: ', wint, '  (expect 1.0)'

      write(*,*)
      write(*,*) 'All tests complete.'

      stop
      end


!==============================================================
      subroutine integrate_W(r, eps, wint)
      implicit none
      real*8 r(15), eps, wint
      real*8 sdme_W
      external sdme_W
      real*8 aPI, dcth, dpd, dpl
      real*8 act, apd, apl
      integer nc, nd, nl, ic, jd, kl

      aPI = acos(-1d0)
      nc = 100; nd = 100; nl = 100

      wint = 0d0
      dcth = 2d0/dble(nc)
      dpd  = 2d0*aPI/dble(nd)
      dpl  = 2d0*aPI/dble(nl)

      do ic = 1, nc
        act = -1d0 + (dble(ic)-0.5d0)*dcth
        do jd = 1, nd
          apd = (dble(jd)-0.5d0)*dpd
          do kl = 1, nl
            apl = (dble(kl)-0.5d0)*dpl
            wint = wint + sdme_W(act, apd, apl, eps, r) &
                        * dcth * dpd * dpl
          enddo
        enddo
      enddo

      return
      end
