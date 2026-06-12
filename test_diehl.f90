!==============================================================
!     test_diehl.f90: Test Diehl W_UU, translation, MC sampling
!     Compile:
!       gfortran -O2 -ffree-line-length-none -std=legacy \
!         -o test_diehl test_diehl.f90 sdme_diehl.f90 \
!         <diffrad_subs.f90>
!     (needs sdme_W, sdme_Wmax, init_sdme, sample_sdme_angles,
!      calc_epsilon, URAND from diffrad_vpk.f90)
!==============================================================
      program test_diehl
      implicit none
      real*4  urand
      real*8  sdme_W, sdme_Wmax, calc_epsilon
      real*8  diehl_W, diehl_Wmax
      external sdme_W, sdme_Wmax, calc_epsilon
      external diehl_W, diehl_Wmax

      real*8 r(15), d(15), r2(15)
      real*8 eps, aPI, wval_sw, wval_di
      real*8 costh, phid, aPhi
      real*8 maxdiff, sumdiff, avgdiff, reldiff, maxrel
      integer i, j, ipt, npt
      integer iy, nevt, nbin
      parameter(nevt=2000000, nbin=80)
      real*8 h_sw(nbin), h_di(nbin)
      real*8 wmax_sw, wmax_di, binval, avgbin
      character(len=60) :: stars
      real*8 wint_di

      aPI = acos(-1d0)
      iy = 77777
      stars = '************************************************************'

!     ══════════════════════════════════════════════════════
!     TEST 1: Round-trip conversion SW -> Diehl -> SW
!     ══════════════════════════════════════════════════════
      write(*,*) '=================================================='
      write(*,*) ' TEST 1: Round-trip SW -> Diehl -> SW'
      write(*,*) '=================================================='

!     Set up a nontrivial SW parameter set
      call init_sdme(r, 0)
      r(1)  = 0.30d0    ! r04_00
      r(2)  = -0.05d0   ! Re r04_10
      r(3)  = -0.02d0   ! r04_1-1
      r(4)  = 0.40d0    ! r1_11
      r(5)  = 0.05d0    ! r1_00
      r(6)  = -0.03d0   ! Re r1_10
      r(7)  = 0.30d0    ! r1_1-1
      r(8)  = 0.01d0    ! Im r2_10
      r(9)  = -0.30d0   ! Im r2_1-1
      r(10) = 0.10d0    ! r5_11
      r(11) = 0.15d0    ! r5_00
      r(12) = -0.04d0   ! Re r5_10
      r(13) = 0.05d0    ! r5_1-1
      r(14) = 0.02d0    ! Im r6_10
      r(15) = -0.05d0   ! Im r6_1-1

      write(*,*) '  Original SW r(15):'
      do i = 1, 15
        write(*,'(A,I2,A,F10.6)') '    r(',i,') = ', r(i)
      enddo

      call sw_to_diehl(r, d)
      write(*,*)
      write(*,*) '  Converted to Diehl d(15):'
      do i = 1, 15
        write(*,'(A,I2,A,F10.6)') '    d(',i,') = ', d(i)
      enddo

      call diehl_to_sw(d, r2)
      write(*,*)
      write(*,*) '  Converted back to SW r2(15):'
      maxdiff = 0d0
      do i = 1, 15
        reldiff = abs(r2(i) - r(i))
        if(reldiff .gt. maxdiff) maxdiff = reldiff
        write(*,'(A,I2,A,F10.6,A,F10.6,A,ES10.2)') &
          '    r2(',i,') = ', r2(i), &
          '  orig = ', r(i), '  diff = ', r2(i)-r(i)
      enddo
      write(*,'(A,ES10.2)') '  Max |diff| = ', maxdiff
      if(maxdiff .lt. 1d-12) then
        write(*,*) '  >> PASS: round-trip exact'
      else
        write(*,*) '  >> FAIL: round-trip error!'
      endif

!     ══════════════════════════════════════════════════════
!     TEST 2: Point-by-point comparison W_SW vs W_Diehl
!     ══════════════════════════════════════════════════════
      write(*,*)
      write(*,*) '=================================================='
      write(*,*) ' TEST 2: W_SW vs W_Diehl at random points'
      write(*,*) '  (W_Diehl should = 2*pi * W_SW)'
      write(*,*) '=================================================='

      eps = 0.85d0

      maxdiff = 0d0
      sumdiff = 0d0
      maxrel  = 0d0
      npt = 10000

      do ipt = 1, npt
        costh = 2d0*dble(urand(iy)) - 1d0
        phid  = 2d0*aPI*dble(urand(iy))
        aPhi  = 2d0*aPI*dble(urand(iy))

        wval_sw = sdme_W(costh, phid, aPhi, eps, r)
        wval_di = diehl_W(costh, phid, aPhi, eps, d)

        reldiff = abs(wval_di - 2d0*aPI*wval_sw)
        sumdiff = sumdiff + reldiff
        if(reldiff .gt. maxdiff) maxdiff = reldiff
        if(wval_sw .gt. 1d-15) then
          reldiff = abs(wval_di/(2d0*aPI*wval_sw) - 1d0)
          if(reldiff .gt. maxrel) maxrel = reldiff
        endif
      enddo

      avgdiff = sumdiff/dble(npt)
      write(*,'(A,I6,A)') '  Compared ', npt, ' random points'
      write(*,'(A,ES12.4)') '  Max |W_D - 2pi*W_SW| = ', maxdiff
      write(*,'(A,ES12.4)') '  Avg |W_D - 2pi*W_SW| = ', avgdiff
      write(*,'(A,ES12.4)') '  Max relative diff     = ', maxrel

      write(*,*)
      write(*,*) '  Sample points:'
      write(*,'(A)') '   costh     phi       Phi      ' // &
        ' W_SW          W_Diehl       2pi*W_SW'
      do ipt = 1, 8
        costh = 2d0*dble(urand(iy)) - 1d0
        phid  = 2d0*aPI*dble(urand(iy))
        aPhi  = 2d0*aPI*dble(urand(iy))
        wval_sw = sdme_W(costh, phid, aPhi, eps, r)
        wval_di = diehl_W(costh, phid, aPhi, eps, d)
        write(*,'(3F10.4,3ES14.6)') costh, phid, aPhi, &
          wval_sw, wval_di, 2d0*aPI*wval_sw
      enddo

      if(maxrel .lt. 1d-10) then
        write(*,*) '  >> PASS: W_Diehl = 2pi * W_SW'
      else
        write(*,*) '  >> FAIL: mismatch!'
      endif

!     ══════════════════════════════════════════════════════
!     TEST 3: Normalization of W_Diehl
!     int dPhi/(2pi) int dphi dcosth W_UU = 1
!     ══════════════════════════════════════════════════════
      write(*,*)
      write(*,*) '=================================================='
      write(*,*) ' TEST 3: Normalization of W_Diehl'
      write(*,*) '  int dPhi/(2pi) int dphi dcosth W_UU = 1'
      write(*,*) '=================================================='

!     Isotropic
      call init_sdme(r, 0)
      call sw_to_diehl(r, d)
      call integrate_diehl(d, eps, wint_di)
      write(*,'(A,F10.6,A)') '  Isotropic:  ', wint_di, &
        '  (expect 1.0)'

!     SCHC+NPE
      call init_sdme(r, 1)
      call sw_to_diehl(r, d)
      call integrate_diehl(d, eps, wint_di)
      write(*,'(A,F10.6,A)') '  SCHC+NPE:   ', wint_di, &
        '  (expect 1.0)'

!     Realistic
      r(1) = 0.3d0; r(2) = -0.05d0; r(3) = -0.02d0
      r(4) = 0.4d0; r(5) = 0.05d0; r(6) = -0.03d0
      r(7) = 0.3d0; r(8) = 0.01d0; r(9) = -0.3d0
      r(10) = 0.1d0; r(11) = 0.15d0; r(12) = -0.04d0
      r(13) = 0.05d0; r(14) = 0.02d0; r(15) = -0.05d0
      call sw_to_diehl(r, d)
      call integrate_diehl(d, eps, wint_di)
      write(*,'(A,F10.6,A)') '  Realistic:  ', wint_di, &
        '  (expect 1.0)'

!     ══════════════════════════════════════════════════════
!     TEST 4: MC sampling comparison SW vs Diehl
!     ══════════════════════════════════════════════════════
      write(*,*)
      write(*,*) '=================================================='
      write(*,*) ' TEST 4: MC sampling SW vs Diehl (costh hists)'
      write(*,*) '  r04_00=0.3, eps=0.85, 2M events each'
      write(*,*) '=================================================='

      eps = 0.85d0
      wmax_sw = sdme_Wmax(eps, r)
      wmax_di = diehl_Wmax(eps, d)

      do i = 1, nbin
        h_sw(i) = 0d0
        h_di(i) = 0d0
      enddo

      iy = 11111
      do i = 1, nevt
        call sample_sdme_angles(costh, phid, aPhi, &
                                eps, r, wmax_sw, iy)
        j = min(nbin, max(1, int((costh+1d0)/2d0*nbin) + 1))
        h_sw(j) = h_sw(j) + 1d0
      enddo

      iy = 22222
      do i = 1, nevt
        call sample_diehl_angles(costh, phid, aPhi, &
                                 eps, d, wmax_di, iy)
        j = min(nbin, max(1, int((costh+1d0)/2d0*nbin) + 1))
        h_di(j) = h_di(j) + 1d0
      enddo

      avgbin = dble(nevt)/dble(nbin)
      write(*,*)
      write(*,*) '  cos(th)   SW_count  Diehl_ct  ' // &
        'SW_hist                      Diehl_hist'
      do i = 1, nbin
        binval = -1d0 + (dble(i)-0.5d0)*2d0/dble(nbin)
        write(*,'(F8.3,2I9,2X,A,2X,A)') binval, &
          nint(h_sw(i)), nint(h_di(i)), &
          stars(1:min(60,max(0,nint(h_sw(i)/avgbin*25d0)))), &
          stars(1:min(60,max(0,nint(h_di(i)/avgbin*25d0))))
      enddo

!     ══════════════════════════════════════════════════════
!     TEST 5: Special cases for translation formulas
!     ══════════════════════════════════════════════════════
      write(*,*)
      write(*,*) '=================================================='
      write(*,*) ' TEST 5: Translation for special SDME sets'
      write(*,*) '=================================================='

!     Case A: SCHC+NPE
      write(*,*) '  --- SCHC+NPE ---'
      call init_sdme(r, 1)
      call sw_to_diehl(r, d)
      write(*,*) '  SW r:   r04_00=0, r1_11=0.5, r1_1-1=0.5,' // &
        ' Im r2_1-1=-0.5'
      write(*,*) '  Diehl d:'
      do i = 1, 15
        if(abs(d(i)) .gt. 1d-15) then
          write(*,'(A,I2,A,F10.6)') '    d(',i,') = ', d(i)
        endif
      enddo

!     Case B: Pure longitudinal r04_00=1
      write(*,*) '  --- Pure longitudinal r04_00=1 ---'
      call init_sdme(r, 0)
      r(1) = 1d0
      call sw_to_diehl(r, d)
      write(*,*) '  Diehl d:'
      do i = 1, 15
        if(abs(d(i)) .gt. 1d-15) then
          write(*,'(A,I2,A,F10.6)') '    d(',i,') = ', d(i)
        endif
      enddo
      write(*,*) '  (d(1)=1 means W^LL=1, sin^2 coeff = 0)'

!     Generate data files for plotting
      write(*,*)
      write(*,*) '=================================================='
      write(*,*) ' Generating data files for plots...'
      write(*,*) '=================================================='

!     Realistic scenario: output both SW and Diehl histograms
      call init_sdme(r, 0)
      r(1)=0.3d0; r(2)=-0.05d0; r(3)=-0.02d0
      r(4)=0.4d0; r(5)=0.05d0; r(6)=-0.03d0
      r(7)=0.3d0; r(8)=0.01d0; r(9)=-0.3d0
      r(10)=0.1d0; r(11)=0.15d0; r(12)=-0.04d0
      r(13)=0.05d0; r(14)=0.02d0; r(15)=-0.05d0
      call sw_to_diehl(r, d)
      eps = 0.85d0

      wmax_sw = sdme_Wmax(eps, r)
      wmax_di = diehl_Wmax(eps, d)

!     3 angle histograms, both formalisms
      call gen_comparison(r, d, eps, wmax_sw, wmax_di, &
                          nevt, nbin)

      write(*,*) 'All tests complete.'
      stop
      end


!==============================================================
!     integrate_diehl: numerical integration
!     int dPhi/(2pi) int dphi dcosth W_UU
!==============================================================
      subroutine integrate_diehl(d, eps, wint)
      implicit none
      real*8 d(15), eps, wint
      real*8 diehl_W
      external diehl_W
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
!           Diehl norm: int dPhi/(2pi) int dphi dcosth W_UU = 1
!           so we integrate W_UU × dcosth × dphi × dPhi/(2pi)
            wint = wint + diehl_W(act, apd, apl, eps, d) &
                        * dcth * dpd * dpl / (2d0*aPI)
          enddo
        enddo
      enddo

      return
      end


!==============================================================
!     gen_comparison: generate histogram data files
!==============================================================
      subroutine gen_comparison(r, d, eps, wmax_sw, wmax_di, &
                                nevt, nbin)
      implicit none
      real*4 urand
      real*8 sdme_W, diehl_W
      external sdme_W, diehl_W
      real*8 r(15), d(15), eps, wmax_sw, wmax_di
      integer nevt, nbin

      real*8 ct_sw(200), ct_di(200)
      real*8 pd_sw(200), pd_di(200)
      real*8 pl_sw(200), pl_di(200)
      real*8 costh, phid, aPhi, aPI, binval
      integer i, ib, iy

      aPI = acos(-1d0)

      do i = 1, nbin
        ct_sw(i) = 0d0; ct_di(i) = 0d0
        pd_sw(i) = 0d0; pd_di(i) = 0d0
        pl_sw(i) = 0d0; pl_di(i) = 0d0
      enddo

!     SW sampling
      iy = 33333
      do i = 1, nevt
        call sample_sdme_angles(costh, phid, aPhi, &
                                eps, r, wmax_sw, iy)
        ib = min(nbin, max(1, int((costh+1d0)/2d0*nbin)+1))
        ct_sw(ib) = ct_sw(ib) + 1d0
        ib = min(nbin, max(1, int(phid/(2d0*aPI)*nbin)+1))
        pd_sw(ib) = pd_sw(ib) + 1d0
        ib = min(nbin, max(1, int(aPhi/(2d0*aPI)*nbin)+1))
        pl_sw(ib) = pl_sw(ib) + 1d0
      enddo

!     Diehl sampling
      iy = 44444
      do i = 1, nevt
        call sample_diehl_angles(costh, phid, aPhi, &
                                 eps, d, wmax_di, iy)
        ib = min(nbin, max(1, int((costh+1d0)/2d0*nbin)+1))
        ct_di(ib) = ct_di(ib) + 1d0
        ib = min(nbin, max(1, int(phid/(2d0*aPI)*nbin)+1))
        pd_di(ib) = pd_di(ib) + 1d0
        ib = min(nbin, max(1, int(aPhi/(2d0*aPI)*nbin)+1))
        pl_di(ib) = pl_di(ib) + 1d0
      enddo

!     Write files
      open(30, file='diehl_cmp_costh.dat')
      open(31, file='diehl_cmp_phi.dat')
      open(32, file='diehl_cmp_Phi.dat')
      do i = 1, nbin
        binval = -1d0 + (dble(i)-0.5d0)*2d0/dble(nbin)
        write(30,'(F10.4,2(1X,F12.1))') binval, ct_sw(i), ct_di(i)
        binval = (dble(i)-0.5d0)*360d0/dble(nbin)
        write(31,'(F10.4,2(1X,F12.1))') binval, pd_sw(i), pd_di(i)
        write(32,'(F10.4,2(1X,F12.1))') binval, pl_sw(i), pl_di(i)
      enddo
      close(30); close(31); close(32)
      write(*,*) '  Written: diehl_cmp_costh/phi/Phi.dat'

      return
      end
