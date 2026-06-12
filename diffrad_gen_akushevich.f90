!==============================================================
!     diffrad_gen.f  -- version 2.0 (fully fixed)
!
!     MC Generator for diffractive vector meson electroproduction
!     with QED radiative corrections (collinear approximation)
!     Based on DIFFRAD by I.Akushevich (1998)
!
!     LUND output format (px py pz E mass):
!       1: beam e-      (initial, shifted if ISR)
!       2: target p     (initial, at rest)
!       3: scattered e- (final, shifted if FSR)
!       4: pi+          (from rho decay)
!       5: pi-          (from rho decay)
!       6: recoil p     (final)
!       7: gamma        (final, hard events only)
!==============================================================

      program diffrad_gen
      implicit real*8(a-h,o-z)
      real*4 urand

      common/cmp/pi,alpha,amp,amp2,ap,ap2,aml2,amc2,amv,barn
      common/bwpar/amv0,gamv
      common/sxy/s,x,sx,sxp,q2,w2,aly,anu,sqly,an,tamin,tamax,xs,ys
      common/phi/phirad,tdif,phidif,tq,vmax,ivec
      common/tpar/tslope
      common/vv1vv2/aa1,aa2,bb,sib
      common/ivv/vcurr,cutv
      common/cuts/wmin2
      common/amf2/taa,atm(8,6),sfm0(8)
      common/sigsig/sigmat0,sigmal0
      common/pri/ipri

      real*8 k1(4),k2(4),ptar(4),ph(4),pp(4),kgam(4),pip(4),pim(4)
      integer*4 iy
      integer*8 ntry, maxtry
      character(len=256) :: input_file, lund_file, stat_file, vdist_file
      character(len=256) :: cl_arg
      integer :: narg_tot, iarg_cur, llen

!     Parse command-line arguments: -input <file>  -lund <file>
      input_file = 'gen_input.dat'
      lund_file  = ''
      narg_tot   = command_argument_count()
      iarg_cur   = 1
      do while(iarg_cur .le. narg_tot)
        call get_command_argument(iarg_cur, cl_arg)
        if(trim(cl_arg) .eq. '-input') then
          iarg_cur = iarg_cur + 1
          call get_command_argument(iarg_cur, input_file)
        elseif(trim(cl_arg) .eq. '-lund') then
          iarg_cur = iarg_cur + 1
          call get_command_argument(iarg_cur, lund_file)
        endif
        iarg_cur = iarg_cur + 1
      enddo

!     Read input
      open(unit=8, file=trim(input_file), status='old')
      read(8,*) bmom
      read(8,*) tmom
      read(8,*) lepton
      read(8,*) ivec
      read(8,*) cutv
      read(8,*) nev
      read(8,*) iy
      read(8,*) q2min
      read(8,*) q2max
      read(8,*) ymin
      read(8,*) ymax
      read(8,*) tmin
      read(8,*) tmax
      read(8,*) tslope      ! exponential t-slope b (GeV^-2)
      read(8,*) iborn       ! 0=full RC,  1=Born only
      read(8,*) wmin        ! minimum W (GeV); events with W < wmin rejected
      close(8)
      wmin2 = wmin*wmin

      call setcon(ivec,lepton)
      s = 2d0*(sqrt(tmom**2+amp**2)*sqrt(bmom**2+aml2)+bmom*tmom)
      ebeam = bmom

      write(*,'(a)') '=============================='
      write(*,'(a,f8.3,a)') ' Ebeam  = ',ebeam,' GeV'
      write(*,'(a,f8.3,a)') ' sqrt(S)= ',sqrt(s),' GeV'
      write(*,'(a,i2)')     ' ivec   = ',ivec
      write(*,'(a,f6.3,a)') ' cutv   = ',cutv,' GeV^2'
      write(*,'(a,f6.3,a)') ' wmin   = ',wmin,' GeV'
      write(*,'(a,i8)')     ' nev    = ',nev
      if(iborn.eq.1) write(*,'(a)') ' Mode: BORN ONLY (no RC)'
      write(*,'(a)') '=============================='

!     Default lund filename if not given on command line
      if(len_trim(lund_file) .eq. 0) then
        if(iborn.eq.1) then
          lund_file = 'born_events.lund'
        else
          lund_file = 'rc_events.lund'
        endif
      endif
!     Derive stat and vdist filenames: strip .lund, append _stat.dat / _vdist.dat
      llen = len_trim(lund_file)
      if(llen .gt. 5 .and. lund_file(llen-4:llen) .eq. '.lund') then
        stat_file  = lund_file(1:llen-5)//'_stat.dat'
        vdist_file = lund_file(1:llen-5)//'_vdist.dat'
      else
        stat_file  = trim(lund_file)//'_stat.dat'
        vdist_file = trim(lund_file)//'_vdist.dat'
      endif
      write(*,'(a)') ' Output: '//trim(lund_file)
      write(*,'(a)') ' Stats:  '//trim(stat_file)
      open(unit=10, file=trim(lund_file))
      open(unit=11, file=trim(stat_file))
      open(unit=12, file=trim(vdist_file))

      ngen=0; nsoft=0; nhard=0; nisr=0; nfsr=0; ntry=0
      nf_born=0; nf_thresh=0; nf_tkin=0; nf_vmax=0
      nf_sib=0; nf_sshxxh=0; nf_sxqv=0; nf_sigtot=0
      wsum=0d0; wsum2=0d0   ! for cross section integration
      vcurr=0d0; ipri=0      ! for qqt/podinl
      wmax=0d0              ! maximum weight for accept/reject

!     ── Warm-up pass to find wmax ──────────────────────────────
      write(*,'(a)') ' Finding wmax (warm-up)...'
      nwarm = 10000
      do iwarm = 1, nwarm
        call sample_born(q2min,q2max,ymin,ymax,tmin,tmax,iy,sg_born,iacc)
        if(iacc.eq.0) cycle
        call conkin(s)
        call sample_bw(amv0,gamv,ivec,iy,sqrt(max(0d0,w2))-amp,amv)
        if(sqrt(w2).lt.amp+amv) cycle
        tt1=w2-q2-amp2; tt2=w2-amp2+amv**2
        disc1=tt1**2+4d0*q2*w2; disc2=tt2**2-4d0*amv**2*w2
        if(disc1.lt.0d0.or.disc2.lt.0d0) cycle
        tdmink=-q2+amv**2-.5d0/w2*(tt1*tt2+sqrt(disc1)*sqrt(disc2))
        tdmaxk=-q2+amv**2-.5d0/w2*(tt1*tt2-sqrt(disc1)*sqrt(disc2))
        if(tdif.lt.tdmink.or.tdif.gt.tdmaxk) cycle
        call bornin(sib)
        if(sib.le.0d0) cycle
!       Warm-up: use sg_born*sib for both Born and RC runs
!       For RC run, multiply by factor 2 as safety margin
!       (exact qqt can exceed Born but typically by <2x)
          wtrial = sg_born * sib
        if(wtrial.gt.wmax) wmax = wtrial
      enddo
      if(iborn.eq.1)then
        wmax = wmax * 1.5d0   ! 50% margin for Born
      else
        wmax = wmax * 3.0d0   ! 3x margin for RC (covers RC enhancement)
      endif
      write(*,'(a,g12.4)') ' wmax = ',wmax
      write(*,'(a)') ' Starting main generation...'

!     ── Main event loop ──────────────────────────────────────────
      do while (ngen .lt. nev)

        ntry = ntry + 1
!       Attempt limit: large for J/psi (low efficiency near threshold),
!       moderate for rho/phi/omega
        maxtry = 1000_8*nev
        if(ivec.eq.4) maxtry = 100000_8*nev
        if(ntry .gt. maxtry) then
          write(*,*) 'ERROR: too many attempts'
          goto 999
        endif

!       Step 1: Sample Born kinematics into common blocks
        call sample_born(q2min,q2max,ymin,ymax,tmin,tmax, &
                         iy,sg_born,iacc)
        if(iacc.eq.0)then; nf_born=nf_born+1; goto 100; endif

!       Step 2: Compute derived kinematics
        call conkin(s)

!       Sample vector meson mass from Breit-Wigner
        call sample_bw(amv0,gamv,ivec,iy,sqrt(max(0d0,w2))-amp,amv)

        if(sqrt(w2).lt.amp+amv)then
          nf_thresh=nf_thresh+1; goto 100
        endif

        tt1 = w2-q2-amp2
        tt2 = w2-amp2+amv**2
        disc1 = tt1**2+4d0*q2*w2
        disc2 = tt2**2-4d0*amv**2*w2
        if(disc1.lt.0d0.or.disc2.lt.0d0)then
          nf_born=nf_born+1; goto 100
        endif
        tdmink = -q2+amv**2-.5d0/w2*(tt1*tt2+sqrt(disc1)*sqrt(disc2))
        tdmaxk = -q2+amv**2-.5d0/w2*(tt1*tt2-sqrt(disc1)*sqrt(disc2))
        if(tdif.lt.tdmink.or.tdif.gt.tdmaxk)then
          nf_tkin=nf_tkin+1; goto 100
        endif

!       Step 3: RC quantities
        sxt = sx+tdif
        tq  = q2+tdif-amv**2
        aa1 = (q2*sxp*sxt-(s*sx+2d0*amp2*q2)*tq)/2d0/aly
        aa2 = (q2*sxp*sxt-(x*sx-2d0*amp2*q2)*tq)/2d0/aly
        sqbb1 = sqrt(max(0d0,q2*sxt**2-sxt*sx*tq-amp2*tq**2-amv**2*aly))
        sqbb2 = sqrt(max(0d0,q2*(s*x-amp2*q2)-aml2*aly))
        bb = sqbb1*sqbb2/aly

        vmax_kin = tt2+.5d0/q2*(-tt1*tq &
          -sqrt(max(0d0,tt1**2+4d0*q2*w2)) &
          *sqrt(max(0d0,tq**2+4d0*amv**2*q2))) &
          - 1d-8
        if(cutv.gt.1d-12)then
          vmax = min(vmax_kin,cutv)
        else
          vmax = vmax_kin
        endif
        if(vmax.le.0d0)then; nf_vmax=nf_vmax+1; goto 100; endif

        call bornin(sib)
        if(sib.le.0d0)then; nf_sib=nf_sib+1; goto 100; endif

        vv1 = aa1/2d0
        vv2 = aa2/2d0
        ssh = x+q2-vv2
        xxh = s-q2-vv1
        if(ssh.le.0d0.or.xxh.le.0d0)then
          nf_sshxxh=nf_sshxxh+1; goto 100
        endif

        dlm = log(q2/aml2)
        delinf_val = (dlm-1d0)*log(vmax**2/ssh/xxh)
        deltavr = (1.5d0*dlm-2d0-.5d0*log(xxh/ssh)**2 &
                   +fspen(1d0-amp2*q2/ssh/xxh)-pi**2/6d0)
        delta_vac = vacpol(q2)
        extai1 = exp(alpha/pi*delinf_val)
        sig_soft = sib*extai1*(1d0+alpha/pi*(deltavr+delta_vac))

!       Exact hard cross section via qqt (integrates podinl)
!       Only needed for RC run; skip for Born-only to avoid slow integration
        if(iborn.eq.0) then
          phidif = phirad
          call difflt(q2,w2,tdif,sigmal0,sigmat0)
          call qqt(sig_hard_exact)
          sig_hard = max(0d0, sig_hard_exact)
        else
          sig_hard = 0d0
        endif
        sig_total = sig_soft + sig_hard
        if(sig_total.le.0d0)then; nf_sigtot=nf_sigtot+1; goto 100; endif

        sg_born_full = sg_born * sib
        if(sg_born_full.le.0d0)then
          nf_sigtot=nf_sigtot+1
          goto 100
        endif

!       ── Accept/reject ─────────────────────────────────────
!       Born run: weight = sg_born * sib
!       RC run:   weight = sg_born * sig_total (includes RC)
        r_ar = dble(urand(iy))
        if(iborn.eq.1)then
          ar_weight = sg_born_full
        else
          ar_weight = sg_born * sig_total
        endif
        if(r_ar .gt. ar_weight/wmax) goto 100

        ngen = ngen + 1
        wsum  = wsum  + ar_weight
        wsum2 = wsum2 + ar_weight**2

!       Step 4: Born only or full RC?
        if(iborn.eq.1)then
!         BORN ONLY - skip all RC
          nsoft = nsoft + 1
          call build_4vectors(ebeam,xs,ys,tdif,phirad,k1,ptar,k2,ph,pp)
          call rotz_event(k2,ph,pp,kgam,pi,dble(urand(iy)))
          call write_lund(10,ngen,k1,ptar,k2,ph,pp,kgam, &
                          .false.,ivec,iy,pip,pim,ebeam,ar_weight)
          goto 100
        endif

!       Full RC: Soft or hard?
        prob_hard = sig_hard/sig_total
        r = dble(urand(iy))

        if(r .gt. prob_hard) then
!         NON-RADIATED EVENT
          nsoft = nsoft + 1
          call build_4vectors(ebeam,xs,ys,tdif,phirad,k1,ptar,k2,ph,pp)
          call rotz_event(k2,ph,pp,kgam,pi,dble(urand(iy)))
          call write_lund(10,ngen,k1,ptar,k2,ph,pp,kgam, &
                          .false.,ivec,iy,pip,pim,ebeam,ar_weight)

        else
!         HARD RADIATED EVENT
          nhard = nhard + 1

!         Sample v by accept/reject weighted by podinl (Bug #2 fix)
          vcut_use = 1d-4
          rvlogmin = log(vcut_use)
          rvlogmax = log(vmax)
          if(rvlogmax.le.rvlogmin)then
            nf_born=nf_born+1
            goto 100
          endif
          call sample_vrad(vcut_use,vmax,tamin,tamax,phidif,iy,vrad, &
                           iok_vrad)
          if(iok_vrad.eq.0)then
            nf_born=nf_born+1
            goto 100
          endif
          write(12,'(f12.6)') vrad

!         ISR vs FSR
          Ee1_val = ebeam
          Ee2_val = ebeam*(1d0-ys)
          p_isr   = (Ee2_val**2)/(Ee1_val**2+Ee2_val**2)

          r2 = dble(urand(iy))
          if(r2 .lt. p_isr) then
!           ISR
            nisr = nisr + 1
            omega  = vrad/(2d0*amp)      ! Bug #1 fix: v=2*Mp*omega
            Ee1_rc = Ee1_val - omega
            if(Ee1_rc.le.0d0)then
              nf_thresh=nf_thresh+1
              goto 100
            endif
            call build_4vectors(Ee1_rc,xs,ys,tdif,phirad, &
                                    k1,ptar,k2,ph,pp)
            kgam(1) = omega
            kgam(2) = 0d0
            kgam(3) = 0d0
            kgam(4) = omega
!           Restore original beam in LUND: pp was computed with Ee1_rc so
!           pp = k1_orig + ptar - kgam - k2 - rho (4-momentum conserved)
            k1(1) = Ee1_val
            k1(2) = 0d0
            k1(3) = 0d0
            k1(4) = sqrt(max(0d0,Ee1_val**2-aml2))
          else
!           FSR
            nfsr = nfsr + 1
            omega  = vrad/(2d0*amp)      ! Bug #1 fix: v=2*Mp*omega
            Ee2_rc = Ee2_val - omega
            if(Ee2_rc.le.0d0)then
              nf_thresh=nf_thresh+1
              goto 100
            endif
            call build_4vectors(ebeam,xs,ys,tdif,phirad, &
                                    k1,ptar,k2,ph,pp)
!           pp from build_4vectors is Born-level: pp = k1+ptar-k2_born-rho
!           Do NOT recompute pp with k2_fsr: k2_fsr+kgam=k2_born exactly,
!           so 4-momentum is conserved with pp_Born.
            ascale = Ee2_rc/Ee2_val
            k2(1) = k2(1)*ascale
            k2(2) = k2(2)*ascale
            k2(3) = k2(3)*ascale
            k2(4) = k2(4)*ascale
            kgam(1) = omega
            kgam(2) = (k2(2)/Ee2_rc)*omega
            kgam(3) = (k2(3)/Ee2_rc)*omega
            kgam(4) = (k2(4)/Ee2_rc)*omega
          endif

          call rotz_event(k2,ph,pp,kgam,pi,dble(urand(iy)))
          call write_lund(10,ngen,k1,ptar,k2,ph,pp,kgam, &
                          .true.,ivec,iy,pip,pim,ebeam,ar_weight)
        endif

  100   continue
      enddo

  999 continue

      write(*,'(a)') ''
      write(*,'(a)') '====== Debug failure counts ======'
      write(*,'(a,i8)') ' sample_born fail : ',nf_born
      write(*,'(a,i8)') ' W threshold fail : ',nf_thresh
      write(*,'(a,i8)') ' t kinemat  fail  : ',nf_tkin
      write(*,'(a,i8)') ' vmax<=0    fail  : ',nf_vmax
      write(*,'(a,i8)') ' sib<=0     fail  : ',nf_sib
      write(*,'(a,i8)') ' ssh/xxh    fail  : ',nf_sshxxh
      write(*,'(a,i8)') ' sx-qv      fail  : ',nf_sxqv
      write(*,'(a,i8)') ' sig_total  fail  : ',nf_sigtot
      write(*,'(a)') '=================================='
      write(*,'(a)') '====== Generator Statistics ======'
      anorm = 1d0/max(1_8,ntry)
!     With accept/reject: sigma = wmax * efficiency
      axsec = wmax * dble(ngen)/max(1_8,ntry)
      axsec_err = axsec/sqrt(dble(max(1,ngen)))
      write(*,'(a)') '====== Cross Section ======'
      write(*,'(a,g12.4,a)') ' sigma_Born   = ',axsec,' nb'
      write(*,'(a,g12.4,a)') ' stat error   = ',axsec_err,' nb'
      write(*,'(a,g12.4)')   ' efficiency   = ',dble(ngen)/max(1_8,ntry)
      write(*,'(a)') '=================================='
      write(*,'(a,i8)') ' Events generated : ',ngen
      write(*,'(a,i12)') ' Total attempts   : ',ntry
      write(*,'(a,i8)') ' Non-radiated     : ',nsoft
      write(*,'(a,i8)') ' Hard radiated    : ',nhard
      write(*,'(a,i8)') '   ISR events     : ',nisr
      write(*,'(a,i8)') '   FSR events     : ',nfsr
      write(*,'(a,f8.4)') ' Hard fraction    : ', &
                           dble(nhard)/max(1,ngen)
      write(*,'(a)') '=================================='

      write(11,'(a,i8)') 'ngen      = ',ngen
      write(11,'(a,i8)') 'nsoft     = ',nsoft
      write(11,'(a,i8)') 'nhard     = ',nhard
      write(11,'(a,f8.4)') 'hard_frac = ',dble(nhard)/max(1,ngen)
      write(11,'(a,g12.4)') 'sigma_nb  = ',axsec
      write(11,'(a,g12.4)') 'sigma_err = ',axsec_err

      close(10); close(11); close(12)
      end


!==============================================================
!     sample_born
!==============================================================
      subroutine sample_born(q2min,q2max,ymin,ymax,tmin,tmax, &
                             iy,sg_born,iacc)
      implicit real*8(a-h,o-z)
      real*4 urand
      common/cmp/pi,alpha,amp,amp2,ap,ap2,aml2,amc2,amv,barn
      common/sxy/s,x,sx,sxp,q2,w2,aly,anu,sqly,an,tamin,tamax,xs,ys
      common/phi/phirad,tdif,phidif,tq,vmax,ivec
      common/tpar/tslope
      common/cuts/wmin2
      integer*4 iy

      iacc = 0

      rq2min = log(q2min)
      rq2max = log(q2max)
      q2 = exp(rq2min + dble(urand(iy))*(rq2max-rq2min))
      wjacq2 = q2*(rq2max-rq2min)

      ys = ymin + dble(urand(iy))*(ymax-ymin)
      wjacy = ymax-ymin

      xs = q2/(s*ys)
      if(xs.ge.1d0.or.xs.le.0d0) return

!     W cuts: W^2 = M_p^2 + s*y - Q^2
      w2loc = amp2 + s*ys - q2
!     (1) hard physics threshold: W > M_p + M_V
      w2thr = (amp + amv)**2
      if(w2loc.lt.w2thr) return
!     (2) user analysis cut: W > wmin
      if(w2loc.lt.wmin2) return

      abslope = tslope
      tdmin_u = -tmax
      tdmax_u = -tmin
      abt1 = exp(abslope*tdmin_u)
      abt2 = exp(abslope*tdmax_u)
      abst = abt2 - abt1
      if(abst.le.0d0) return
      ranexp = abt1 + dble(urand(iy))*abst
      if(ranexp.le.0d0) return
      tdif  = log(ranexp)/abslope
      wjact = abst/abslope/exp(abslope*tdif)

      phirad = 2d0*pi*dble(urand(iy))

      sg_born = wjacq2 * wjacy * wjact * 2d0*pi
      iacc = 1
      return
      end


!==============================================================
!     build_4vectors
!     NOTE: all mass variables use 'am' prefix to avoid implicit integer
!     (m,M,k,j,l,n,i are integer by default in Fortran)
!==============================================================
      subroutine build_4vectors(ebeam,xs,ys,tdif,phirad, &
                                k1,ptar,k2,ph,pp)
      implicit real*8(a-h,o-z)
      common/cmp/pi,alpha,amp,amp2,ap,ap2,aml2,amc2,amv,barn
      common/sxy/s,x,sx,sxp,q2,w2,aly,anu,sqly,an,tamin,tamax,xs_c,ys_c
      real*8 k1(4),ptar(4),k2(4),ph(4),pp(4)
      real*8 qvec(3),qhat(3),e1v(3),e2v(3),ph3(3),tmp(3)

!     Beam electron (on-shell with mass)
      ap1_e = sqrt(max(0d0,ebeam**2-aml2))
      k1(1) = ebeam
      k1(2) = 0d0
      k1(3) = 0d0
      k1(4) = ap1_e

!     Target proton at rest
      ptar(1) = amp
      ptar(2) = 0d0
      ptar(3) = 0d0
      ptar(4) = 0d0

!     Scattered electron (on-shell with mass)
      aEe2 = ebeam*(1d0-ys)
      if(aEe2.le.0d0) return
      ap2_e = sqrt(max(0d0,aEe2**2-aml2))
      if(ap2_e.le.0d0) return
!     Q2 = 2*E1*E2 - 2*me^2 - 2*|p1|*|p2|*cos(theta)
      aq2_e = s*xs*ys
      acosthe = (2d0*ebeam*aEe2 - 2d0*aml2 - aq2_e) &
                / (2d0*ap1_e*ap2_e)
      if(abs(acosthe).gt.1d0) return
      asinthe = sqrt(max(0d0,1d0-acosthe**2))

      k2(1) = aEe2
      k2(2) = ap2_e*asinthe
      k2(3) = 0d0
      k2(4) = ap2_e*acosthe

!     Virtual photon
      aEq     = k1(1)-k2(1)
      qvec(1) = k1(2)-k2(2)
      qvec(2) = k1(3)-k2(3)
      qvec(3) = k1(4)-k2(4)
      aqmag   = sqrt(qvec(1)**2+qvec(2)**2+qvec(3)**2)

!     Rho meson 4-vector
!     Use 'amrho2' not 'Mrho2' -- M->m is IMPLICIT INTEGER!
      amrho2 = amv**2
      anu_loc = ebeam - aEe2
      aq2_loc = s*xs*ys

!     q.ph = (-Q2 + Mrho2 - t) / 2
      aqdotph = (-aq2_loc + amrho2 - tdif)/2d0

!     Erho from 4-momentum conservation
      aErho = (2d0*amp*anu_loc - aq2_loc + amrho2 - 2d0*aqdotph) &
              /(2d0*amp)
      if(aErho.lt.amv) return

      aphrho = sqrt(max(0d0,aErho**2-amrho2))

!     Angle of rho wrt virtual photon
      if(aqmag*aphrho.le.0d0) return
      acosalpha = (aEq*aErho - aqdotph)/(aqmag*aphrho)
      if(abs(acosalpha).gt.1d0) acosalpha=sign(1d0,acosalpha)
      asinalpha = sqrt(max(0d0,1d0-acosalpha**2))

!     Orthonormal basis around q
      qhat(1) = qvec(1)/aqmag
      qhat(2) = qvec(2)/aqmag
      qhat(3) = qvec(3)/aqmag

      tmp(1) = 0d0; tmp(2) = 1d0; tmp(3) = 0d0
      call cross3(qhat,tmp,e1v)
      aemag = sqrt(e1v(1)**2+e1v(2)**2+e1v(3)**2)
      if(aemag.lt.1d-10) then
        tmp(1)=0d0; tmp(2)=0d0; tmp(3)=1d0
        call cross3(qhat,tmp,e1v)
        aemag=sqrt(e1v(1)**2+e1v(2)**2+e1v(3)**2)
      endif
      e1v(1)=e1v(1)/aemag; e1v(2)=e1v(2)/aemag; e1v(3)=e1v(3)/aemag
      call cross3(qhat,e1v,e2v)

!     Rho 3-momentum direction
      ph3(1) = aphrho*(acosalpha*qhat(1) &
              + asinalpha*cos(phirad)*e1v(1) &
              + asinalpha*sin(phirad)*e2v(1))
      ph3(2) = aphrho*(acosalpha*qhat(2) &
              + asinalpha*cos(phirad)*e1v(2) &
              + asinalpha*sin(phirad)*e2v(2))
      ph3(3) = aphrho*(acosalpha*qhat(3) &
              + asinalpha*cos(phirad)*e1v(3) &
              + asinalpha*sin(phirad)*e2v(3))

!     Put rho on mass shell: E = sqrt(|p|^2 + Mrho^2)
      ph(2) = ph3(1)
      ph(3) = ph3(2)
      ph(4) = ph3(3)
      ph(1) = sqrt(ph(2)**2+ph(3)**2+ph(4)**2+amrho2)

!     Recoil proton from 4-momentum conservation
      pp(1) = k1(1)+ptar(1)-k2(1)-ph(1)
      pp(2) = k1(2)+ptar(2)-k2(2)-ph(2)
      pp(3) = k1(3)+ptar(3)-k2(3)-ph(3)
      pp(4) = k1(4)+ptar(4)-k2(4)-ph(4)

      return
      end


!==============================================================
!     cross3
!==============================================================
      subroutine cross3(a,b,c)
      implicit real*8(a-h,o-z)
      dimension a(3),b(3),c(3)
      c(1) = a(2)*b(3)-a(3)*b(2)
      c(2) = a(3)*b(1)-a(1)*b(3)
      c(3) = a(1)*b(2)-a(2)*b(1)
      return
      end


!==============================================================
!     decay_rho: vector meson -> h+ h- isotropic in rest frame
!       ivec=1 (rho): pi+ pi-  (M_pi = 0.13957 GeV)
!       ivec=3 (phi): K+  K-   (M_K  = 0.49368 GeV)
!       others: pi+ pi- (default)
!==============================================================
      subroutine decay_rho(ph,pip,pim,iy,ivec)
      implicit real*8(a-h,o-z)
      real*4 urand
      common/cmp/pi,alpha,amp,amp2,ap,ap2,aml2,amc2,amv,barn
      real*8 ph(4),pip(4),pim(4)
      real*8 abeta(3),pip_rf(4),pim_rf(4)
      integer*4 iy, ivec

      if(ivec.eq.3) then
        ampi = 0.493677d0   ! K+/K- mass
      elseif(ivec.eq.4) then
        ampi = 0.000511d0   ! e+/e- mass
      else
        ampi = 0.13957d0    ! pi+/pi- mass
      endif
      amrho = amv

!     Pion momentum in rho rest frame
      appion = sqrt(max(0d0,(amrho/2d0)**2-ampi**2))

!     Random decay direction
      acosth = 2d0*dble(urand(iy))-1d0
      asinth = sqrt(max(0d0,1d0-acosth**2))
      aphid  = 2d0*pi*dble(urand(iy))

      pip_rf(1) = amrho/2d0
      pip_rf(2) = appion*asinth*cos(aphid)
      pip_rf(3) = appion*asinth*sin(aphid)
      pip_rf(4) = appion*acosth

      pim_rf(1) = amrho/2d0
      pim_rf(2) = -pip_rf(2)
      pim_rf(3) = -pip_rf(3)
      pim_rf(4) = -pip_rf(4)

!     Boost to lab (rho velocity)
      aErho   = ph(1)
      abeta(1) = ph(2)/aErho
      abeta(2) = ph(3)/aErho
      abeta(3) = ph(4)/aErho
      abeta2   = abeta(1)**2+abeta(2)**2+abeta(3)**2

!     Use actual rho mass (ph is on mass shell from build_4vectors)
      agamma = aErho/amrho

      call lorentz_boost(pip_rf,abeta,agamma,abeta2,pip)
      call lorentz_boost(pim_rf,abeta,agamma,abeta2,pim)

      return
      end


!==============================================================
!     lorentz_boost: p_out = boost(p_in) by velocity abeta
!     Standard formula: E'=g(E+b.p), p'=p+(g-1)/b2*(b.p)*b+g*E*b
!==============================================================
      subroutine lorentz_boost(p_in,abeta,agamma,abeta2,p_out)
      implicit real*8(a-h,o-z)
      real*8 p_in(4),abeta(3),p_out(4)

      if(abeta2.lt.1d-20) then
        p_out(1)=p_in(1); p_out(2)=p_in(2)
        p_out(3)=p_in(3); p_out(4)=p_in(4)
        return
      endif

      abdotp = abeta(1)*p_in(2)+abeta(2)*p_in(3)+abeta(3)*p_in(4)
      acoef  = (agamma-1d0)/abeta2*abdotp + agamma*p_in(1)

      p_out(1) = agamma*(p_in(1)+abdotp)
      p_out(2) = p_in(2) + acoef*abeta(1)
      p_out(3) = p_in(3) + acoef*abeta(2)
      p_out(4) = p_in(4) + acoef*abeta(3)

      return
      end


!==============================================================
!     rotz_event  --  Random rotation around beam (z) axis
!     Rotates k2, ph, pp, kgam by uniform random angle in [0,2pi]
!     k1 (beam along z) is invariant under this rotation
!==============================================================
      subroutine rotz_event(k2,ph,pp,kgam,pi,rnd)
      implicit real*8(a-h,o-z)
      real*8 k2(4),ph(4),pp(4),kgam(4)
      real*8 cphi,sphi,px_tmp
      phi_rot = 2d0*pi*rnd
      cphi = cos(phi_rot)
      sphi = sin(phi_rot)
!     Rotate each 4-vector: (E unchanged, px'=px*c-py*s, py'=px*s+py*c)
      px_tmp = k2(2)*cphi - k2(3)*sphi
      k2(3)  = k2(2)*sphi + k2(3)*cphi
      k2(2)  = px_tmp
      px_tmp = ph(2)*cphi - ph(3)*sphi
      ph(3)  = ph(2)*sphi + ph(3)*cphi
      ph(2)  = px_tmp
      px_tmp = pp(2)*cphi - pp(3)*sphi
      pp(3)  = pp(2)*sphi + pp(3)*cphi
      pp(2)  = px_tmp
      px_tmp   = kgam(2)*cphi - kgam(3)*sphi
      kgam(3)  = kgam(2)*sphi + kgam(3)*cphi
      kgam(2)  = px_tmp
      return
      end

!==============================================================
!     write_lund  --  CLAS12 GEMC LUND format
!     Header (10 mandatory + 5 optional fields):
!       npart  1  1  0  0  11  Ebeam  2212  1  sigma  xB  y  W2  Q2  nu
!     Body (14 fields per particle):
!       idx  lifetime  type  pid  parent  daughter  px py pz E mass  vx vy vz
!     type=1: sent to Geant4;  type=0: not sent (beam e-, vector meson)
!     Particle order: 1=beam e-, 2=scattered e-, 3=recoil p,
!                     4=vector meson, 5=h+, 6=h-, [7=photon]
!==============================================================
      subroutine write_lund(lun,iev,k1,ptar,k2,ph,pp,kgam, &
                            has_photon,ivec,iy,pip,pim,ebeam,sigma_ev)
      implicit real*8(a-h,o-z)
      common/cmp/pi,alpha,amp,amp2,ap,ap2,aml2,amc2,amv,barn
      real*8 k1(4),ptar(4),k2(4),ph(4),pp(4),kgam(4),pip(4),pim(4)
      logical has_photon
      integer*4 iy
      integer pidp, pidm, pidmeson

      call decay_rho(ph,pip,pim,iy,ivec)

!     npart: beam e-, scattered e-, recoil p, meson, h+, h-, [photon]
      npart = 6
      if(has_photon) npart = 7

!     Kinematics for optional header fields
      anu_ev = k1(1) - k2(1)
      q2_ev  = 2d0*(k1(1)*k2(1) - k1(2)*k2(2) &
                   -k1(3)*k2(3) - k1(4)*k2(4)) - 2d0*aml2
      w2_ev  = amp2 + 2d0*amp*anu_ev - q2_ev
      yy_ev  = anu_ev / ebeam
      xb_ev  = q2_ev / (2d0*amp*anu_ev)

!     Header: npart 1 1 0 0 11 Ebeam 2212 1 sigma xB y W2 Q2 nu
      write(lun,'(i5,4i3,i5,f8.3,i6,i3,e14.6,5f9.4)') &
            npart, 1, 1, 0, 0, &
            11, ebeam, 2212, 1, sigma_ev, &
            xb_ev, yy_ev, w2_ev, q2_ev, anu_ev

!     Meson PID and daughter masses/PIDs
      if(ivec.eq.3) then
!       phi -> K+ K-
        pidmeson = 333
        amh  = 0.493677d0
        pidp =  321
        pidm = -321
      elseif(ivec.eq.4) then
!       J/psi -> e+ e-
        pidmeson = 443
        amh  = 0.000511d0
        pidp = -11
        pidm =  11
      else
!       rho/omega -> pi+ pi-
        pidmeson = 113
        amh  = 0.13957d0
        pidp =  211
        pidm = -211
      endif

      aml_ev = sqrt(aml2)   ! electron mass

!     Body: idx  lifetime  type  pid  parent  daughter  px py pz E mass  vx vy vz
!     1: beam electron (type=0, not to Geant4)
      write(lun,'(i4,f6.1,i3,i6,2i3,4f10.4,f12.6,3f8.3)') &
            1, 0.0d0, 0, 11, 0, 0, &
            k1(2),k1(3),k1(4),k1(1), aml_ev, 0d0,0d0,0d0
!     2: scattered electron (type=1)
      write(lun,'(i4,f6.1,i3,i6,2i3,4f10.4,f12.6,3f8.3)') &
            2, 0.0d0, 1, 11, 0, 0, &
            k2(2),k2(3),k2(4),k2(1), aml_ev, 0d0,0d0,0d0
!     3: recoil proton (type=1)
      write(lun,'(i4,f6.1,i3,i6,2i3,4f10.4,f12.6,3f8.3)') &
            3, 0.0d0, 1, 2212, 0, 0, &
            pp(2),pp(3),pp(4),pp(1), amp, 0d0,0d0,0d0
!     4: vector meson (type=0, not to Geant4, daughter=5)
      write(lun,'(i4,f6.1,i3,i6,2i3,4f10.4,f12.6,3f8.3)') &
            4, 0.0d0, 0, pidmeson, 0, 5, &
            ph(2),ph(3),ph(4),ph(1), amv, 0d0,0d0,0d0
!     5: positive daughter (type=1, parent=4)
      write(lun,'(i4,f6.1,i3,i6,2i3,4f10.4,f12.6,3f8.3)') &
            5, 0.0d0, 1, pidp, 4, 0, &
            pip(2),pip(3),pip(4),pip(1), amh, 0d0,0d0,0d0
!     6: negative daughter (type=1, parent=4)
      write(lun,'(i4,f6.1,i3,i6,2i3,4f10.4,f12.6,3f8.3)') &
            6, 0.0d0, 1, pidm, 4, 0, &
            pim(2),pim(3),pim(4),pim(1), amh, 0d0,0d0,0d0
!     7: radiated photon (type=1, RC events only)
      if(has_photon) then
        write(lun,'(i4,f6.1,i3,i6,2i3,4f10.4,f12.6,3f8.3)') &
              7, 0.0d0, 1, 22, 0, 0, &
              kgam(2),kgam(3),kgam(4),kgam(1), 0d0, 0d0,0d0,0d0
      endif

      return
      end


!==============================================================
!     conkin
!==============================================================
      subroutine conkin(snuc)
      implicit real*8(a-h,o-z)
      common/cmp/pi,alpha,amp,amp2,ap,ap2,aml2,amc2,amv,barn
      common/sxy/s,x,sx,sxp,q2,w2,aly,anu,sqly,an,tamin,tamax,xs,ys
      s=snuc; x=s*(1d0-ys); sx=s-x; sxp=s+x
      w2=amp2+s-q2-x; aly=sx**2+4d0*amp2*q2; sqly=sqrt(aly)
      anu=sx/ap; an=alpha*ys/(8d0*pi**2)*barn
      tamax=(sx+sqly)/ap2; tamin=-q2/amp2/tamax
      return
      end


!==============================================================
!     sample_bw: sample vector meson mass from Breit-Wigner
!       amv0   -- nominal mass (GeV)
!       gamv   -- full width (GeV)
!       iy     -- random seed
!       amwmax -- kinematic upper limit: sqrt(W2) - Mp
!       amv    -- output: sampled mass (GeV)
!
!     Uses: M^2 = M0^2 + M0*Gamma*tan(theta),
!           theta uniform in [theta_min, theta_max]
!           with M_min = 2*M_pi, M_max = min(amwmax, M0+5*Gamma)
!==============================================================
      subroutine sample_bw(amv0,gamv,ivec,iy,amwmax,amv)
      implicit real*8(a-h,o-z)
      real*4 urand
      integer*4 iy, ivec
!     Sample M from Breit-Wigner with running width (rho only)
!     For rho (ivec=1): Gamma(M) = Gamma0*(p(M)/p(M0))^3*(M0/M)
!     For others: constant width
      if(ivec.eq.3) then
        ampi = 0.493677d0   ! phi -> K+K-
      elseif(ivec.eq.4) then
        ampi = 0.000511d0   ! J/psi -> e+e-
      else
        ampi = 0.13957d0    ! rho/omega -> pi+pi-
      endif
      ammin = 2d0*ampi
      ammax = min(amwmax, amv0+5d0*gamv)
!     Treat as delta function if width < 1 MeV (catches J/psi: 92.9 keV)
      if(ammax.le.ammin .or. gamv.lt.1d-3) then
        amv = amv0
        return
      endif
!     For rho: integrate using constant-width BW as envelope,
!     then apply accept/reject for running width correction
      p0 = sqrt(max(0d0, amv0**2/4d0 - ampi**2))
      thmin = atan((ammin**2-amv0**2)/(amv0*gamv))
      thmax = atan((ammax**2-amv0**2)/(amv0*gamv))
      ibw_try = 0
  10  continue
      ibw_try = ibw_try + 1
      if(ibw_try .gt. 10000) then
!       Give up: return nominal mass so the caller can reject the event
        amv = amv0
        return
      endif
      theta = thmin + dble(urand(iy))*(thmax-thmin)
      amv2  = amv0**2 + amv0*gamv*tan(theta)
      amv   = sqrt(max(ammin**2, amv2))
      if(ivec.eq.1) then
!       Running width for rho: Gamma(M) = Gamma0*(p/p0)^3*(M0/M)
        pm = sqrt(max(0d0, amv**2/4d0 - ampi**2))
        if(p0.gt.0d0) then
          gamrun = gamv*(pm/p0)**3*(amv0/amv)
        else
          gamrun = gamv
        endif
!       Accept/reject: w(M) = BW(M,Gamma_run) / BW(M,Gamma0)
        bw_run  = 1d0/((amv2-amv0**2)**2 + amv0**2*gamrun**2)
        bw_const= 1d0/((amv2-amv0**2)**2 + amv0**2*gamv**2)
        ratio   = bw_run / bw_const * (gamrun/gamv)**2
        if(dble(urand(iy)).gt.ratio) goto 10
      endif
      end

!==============================================================
!     Akushevich physics routines (unchanged)
!==============================================================
      subroutine setcon(ivec,lepton)
      implicit real*8(a-h,o-z)
      common/cmp/pi,alpha,amp,amp2,ap,ap2,aml2,amc2,amv,barn
      common/bwpar/amv0,gamv
      dimension amhad(4),gamhad(4)
      data amhad /0.7683d0,   0.78195d0,  1.019412d0, 3.0969d0  /
      data gamhad/0.1502d0,   0.00849d0,  0.004266d0, 0.0000929d0/
      if(lepton.eq.1)aml2=.261112d-6
      if(lepton.eq.2)aml2=.111637d-1
      pi=3.1415926d0; alpha=.729735d-2; barn=.389379d6
      amv0=amhad(ivec); gamv=gamhad(ivec)
      amv=amv0; amp=.938272d0; amp2=amp**2
      ap=2d0*amp; ap2=2d0*amp2; amc2=amp2
      end

      subroutine bornin(sibor)
      implicit real*8(a-h,o-z)
      common/cmp/pi,alpha,amp,amp2,ap,ap2,aml2,amc2,amv,barn
      common/sxy/s,x,sx,sxp,q2,w2,aly,anu,sqly,an,tamin,tamax,xs,ys
      common/phi/phirad,tdif,phidif,tq,vmax,ivec
      common/sigsig/sigmat,sigmal
      common/pri/ipri
      aga2=q2/anu**2
      ipri=0
      call difflt(q2,w2,tdif,sigmal,sigmat)
      sibor=2d0*an/(xs*ys**2)*(ys**2*sigmat &
            +2d0*(1d0-ys-.25d0*ys**2*aga2)*(sigmal+sigmat))
      end

      subroutine difflt(q2,w2,t,sigl,sigt)
      implicit real*8(a-h,o-z)
      common/cmp/pi,alpha,amp,amp2,ap,ap2,aml2,amc2,amv,barn
      common/phi/phirad,tdif,phidif,tq,vmax,ivec
      common/tpar/tslope
      common/pri/ipri
      dimension ggam(4)
      data p02/.5d0/al_s/.25d0/
      data ggam/6.77d-6,0.6d-6,1.37d-6,5.36d-6/
      asx=w2+q2-amp2; aanu=asx/ap; asxt=asx+t; aeh=asxt/2d0/amp
      amv2=amv**2; aeta=1d0; aff2=exp(tslope*t)
      atqt=t+q2-amv2
      apt2=(-(4d0*(aanu**2+q2)*amv2+4d0*aanu*aeh*atqt &
           -4d0*aeh**2*q2+atqt**2))/(4d0*(aanu**2+q2))
      if(apt2.lt.0d0)apt2=0d0
      if(w2.lt.(amp+amv)**2)then; sigl=0d0; sigt=0d0; return; endif
      axsb=(q2+amv2+apt2)/w2; aq2b=(q2+amv2+apt2)/4d0
      if(apt2.le.p02)then
        afm=log((4d0*aq2b-apt2+p02)/(apt2+p02))
      else
        afm=log((apt2+p02)/(4d0*aq2b-apt2+p02)*4d0*aq2b**2/apt2**2)
      endif
      axsbgm=3d0*(1d0-axsb)**5
      ask=axsbgm*afm/(2d0*aq2b*(2d0*aq2b-apt2)*log(8d0*aq2b/p02))
      sigt=al_s**2*ggam(ivec)*amv**3/3d0/alpha*pi**3*ask**2*aff2*aeta**2
      sigl=q2/amv2*sigt
      end

      double precision function vacpol(t)
      implicit real*8(a-h,o-z)
      common/cmp/pi,alpha,amp,amp2,ap,ap2,aml2,amc2,amv,barn
      dimension am2(3)
      data am2/.26110d-6,.111637d-1,3.18301d0/
      asuml=0d0
      do 10 i=1,3
        aa2=2d0*am2(i)
        asqlmi=sqrt(t*t+2d0*aa2*t)
        aallmi=log((asqlmi+t)/(asqlmi-t))/asqlmi
  10  asuml=asuml+2d0*(t+aa2)*aallmi/3d0-10d0/9d0 &
                  +4d0*aa2*(1d0-aa2*aallmi)/3d0/t
      if(t.lt.1d0)then
        aaaa=-1.345d-9; abbb=-2.302d-3; accc=4.091d0
      elseif(t.lt.64d0)then
        aaaa=-1.512d-3; abbb=-2.822d-3; accc=1.218d0
      else
        aaaa=-1.1344d-3; abbb=-3.0680d-3; accc=9.9992d-1
      endif
      asumh=-(aaaa+abbb*log(1d0+accc*t))*2d0*pi/alpha
      vacpol=asuml+asumh
      end

      double precision function fspens(x)
      implicit real*8(a-h,o-z)
      af=0d0; aa=1d0; aan=0d0; atch=1d-16
  1   aan=aan+1d0; aa=aa*x; ab=aa/aan**2; af=af+ab
      if(ab-atch)2,2,1
  2   fspens=af
      return
      end

      double precision function fspen(x)
      implicit real*8(a-h,o-z)
      data af1/1.644934d0/
      if(x)8,1,1
  1   if(x-.5d0)2,2,3
  2   fspen=fspens(x); return
  3   if(x-1d0)4,4,5
  4   fspen=af1-log(x)*log(1d0-x+1d-10)-fspens(1d0-x); return
  5   if(x-2d0)6,6,7
  6   fspen=af1-.5d0*log(x)*log((x-1d0)**2/x)+fspens(1d0-1d0/x)
      return
  7   fspen=2d0*af1-.5d0*log(x)**2-fspens(1d0/x); return
  8   if(x+1d0)10,9,9
  9   fspen=-.5d0*log(1d0-x)**2-fspens(x/(x-1d0)); return
 10   fspen=-.5d0*log(1d0-x)*log(x**2/(1d0-x)) &
            -af1+fspens(1d0/(1d0-x))
      return
      end


!==============================================================
!  Exact RC subroutines from Akushevich idiffrad.f
!  qqt -> qqtphi -> rv2ln -> podinl
!  These compute sigma_hard exactly (eq.15 of the paper)
!==============================================================

      subroutine qqt(tai)
      implicit real*8(a-h,o-z)
      external qqtphi
      real*8 qqtphi
      common/cmp/pi,alpha,amp,amp2,ap,ap2,aml2,amc2,amv,barn
!     Use 32-point Gauss quadrature for smooth phi integrand
      call dqg32(0d0,2.d0*pi,qqtphi,tai)
      tai=tai/2.d0
      end

      double precision function qqtphi(phi)
      implicit real*8(a-h,o-z)
      external rv2ln
      real*8 rv2ln
      common/cmp/pi,alpha,amp,amp2,ap,ap2,aml2,amc2,amv,barn
      common/sxy/s,x,sx,sxp,q2,w2,aly,anu,sqly,an,tamin,tamax,xs,ys
      common/phi/phirad,tdif,phidif,tq,vmax,ivec
      common/ivv/vcurr,cutv
      dimension tlm(4)
      phirad=phi
      tlm(1)=log(xs+tamin)
      tlm(4)=log(xs+tamax)
      tlm(2)=log(xs-q2/s)
      tlm(3)=log(xs+q2/x)
      res=0d0
      do ii=1,3
        ep=1d-10
        call simptx(tlm(ii)+ep,tlm(ii+1)-ep,10,1d-2,rv2ln,re)
        tai=an*alpha/pi*re
        res=res+tai
      enddo
      qqtphi=res
      end

      double precision function rv2ln(taln)
      implicit real*8(a-h,o-z)
      external podinl
      real*8 podinl
      common/cmp/pi,alpha,amp,amp2,ap,ap2,aml2,amc2,amv,barn
      common/sxy/s,x,sx,sxp,q2,w2,aly,anu,sqly,an,tamin,tamax,xs,ys
      common/amf2/taa,atm(8,6),sfm0(8)
      common/phi/phirad,tdif,phidif,tq,vmax,ivec
      common/ivv/vcurr,cutv
      ta=exp(taln)-q2/sx
      taa=ta
      sqrtmb=sqrt((ta-tamin)*(tamax-ta)*(s*x*q2-q2**2*amp2-aml2*aly))
      z1=(q2*sxp+ta*(s*sx+ap2*q2)-ap*cos(phirad)*sqrtmb)/aly
      z2=(q2*sxp+ta*(x*sx-ap2*q2)-ap*cos(phirad)*sqrtmb)/aly
      abb=1.d0/sqly/pi
      bi12=abb/(z1*z2)
      bi1pi2=abb/z2+abb/z1
      abis=abb/z2**2+abb/z1**2
      abir=abb/z2**2-abb/z1**2
      ahi2=aml2*abis-q2*bi12
      atm(1,1)=4.d0*q2*ahi2
      atm(1,2)=4.d0*ahi2*ta
      atm(1,3)=-2.d0*(2.d0*abb+bi12*ta**2)
      atm(2,1)=2d0*(s*x-amp2*q2)*ahi2/amp2
      atm(2,2)=(2.d0*aml2*abir*sxp-4.d0*amp2*ahi2*ta-bi12*sxp**2*ta+ &
                bi1pi2*sxp*sx+2.d0*ahi2*sx)/(2.d0*amp2)
      atm(2,3)=(2.d0*(2.d0*abb+bi12*ta**2)*amp2-bi12* &
                sx*ta-bi1pi2*sxp)/(2.d0*amp2)
      vmin=1d-4
      call simpux(vmin,vmax,10,5d-3,podinl,res)
      rv2ln=res*(q2/sx+ta)
      end

      double precision function podinl(v)
      implicit real*8(a-h,o-z)
      common/cmp/pi,alpha,amp,amp2,ap,ap2,aml2,amc2,amv,barn
      common/sxy/s,x,sx,sxp,q2,w2,aly,anu,sqly,an,tamin,tamax,xs,ys
      common/amf2/ta,atm(8,6),sfm0(8)
      common/phi/phirad,tdif,phidif,tq,vmax,ivec
      common/sigsig/sigmat0,sigmal0
      common/ivv/vcurr,cutv
      dimension sfm(8)
      sxtm=sx+tdif-v
      aak=((2d0*q2+ta*sx)*sxtm-(sx-2d0*ta*amp2)*tq)/aly/2d0
      sqbbk1=sqrt(max(0d0,q2*sxtm**2-sxtm*sx*tq-amp2*tq**2-amv**2*aly))
      sqbbk2=amp*sqrt(max(0d0,(tamax-ta)*(ta-tamin)))
      abbk=sqbbk1*sqbbk2/aly
      d2kvir=2d0*(aak+abbk*cos(phirad-phidif))
      factor=1d0+ta-d2kvir
      if(abs(factor).lt.1d-10)then; podinl=0d0; return; endif
      r=v/factor
      tldq2=q2+r*ta
      tldw2=w2-r*(1.d0+ta)
      tldtd=tdif-r*(ta-d2kvir)
      call difflt(tldq2,tldw2,tldtd,sigmal,sigmat)
      sfm(1)=(sx-r)*sigmat
      sfm(2)=2.d0*ap2/(sx-r)*tldq2*(sigmat+sigmal)
      sfm0(1)=sx*sigmat0
      sfm0(2)=2.d0*ap2/sx*q2*(sigmat0+sigmal0)
      podinl=0.d0
      do isf=1,2
        do irr=1,3
          app=sfm(isf)
          if(irr.eq.1)app=app-sfm0(isf)*(1.d0+r*ta/q2)**2
          pres=app*r**(irr-2)/(q2+r*ta)**2/2.d0
          podinl=podinl-atm(isf,irr)*pres
        enddo
      enddo
      podinl=podinl/factor
      end

!--------------------------------------------------------------
!  sample_vrad: sample vrad by accept/reject weighted by podinl.
!  Bug #2 fix: replaces log-flat sampling with physical distribution.
!  podinl is declared external here, NOT in the main program,
!  to avoid interfering with simpux/rv2ln function-argument passing.
!--------------------------------------------------------------
      subroutine sample_vrad(vmin_in,vmax_in,tamin_in,tamax_in, &
                             phidif_in,iy,vrad_out,iok)
      implicit real*8(a-h,o-z)
      external podinl
      real*8 podinl
      real*4 urand
      integer*4 iy
      common/cmp/pi,alpha,amp,amp2,ap,ap2,aml2,amc2,amv,barn
      common/sxy/s,x,sx,sxp,q2,w2,aly,anu,sqly,an,tamin,tamax,xs,ys
      common/amf2/taa,atm(8,6),sfm0(8)
      common/phi/phirad,tdif,phidif,tq,vmax,ivec

!     Set up atm at near-collinear ta with phirad = phidif
      phirad = phidif_in
      ta_c = tamin_in + 1d-3*(tamax_in-tamin_in)
      if(ta_c.gt.tamax_in) ta_c = 0.5d0*(tamin_in+tamax_in)

!     Populate atm for ta_c (mirrors rv2ln logic)
      ta   = ta_c
      taa  = ta
      sqrtmb = sqrt(max(0d0,(ta-tamin_in)*(tamax_in-ta) &
                    *(s*x*q2-q2**2*amp2-aml2*aly)))
      z1 = (q2*sxp+ta*(s*sx+ap2*q2)-ap*cos(phirad)*sqrtmb)/aly
      z2 = (q2*sxp+ta*(x*sx-ap2*q2)-ap*cos(phirad)*sqrtmb)/aly
      if(abs(z1).lt.1d-30.or.abs(z2).lt.1d-30)then
        iok = 0
        return
      endif
      abb    = 1.d0/sqly/pi
      bi12   = abb/(z1*z2)
      bi1pi2 = abb/z2+abb/z1
      abis   = abb/z2**2+abb/z1**2
      abir   = abb/z2**2-abb/z1**2
      ahi2   = aml2*abis-q2*bi12
      atm(1,1) = 4.d0*q2*ahi2
      atm(1,2) = 4.d0*ahi2*ta
      atm(1,3) = -2.d0*(2.d0*abb+bi12*ta**2)
      atm(2,1) = 2d0*(s*x-amp2*q2)*ahi2/amp2
      atm(2,2) = (2.d0*aml2*abir*sxp-4.d0*amp2*ahi2*ta &
                  -bi12*sxp**2*ta+bi1pi2*sxp*sx &
                  +2.d0*ahi2*sx)/(2.d0*amp2)
      atm(2,3) = (2.d0*(2.d0*abb+bi12*ta**2)*amp2 &
                  -bi12*sx*ta-bi1pi2*sxp)/(2.d0*amp2)

!     Scan podinl over v to find local maximum
      rvlogmin = log(vmin_in)
      rvlogmax = log(vmax_in)
      pmax = 0d0
      do iscan = 1, 50
        frac  = dble(iscan-1)/49d0
        vscan = exp(rvlogmin + frac*(rvlogmax-rvlogmin))
        pval  = podinl(vscan)
        if(pval.gt.pmax) pmax = pval
      enddo
      pmax = pmax*2d0   ! safety margin

      if(pmax.le.0d0)then
        iok = 0
        return
      endif

!     Accept/reject loop
      vrad_out = exp(rvlogmin+dble(urand(iy))*(rvlogmax-rvlogmin))
      do itry = 1, 10000
        vtry = exp(rvlogmin+dble(urand(iy))*(rvlogmax-rvlogmin))
        pval = max(0d0, podinl(vtry))
        if(dble(urand(iy)).lt.pval/pmax)then
          vrad_out = vtry
          exit
        endif
      enddo
      iok = 1
      end

      subroutine simpsx(a,b,np,ep,func,res)
      implicit real*8(a-h,o-z)
      external func
      step=(b-a)/np
      call simps(a,b,step,ep,1d-18,func,ra,res,r2,r3)
      end

      subroutine simptx(a,b,np,ep,func,res)
      implicit real*8(a-h,o-z)
      external func
      step=(b-a)/np
      call simpt(a,b,step,ep,1d-18,func,ra,res,r2,r3)
      end

      subroutine simpux(a,b,np,ep,func,res)
      implicit real*8(a-h,o-z)
      external func
      step=(b-a)/np
      call simpu(a,b,step,ep,1d-18,func,ra,res,r2,r3)
      end

      subroutine dqg32(xl,xu,fct,y)
      double precision xl,xu,y,a,b,c,fct
      a=.5d0*(xu+xl); b=xu-xl
      c=.49863193092474078d0*b
      y=.35093050047350483d-2*(fct(a+c)+fct(a-c))
      c=.49280575577263417d0*b
      y=y+.8137197365452835d-2*(fct(a+c)+fct(a-c))
      c=.48238112779375322d0*b
      y=y+.12696032654631030d-1*(fct(a+c)+fct(a-c))
      c=.46745303796886984d0*b
      y=y+.17136931456510717d-1*(fct(a+c)+fct(a-c))
      c=.44816057788302606d0*b
      y=y+.21417949011113340d-1*(fct(a+c)+fct(a-c))
      c=.42468380686628499d0*b
      y=y+.25499029631188088d-1*(fct(a+c)+fct(a-c))
      c=.39724189798397120d0*b
      y=y+.29342046739267774d-1*(fct(a+c)+fct(a-c))
      c=.36609105937014484d0*b
      y=y+.32911111388180923d-1*(fct(a+c)+fct(a-c))
      c=.33152213346510760d0*b
      y=y+.36172897054424253d-1*(fct(a+c)+fct(a-c))
      c=.29385787862038116d0*b
      y=y+.39096947893535153d-1*(fct(a+c)+fct(a-c))
      c=.25344995446611470d0*b
      y=y+.41655962113473378d-1*(fct(a+c)+fct(a-c))
      c=.21067563806531767d0*b
      y=y+.43826046502201906d-1*(fct(a+c)+fct(a-c))
      c=.16593430114106382d0*b
      y=y+.45586939347881942d-1*(fct(a+c)+fct(a-c))
      c=.11964368112606854d0*b
      y=y+.46922199540402283d-1*(fct(a+c)+fct(a-c))
      c=.7223598079139825d-1*b
      y=y+.47819360039637430d-1*(fct(a+c)+fct(a-c))
      c=.24153832843869158d-1*b
      y=b*(y+.48270044257363900d-1*(fct(a+c)+fct(a-c)))
      end

      subroutine simps(a1,b1,h1,reps1,aeps1,funct,x,ai,aih,aiabs)
      implicit real*8(a-h,o-z)
      dimension f(7),p(5)
      h=dsign(h1,b1-a1); s=dsign(1.d0,h)
      a=a1; b=b1; ai=0.d0; aih=0.d0; aiabs=0.d0
      p(2)=4.d0; p(4)=4.d0; p(3)=2.d0; p(5)=1.d0
      if(b-a) 1,2,1
    1 reps=dabs(reps1); aeps=dabs(aeps1)
      do 3 k=1,7
  3   f(k)=10.d16
      x=a; c=0.d0; f(1)=funct(x)/3.d0
    4 x0=x
      if((x0+4.*h-b)*s) 5,5,6
    6 h=(b-x0)/4.; if(h) 7,2,7
    7 do 8 k=2,7
  8   f(k)=10.d16
      c=1.d0
    5 di2=f(1); di3=dabs(f(1))
      do 9 k=2,5
      x=x+h
      if((x-b)*s) 23,24,24
   24 x=b
   23 if(f(k)-10.d16) 10,11,10
   11 f(k)=funct(x)/3.
   10 di2=di2+p(k)*f(k)
    9 di3=di3+p(k)*abs(f(k))
      di1=(f(1)+4.*f(3)+f(5))*2.*h
      di2=di2*h; di3=di3*h
      if(reps) 12,13,12
   13 if(aeps) 12,14,12
   12 eps=dabs((aiabs+di3)*reps)
      if(eps-aeps) 15,16,16
   15 eps=aeps
   16 delta=dabs(di2-di1)
      if(delta-eps) 20,21,21
   20 if(delta-eps/8.) 17,14,14
   17 h=2.*h; f(1)=f(5); f(2)=f(6); f(3)=f(7)
      do 19 k=4,7
  19  f(k)=10.d16
      go to 18
   14 f(1)=f(5); f(3)=f(6); f(5)=f(7)
      f(2)=10.d16; f(4)=10.d16; f(6)=10.d16; f(7)=10.d16
   18 di1=di2+(di2-di1)/15.
      ai=ai+di1; aih=aih+di2; aiabs=aiabs+di3
      go to 22
   21 h=h/2.; f(7)=f(5); f(6)=f(4); f(5)=f(3)
      f(3)=f(2); f(2)=10.d16; f(4)=10.d16
      x=x0; c=0.d0
      go to 5
   22 if(c) 2,4,2
    2 return
      end

      subroutine simpt(a1,b1,h1,reps1,aeps1,funct,x,ai,aih,aiabs)
      implicit real*8(a-h,o-z)
      dimension f(7),p(5)
      h=dsign(h1,b1-a1); s=dsign(1.d0,h)
      a=a1; b=b1; ai=0.d0; aih=0.d0; aiabs=0.d0
      p(2)=4.d0; p(4)=4.d0; p(3)=2.d0; p(5)=1.d0
      if(b-a) 1,2,1
    1 reps=dabs(reps1); aeps=dabs(aeps1)
      do 3 k=1,7
  3   f(k)=10.d16
      x=a; c=0.d0; f(1)=funct(x)/3.
    4 x0=x
      if((x0+4.*h-b)*s) 5,5,6
    6 h=(b-x0)/4.; if(h) 7,2,7
    7 do 8 k=2,7
  8   f(k)=10.d16
      c=1.d0
    5 di2=f(1); di3=dabs(f(1))
      do 9 k=2,5
      x=x+h
      if((x-b)*s) 23,24,24
   24 x=b
   23 if(f(k)-10.d16) 10,11,10
   11 f(k)=funct(x)/3.
   10 di2=di2+p(k)*f(k)
    9 di3=di3+p(k)*abs(f(k))
      di1=(f(1)+4.*f(3)+f(5))*2.*h
      di2=di2*h; di3=di3*h
      if(reps) 12,13,12
   13 if(aeps) 12,14,12
   12 eps=dabs((aiabs+di3)*reps)
      if(eps-aeps) 15,16,16
   15 eps=aeps
   16 delta=dabs(di2-di1)
      if(delta-eps) 20,21,21
   20 if(delta-eps/8.) 17,14,14
   17 h=2.*h; f(1)=f(5); f(2)=f(6); f(3)=f(7)
      do 19 k=4,7
  19  f(k)=10.d16
      go to 18
   14 f(1)=f(5); f(3)=f(6); f(5)=f(7)
      f(2)=10.d16; f(4)=10.d16; f(6)=10.d16; f(7)=10.d16
   18 di1=di2+(di2-di1)/15.
      ai=ai+di1; aih=aih+di2; aiabs=aiabs+di3
      go to 22
   21 h=h/2.; f(7)=f(5); f(6)=f(4); f(5)=f(3)
      f(3)=f(2); f(2)=10.d16; f(4)=10.d16
      x=x0; c=0.d0
      go to 5
   22 if(c) 2,4,2
    2 return
      end

      subroutine simpu(a1,b1,h1,reps1,aeps1,funct,x,ai,aih,aiabs)
      implicit real*8(a-h,o-z)
      dimension f(7),p(5)
      h=dsign(h1,b1-a1); s=dsign(1.d0,h)
      a=a1; b=b1; ai=0.d0; aih=0.d0; aiabs=0.d0
      p(2)=4.d0; p(4)=4.d0; p(3)=2.d0; p(5)=1.d0
      if(b-a) 1,2,1
    1 reps=dabs(reps1); aeps=dabs(aeps1)
      do 3 k=1,7
  3   f(k)=10.d16
      x=a; c=0.d0; f(1)=funct(x)/3.
    4 x0=x
      if((x0+4.*h-b)*s) 5,5,6
    6 h=(b-x0)/4.; if(h) 7,2,7
    7 do 8 k=2,7
  8   f(k)=10.d16
      c=1.d0
    5 di2=f(1); di3=dabs(f(1))
      do 9 k=2,5
      x=x+h
      if((x-b)*s) 23,24,24
   24 x=b
   23 if(f(k)-10.d16) 10,11,10
   11 f(k)=funct(x)/3.
   10 di2=di2+p(k)*f(k)
    9 di3=di3+p(k)*abs(f(k))
      di1=(f(1)+4.*f(3)+f(5))*2.*h
      di2=di2*h; di3=di3*h
      if(reps) 12,13,12
   13 if(aeps) 12,14,12
   12 eps=dabs((aiabs+di3)*reps)
      if(eps-aeps) 15,16,16
   15 eps=aeps
   16 delta=dabs(di2-di1)
      if(delta-eps) 20,21,21
   20 if(delta-eps/8.) 17,14,14
   17 h=2.*h; f(1)=f(5); f(2)=f(6); f(3)=f(7)
      do 19 k=4,7
  19  f(k)=10.d16
      go to 18
   14 f(1)=f(5); f(3)=f(6); f(5)=f(7)
      f(2)=10.d16; f(4)=10.d16; f(6)=10.d16; f(7)=10.d16
   18 di1=di2+(di2-di1)/15.
      ai=ai+di1; aih=aih+di2; aiabs=aiabs+di3
      go to 22
   21 h=h/2.; f(7)=f(5); f(6)=f(4); f(5)=f(3)
      f(3)=f(2); f(2)=10.d16; f(4)=10.d16
      x=x0; c=0.d0
      go to 5
   22 if(c) 2,4,2
    2 return
      end

      FUNCTION URAND(IY)
      INTEGER*4 IY,M2,IA,IC
      DATA S,M2,IA,IC/.46566128E-9,1073741824,843314861,453816693/
      IY=IY*IA+IC
      IF(IY.LT.0)IY=(IY+M2)+M2
      URAND=FLOAT(IY)*S
      END
