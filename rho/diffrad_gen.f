c==============================================================
c     diffrad_gen.f  -- version 2.0 (fully fixed)
c
c     MC Generator for diffractive vector meson electroproduction
c     with QED radiative corrections (collinear approximation)
c     Based on DIFFRAD by I.Akushevich (1998)
c
c     LUND output format (px py pz E mass):
c       1: beam e-      (initial, shifted if ISR)
c       2: target p     (initial, at rest)
c       3: scattered e- (final, shifted if FSR)
c       4: pi+          (from rho decay)
c       5: pi-          (from rho decay)
c       6: recoil p     (final)
c       7: gamma        (final, hard events only)
c==============================================================

      program diffrad_gen
      implicit real*8(a-h,o-z)
      real*4 urand

      common/cmp/pi,alpha,amp,amp2,ap,ap2,aml2,amc2,amv,barn
      common/sxy/s,x,sx,sxp,q2,w2,aly,anu,sqly,an,tamin,tamax,xs,ys
      common/phi/phirad,tdif,phidif,tq,vmax,ivec
      common/vv1vv2/aa1,aa2,bb,sib
      common/ivv/vcurr,cutv

      real*8 k1(4),k2(4),ptar(4),ph(4),pp(4),kgam(4),pip(4),pim(4)
      integer*4 iy

c     Read input
      open(unit=8, file='gen_input.dat', status='old')
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
      close(8)

      call setcon(ivec,lepton)
      s = 2d0*(sqrt(tmom**2+amp**2)*sqrt(bmom**2+aml2)+bmom*tmom)
      ebeam = bmom

      write(*,'(a)') '=============================='
      write(*,'(a,f8.3,a)') ' Ebeam  = ',ebeam,' GeV'
      write(*,'(a,f8.3,a)') ' sqrt(S)= ',sqrt(s),' GeV'
      write(*,'(a,i2)')     ' ivec   = ',ivec
      write(*,'(a,f6.3,a)') ' cutv   = ',cutv,' GeV^2'
      write(*,'(a,i8)')     ' nev    = ',nev
      write(*,'(a)') '=============================='

      open(unit=10, file='events.lund')
      open(unit=11, file='gen_stat.dat')
      open(unit=12, file='v_dist.dat')

      ngen=0; nsoft=0; nhard=0; nisr=0; nfsr=0; ntry=0
      nf_born=0; nf_thresh=0; nf_tkin=0; nf_vmax=0
      nf_sib=0; nf_sshxxh=0; nf_sxqv=0; nf_sigtot=0

      do while (ngen .lt. nev)

        ntry = ntry + 1
        if(ntry .gt. 100*nev) then
          write(*,*) 'ERROR: too many attempts'
          goto 999
        endif

c       Step 1: Sample Born kinematics into common blocks
        call sample_born(q2min,q2max,ymin,ymax,tmin,tmax,
     .                   iy,sg_born,iacc)
        if(iacc.eq.0)then; nf_born=nf_born+1; goto 100; endif

c       Step 2: Compute derived kinematics
        call conkin(s)

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

c       Step 3: RC quantities
        sxt = sx+tdif
        tq  = q2+tdif-amv**2
        aa1 = (q2*sxp*sxt-(s*sx+2d0*amp2*q2)*tq)/2d0/aly
        aa2 = (q2*sxp*sxt-(x*sx-2d0*amp2*q2)*tq)/2d0/aly
        sqbb1 = sqrt(max(0d0,q2*sxt**2-sxt*sx*tq-amp2*tq**2-amv**2*aly))
        sqbb2 = sqrt(max(0d0,q2*(s*x-amp2*q2)-aml2*aly))
        bb = sqbb1*sqbb2/aly

        vmax_kin = tt2+.5d0/q2*(-tt1*tq
     .    -sqrt(max(0d0,tt1**2+4d0*q2*w2))
     .    *sqrt(max(0d0,tq**2+4d0*amv**2*q2)))
     .    - 1d-8
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
        deltavr = (1.5d0*dlm-2d0-.5d0*log(xxh/ssh)**2
     .             +fspen(1d0-amp2*q2/ssh/xxh)-pi**2/6d0)
        delta_vac = vacpol(q2)
        extai1 = exp(alpha/pi*delinf_val)
        sig_soft = sib*extai1*(1d0+alpha/pi*(deltavr+delta_vac))

        bslope = 5d0
        qv = q2+amv**2
        if(sx-qv.le.0d0)then; nf_sxqv=nf_sxqv+1; goto 100; endif
        del_rad  = 2d0*bslope*vmax*(dlm-1d0)*qv/(sx-qv)
        sig_hard = alpha/pi*del_rad*sib
        sig_total = sig_soft+sig_hard
        if(sig_total.le.0d0)then; nf_sigtot=nf_sigtot+1; goto 100; endif

        sg_born_full = sg_born * sib
        if(sg_born_full.le.0d0)then
          nf_sigtot=nf_sigtot+1
          goto 100
        endif

        ngen = ngen + 1

c       Step 4: Soft or hard?
        prob_hard = sig_hard/sig_total
        r = dble(urand(iy))

        if(r .gt. prob_hard) then
c         NON-RADIATED EVENT
          nsoft = nsoft + 1
          call build_4vectors(ebeam,xs,ys,tdif,phirad,k1,ptar,k2,ph,pp)
          call write_lund(10,ngen,k1,ptar,k2,ph,pp,kgam,
     .                    .false.,ivec,iy,pip,pim)

        else
c         HARD RADIATED EVENT
          nhard = nhard + 1

c         Sample v log-flat
          vcut_use = 1d-4
          rvlogmin = log(vcut_use)
          rvlogmax = log(vmax)
          if(rvlogmax.le.rvlogmin)then
            nf_born=nf_born+1
            goto 100
          endif
          rlogv = rvlogmin + dble(urand(iy))*(rvlogmax-rvlogmin)
          vrad  = exp(rlogv)
          write(12,'(f12.6)') vrad

c         ISR vs FSR
          Ee1_val = ebeam
          Ee2_val = ebeam*(1d0-ys)
          p_isr   = (Ee2_val**2)/(Ee1_val**2+Ee2_val**2)

          r2 = dble(urand(iy))
          if(r2 .lt. p_isr) then
c           ISR
            nisr = nisr + 1
            omega  = vrad/(2d0*Ee1_val)
            Ee1_rc = Ee1_val - omega
            if(Ee1_rc.le.0d0)then
              nf_thresh=nf_thresh+1
              goto 100
            endif
            call build_4vectors(Ee1_rc,xs,ys,tdif,phirad,
     .                              k1,ptar,k2,ph,pp)
            kgam(1) = omega
            kgam(2) = 0d0
            kgam(3) = 0d0
            kgam(4) = omega
          else
c           FSR
            nfsr = nfsr + 1
            omega  = vrad/(2d0*Ee2_val)
            Ee2_rc = Ee2_val - omega
            if(Ee2_rc.le.0d0)then
              nf_thresh=nf_thresh+1
              goto 100
            endif
            call build_4vectors(ebeam,xs,ys,tdif,phirad,
     .                              k1,ptar,k2,ph,pp)
            ascale = Ee2_rc/Ee2_val
            k2(1) = k2(1)*ascale
            k2(2) = k2(2)*ascale
            k2(3) = k2(3)*ascale
            k2(4) = k2(4)*ascale
c           Recompute recoil proton with shifted k2
            pp(1) = k1(1)+ptar(1)-k2(1)-ph(1)
            pp(2) = k1(2)+ptar(2)-k2(2)-ph(2)
            pp(3) = k1(3)+ptar(3)-k2(3)-ph(3)
            pp(4) = k1(4)+ptar(4)-k2(4)-ph(4)
            kgam(1) = omega
            kgam(2) = (k2(2)/Ee2_rc)*omega
            kgam(3) = (k2(3)/Ee2_rc)*omega
            kgam(4) = (k2(4)/Ee2_rc)*omega
          endif

          call write_lund(10,ngen,k1,ptar,k2,ph,pp,kgam,
     .                    .true.,ivec,iy,pip,pim)
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
      write(*,'(a,i8)') ' Events generated : ',ngen
      write(*,'(a,i8)') ' Total attempts   : ',ntry
      write(*,'(a,i8)') ' Non-radiated     : ',nsoft
      write(*,'(a,i8)') ' Hard radiated    : ',nhard
      write(*,'(a,i8)') '   ISR events     : ',nisr
      write(*,'(a,i8)') '   FSR events     : ',nfsr
      write(*,'(a,f8.4)') ' Hard fraction    : ',
     .                     dble(nhard)/max(1,ngen)
      write(*,'(a)') '=================================='

      write(11,'(a,i8)') 'ngen      = ',ngen
      write(11,'(a,i8)') 'nsoft     = ',nsoft
      write(11,'(a,i8)') 'nhard     = ',nhard
      write(11,'(a,f8.4)') 'hard_frac = ',dble(nhard)/max(1,ngen)

      close(10); close(11); close(12)
      end


c==============================================================
c     sample_born
c==============================================================
      subroutine sample_born(q2min,q2max,ymin,ymax,tmin,tmax,
     .                       iy,sg_born,iacc)
      implicit real*8(a-h,o-z)
      real*4 urand
      common/cmp/pi,alpha,amp,amp2,ap,ap2,aml2,amc2,amv,barn
      common/sxy/s,x,sx,sxp,q2,w2,aly,anu,sqly,an,tamin,tamax,xs,ys
      common/phi/phirad,tdif,phidif,tq,vmax,ivec
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

      abslope = 5d0
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


c==============================================================
c     build_4vectors
c     NOTE: all mass variables use 'am' prefix to avoid implicit integer
c     (m,M,k,j,l,n,i are integer by default in Fortran)
c==============================================================
      subroutine build_4vectors(ebeam,xs,ys,tdif,phirad,
     .                          k1,ptar,k2,ph,pp)
      implicit real*8(a-h,o-z)
      common/cmp/pi,alpha,amp,amp2,ap,ap2,aml2,amc2,amv,barn
      real*8 k1(4),ptar(4),k2(4),ph(4),pp(4)
      real*8 qvec(3),qhat(3),e1v(3),e2v(3),ph3(3),tmp(3)

c     Beam electron
      k1(1) = ebeam
      k1(2) = 0d0
      k1(3) = 0d0
      k1(4) = ebeam

c     Target proton at rest
      ptar(1) = amp
      ptar(2) = 0d0
      ptar(3) = 0d0
      ptar(4) = 0d0

c     Scattered electron (ultrarelativistic)
      aEe2 = ebeam*(1d0-ys)
      if(aEe2.le.0d0) return
      acosthe = 1d0 - (ebeam*ebeam*ys*ys*xs*2d0*amp/
     .          ((ebeam-ebeam*ys)*ebeam*2d0))
c     Simpler: Q2 = 2*E1*E2*(1-cos)
      acosthe = 1d0 - (s*xs*ys)/(2d0*ebeam*aEe2)
      if(abs(acosthe).gt.1d0) return
      asinthe = sqrt(max(0d0,1d0-acosthe**2))

      k2(1) = aEe2
      k2(2) = aEe2*asinthe
      k2(3) = 0d0
      k2(4) = aEe2*acosthe

c     Virtual photon
      aEq     = k1(1)-k2(1)
      qvec(1) = k1(2)-k2(2)
      qvec(2) = k1(3)-k2(3)
      qvec(3) = k1(4)-k2(4)
      aqmag   = sqrt(qvec(1)**2+qvec(2)**2+qvec(3)**2)

c     Rho meson 4-vector
c     Use 'amrho2' not 'Mrho2' -- M->m is IMPLICIT INTEGER!
      amrho2 = amv**2
      anu_loc = ebeam - aEe2
      aq2_loc = s*xs*ys

c     q.ph = (-Q2 + Mrho2 - t) / 2
      aqdotph = (-aq2_loc + amrho2 - tdif)/2d0

c     Erho from 4-momentum conservation
      aErho = (2d0*amp*anu_loc - aq2_loc + amrho2 - 2d0*aqdotph)
     .        /(2d0*amp)
      if(aErho.lt.amv) return

      aphrho = sqrt(max(0d0,aErho**2-amrho2))

c     Angle of rho wrt virtual photon
      if(aqmag*aphrho.le.0d0) return
      acosalpha = (aEq*aErho - aqdotph)/(aqmag*aphrho)
      if(abs(acosalpha).gt.1d0) acosalpha=sign(1d0,acosalpha)
      asinalpha = sqrt(max(0d0,1d0-acosalpha**2))

c     Orthonormal basis around q
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

c     Rho 3-momentum direction
      ph3(1) = aphrho*(acosalpha*qhat(1)
     .        + asinalpha*cos(phirad)*e1v(1)
     .        + asinalpha*sin(phirad)*e2v(1))
      ph3(2) = aphrho*(acosalpha*qhat(2)
     .        + asinalpha*cos(phirad)*e1v(2)
     .        + asinalpha*sin(phirad)*e2v(2))
      ph3(3) = aphrho*(acosalpha*qhat(3)
     .        + asinalpha*cos(phirad)*e1v(3)
     .        + asinalpha*sin(phirad)*e2v(3))

c     Put rho on mass shell: E = sqrt(|p|^2 + Mrho^2)
      ph(2) = ph3(1)
      ph(3) = ph3(2)
      ph(4) = ph3(3)
      ph(1) = sqrt(ph(2)**2+ph(3)**2+ph(4)**2+amrho2)

c     Recoil proton from 4-momentum conservation
      pp(1) = k1(1)+ptar(1)-k2(1)-ph(1)
      pp(2) = k1(2)+ptar(2)-k2(2)-ph(2)
      pp(3) = k1(3)+ptar(3)-k2(3)-ph(3)
      pp(4) = k1(4)+ptar(4)-k2(4)-ph(4)

      return
      end


c==============================================================
c     cross3
c==============================================================
      subroutine cross3(a,b,c)
      implicit real*8(a-h,o-z)
      dimension a(3),b(3),c(3)
      c(1) = a(2)*b(3)-a(3)*b(2)
      c(2) = a(3)*b(1)-a(1)*b(3)
      c(3) = a(1)*b(2)-a(2)*b(1)
      return
      end


c==============================================================
c     decay_rho: rho -> pi+ pi- isotropic in rho rest frame
c     KEY FIX: use 'ampi','amrho' NOT 'mpi','Mrho' (m is implicit integer!)
c==============================================================
      subroutine decay_rho(ph,pip,pim,iy)
      implicit real*8(a-h,o-z)
      real*4 urand
      common/cmp/pi,alpha,amp,amp2,ap,ap2,aml2,amc2,amv,barn
      real*8 ph(4),pip(4),pim(4)
      real*8 abeta(3),pip_rf(4),pim_rf(4)
      integer*4 iy

c     Use 'am' prefix -- avoids implicit integer for m,M letters
      ampi  = 0.13957d0
      amrho = amv

c     Pion momentum in rho rest frame
      appion = sqrt(max(0d0,(amrho/2d0)**2-ampi**2))

c     Random decay direction
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

c     Boost to lab (rho velocity)
      aErho   = ph(1)
      abeta(1) = ph(2)/aErho
      abeta(2) = ph(3)/aErho
      abeta(3) = ph(4)/aErho
      abeta2   = abeta(1)**2+abeta(2)**2+abeta(3)**2

c     Use actual rho mass (ph is on mass shell from build_4vectors)
      agamma = aErho/amrho

      call lorentz_boost(pip_rf,abeta,agamma,abeta2,pip)
      call lorentz_boost(pim_rf,abeta,agamma,abeta2,pim)

      return
      end


c==============================================================
c     lorentz_boost: p_out = boost(p_in) by velocity abeta
c     Standard formula: E'=g(E+b.p), p'=p+(g-1)/b2*(b.p)*b+g*E*b
c==============================================================
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


c==============================================================
c     write_lund
c==============================================================
      subroutine write_lund(lun,iev,k1,ptar,k2,ph,pp,kgam,
     .                      has_photon,ivec,iy,pip,pim)
      implicit real*8(a-h,o-z)
      common/cmp/pi,alpha,amp,amp2,ap,ap2,aml2,amc2,amv,barn
      real*8 k1(4),ptar(4),k2(4),ph(4),pp(4),kgam(4),pip(4),pim(4)
      logical has_photon
      integer*4 iy

      call decay_rho(ph,pip,pim,iy)

      npart = 6
      if(has_photon) npart = 7

      write(lun,'(i5,2(5x,i1),6(5x,f8.4),5x,i8)')
     .      npart,1,1,0d0,0d0,0d0,0d0,0d0,0d0,iev

      write(lun,'(i4,i6,5f12.6)') 1, 11,
     .      k1(2),k1(3),k1(4),k1(1),0.000511d0
      write(lun,'(i4,i6,5f12.6)') 2, 2212,
     .      ptar(2),ptar(3),ptar(4),ptar(1),amp
      write(lun,'(i4,i6,5f12.6)') 3, 11,
     .      k2(2),k2(3),k2(4),k2(1),0.000511d0
      write(lun,'(i4,i6,5f12.6)') 4, 211,
     .      pip(2),pip(3),pip(4),pip(1),0.13957d0
      write(lun,'(i4,i6,5f12.6)') 5,-211,
     .      pim(2),pim(3),pim(4),pim(1),0.13957d0
      write(lun,'(i4,i6,5f12.6)') 6, 2212,
     .      pp(2),pp(3),pp(4),pp(1),amp
      if(has_photon)then
        write(lun,'(i4,i6,5f12.6)') 7, 22,
     .        kgam(2),kgam(3),kgam(4),kgam(1),0d0
      endif

      return
      end


c==============================================================
c     conkin
c==============================================================
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


c==============================================================
c     Akushevich physics routines (unchanged)
c==============================================================
      subroutine setcon(ivec,lepton)
      implicit real*8(a-h,o-z)
      common/cmp/pi,alpha,amp,amp2,ap,ap2,aml2,amc2,amv,barn
      dimension amhad(4)
      data amhad/0.7683d0,0.78195d0,1.019412d0,3.0969d0/
      if(lepton.eq.1)aml2=.261112d-6
      if(lepton.eq.2)aml2=.111637d-1
      pi=3.1415926d0; alpha=.729735d-2; barn=.389379d6
      amv=amhad(ivec); amp=.938272d0; amp2=amp**2
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
      sibor=2d0*an/(xs*ys**2)*(ys**2*sigmat
     .      +2d0*(1d0-ys-.25d0*ys**2*aga2)*(sigmal+sigmat))
      end

      subroutine difflt(q2,w2,t,sigl,sigt)
      implicit real*8(a-h,o-z)
      common/cmp/pi,alpha,amp,amp2,ap,ap2,aml2,amc2,amv,barn
      common/phi/phirad,tdif,phidif,tq,vmax,ivec
      common/pri/ipri
      dimension ggam(4)
      data p02/.5d0/al_s/.25d0/
      data ggam/6.77d-6,0.6d-6,1.37d-6,5.36d-6/
      asx=w2+q2-amp2; aanu=asx/ap; asxt=asx+t; aeh=asxt/2d0/amp
      amv2=amv**2; aeta=1d0; aff2=exp(5d0*t)
      atqt=t+q2-amv2
      apt2=(-(4d0*(aanu**2+q2)*amv2+4d0*aanu*aeh*atqt
     .     -4d0*aeh**2*q2+atqt**2))/(4d0*(aanu**2+q2))
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
  10  asuml=asuml+2d0*(t+aa2)*aallmi/3d0-10d0/9d0
     .            +4d0*aa2*(1d0-aa2*aallmi)/3d0/t
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
 10   fspen=-.5d0*log(1d0-x)*log(x**2/(1d0-x))
     .      -af1+fspens(1d0/(1d0-x))
      return
      end

      FUNCTION URAND(IY)
      INTEGER*4 IY,M2,IA,IC
      DATA S,M2,IA,IC/.46566128E-9,1073741824,843314861,453816693/
      IY=IY*IA+IC
      IF(IY.LT.0)IY=(IY+M2)+M2
      URAND=FLOAT(IY)*S
      END
