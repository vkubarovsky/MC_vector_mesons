!
!     Works for all 1s, 2s and d waves
!     Gluon density tabulated
!
!    Classic DIVE 2003, checked with various data from Review.
!
      program DIVE
      implicit real*8(a-h,k-m,o-z)
      external diveinit,dsidtf
      common/maincom/alpha,pi,anc,amvm(5,4),mq,cv,gamm(4,3),skew
      common/keycom/ima,nout,nvm,nst,nwf
      common/poicom/xeff,q2,w2,delta,ep
      parameter(npoimax=100)
      dimension amqm(5)
      dimension arx(npoimax),aq2(npoimax),ade(npoimax)
      data amqm/.22d0,.37d0,1.50d0,5.00d0,5.00/


      open(2,file='input.dat',status='old')
      open(3,file='data2/file.dat')
      open(4,file='data2/file1.dat')


      read(2,*)nout  ! -- output  
      read(2,*)nvm   ! -- vector meson
      read(2,*)nst   ! -- vector meson state
      read(2,*)nwf   ! -- wave function ansatz   

c      iglu = 1 ! -- gluon density
      alpha=.729735d-2
      pi=atan(1d0)*4d0
      anc=3.d0             ! -- N_c
      ss=4d0*30d0*920d0 ! -- HERA only (NEW)
c      ss=4d0*27.57d0*810d0 ! -- HERA only
c      ss=4d0*1000d0*25000d0 ! -- sqrt(s) = 10 TeV
      if(nvm.eq.1)cv=1d0/sqrt(2d0)
      if(nvm.eq.2)cv=1d0/3d0
      if(nvm.eq.3)cv=2d0/3d0
      if(nvm.eq.4)cv=1d0/3d0
      if(nvm.eq.5)cv=1d0/3d0
      mq = amqm(nvm)
      skew = 0.41d0

c     outputting input
      write(*,*)'nvm =',nvm,',nst =',nst,'nwf =',nwf
c      write(10,*)'nvm =',nvm,',nst =',nst,'nwf =',nwf

      call wfinit ! -- initializing wave functions
      IF (nout.EQ.0) STOP
c      wave function output is done inside the subroutine

c     Reading kinematic data
      read(2,*)npoi
      read(2,*)(arx(i),i=1,npoi) ! -- reading W/x data
      read(2,*)(aq2(i),i=1,npoi) ! -- reading Q^2 data
      read(2,*)(ade(i),i=1,npoi) ! -- reading |t| data

c    *** Starting the main kinematic loop ***
      DO ipoi=1,npoi

	 q2=aq2(ipoi)
	 if(arx(ipoi).gt.0d0)then
	   xeff=arx(ipoi)
	   w2= -q2+(q2+amvm(nvm,nst)**2)/(xeff/skew)
	 else
	   w2=arx(ipoi)**2
         xeff=skew*(q2+amvm(nvm,nst)**2)/(q2+w2)
	 endif
	 y=(w2+q2-0.938272**2)/ss
	 ep=(1d0-y)/(1d0-y+y**2/2d0)

c    *** Total cross section ***
         IF (nout.EQ.1) then 
            A=0.35d0
            B=1d0
            N=3
        ! Integration is done with respect to xi = exp(-t)
        ! 2*N+1 = 7 points are enough for 5% accuracy 
            step=(b-a)/2/n
            r_l=0
            r_t=0
            do i=0,2*n
               if(mod(i,2).eq.0)then
                  ki=2
               else
                  ki=4
               endif
               if (i.EQ.0.or.i.EQ.2*n) ki = 1
               xi = a + step*i
               t = -log(xi)
               delta = sqrt(t)
               call dsidtf(dsig_l,dsig_t)
               value_l = dsig_l
               value_t = dsig_t
               r_l = r_l + ki*value_l/xi
               r_t = r_t + ki*value_t/xi
       write (*,*)sngl(t),sngl(dsig_l),sngl(dsig_t)
            enddo
            sigma_l = r_l*step/3d0
            sigma_t = r_t*step/3d0
            sigma_tot = ep*sigma_l+sigma_t 
            b_l = value_l/sigma_l
            b_t = value_t/sigma_t
            b_tot = (ep*value_l+value_t)/sigma_tot
       write (*,*)sngl(q2), sngl(dsqrt(W2)),sngl(sigma_l),sngl(sigma_t)
       write (*,*)sngl(b_l),sngl(b_t),sngl(b_tot)
         ENDIF   
         
c    *** Differential cross section or distribution *** 
         IF (nout.GT.1) then 
	    delta=sqrt(ade(ipoi))
            call dsidtf(dsig_l,dsig_t)
            dsigdt = ep*dsig_l+dsig_t
c            qb2l = 0.3*q2+0.33
c            modif = dsig_l/q2*(qb2l)**4
            write (*,*)sngl(q2),sngl(dsig_l),sngl(dsig_t),sngl(dsigdt)
c            write (3,*)sngl(q2),sngl(dsig_l),sngl(dsig_t),sngl(dsigdt)
         ENDIF   
         

c         pause
      enddo ! -- kinematic loop    

      close(3)
      close(4)
	write(*,*)'End'
c	PAUSE
      END ! - program  


! ============  END OF THE MAIN PROGRAM  ======================== !

      Block data diveinit
      implicit real*8(a-h,k-m,o-z)
      common/maincom/alpha,pi,anc,amvm(5,4),mq,cv,gamm(4,3),skew
      common/keycom/ima,nout,nvm,nst,nwf
      common/poicom/xeff,q2,w2,delta,ep
      data amvm/0.77d0,1.019d0,3.0969d0,9.46d0, 10d0    ! 1s
     .	       ,1.45d0,1.68d0 ,3.6860d0,10.02d0, 10d0    ! 2s
c     .         ,0.77d0,1.019d0,3.0969d0,9.46d0, 10d0 ! 2s false
     .	       ,1.70d0,1.68d0 ,3.77d0  ,10.02d0, 10d0    !  d
     .	       ,0.77d0,1.019d0 ,3.0969d0  ,9.46d0, 10d0/   !  mix as 1s
      data gamm/6.77d-6,1.37d-6,5.36d-6,1.34d-6    ! 1s
     .	 ,1.90d-6,0.53d-6,2.14d-6,0.52d-6    ! 2s
     .	 ,0.14d-6,0.03d-6,0.26d-6,0.01d-6/   !  d

      end

! ================== WAVE FUNCTIONS =================== !

      subroutine wfinit
      implicit real*8(a-h,k-m,o-z)
      external fapar1s,fapar2s,fapard,psidf,psif,psin,psiort,fmix1s
      common/maincom/alpha,pi,anc,amvm(5,4),mq,cv,gamm(4,3),skew
      common/keycom/ima,nout,nvm,nst,nwf
      common/poicom/xeff,q2,w2,delta,ep
      common/psicom/a1psi,a2psi,adpsi,cc1,cc2,ccd,anode
      common/mix/phimix
c      data rvm/2.0d0/ ! -- GeV**-2


         nsttrue = nst 
 
 	 nst = 1

         aaamin=1d-1
	 aaamax=100d0
         write(*,*)'starting 1S...'
	 call mhord(2,aaamin,aaamax,1d-6,fapar1s,apar1s,ikey)
	 if(ikey.ne.0)then ! appears if something was wrong
	   nhor=10
	   do ihor=1,nhor
	     aaa=aaamin+(ihor-1)*(aaamax-aaamin)/(nhor-1)
	     ffun=fapar1s(aaa)
             write(*,'(2g12.4)')aaa,ffun
	   enddo
	   stop 'ikey1s'
	 endif

         call simpsx(0d0,1d2,10000,1d-5,psin,ani)
         an=anc/(2d0*pi)**3*ani
	 write(*,*)' Normalization 1s =',an
c	 write(10,*)' Normalization 1s =',an

	 call simpsx(0d0,1d4,10000,1d-5,psif,fsi)
	 fs1ca=anc/(2d0*pi)**3*fsi
         gam1s=4d0*pi*alpha**2*cv**2/(3d0*amvm(nvm,1)**3)*fs1ca**2
	 write(*,*)' 1S decay width =',gam1s
c	 write(10,*)' 1S decay constant =',fs1ca
         
         nst = 2

	 aaamin=1d-1
	 aaamax=100d0
         write(*,*)'starting 2S...'
	 call mhord(2,aaamin,aaamax,1d-5,fapar2s,apar2s,ikey)
	 if(ikey.ne.0)then ! appears if something was wrong
	   nhor=10
	   do ihor=1,nhor
	     aaa=aaamin+(ihor-1)*(aaamax-aaamin)/(nhor-1)
	     ffun=fapar2s(aaa)
             write(*,'(2g12.4)')aaa,ffun
	   enddo
	   stop 'ikey2s'
	 endif
	
         call simpsx(0d0,1d2,10000,1d-5,psin,ani)
         an=anc/(2d0*pi)**3*ani
	 write(*,*)' Normalization 2s =',an
c	 write(10,*)' Normalization 2s =',an

         call simpsx(0d0,1d2,10000,1d-5,psiort,ani)
         an=anc/(2d0*pi)**3*ani
	 write(*,*)' Ortogonality check  =',an
c	 write(10,*)' Normalization 2s =',an

	 call simpsx(0d0,1d4,10000,1d-5,psif,fsi)
	 fs2ca=anc/(2d0*pi)**3*fsi
         gam2s=4d0*pi*alpha**2*cv**2/(3d0*amvm(nvm,2)**3)*fs2ca**2
	 write(*,*)' 2S decay width, anode =',gam2s,sngl(anode)
c	 write(10,*)' 2S decay constant =',fs2ca
         
 	 nst = 3

         aaamin=1d-1
	 aaamax=100d0
         write(*,*)'starting D wave ...'
	 call mhord(2,aaamin,aaamax,1d-6,fapard,apard,ikey)
	 if(ikey.ne.0)then ! appears if something was wrong
	   nhor=10
	   do ihor=1,nhor
	     aaa=aaamin+(ihor-1)*(aaamax-aaamin)/(nhor-1)
	     ffun=fapard(aaa)
             write(*,'(2g12.4)')aaa,ffun
	   enddo
	   stop 'ikey1s'
 	 endif

         call simpsx(0d0,1d4,10000,1d-5,psin,ani)
         an=anc/(2d0*pi)**3*ani
	 write(*,*)' Normalization D =',an
c	 write(10,*)' Normalization D =',an

	 call simpsx(0d0,1d5,100000,1d-5,psidf,fdi)
	 fdca=anc/(2d0*pi)**3*fdi
         gamd=4d0*pi*alpha**2*cv**2/(3d0*amvm(nvm,3)**3)*fdca**2
	 write(*,*)' D wave decay width =',gamd
c	 write(10,*)' D wave decay constant =',fdca
         


c	 if(nwf.eq.1)then
c	   ccc=sqrt(80d0/9d0)
c	   fdfs=ccc*sqrt(a1psi**3/adpsi**3)/adpsi**2/amvm(nvm,3)**2
c	 elseif(nwf.eq.2)then
c	   ccc=sqrt(160d0/3d0)
c	   fdfs=ccc*sqrt(a1psi**3/adpsi**3)/adpsi**2/amvm(nvm,3)**2
c	 endif
c           fs1camix=fs1ca*cos(phimix)+fdca*sin(phimix)
c           fdcamix=-fs1ca*sin(phimix)+fdca*cos(phimix)
c           fs1ca=fs1camix
c           fdca=fdcamix
c         gam1s=4d0*pi*alpha**2*cv**2/(3d0*amvm(nvm,1)**3)*fs1ca**2
c         gamd =4d0*pi*alpha**2*cv**2/(3d0*amvm(nvm,3)**3)*fdca**2

         nst = nsttrue
        write(*,'(i4,6g12.4)')nvm,mq,fs1ca,fs2ca,fdca
	 write(*,'(16x,6g12.4)')      gam1s,gam2s,gamd
	 write(*,'(a10,5x,3f12.4)')' a1,a2,ad:', a1psi,a2psi,adpsi
c        write(10,'(i4,6g12.4)')nvm,mq,fs1ca,fs2ca,fdca
c	 write(10,'(16x,6g12.4)')      gam1s,gam2s,gamd
c	 write(10,'(a10,5x,3f12.4)')' a1,a2,ad:', a1psi,a2psi,adpsi
c	 write(*,'(a5,10x,3f12.4)')' m_v ',(amvm(ivm,iii),iii=1,3)
c        write(*,'(2g12.4,f7.3)')fdca/fs1ca,fdfs,fdca/fs1ca/fdfs


      end


      real*8 function fapar1s(a)
      implicit real*8(a-h,k-m,o-z)
      external psif
      common/psicom/a1,a2,ad,cc1,cc2,ccd,anode
      common/maincom/alpha,pi,anc,amvm(5,4),mq,cv,gamm(4,3),skew
      common/keycom/ima,nout,nvm,nst,nwf
      common/poicom/xeff,q2,w2,delta,ep

      a1=a
            AA=0.0
            BB=20d0
            N=1000
            step=(bb-aa)/2/n
            r=0
            do i=0,2*n
               if(mod(i,2).eq.0)then
                  ki=2
               else
                  ki=4
               endif
               if (i.EQ.0.or.i.EQ.2*n) ki = 1
               p = aa + step*i
               amm=2d0*sqrt(mq**2+p**2)
               if(nwf.eq.1)then ! --- suppressed coulomb 
                 radial=1d0/sqrt(amm)/(1d0+a1**2*p**2)**2
               elseif(nwf.eq.2)then ! --- suppressed oscillator
	         radial=1d0/sqrt(amm)*exp(-p**2*a1**2/2d0)
               elseif(nwf.eq.3)then ! --- pure oscillator
	         radial=exp(-p**2*a1**2/2d0)
               endif
               value = p*p*DSQRT(mq**2+p**2)*radial**2
               r = r + ki*value
            enddo
            resint = r*step/3d0
            cc1 = DSQRT(pi*pi/4d0/anc/resint)
            
      call simpsx(0d0,1.d4,10000,1d-5,psif,fsi)
      fs=sqrt((3d0*amvm(nvm,1)**3*gamm(nvm,1))/(4d0*pi*alpha**2*cv**2))
      fsca=anc/(2d0*pi)**3*fsi
      fapar1s=fsca/fs-1d0
c      write(*,'(6f11.4)')a,fs,fsca,fapar1s
c      write(10,'(6f11.4)')a,fs,fsca,fapar1s
      end



      real*8 function fapar2s(a)
      implicit real*8(a-h,k-m,o-z)
      external psif
      common/psicom/a1,a2,ad,cc1,cc2,ccd,anode
      common/maincom/alpha,pi,anc,amvm(5,4),mq,cv,gamm(4,3),skew
      common/keycom/ima,nout,nvm,nst,nwf
      common/poicom/xeff,q2,w2,delta,ep

            a2=a
      ! -- ortogonality integral
            AA=0.0
            BB=20d0
            N=1000
            step=(bb-aa)/2/n
            r1=0
            r2=0

      if(nwf.eq.1)then ! --- suppressed coulomb 
            do i=0,2*n
               if(mod(i,2).eq.0)then
                  ki=2
               else
                  ki=4
               endif
               if (i.EQ.0.or.i.EQ.2*n) ki = 1
               p = aa + step*i
               amm=2d0*sqrt(mq**2+p**2)
               f1=1d0/amm/(1d0+a1**2*p**2)**2
     .     *p*p*a2*a2/(1d0+a2**2*p**2)**3
               f2=1d0/amm/(1d0+a1**2*p**2)**2
     .     /(1d0+a2**2*p**2)**3
               value1 = p*p*DSQRT(mq**2+p**2)*f1
               value2 = p*p*DSQRT(mq**2+p**2)*f2
               r1 = r1 + ki*value1
               r2 = r2 + ki*value2
            enddo

       elseif(nwf.eq.2)then ! --- suppressed oscillator
            do i=0,2*n
               if(mod(i,2).eq.0)then
                  ki=2
               else
                  ki=4
               endif
               if (i.EQ.0.or.i.EQ.2*n) ki = 1
               p = aa + step*i
               amm=2d0*sqrt(mq**2+p**2)

	       f1=1d0/amm*exp(-p**2*a1**2/2d0)
     .     *p*p*a2*a2*exp(-p**2*a2**2/2d0)
	       f2=1d0/amm*exp(-p**2*a1**2/2d0)
     .     *exp(-p**2*a2**2/2d0)

               value1 = p*p*DSQRT(mq**2+p**2)*f1
               value2 = p*p*DSQRT(mq**2+p**2)*f2
               r1 = r1 + ki*value1
               r2 = r2 + ki*value2
            enddo

       elseif(nwf.eq.3)then ! --- pure oscillator
            do i=0,2*n
               if(mod(i,2).eq.0)then
                  ki=2
               else
                  ki=4
               endif
               if (i.EQ.0.or.i.EQ.2*n) ki = 1
               p = aa + step*i
               amm=2d0*sqrt(mq**2+p**2)

               f1=exp(-p**2*a1**2/2d0)*exp(-p**2*a2**2/2d0)*p*p*a2*a2
	       f2=exp(-p**2*a1**2/2d0)*exp(-p**2*a2**2/2d0)

               value1 = p*p*DSQRT(mq**2+p**2)*f1
               value2 = p*p*DSQRT(mq**2+p**2)*f2
               r1 = r1 + ki*value1
               r2 = r2 + ki*value2
            enddo

       endif

            resint1 = r1*step/3d0
            resint2 = r2*step/3d0
            anode = resint1/resint2      

            AA=0.0
            BB=20d0
            N=1000
            step=(bb-aa)/2/n
            r=0
            do i=0,2*n
               if(mod(i,2).eq.0)then
                  ki=2
               else
                  ki=4
               endif
               if (i.EQ.0.or.i.EQ.2*n) ki = 1
               p = aa + step*i
               amm=2d0*sqrt(mq**2+p**2)
       if(nwf.eq.1)then ! --- suppressed coulomb
         radial=1d0*(p*p*a2*a2-anode)/sqrt(amm)/(1d0+a2**2*p**2)**3
       elseif(nwf.eq.2)then ! --- suppressed oscillator
	 radial=1d0/sqrt(amm)*exp(-p**2*a2**2/2d0)*(anode-p*p*a2*a2)
       elseif(nwf.eq.3)then ! --- pure oscillator
	 radial=1d0*exp(-p**2*a2**2/2d0)*(anode-p*p*a2*a2)
       endif
               value = p*p*DSQRT(mq**2+p**2)*radial**2
               r = r + ki*value
            enddo
            resint = r*step/3d0
            cc2 = DSQRT(pi*pi/4d0/anc/resint)
            
      call simpsx(0d0,1.d4,100000,1d-5,psif,fsi)
      fs=sqrt((3d0*amvm(nvm,2)**3*gamm(nvm,2))/(4d0*pi*alpha**2*cv**2))
      fsca=anc/(2d0*pi)**3*fsi
      fapar2s=abs(fsca/fs)-1d0
c      write(*,*)'My: =',sngl(a2),sngl(anode),sngl(fs),sngl(fsi)
c      write(*,*)'My: =',sngl(a),sngl(anode),sngl(fs),sngl(fsca)
c      write(*,'(6f11.4)')a,fs,fsca,fapar2s,anode

      end


      real*8 function fapard(a)
      implicit real*8(a-h,k-m,o-z)
      external psidf
      common/psicom/a1,a2,ad,cc1,cc2,ccd,anode
      common/maincom/alpha,pi,anc,amvm(5,4),mq,cv,gamm(4,3),skew
      common/keycom/ima,nout,nvm,nst,nwf
      common/poicom/xeff,q2,w2,delta,ep

      ad=a	  
            AA=0.0
            BB=20d0
            N=1000
            step=(bb-aa)/2/n
            r=0
            do i=0,2*n
               if(mod(i,2).eq.0)then
                  ki=2
               else
                  ki=4
               endif
               if (i.EQ.0.or.i.EQ.2*n) ki = 1
               p = aa + step*i
               amm=2d0*sqrt(mq**2+p**2)
               if(nwf.eq.1)then ! --- suppressed coulomb 
                 radial=1d0/sqrt(amm)/(1d0+ad**2*p**2)**4
               elseif(nwf.eq.2)then ! --- suppressed oscillator
	         radial=1d0/sqrt(amm)*exp(-p**2*ad**2/2d0)
               elseif(nwf.eq.3)then ! --- pure oscillator
	         radial=exp(-p**2*ad**2/2d0)
               endif
               value = p**6*DSQRT(mq**2+p**2)*radial**2
               r = r + ki*value
            enddo
            resint = r*step/3d0
            ccd = DSQRT(pi*pi/8d0/anc/resint)
            
      call simpsx(0d0,1.d4,10000,1d-5,psidf,fdi)
      fd=sqrt((3d0*amvm(nvm,3)**3*gamm(nvm,3))/(4d0*pi*alpha**2*cv**2))
      fdca=anc/(2d0*pi)**3*fdi
      fapard=fdca/fd-1d0
c      write(*,'(6f11.4)')a,fs,fsca,fapar1s
c      write(10,'(6f11.4)')a,fs,fsca,fapar1s
      end


      real*8 function phi1s(p)
      implicit real*8(a-h,k-m,o-z)
	common/maincom/alpha,pi,anc,amvm(5,4),mq,cv,gamm(4,3),skew
	common/keycom/ima,nout,nvm,nst,nwf
	common/poicom/xeff,q2,w2,delta,ep
	common/psicom/a1,a2,ad,cc1,cc2,ccd,anode
	
      
	amm=2d0*sqrt(mq**2+p**2)
	if(nwf.eq.1)then ! --- suppressed coulomb
	phi1s=cc1/sqrt(amm)/(1d0+a1**2*p**2)**2
	elseif(nwf.eq.2)then ! --- suppressed oscillator
	phi1s=cc1/sqrt(amm)*exp(-p**2*a1**2/2d0)
	elseif(nwf.eq.3)then ! --- pure oscillator
	phi1s=cc1*exp(-p**2*a1**2/2d0)
	endif
      end

      real*8 function phi2s(p)
      implicit real*8(a-h,k-m,o-z)
      common/maincom/alpha,pi,anc,amvm(5,4),mq,cv,gamm(4,3),skew
      common/keycom/ima,nout,nvm,nst,nwf
      common/poicom/xeff,q2,w2,delta,ep
      common/psicom/a1,a2,ad,cc1,cc2,ccd,anode

      amm=2d0*sqrt(mq**2+p**2)
      if(nwf.eq.1)then ! --- suppressed coulomb
	phi2s=cc2*(p*p*a2*a2-anode)/sqrt(amm)/(1d0+a2**2*p**2)**3
      elseif(nwf.eq.2)then ! --- suppressed oscillator
	phi2s=cc2/sqrt(amm)*exp(-p**2*a2**2/2d0)
     .     *(anode-p*p*a2*a2)
      elseif(nwf.eq.3)then ! --- pure oscillator
	phi2s=cc2*exp(-p**2*a2**2/2d0)*(anode-p*p*a2*a2)
      endif

c      write(*,'(4g11.4)')p,cc2,a2,phi2s

      end


      real*8 function phid(p)
      implicit real*8(a-h,k-m,o-z)
      common/maincom/alpha,pi,anc,amvm(5,4),mq,cv,gamm(4,3),skew
      common/keycom/ima,nout,nvm,nst,nwf
      common/poicom/xeff,q2,w2,delta,ep
      common/psicom/a1,a2,ad,cc1,cc2,ccd,anode

      amm=2d0*sqrt(mq**2+p**2)
      if(nwf.eq.1)then ! --- suppressed coulomb
	phid=ccd/sqrt(amm)/(1d0+ad**2*p**2)**4
      elseif(nwf.eq.2)then ! --- suppressed oscillator
	phid=ccd/sqrt(amm)*exp(-p**2*ad**2/2d0)
      elseif(nwf.eq.3)then ! --- pure oscillator
	phid=ccd*exp(-p**2*ad**2/2d0)
      endif

c      write(*,'(4g11.4)')p,cc1,a1,phi1s

      end



      real*8 function phi(p)
      implicit real*8(a-h,k-m,o-z)
      external phi1s,phi2s,phid
      common/maincom/alpha,pi,anc,amvm(5,4),mq,cv,gamm(4,3),skew
      common/keycom/ima,nout,nvm,nst,nwf
      common/poicom/xeff,q2,w2,delta,ep
       if(nst.EQ.1)then
	phi=phi1s(p)
       elseif(nst.eq.2)then
	phi=phi2s(p)
       elseif(nst.eq.3)then
	phi=phid(p)
       endif
      end

      real*8 function psin(p)  ! - only for normalization check
      implicit real*8(a-h,l,m,o-z)
      external phi
      common/maincom/alpha,pi,anc,amvm(5,4),mq,cv,gamm(4,3),skew
      common/keycom/ima,nout,nvm,nst,nwf
      common/poicom/xeff,q2,w2,delta,ep

      amm=2d0*sqrt(mq**2+p**2)
      psiics=phi(p)
      if(nst.eq.3)then
         psin=4d0*pi*p**6 * 8d0*amm * psiics**2    ! d
      else
         psin=4d0*pi*p**2 * 4d0*amm * psiics**2    ! 1s,2s
      endif
      end


      real*8 function psiort(p)  ! - ortogonality check
      implicit real*8(a-h,l,m,o-z)
      external phi1s,phi2s
      common/maincom/alpha,pi,anc,amvm(5,4),mq,cv,gamm(4,3),skew
      common/keycom/ima,nout,nvm,nst,nwf
      common/poicom/xeff,q2,w2,delta,ep

      amm=2d0*sqrt(mq**2+p**2)
      pss=phi1s(p)*phi2s(p)
      psiort=4d0*pi*p**2 * 4d0*amm * pss     ! 1s,2s
      end



      real*8 function psif(p)
      implicit real*8(a-h,k-m,o-z)
      external phi
      common/maincom/ alpha,pi,anc,amvm(5,4),mq,cv,gamm(4,3)

      amm=2d0*sqrt(mq**2+p**2)
      psiics=phi(p)
      psif=4d0*pi*p**2 * 8d0/3d0*(amm+mq) * psiics
c     print *,p,psif1
      end

      real*8 function psidf(p)
      implicit real*8(a-h,k-m,o-z)
      external phi
      common/maincom/ alpha,pi,anc,amvm(5,4),mq,cv,gamm(4,3)

      amm=2d0*sqrt(mq**2+p**2)
      psiics=phi(p)
      psidf=4d0*pi*p**2 * 32d0/3d0*p**4/(amm+2d0*mq) * psiics
c     print *,p,psidf
      end



! -------------------------------------------------------- !

! =================== AMPLITUDES ========================= !

      subroutine dsidtf(dsig_l,dsig_t)
      implicit real*8(a-h,k-m,o-z)
      external resz,resk,reskap,fnul
      common/int/alz,blz,alk2,blk2,alkap2,blkap2,nlz,nlk2,nlkap2
      common/int2/lz,lk2,lkap2
      common/maincom/alpha,pi,anc,amvm(5,4),mq,cv,gamm(4,3),skew
      common/keycom/ima,nout,nvm,nst,nwf
      common/poicom/xeff,q2,w2,delta,ep
      common/psicom/a1psi,a2psi,adpsi,cc1
      common/tab/dgdtab
      dimension dgdtab(0:100,0:100)
      dimension matre(10)
      dimension at(0:15)
      character*11 char(0:15)
      data char/'  R=NL/NT  '   ! 0
     .	       ,'    r04_00 '   ! 1
     .	       ,' Re r04_10 '   ! 2
     .	       ,'    r04_1-1'   ! 3
     .	       ,'    r^1_11 '   ! 4
     .	       ,' Re r^1_10 '   ! 5
     .	       ,'    r^1_00 '   ! 6
     .	       ,'    r^1_1-1'   ! 7
     .	       ,' Im r^2_10 '   ! 8
     .	       ,' Im r^2_1-1'   ! 9
     .	       ,'    r^5_11 '   ! 10
     .	       ,' Re r^5_10 '   ! 11
     .	       ,'    r^5_00 '   ! 12
     .	       ,'    r^5_1-1'   ! 13
     .	       ,' Im r^6_10 '   ! 14
     .	       ,' Im r^6_1-1'/  ! 15

      coef=cv*sqrt(4d0*pi*alpha)
      barn=.389379d6       ! -- transition of GeV**-2 to mkbn
      mv=amvm(nvm,nst)

!------------- starting tabulation of DGD -------------------------!

      lxga = 0d0
      lxgb = 8d0
      Nxg = 100
      stepxg = (lxgb - lxga)/Nxg

      kap02 = 1d-4  ! kappa**2 min = 0.0001 GeV**2
      lka = 0d0
      lkb = 9d0    ! kappa**2 max = 10**6 GeV**2
      Nk = 100
      stepk = (lkb - lka)/Nk

      DO ixg = 0,Nxg
         lxg = lxga + ixg*stepxg
         xg = DEXP(-DLOG(10d0)*lxg) 

         DO ik = 0,Nk
            lk = lka + ik*stepk
            kap2 = kap02*DEXP(DLOG(10d0)*lk) 
            
c	    write(*,*)xg,kap2
            dgdtab(ixg,ik) = dgd(xg,kap2,delta**2,1)
c            dgdtab(ixg,ik) = dgd(xg,kap2,0d0,1)
c	    write(*,*)ixg,ik,dgdtab(ixg,ik)
            
         ENDDO

      ENDDO

c      STOP

!------------- end of tabulation of DGD ----------------------------!


      alz=0d0
c      blz=log(10d0*(q2/mv**2+1d0) )
      blz=7d0
      alk2=-10.d0
      blk2=5.d0
      alkap2=-10.d0-log(q2+mv**2)
c      alkap2=-10.d0
      blkap2=7d0
      nlz=10
      nlk2=10
      nlkap2=16

      do ima=1,5

         ! ==== Imaginary parts ==== !
c      IF (nout.LE.3) res_im=resz(resk,reskap)
      IF (nout.LE.3) res_im=reskap(resz,resk)
c      IF (nout.EQ.4) res_im=resk(resz,reskap)
      IF (nout.EQ.5) then
         z=0.5d0
         lz=0d0
         lk2=-log(100d0)
         k2=exp(lk2)*1d0        ! GeV^2
         res_im=reskap(fnul,fnul)
      END IF   
c      IF (nout.EQ.6) then       ! without kappa integration
c         lkap2=-log(1d0)
c         kap2=exp(lkap2)*1d0        ! GeV^2
c         res_im=resz(resk,fnul)
c      END IF   


         ! ==== Real parts ==== ! 
         !!! tested only for nout = 1 !!!

         factorW2 = 1.1d0
         W2old = W2

         W2 = W2old*factorW2
         xeff=skew*(q2+amvm(nvm,nst)**2)/(q2+w2)
         res1=resz(resk,reskap)
         W2 = W2old/factorW2
         xeff=skew*(q2+amvm(nvm,nst)**2)/(q2+w2)
         res2=resz(resk,reskap)

         IF(abs(res1).GT.1d-15)then
            deriv = DLOG(abs(res1/res2))/(2d0*DLOG(factorW2))
         ELSE
            deriv = 0.0d0
         ENDIF
         deriv = min(deriv,1d0)
         deriv = max(-1d0,deriv)
c         write(*,*)'lam = ',SNGL(deriv)


         res_re = -Pi/2d0*deriv*res_im ! - re

         W2 = W2old
         xeff=skew*(q2+amvm(nvm,nst)**2)/(q2+w2)


c         Rfactor = 2d0**(deriv*2+3)/DSQRT(pi)
c     .        *GammaF(deriv+2.5d0)/GammaF(deriv+4d0)
c         res_im = res_im*Rfactor  ! -im

c        write(*,'(a6,i4,g12.4)')'res: ',ima,res_im
c         write(10,'(a6,i4,g12.4)')'res: ',ima,res_im
c         write(*,*)' res: ',ima+5,res_re
c         write(10,*)' res: ',ima+5,res_re
         
	 matre(ima)=coef*res_im
	 matre(ima+5)=coef*res_re

c         write(3,*)SNGL(Q2),SNGL(res_im),SNGL(res_re),
c     .     SNGL(dsqrt(res_im**2+res_re**2)),SNGL(DATAN(res_im/res_re))

      enddo ! ima 
      

	ml0im=matre(1)
	mttim=matre(2)
	mdfim=matre(3)
	mltim=matre(4)
	mt0im=matre(5)
	ml0re=matre(6)
	mttre=matre(7)
	mdfre=matre(8)
	mltre=matre(9)
	mt0re=matre(10)

      ant= mttim**2+mttre**2
     .	  +mdfim**2+mdfre**2
     .	  +mltim**2+mltre**2
      anl= ml0im**2+ml0re**2
     .	  +2d0*(mt0im**2+mt0re**2)

c      write(*,*)sngl(q2),sngl(ml0im**2+ml0re**2),
c     . sngl(mttim**2+mttre**2),sngl(mdfim**2+mdfre**2),
c     . sngl(mltim**2+mltre**2),sngl(mt0im**2+mt0re**2)
c      write(3,*)sngl(q2),sngl(dsqrt(ml0im**2+ml0re**2)),
c     . sngl(dsqrt(mttim**2+mttre**2)),sngl(dsqrt(mdfim**2+mdfre**2)),
c     . sngl(dsqrt(mltim**2+mltre**2)),sngl(dsqrt(mt0im**2+mt0re**2))
      dsig_t = ant/16./pi*barn
      dsig_l = anl/16./pi*barn
      dsigmadt = dsig_t + ep*dsig_l
!      dsigmadt = dsig_t

c      write(*,'(7g12.6)')sngl(q2),sngl(sqrt(w2)),sngl(delta**2),dsigmadt
c      write(10,*)sngl(q2),sngl(sqrt(w2)),sngl(delta**2),
c     . sngl(dsig_l),sngl(dsig_t)

      IF (ep.LT.0.9d0) write(*,*)'Warning: small ep!',ep


      ! Igor's old...
      bbb=1d0/(ant+ep*anl)
      AT(0)=anl/ant
      AT(1)=((ML0IM**2+ML0RE**2)*EP + MltIM**2+MltRE**2)*BBB
      AT(2)=(2.*(ML0IM*Mt0IM+ML0RE*Mt0RE)*EP + MDFIM*MltIM+MDFRE*
     . MltRE - MltIM*MTTIM-MltRE*MTTRE)*BBB/2.  ! %%%%% Why minus here???
      AT(3)=-((Mt0IM**2+Mt0RE**2)*EP-MDFIM*MTTIM-MDFRE*MTTRE)*
     . BBB
      AT(4)=(MDFIM*MTTIM+MDFRE*MTTRE)*BBB
      AT(5)=(-MDFIM*MltIM-MDFRE*MltRE + MltIM*MTTIM+MltRE*MTTRE)*
     . BBB/2. ! %%%%% why minus here???
      AT(6)=-(MltIM**2+MltRE**2)*BBB
      AT(7)=(MDFIM**2+MDFRE**2 + MTTIM**2+MTTRE**2)*BBB/2.
      AT(8)=-(MDFIM*MltIM+MDFRE*MltRE + MltIM*MTTIM+MltRE*MTTRE)*
     . BBB/2. ! %%%%% why minus in the first term?
      AT(9)=(MDFIM**2+MDFRE**2-MTTIM**2-MTTRE**2)*BBB/2.
      AT(10)=-(MDFIM*Mt0IM+MDFRE*Mt0RE - Mt0IM*MTTIM-Mt0RE*MTTRE)
     . *BBB/sqrt(2.0)
      AT(11)=-(MDFIM*ML0IM+MDFRE*ML0RE - ML0IM*MTTIM-ML0RE*MTTRE + 
     . 2.*MLTIM*MT0IM+2.*MLTRE*MT0RE)*BBB/(2.*sqrt(2.0)) 
                ! %%%%% why plus in the last?
      AT(12)=-(ML0IM*MltIM+ML0RE*MltRE)*BBB*sqrt(2.0) ! %%% why minus??
      AT(13)=(MDFIM*Mt0IM+MDFRE*Mt0RE - Mt0IM*MTTIM-Mt0RE*MTTRE)*
     . BBB/sqrt(2.0)  ! %%%%% why minus in the second?
      AT(14)=-(MDFIM*ML0IM+MDFRE*ML0RE + ML0IM*MTTIM+ML0RE*MTTRE)
     . *BBB/(2.*sqrt(2.0))
      AT(15)=((MDFIM*Mt0IM+MDFRE*Mt0RE + Mt0IM*MTTIM+Mt0RE*MTTRE)*
     . BBB)/sqrt(2.0)  ! %%%%% why plus in the second?
  

c      AT(0)=anl/ant
c      AT(1)=((ML0IM**2+ML0RE**2)*EP + MltIM**2+MltRE**2)*BBB
c      AT(2)=(2.*(ML0IM*Mt0IM+ML0RE*Mt0RE)*EP - MDFIM*MltIM-MDFRE*
c     . MltRE + MltIM*MTTIM+MltRE*MTTRE)*BBB/2.  ! %%%%% Why minus here???
c      AT(3)=(-(Mt0IM**2+Mt0RE**2)*EP+MDFIM*MTTIM+MDFRE*MTTRE)*
c     . BBB
c      AT(4)=(MDFIM*MTTIM+MDFRE*MTTRE)*BBB
c      AT(5)=(MDFIM*MltIM+MDFRE*MltRE - MltIM*MTTIM-MltRE*MTTRE)*
c     . BBB/2. ! %%%%% why minus here???
c      AT(6)=-(MltIM**2+MltRE**2)*BBB
c      AT(7)=(MDFIM**2+MDFRE**2 + MTTIM**2+MTTRE**2)*BBB/2.
c      AT(8)=(MDFIM*MltIM+MDFRE*MltRE + MltIM*MTTIM+MltRE*MTTRE)*
c     . BBB/2. ! %%%%% why minus in the first term?
c      AT(9)=(MDFIM**2+MDFRE**2-MTTIM**2-MTTRE**2)*BBB/2.
c      AT(10)=-(MDFIM*Mt0IM+MDFRE*Mt0RE - Mt0IM*MTTIM-Mt0RE*MTTRE)
c     . *BBB/sqrt(2.0)
c      AT(11)=(-MDFIM*ML0IM-MDFRE*ML0RE + ML0IM*MTTIM+ML0RE*MTTRE + 
c     . 2.*MLTIM*MT0IM+2.*MLTRE*MT0RE)*BBB/(2.*sqrt(2.0)) 
c                ! %%%%% why plus in the last?
c      AT(12)=(ML0IM*MltIM+ML0RE*MltRE)*BBB*sqrt(2.0) ! %%% why minus??
c      AT(13)=(MDFIM*Mt0IM+MDFRE*Mt0RE - Mt0IM*MTTIM-Mt0RE*MTTRE)*
c     . BBB/sqrt(2.0)  ! %%%%% why minus in the second?
c      AT(14)=-(MDFIM*ML0IM+MDFRE*ML0RE + ML0IM*MTTIM+ML0RE*MTTRE)
c     . *BBB/(2.*sqrt(2.0))
c      AT(15)=((MDFIM*Mt0IM+MDFRE*Mt0RE + Mt0IM*MTTIM+Mt0RE*MTTRE)*
c     . BBB)/sqrt(2.0)  ! %%%%% why plus in the second?


       do i15=0,15
	 write(*,'(a11,2f10.4)')char(i15),(at(i15),iat=1,1)
c	 write(*,*)sngl(q2),sngl(at)
c	 write(10,*)sngl(q2),sngl(at)
       enddo

c      write(*,*)q2,at
c      write(3,*)sngl(q2),sngl(at)



      end


      real*8 function resz(f1,f2)
      implicit real*8(a-h,k-m,o-z)
      COMMON/INT/AZ,BZ,AK,BK,AKAP,BKAP,NZ,NK,NKAP
      common/int2/lz,lk2,lkap2
      common/maincom/alpha,pi,anc,amvm(5,4),mq,cv,gamm(4,3),skew
      common/keycom/ima,nout,nvm,nst,nwf
      common/poicom/xeff,q2,w2,delta,ep
      EXTERNAL F1,F2,FNUL
       A=AZ
       B=BZ
       N=NZ

       step=(b-a)/2/n
       r=0
       do i=0,2*n
         if(mod(i,2).eq.0)then
           ki=2
         else
           ki=4
         endif
         if (i.EQ.0.or.i.EQ.2*n) ki = 1
         lz=a+step*i
         z = 0.5d0*exp(-lz)
         value = f1(f2,fnul)
         r=r+ki*value
         IF (nout.EQ.3) write (*,*) SNGL(z),SNGL(value)
       enddo
       resz=r*step/3d0
      end



      real*8 function resk(f1,f2)
      implicit real*8(a-h,k-m,o-z)
      COMMON/INT/AZ,BZ,AK,BK,AKAP,BKAP,NZ,NK,NKAP
      common/int2/lz,lk2,lkap2
      common/maincom/alpha,pi,anc,amvm(5,4),mq,cv,gamm(4,3),skew
      common/keycom/ima,nout,nvm,nst,nwf
      common/poicom/xeff,q2,w2,delta,ep
      EXTERNAL F1,F2,FNUL
       A=Ak
       B=Bk
       N=Nk

       step=(b-a)/2/n
       r=0
       do i=0,2*n
         if(mod(i,2).eq.0)then
           ki=2
         else
           ki=4
         endif
         if (i.EQ.0.or.i.EQ.2*n) ki = 1
         lk2=a+step*i
         k2=exp(lk2)*1d0        ! GeV^2
         value = f1(f2,fnul)
         r=r+ki*value
         IF (nout.EQ.4) write (*,*) SNGL(sqrt(k2)),SNGL(value)
       enddo
       resk=r*step/3d0
      end



      real*8 function reskap(f1,f2)
      implicit real*8(a-h,k-m,o-z)
      COMMON/INT/AZ,BZ,AK,BK,AKAP,BKAP,NZ,NK,NKAP
      common/int2/lz,lk2,lkap2
      common/maincom/alpha,pi,anc,amvm(5,4),mq,cv,gamm(4,3),skew
      common/keycom/ima,nout,nvm,nst,nwf
      common/poicom/xeff,q2,w2,delta,ep
      EXTERNAL F1,F2,FNUL
       A=Akap
       B=Bkap
       N=Nkap

       step=(b-a)/2/n
       r=0
       do i=0,2*n
         if(mod(i,2).eq.0)then
           ki=2
         else
           ki=4
         endif
         if (i.EQ.0.or.i.EQ.2*n) ki = 1
         lkap2=a+step*i
         kap2=(q2+amvm(nvm,nst)**2)*exp(lkap2)
c         kap2=exp(lkap2)
         value = f1(f2,fnul) !*dgd(xeff,kap2,delta**2,1)
c         value = dgd(xeff,kap2,delta**2,1)
c         qb2 = (q2+amvm(nvm,nst)**2)/4d0
         r=r+ki*value
c          write (*,*) SNGL(kap2),SNGL(value)
c         write (*,*) SNGL(xeff),SNGL(kap2),SNGL(value),SNGL(r*step/3d0)
       enddo
       reskap=r*step/3d0
c         IF (nout.EQ.2)write(*,*) 'Q2 = ',SNGL(q2),' G = ',SNGL(reskap)
      end


      real*8 function fnul(fun1,fun2)
      implicit real*8(a-h,k-m,o-z)
      external amplit_f4,fun1,fun2
      common/int2/lz,lk2,lkap2

      common/maincom/alpha,pi,anc,amvm(5,4),mq,cv,gamm(4,3),skew
      common/keycom/ima,nout,nvm,nst,nwf
      common/poicom/xeff,q2,w2,delta,ep
      common/intcom/z,k,kap,qb2,iint
      common/internal/amm,mq2,kap2,k2,de2,sq2,mqim,mmim,q,kz2,r2
      common/mix/phimix
      common/tab/dgdtab
      dimension dgdtab(0:100,0:100)


      z=.5*exp(-lz)
      k2=exp(lk2)*1d0  ! GeV^2
      kap2=(q2+amvm(nvm,nst)**2)*exp(lkap2)
c      kap2=exp(lkap2)

      k=sqrt(k2)
      kap=sqrt(kap2)

      ako=k2*kap2*z/2d0 ! Jacobian

      kappa=kap
      sq2=sqrt(2d0)
      de2=delta**2
      mq2=mq**2
      z1=1d0-z
      z21=2d0*z-1d0
      amm=sqrt((mq2+k2)/z/z1)
      r2=.25d0*amm**2-mq2
      kz2=r2-k2
      qb2=mq2+z*z1*q2
      q=sqrt(q2)
      mmim=amm/(amm+2d0*mq)
      mqim=mq/(amm+2d0*mq)

        phimin=0d0
        phimax=pi
	accu=5d-2
        n200=5

        IF(nst.LE.3)then  ! --- no mixing
           
           call simpsx(phimin,phimax,n200,accu,amplit_f4,tiima1)
           tiima=tiima1
c           tiima = 4d0*sqrt(Q2)*amm*z*z*z1*z1*2d0*
c     .    (1d0/dsqrt((qb2+kap2+k2)**2-4d0*k2*kap2) - 1/(qb2 + k2))

           sqr2=sqrt(r2)
           psi=phi(sqr2)        ! --- true wave function
           tiimapsi=tiima*psi

         ELSE             ! --- mixing

           nst = 1   ! - first 1S wave part

c           z=1-z
           call simpsx(phimin,phimax,n200,accu,amplit_f4,tiima1)
c           z=1-z
c           call simpsx(phimin,phimax,n200,accu,amplit_f4,tiima2)
           tiima_1s=tiima1
           sqr2=sqrt(r2)
           psi_1s=phi(sqr2)        ! --- true wave function

           nst = 3   ! - then D wave part

c           z=1-z
           call simpsx(phimin,phimax,n200,accu,amplit_f4,tiima1)
c           z=1-z
c           call simpsx(phimin,phimax,n200,accu,amplit_f4,tiima2)
c           tiima_d=.5d0*(tiima1+tiima2)
           tiima_d=tiima1
           sqr2=sqrt(r2)
           psi_d=phi(sqr2)        ! --- true wave function

           tiimapsi = tiima_1s*psi_1s*cos(phimix) 
     .          + tiima_d*psi_d*sin(phimix) ! - together
           nst = 4   ! - back
c           write(*,*)
            
         ENDIF

c       xeff = 1e-4
c      dgsf = dgd(xeff,kap**2,delta**2,1)

           kap02 = 1d-4         ! kappa**2 min = 0.0001 GeV**2
           lxgtest = -DLOG(xeff)/DLOG(10d0)
           lkap2test = DLOG(kap2/kap02)/DLOG(10d0)
           
           lxga = 0d0
           lxgb = 8d0
           Nxg = 100
           stepxg = (lxgb - lxga)/Nxg

           lka = 0d0
           lkb = 9d0            ! kappa**2 max = 10**6 GeV**2
           Nk = 100
           stepk = (lkb - lka)/Nk

           ix = INT(lxgtest/stepxg)
           delx = lxgtest/stepxg - ix
           ik = INT(lkap2test/stepk)
           delk = lkap2test/stepk - ik
           dgsf = dgdtab(ix,ik) + (dgdtab(ix+1,ik)-dgdtab(ix,ik))*delx
     .          + (dgdtab(ix,ik+1)-dgdtab(ix,ik))*delk
     .          + (dgdtab(ix+1,ik+1)-dgdtab(ix+1,ik)
     .          - dgdtab(ix,ik+1)+dgdtab(ix,ik))*delx*delk
c        write(*,*)ix,ik,lxgtest,lkap2test,kap21,dgsf

      qeff2 = max(qb2+k**2,kappa**2)

      fnul=ako*tiimapsi/z/z1*alpha_s(qeff2)*dgsf/kap2**2  ! -- exact formula
c      fnul=ako*tiimapsi/z/z1*alpha_s(qeff2)/kap2**2  ! -- without alphas and dgd
c      fnul=ako/zn*psi/kap**2  ! -- only psi
c      fnul=ako*dgsf/kap**4  ! -- gluon density

      end

      double precision function amplit_f4(phi)
      implicit real*8(a-h,k-m,o-z)
      common/maincom/alpha,pi,anc,amvm(5,4),mq,cv,gamm(4,3),skew
      common/keycom/ima,nout,nvm,nst,nwf
      common/poicom/xeff,q2,w2,delta,ep
      common/intcom/z,k,kappa,qb2,iint
      common/internal/amm,mq2,kap2,k2,de2,sq2,
     .   mqim,mmim,q,kz2,r2


c      delta1 = delta
c      de21 = de2
c      delta = 0
c      de2 = 0

      z1=1d0-z
      z21=2d0*z-1d0

      cf=cos(phi)  ! k-delta azimutal angle
      c2f=cos(2d0*phi)
      rr2=k2+.25d0*z21**2*de2+z21*k*delta*cf
      sqf=sqrt((rr2+kap2+qb2)**2-4d0*rr2*kap2)
      f0=1d0/sqf
      fnon=1d0/sqf*(1d0-2d0*kap2/(sqf+rr2+kap2+qb2))
      f1=1d0/(k2+z1**2*de2-2d0*k*delta*z1*cf+qb2)
      f2=1d0/(k2+z**2*de2+2d0*k*delta*z*cf+qb2)

      fi2=-2d0*f0+f1+f2
      fi1k=-2d0*(k2+.5d0*z21*k*delta*cf)*fnon
     .  +(k2+z*k*delta*cf)*f2
     .  +(k2-z1*k*delta*cf)*f1
      fi1ekv=-.5d0*(-2d0*(k2*c2f+.5d0*k*delta*z21*cf)*fnon
     .  +(k2*c2f+z*delta*k*cf)*F2
     .  +(k2*c2f-z1*delta*k*cf)*F1)
      fi2ek=-1d0/sq2*k*cf*fi2
      fi2vk=-fi2ek
      fi1e=-1d0/sq2*(-2d0*(k*cf+.5d0*delta*z21)*fnon
     .  +(k*cf+z*delta)*F2
     .  +(k*cf-z1*delta)*F1)


      if(nst.le.2)then

	 sko1=1d0+.5d0*z21**2*mqim/z/z1

	 if(ima.eq.1.or.ima.eq.6)then
	   tiima=-4d0*q*amm*z**2*z1**2*fi2*sko1
	 elseif(ima.eq.2.or.ima.eq.7)then
	   tiima=mq2*fi2+fi1k + z21**2*0.5d0*fi1k*mmim
     .	    -.5d0*fi1k+mqim*k2*fi2
	 elseif(ima.eq.3.or.ima.eq.8)then
	   tiima= z21**2*fi1ekv*mmim
     .	    -fi1ekv-mqim*k2*fi2*c2f
	 elseif(ima.eq.4.or.ima.eq.9)then
	   tiima=-2d0*z*z1*z21*amm*fi1e*sko1
     .	     +mq*mmim*z21*fi2ek
	 elseif(ima.eq.5.or.ima.eq.10)then
	   tiima=2d0*q*z*z1*z21*mmim*fi2vk
	 endif

      endif

      if(nst.eq.3)then

c	 kz2=.25d0*amm**2-mq2-k2
	 sk2=k2-4d0*mq/amm*kz2

	 if(ima.eq.1.or.ima.eq.6)then
	   tiima=-q*amm*z*z1*sk2*fi2
	 elseif(ima.eq.2.or.ima.eq.7)then
	   tiima=r2*(mq2*fi2+fi1k)
     .	    + z21**2*(r2+mq2+amm*mq)*0.5d0*fi1k
     .	    -r2*0.5d0*fi1k-mq*(amm+mq)*0.5d0*k2*fi2
	 elseif(ima.eq.3.or.ima.eq.8)then
	   tiima= z21**2*(r2+mq2+amm*mq)*fi1ekv
     .	    -r2*fi1ekv+mq*(amm+mq)*0.5d0*k2*fi2*c2f
	 elseif(ima.eq.4.or.ima.eq.9)then
	   tiima=-.5d0*z21*amm*(fi1e*sk2
     .	     +mq*(amm+mq)*fi2ek)
	 elseif(ima.eq.5.or.ima.eq.10)then
	   tiima=2d0*q*z*z1*z21*fi2vk*(r2+mq2+amm*mq)
	 endif

      endif

      if(ima.eq.1.or.ima.eq.4.or.
     .   ima.eq.6.or.ima.eq.9)tiima=-tiima
c      write(*,*)tiima

      amplit_f4=tiima/pi

c      delta = delta1
c      de2 = de21

      end

! ===============  GLUON DENSITY BLOCK ======================== !

c      include 'dgd_2.2.f'
      include 'dgd.f'

! ===============  DENSITY MATRIX BLOCK ======================== !




! ========================  MATH  ============================= ! 

      subroutine mhord(ime,a1,b1,eps,f,r,ikey)
      implicit real*8(a-h,k-m,o-z)
      external f
      ikey=0
      a=a1
      b=b1
      f1=f(a)
      f2=f(b)
      if(f1*f2.gt.0d0)then
	 ikey=1
	 return
      endif
22    continue
      if(ime.eq.1)x1=b-(b-a)/(f2-f1)*f2
      if(ime.eq.2)x1=(a+b)/2d0
      ff=f(x1)
      if(ff*f1.lt.0d0)then
	 rat=b-x1
	 b=x1
	 f2=ff
	 if(rat.gt.eps)goto 22
	 r=b+rat/2d0
      else
	 rat=x1-a
	 a=x1
	 f1=ff
	 if(rat.gt.eps)goto 22
	 r=a-rat/2d0
      endif
      end


      subroutine simps(a1,b1,h1,reps1,aeps1,funct,x,ai,aih,aiabs)
      implicit real*8(a-h,k-m,o-z)
      dimension f(7),p(5)
      h=dsign(h1,b1-a1)
      s=dsign(1.d0,h)
      a=a1
      b=b1
      ai=0.d0
      aih=0.d0
      aiabs=0.d0
      p(2)=4.d0
      p(4)=4.d0
      p(3)=2.d0
      p(5)=1.d0
      if(b-a) 1,2,1
    1 reps=dabs(reps1)
      aeps=dabs(aeps1)
      do 3 k=1,7
  3   f(k)=10.d16
      x=a
      c=0.d0
      f(1)=funct(x)/3.
    4 x0=x
      if((x0+4.*h-b)*s) 5,5,6
    6 h=(b-x0)/4.
      if(h) 7,2,7
    7 do 8 k=2,7
  8   f(k)=10.d16
      c=1.d0
    5 di2=f(1)
      di3=dabs(f(1))
      do 9 k=2,5
      x=x+h
      if((x-b)*s) 23,24,24
   24 x=b
   23 if(f(k)-10.d16) 10,11,10
   11 f(k)=funct(x)/3.
   10 di2=di2+p(k)*f(k)
    9 di3=di3+p(k)*abs(f(k))
      di1=(f(1)+4.*f(3)+f(5))*2.*h
      di2=di2*h
      di3=di3*h
      if(reps) 12,13,12
   13 if(aeps) 12,14,12
   12 eps=dabs((aiabs+di3)*reps)
      if(eps-aeps) 15,16,16
   15 eps=aeps
   16 delta=dabs(di2-di1)
      if(delta-eps) 20,21,21
   20 if(delta-eps/8.) 17,14,14
   17 h=2.*h
      f(1)=f(5)
      f(2)=f(6)
      f(3)=f(7)
      do 19 k=4,7
  19  f(k)=10.d16
      go to 18
   14 f(1)=f(5)
      f(3)=f(6)
      f(5)=f(7)
      f(2)=10.d16
      f(4)=10.d16
      f(6)=10.d16
      f(7)=10.d16
   18 di1=di2+(di2-di1)/15.
      ai=ai+di1
      aih=aih+di2
      aiabs=aiabs+di3
      go to 22
   21 h=h/2.
      f(7)=f(5)
      f(6)=f(4)
      f(5)=f(3)
      f(3)=f(2)
      f(2)=10.d16
      f(4)=10.d16
      x=x0
      c=0.d0
      go to 5
   22 if(c) 2,4,2
    2 return
      end

      subroutine simpsx(a,b,np,ep,func,res)
      implicit real*8 (a-h,k-m,o-z)
      external func
      step=(b-a)/np
      call simps(a,b,step,ep,1d-3,func,ra,res,r2,r3)
      end


      Double precision Function GammaF(xf)
      implicit real*8 (a-h,k-m,o-z)
      dimension ci(0:8)
      data ci/1.0, -0.577216, 0.989056, -0.907479, 0.981728,
     .   -0.981995, 0.993149, -0.996002, 0.998106/

      IF (xf.LT.0.5)write(*,*)'Warning!!!',xf

      b = 1.0
      num = int(xf-0.5)
      IF (num.GT.0)then
         DO i=1,num
            b = b*(xf-1)
            xf = xf-1
         END DO

      ENDIF
   
      z = xf-1
      a = 0.0
      DO i=0,8
         a = a + ci(i)*z**i
      END DO   
      GammaF = a*b
c      write(*,*)a*b
      
      end
