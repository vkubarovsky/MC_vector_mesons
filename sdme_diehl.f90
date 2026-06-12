!==============================================================
!     sdme_diehl.f90
!     Diehl formalism for W_UU and translation to/from
!     Schilling-Wolf SDMEs
!
!     Diehl decomposes W_UU by meson helicity:
!       W_UU = (3/4pi)[cos^2(th) W^LL(Phi)
!              + sqrt(2) cos(th) sin(th) W^LT(Phi,phi)
!              + sin^2(th) W^TT(Phi,phi)]
!
!     15 effective Diehl parameters d(1..15):
!       d(1)  = u^{00}_{++} + eps*u^{00}_{00}    WLL const
!       d(2)  = Re u^{00}_{0+}                    WLL cosΦ
!       d(3)  = u^{00}_{-+}                       WLL cos2Φ
!       d(4)  = Re(u^{0+}_{++} - u^{-0}_{++} + 2eps u^{0+}_{00})
!                                                  WLT cos(phi)
!       d(5)  = Re(u^{0+}_{0+} - u^{-0}_{0+})    WLT cos(Φ+phi)
!       d(6)  = Re u^{0+}_{-+}                    WLT cos(2Φ+phi)
!       d(7)  = Re(u^{0-}_{0+} - u^{+0}_{0+})    WLT cos(Φ-phi)
!       d(8)  = Re u^{+0}_{-+}                    WLT cos(2Φ-phi)
!       d(9)  = Re(u^{-+}_{++} + eps u^{-+}_{00}) WTT cos(2phi)
!       d(10) = Re(u^{++}_{0+} + u^{--}_{0+})     WTT cosΦ
!       d(11) = Re u^{-+}_{0+}                    WTT cos(Φ+2phi)
!       d(12) = Re u^{++}_{-+}                    WTT cos2Φ
!       d(13) = Re u^{+-}_{0+}                    WTT cos(Φ-2phi)
!       d(14) = u^{-+}_{-+}                       WTT cos(2Φ+2phi)
!       d(15) = u^{+-}_{-+}                       WTT cos(2Φ-2phi)
!
!     Normalization: int dPhi/(2pi) int dphi d(costh) W_UU = 1
!       so W_UU = 2*pi * W_SW
!
!     Convention: Phi = Diehl's phi (production vs lepton plane)
!                 phi = Diehl's varphi (decay vs production plane)
!                 th  = Diehl's vartheta (decay polar angle)
!==============================================================


!==============================================================
!     diehl_W: W_UU in Diehl decomposition
!     Returns W_UU(costh, phi, Phi; eps, d)
!==============================================================
      function diehl_W(costh, phid, aPhi, eps, d)
      implicit real*8(a-h,o-z)
      real*8 d(15), diehl_W

      aPI = acos(-1d0)
      sinth2 = 1d0 - costh*costh
      sinth  = sqrt(max(0d0, sinth2))
      s2th   = 2d0*sinth*costh         ! sin(2*th)

      eTL = sqrt(max(0d0, eps*(1d0+eps)))

!     ── W^LL (multiplies cos^2 theta) ────────────────────
      wLL = d(1) &
          - 2d0*cos(aPhi)*eTL*d(2) &
          - cos(2d0*aPhi)*eps*d(3)

!     ── W^LT (multiplies sqrt(2)*cos(th)*sin(th)) ───────
      wLT = sqrt(eps*(1d0+eps))*cos(aPhi+phid)*d(5) &
           - cos(phid)*d(4) &
           + eps*cos(2d0*aPhi+phid)*d(6) &
           - sqrt(eps*(1d0+eps))*cos(aPhi-phid)*d(7) &
           + eps*cos(2d0*aPhi-phid)*d(8)

!     ── W^TT (multiplies sin^2 theta) ───────────────────
      wTT = 0.5d0*(1d0 - d(1)) &
          + 0.5d0*cos(2d0*aPhi+2d0*phid)*eps*d(14) &
          - cos(aPhi)*eTL*d(10) &
          + cos(aPhi+2d0*phid)*eTL*d(11) &
          - cos(2d0*phid)*d(9) &
          - cos(2d0*aPhi)*eps*d(12) &
          + cos(aPhi-2d0*phid)*eTL*d(13) &
          + 0.5d0*cos(2d0*aPhi-2d0*phid)*eps*d(15)

!     ── Full W_UU ────────────────────────────────────────
      diehl_W = (3d0/(4d0*aPI)) * ( &
                costh*costh * wLL &
              + s2th/sqrt(2d0) * wLT &
              + sinth2 * wTT )

      if(diehl_W .lt. 0d0) diehl_W = 0d0

      return
      end


!==============================================================
!     sw_to_diehl: Convert Schilling-Wolf r(15) to Diehl d(15)
!     Pure algebraic mapping, epsilon-independent
!
!     SW r(15):                      Diehl d(15):
!     r(1) = r04_00                  d(1)  = r04_00
!     r(2) = Re r04_10              d(2)  = -r5_00/sqrt(2)
!     r(3) = r04_1-1                d(3)  = r1_00
!     r(4) = r1_11                  d(4)  = 2*Re r04_10
!     r(5) = r1_00                  d(5)  = -sqrt(2)*(Re r5_10 + Im r6_10)
!     r(6) = Re r1_10               d(6)  = Re r1_10 + Im r2_10
!     r(7) = r1_1-1                 d(7)  = sqrt(2)*(Re r5_10 - Im r6_10)
!     r(8) = Im r2_10               d(8)  = Re r1_10 - Im r2_10
!     r(9) = Im r2_1-1              d(9)  = r04_1-1
!     r(10)= r5_11                  d(10) = -sqrt(2)*r5_11
!     r(11)= r5_00                  d(11) = -(r5_1-1 + Im r6_1-1)/sqrt(2)
!     r(12)= Re r5_10               d(12) = r1_11
!     r(13)= r5_1-1                 d(13) = -(r5_1-1 - Im r6_1-1)/sqrt(2)
!     r(14)= Im r6_10               d(14) = r1_1-1 + Im r2_1-1
!     r(15)= Im r6_1-1              d(15) = r1_1-1 - Im r2_1-1
!==============================================================
      subroutine sw_to_diehl(r, d)
      implicit real*8(a-h,o-z)
      real*8 r(15), d(15)

      sq2 = sqrt(2d0)

      d(1)  = r(1)                          ! r04_00
      d(2)  = -r(11)/sq2                    ! -r5_00/sqrt(2)
      d(3)  = r(5)                          ! r1_00
      d(4)  = 2d0*r(2)                      ! 2*Re r04_10
      d(5)  = -sq2*(r(12) + r(14))          ! -sqrt(2)*(Re r5_10 + Im r6_10)
      d(6)  = r(6) + r(8)                   ! Re r1_10 + Im r2_10
      d(7)  = sq2*(r(12) - r(14))           ! sqrt(2)*(Re r5_10 - Im r6_10)
      d(8)  = r(6) - r(8)                   ! Re r1_10 - Im r2_10
      d(9)  = r(3)                          ! r04_1-1
      d(10) = -sq2*r(10)                    ! -sqrt(2)*r5_11
      d(11) = -(r(13) + r(15))/sq2          ! -(r5_1-1 + Im r6_1-1)/sqrt(2)
      d(12) = r(4)                          ! r1_11
      d(13) = -(r(13) - r(15))/sq2          ! -(r5_1-1 - Im r6_1-1)/sqrt(2)
      d(14) = r(7) + r(9)                   ! r1_1-1 + Im r2_1-1
      d(15) = r(7) - r(9)                   ! r1_1-1 - Im r2_1-1

      return
      end


!==============================================================
!     diehl_to_sw: Convert Diehl d(15) to Schilling-Wolf r(15)
!     Inverse of sw_to_diehl
!==============================================================
      subroutine diehl_to_sw(d, r)
      implicit real*8(a-h,o-z)
      real*8 r(15), d(15)

      sq2 = sqrt(2d0)

      r(1)  = d(1)                          ! r04_00
      r(2)  = d(4)/2d0                      ! Re r04_10
      r(3)  = d(9)                          ! r04_1-1
      r(4)  = d(12)                         ! r1_11
      r(5)  = d(3)                          ! r1_00
      r(6)  = (d(6) + d(8))/2d0             ! Re r1_10
      r(7)  = (d(14) + d(15))/2d0           ! r1_1-1
      r(8)  = (d(6) - d(8))/2d0             ! Im r2_10
      r(9)  = (d(14) - d(15))/2d0           ! Im r2_1-1
      r(10) = -d(10)/sq2                    ! r5_11
      r(11) = -sq2*d(2)                     ! r5_00
      r(12) = -(d(5) - d(7))/(2d0*sq2)     ! Re r5_10
      r(13) = -(d(11) + d(13))*sq2/2d0      ! r5_1-1 = -(d11+d13)/sqrt(2)
      r(14) = -(d(5) + d(7))/(2d0*sq2)     ! Im r6_10
      r(15) = (-d(11) + d(13))*sq2/2d0      ! Im r6_1-1 = (-d11+d13)/sqrt(2)

      return
      end


!==============================================================
!     diehl_Wmax: upper bound for accept/reject on diehl_W
!==============================================================
      function diehl_Wmax(eps, d)
      implicit real*8(a-h,o-z)
      real*8 d(15), diehl_Wmax, diehl_W
      external diehl_W

      aPI = acos(-1d0)
      wmax_val = 0d0
      ng = 40

      do i = 0, ng
        acosth = -1d0 + 2d0*dble(i)/dble(ng)
        do j = 0, ng
          aphid = 2d0*aPI*dble(j)/dble(ng)
          do k = 0, ng
            aPhiV = 2d0*aPI*dble(k)/dble(ng)
            wval = diehl_W(acosth, aphid, aPhiV, eps, d)
            if(wval .gt. wmax_val) wmax_val = wval
          enddo
        enddo
      enddo

      diehl_Wmax = 1.2d0 * wmax_val
      if(diehl_Wmax .lt. 1d-20) diehl_Wmax = 1d0

      return
      end


!==============================================================
!     sample_diehl_angles: Accept/reject on diehl_W
!==============================================================
      subroutine sample_diehl_angles(costh, phid, aPhiOut, &
                                     eps, d, wmax, iy)
      implicit real*8(a-h,o-z)
      real*4 urand
      real*8 d(15), diehl_W
      external diehl_W
      integer*4 iy

      aPI = acos(-1d0)

 10   continue
        costh   = 2d0*dble(urand(iy)) - 1d0
        phid    = 2d0*aPI*dble(urand(iy))
        aPhiOut = 2d0*aPI*dble(urand(iy))
        rtest   = dble(urand(iy))

        wval = diehl_W(costh, phid, aPhiOut, eps, d)

        if(rtest*wmax .gt. wval) goto 10

      return
      end
