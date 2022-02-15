cc Copyright (C) 2004-2009: Leslie Greengard and June-Yub Lee 
cc Contact: greengard@cims.nyu.edu
cc 
cc This software is being released under a FreeBSD license
cc (see license.txt in this directory). 
cc
************************************************************************
      subroutine dirft4d1(nj,xj,yj,zj,pj,cj,iflag,ms,mt,mu,mv,fk)
      implicit none
      integer nj, iflag, ms, mt, mu, mv
      real*8 xj(nj), yj(nj), zj(nj), pj(nj) 
      complex*16 cj(nj),fk(-ms/2:(ms-1)/2,-mt/2:(mt-1)/2,-mu/2:(mu-1)/2,
     $           -mv/2:(mv-1)/2)
c ----------------------------------------------------------------------
c     direct computation of nonuniform FFT
c
c                        nj
c     fk(k1,k2,k3,k4) =  SUM cj(j) exp(+/-i k1 xj(j)) *
c                        j=1       exp(+/-i k2 yj(j)) *
c                                  exp(+/-i k3 zj(j)) *
c                                  exp(+/-i k4 pj(j))
c
c     for -ms/2 <= k1 <= (ms-1)/2, 
c         -mt/2 <= k2 <= (mt-1)/2
c         -mu/2 <= k3 <= (mu-1)/2
c         -mv/2 <= k4 <= (mv-1)/2
c
c     If (iflag .ge.0) the + sign is used in the exponential.
c     If (iflag .lt.0) the - sign is used in the exponential.
c
************************************************************************
      integer j, k1, k2, k3, k4
      complex*16 zf,cm1,cm2,cm3,z1n(-ms/2:(ms-1)/2), 
     $          z2n(-mt/2:(mt-1)/2),z3n(-mu/2:(mu-1)/2)
c
      do k4 = -mv/2, (mv-1)/2
         do k3 = -mu/2, (mu-1)/2
            do k2 = -mt/2, (mt-1)/2
               do k1 = -ms/2, (ms-1)/2
                  fk(k1,k2,k3,k4) = dcmplx(0d0,0d0)
               enddo
            enddo
         enddo
      enddo
c
      do j = 1, nj
c
c     ----------------------------------------------------------
c     Precompute exponential for exp(+/-i k1 xj)
c     ----------------------------------------------------------
c
         if (iflag .ge. 0) then
            zf = dcmplx(dcos(xj(j)),+dsin(xj(j)))
         else
            zf = dcmplx(dcos(xj(j)),-dsin(xj(j)))
         endif
         z1n(0) = (1d0,0d0)
         do k1 = 1, (ms-1)/2
            z1n(k1) = zf*z1n(k1-1)
            z1n(-k1)= dconjg(z1n(k1))
         enddo
         if (ms/2*2.eq.ms) z1n(-ms/2) = dconjg(zf*z1n(ms/2-1))
c
c     ----------------------------------------------------------
c     Precompute exponential for exp(+/-i k2 yj)
c     ----------------------------------------------------------
         if (iflag .ge. 0) then
            zf = dcmplx(dcos(yj(j)),+dsin(yj(j)))
         else
            zf = dcmplx(dcos(yj(j)),-dsin(yj(j)))
         endif
         z2n(0) = (1d0,0d0)
         do k1 = 1, (mt-1)/2
            z2n(k1) = zf*z2n(k1-1)
            z2n(-k1)= dconjg(z2n(k1))
         enddo
         if (mt/2*2.eq.mt) z2n(-mt/2) = dconjg(zf*z2n(mt/2-1))
c
c     ----------------------------------------------------------
c     Precompute exponential for exp(+/-i k3 zj)
c     ----------------------------------------------------------
         if (iflag .ge. 0) then
            zf = dcmplx(dcos(zj(j)),+dsin(zj(j)))
         else
            zf = dcmplx(dcos(zj(j)),-dsin(zj(j)))
         endif
         z3n(0) = (1d0,0d0)
         do k1 = 1, (mu-1)/2
            z3n(k1) = zf*z3n(k1-1)
            z3n(-k1)= dconjg(z3n(k1))
         enddo
         if (mu/2*2.eq.mu) z3n(-mu/2) = dconjg(zf*z3n(mu/2-1))
c
c     ----------------------------------------------------------
c     Loop over k4 for pj
c     ----------------------------------------------------------
c
         if (iflag .ge. 0) then
            zf = dcmplx(dcos(pj(j)),+dsin(pj(j)))
         else
            zf = dcmplx(dcos(pj(j)),-dsin(pj(j)))
         endif
c
         cm3 = cj(j)
         do k4 = 0, (mv-1)/2
            do k3 = -mu/2, (mu-1)/2
               cm2 = cm3 * z3n(k3)
               do k2 = -mt/2, (mt-1)/2
                  cm1 = cm2 * z2n(k2)
                  do k1 = -ms/2, (ms-1)/2
                     fk(k1,k2,k3,k4) = fk(k1,k2,k3,k4)
     &                               + cm1 * z1n(k1)
                  enddo
               enddo
            enddo
            cm3 = zf*cm3
         enddo
c
         zf = dconjg(zf)
         cm3 = cj(j)
         do k4 = -1, -mv/2, -1
            cm3 = zf*cm3
            do k3 = -mu/2, (mu-1)/2
                cm2 = cm3 * z3n(k3)
                do k2 = -mt/2, (mt-1)/2
                    cm1 = cm2 * z2n(k2)
                    do k1 = -ms/2, (ms-1)/2
                        fk(k1,k2,k3,k4) = fk(k1,k2,k3,k4) 
     &                                  + cm1 * z1n(k1)
                    enddo
                enddo
            enddo
        enddo

      enddo
      end
c
c
c
c
c
************************************************************************
      subroutine dirft4d2(nj,xj,yj,zj,pj,cj,iflag,ms,mt,mu,mv,fk)
      implicit none
      integer nj, iflag, ms, mt, mu, mv
      real*8 xj(nj), yj(nj), zj(nj), pj(nj)
      complex*16 cj(nj),fk(-ms/2:(ms-1)/2,-mt/2:(mt-1)/2,-mu/2:(mu-1)/2,
     $           -mv/2:(mv-1)/2)
c ----------------------------------------------------------------------
c     direct computation of nonuniform FFT
c
c
c     cj(j) = SUM SUM SUM SUM  fk(k1,k2,k3,k4) exp(+/-i k1 xj(j)) * 
c             k1  k2  k3  k4                   exp(+/-i k2 yj(j)) *
c                                              exp(+/-i k3 zj(j)) *
c                                              exp(+/-i k4 pj(j))
c
c                            for j = 1,...,nj
c
c     where -ms/2 <= k1 <= (ms-1)/2 
c           -mt/2 <= k2 <= (mt-1)/2
c           -mu/2 <= k3 <= (mu-1)/2
c           -mv/2 <= k4 <= (mv-1)/2
c
c
c     If (iflag .ge.0) the + sign is used in the exponential.
c     If (iflag .lt.0) the - sign is used in the exponential.
c
************************************************************************
      integer j, k1, k2, k3, k4
      complex*16 zf, cm1,cm2,cm3,cm4,z1n(-ms/2:(ms-1)/2),
     $               z2n(-mt/2:(mt-1)/2),z3n(-mu/2:(mu-1)/2)
c
      do j = 1, nj
c
c     ----------------------------------------------------------
c     Precompute exponential for exp(+/-i k1 xj)
c     ----------------------------------------------------------
         if (iflag .ge. 0) then
            zf = dcmplx(dcos(xj(j)),+dsin(xj(j)))
         else
            zf = dcmplx(dcos(xj(j)),-dsin(xj(j)))
         endif
         z1n(0) = (1d0,0d0)
         do k1 = 1, (ms-1)/2
            z1n(k1) = zf*z1n(k1-1)
            z1n(-k1)= dconjg(z1n(k1))
         enddo
         if (ms/2*2.eq.ms) z1n(-ms/2) = dconjg(zf*z1n(ms/2-1))
c
c     ----------------------------------------------------------
c     Precompute exponential for exp(+/-i k2 yj)
c     ----------------------------------------------------------
c
         if (iflag .ge. 0) then
            zf = dcmplx(dcos(yj(j)),+dsin(yj(j)))
         else
            zf = dcmplx(dcos(yj(j)),-dsin(yj(j)))
         endif
         z2n(0) = (1d0,0d0)
         do k2 = 1, (mt-1)/2
            z2n(k2) = zf*z2n(k2-1)
            z2n(-k2)= dconjg(z2n(k2))
         enddo
         if (mt/2*2.eq.mt) z2n(-mt/2) = dconjg(zf*z2n(mt/2-1))
c
c     ----------------------------------------------------------
c     Precompute exponential for exp(+/-i k3 zj)
c     ----------------------------------------------------------
c
         if (iflag .ge. 0) then
            zf = dcmplx(dcos(zj(j)),+dsin(zj(j)))
         else
            zf = dcmplx(dcos(zj(j)),-dsin(zj(j)))
         endif
         z3n(0) = (1d0,0d0)
         do k3 = 1, (mu-1)/2
            z3n(k3) = zf*z3n(k3-1)
            z3n(-k3)= dconjg(z3n(k3))
         enddo
         if (mu/2*2.eq.mu) z3n(-mu/2) = dconjg(zf*z3n(mu/2-1))
c
c     ----------------------------------------------------------
c     Loop over k4 for pj
c     ----------------------------------------------------------
         if (iflag .ge. 0) then
            zf = dcmplx(dcos(pj(j)),+dsin(pj(j)))
         else
            zf = dcmplx(dcos(pj(j)),-dsin(pj(j)))
         endif
c
         cm3 = (0d0, 0d0)
         do k3 = -mu/2, (mu-1)/2
            cm2 = (0d0, 0d0)
            do k2 = -mt/2, (mt-1)/2
                cm1 = (0d0, 0d0)
                do k1 = -ms/2, (ms-1)/2
                    cm1 = cm1 + z1n(k1) * fk(k1,k2,k3,0)
                enddo
                cm2 = cm2 + z2n(k2) * cm1
            enddo
            cm3 = cm3 + z3n(k3) * cm2
         enddo
         cj(j) = cm3
c
         cm4 = zf
         do k4 = 1, (mv-1)/2
            cm3 = (0d0, 0d0)
            do k3 = -mu/2, (mu-1)/2
                cm2 = (0d0, 0d0)
                do k2 = -mt/2, (mt-1)/2
                    cm1 = (0d0, 0d0)
                    do k1 = -ms/2, (ms-1)/2
                        cm1 = cm1 + z1n(k1) * fk(k1,k2,k3,k4)
                    enddo
                    cm2 = cm2 + z2n(k2) * cm1
                enddo
                cm3 = cm3 + z3n(k3) * cm2
            enddo
            cj(j) = cj(j) + cm4 * cm3

            cm3 = (0d0, 0d0)
            do k3 = -mu/2, (mu-1)/2
                cm2 = (0d0, 0d0)
                do k2 = -mt/2, (mt-1)/2
                    cm1 = (0d0, 0d0)
                    do k1 = -ms/2, (ms-1)/2
                        cm1 = cm1 + z1n(k1) * fk(k1,k2,k3,-k4)
                    enddo
                    cm2 = cm2 + z2n(k2) * cm1
                enddo
                cm3 = cm3 + z3n(k3) * cm2
            enddo
            cj(j) = cj(j) + dconjg(cm4) * cm3
            cm4 = cm4*zf
        enddo
c
         if (mv/2*2.eq.mv) then
            cm3 = (0d0, 0d0)
            do k3 = -mu/2, (mu-1)/2
               cm2 = (0d0, 0d0)
                  do k2 = -mt/2, (mt-1)/2
                    cm1 = (0d0, 0d0)
                    do k1 = -ms/2, (ms-1)/2
                        cm1 = cm1 + z1n(k1) * fk(k1,k2,k3,-mv/2)
                    enddo
                    cm2 = cm2 + z2n(k2) * cm1
                enddo
                cm3 = cm3 + z3n(k3) * cm2
            enddo
            cj(j) = cj(j) + dconjg(cm4) * cm3
         endif
      enddo
      end
c
c
c
c
c
************************************************************************
      subroutine dirft4d3(nj,xj,yj,zj,pj,cj,iflag,nk,sk,tk,uk,vk,fk)
      implicit none
      integer nj,iflag,nk
      real*8 xj(nj),yj(nj),zj(nj),pj(nj),sk(nk),tk(nk),uk(nk),vk(nk)
      complex*16 cj(nj),fk(nk)
c ----------------------------------------------------------------------
c     direct computation of nonuniform FFT
c
c              nj
c     fk(k) = SUM cj(j) exp(+/-i s(k) xj(j)) *
c             j=1       exp(+/-i t(k) yj(j)) *
c                       exp(+/-i u(k) zj(j)) * 
c                       exp(+/-i v(k) pj(j))
c
c                    for k = 1, ..., nk
c
c     If (iflag .ge.0) the + sign is used in the exponential.
c     If (iflag .lt.0) the - sign is used in the exponential.
c
************************************************************************
      integer k, j
      real*8 ssk, stk, suk, svk
c
      do k = 1, nk
         if (iflag .ge. 0) then
            ssk =  sk(k)
            stk =  tk(k)
            suk =  uk(k)
            svk =  vk(k)
         else
            ssk =  -sk(k)
            stk =  -tk(k)
            suk =  -uk(k)
            svk =  -vk(k)
         endif
c
         fk(k) = dcmplx(0d0,0d0)
         do j = 1, nj
            fk(k) = fk(k) + cj(j) * dcmplx
     &        ( dcos(ssk*xj(j)+stk*yj(j)+suk*zj(j)+svk*pj(j)),
     &          dsin(ssk*xj(j)+stk*yj(j)+suk*zj(j)+svk*pj(j)) )
         enddo
      enddo
      end

