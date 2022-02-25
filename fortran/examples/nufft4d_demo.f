c     Demo using FINUFFT for double-precision 4d transforms in legacy fortran.
c     Does types 1,2,3, including math test against direct summation.
c     Default opts only (see simple1d1 for how to change opts).
c
c     A slight modification of drivers from the CMCL NUFFT, (C) 2004-2009,
c     Leslie Greengard and June-Yub Lee. See: cmcl_license.txt.
c
c     Tweaked by Alex Barnett to call FINUFFT 2/17/17.
c     dyn malloc; type 2 uses same input data fk0, 3/8/17
c     Also see: ../README.
c
c     Compile with, eg (GCC, multithreaded, static, paste to a single line):
c
c     gfortran nufft4d_demo.f ../directft/dirft4d.f -o nufft4d_demo
c     ../../lib-static/libfinufft.a -lstdc++ -lfftw3 -lfftw3_omp -lm -fopenmp
c
      program nufft4d_demo
      implicit none

c     our fortran-header, always needed
      include 'finufft.fh'
c
      integer i,ier,iflag,j,k1,k2,k3,k4,mx,n1,n2,n3,n4
      integer*8 ms,mt,mu,mv,nj,nk
      real*8, allocatable :: xj(:),yj(:),zj(:),pj(:),sk(:),tk(:),
     $                       uk(:),vk(:)
      real*8 err,pi,eps,salg,ealg
      parameter (pi=3.141592653589793238462643383279502884197d0)
      complex*16, allocatable :: cj(:),cj0(:),cj1(:),fk0(:),fk1(:)
c     for default opts, make a null pointer...
      type(nufft_opts), pointer :: defopts => null()
c
c     --------------------------------------------------
c     create some test data
c     --------------------------------------------------
c
      ms = 24/2
      mt = 16/2
      mu = 18/2
      mv = 14/2
      n1 = 16/2
      n2 = 18/2
      n3 = 24/2
      n4 = 14/2
      nj = n1*n2*n3*n4
      nk = ms*mt*mu*mv
c     first alloc everything
      allocate(fk0(nk))
      allocate(fk1(nk))
      allocate(sk(nk))
      allocate(tk(nk))
      allocate(uk(nk))
      allocate(vk(nk))
      allocate(xj(nj))
      allocate(yj(nj))
      allocate(zj(nj))
      allocate(pj(nj))
      allocate(cj(nj))
      allocate(cj0(nj))
      allocate(cj1(nj))
      do k4 = -n4/2, (n4-1)/2
            do k3 = -n3/2, (n3-1)/2
                  do k2 = -n2/2, (n2-1)/2
                        do k1 = -n1/2, (n1-1)/2
                        j =  (k1+n1/2+1) + (k2+n2/2)*n1 
     &                    + (k3+n3/2)*n1*n2 + (k4+n4/2)*n1*n2*n3
                        xj(j) = pi*dcos(-pi*k1/n1)
                        yj(j) = pi*dcos(-pi*k2/n2)
                        zj(j) = pi*dcos(-pi*k3/n3)
                        pj(j) = pi*dcos(-pi*k4/n4)
                        cj(j) = dcmplx(dsin(pi*j/n1),dcos(pi*j/n2))
                        enddo
                  enddo
            enddo
      enddo
c
c     -----------------------
c     start tests
c     -----------------------
c
      iflag = 1
      print*,'Starting 4D testing: ','nj =',nj,'ms,mt,mu,mv =',
     &        ms,mt,mu,mv
      do i = 1,4
         if (i.eq.1) eps=1d-4
         if (i.eq.2) eps=1d-8
         if (i.eq.3) eps=1d-12
         if (i.eq.4) eps=1d-16
c extended/quad precision tests
         if (i.eq.5) eps=1d-20
         if (i.eq.6) eps=1d-24
         if (i.eq.7) eps=1d-28
         if (i.eq.8) eps=1d-32
	 print*,' '
	 print*,' Requested precision eps =',eps
	 print*,' '
c
c     -----------------------
c     call 4D Type 1 method
c     -----------------------
c
         call dirft4d1(nj,xj,yj,zj,pj,cj,iflag,ms,mt,mu,mv,fk0)
         call finufft4d1(nj,xj,yj,zj,pj,cj,iflag,eps,ms,mt,mu,mv,
     &                  fk1,defopts,ier)
         print *, ' ier = ',ier
         call errcomp(fk0,fk1,nk,err)
         print *, ' type 1 error = ',err
c
c     -----------------------
c      call 4D Type 2 method
c     -----------------------
         call dirft4d2(nj,xj,yj,zj,pj,cj0,iflag,ms,mt,mu,mv,fk0)
         call finufft4d2(nj,xj,yj,zj,pj,cj1,iflag,eps,ms,mt,mu,mv,fk0,
     &        defopts,ier)
         print *, ' ier = ',ier
         call errcomp(cj0,cj1,nj,err)
         print *, ' type 2 error = ',err
c
c     -----------------------
c      call 4D Type3 method
c     -----------------------
         do k1 = 1, nk
            sk(k1) = 12*(dcos(k1*pi/nk))
            tk(k1) = 8*(dsin(-pi/2+k1*pi/nk))
            uk(k1) = 10*(dcos(k1*pi/nk))
            vk(k1) = 6*(dcos(-pi/2+k1*pi/nk))
      enddo

         call dirft4d3(nj,xj,yj,zj,pj,cj,iflag,nk,sk,tk,uk,vk,fk0)
         call finufft4d3(nj,xj,yj,zj,pj,cj,iflag,eps,nk,sk,tk,uk,vk,
     &                   fk1,defopts,ier)
         print *, ' ier = ',ier
         call errcomp(fk0,fk1,nk,err)
         print *, ' type 3 error = ',err
      enddo 
      stop
      end
c
c
c
c
c
      subroutine errcomp(fk0,fk1,n,err)
      implicit none
      integer*8 k,n
      complex*16 fk0(n), fk1(n)
      real *8 salg,ealg,err
c
      ealg = 0d0
      salg = 0d0
      do k = 1, n
         ealg = ealg + cdabs(fk1(k)-fk0(k))**2
         salg = salg + cdabs(fk0(k))**2
      enddo
      err = sqrt(ealg/salg)
      return
      end
