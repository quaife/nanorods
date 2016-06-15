      subroutine laplaceDLP(q,xs,ys,ns,dir1,dir2,pot)
c     Evaluates the double-layer potential for Laplace equation
c     (s1,s2) is the source strength of length ns
c     xs and ys are the source and target locations of length ns
c     (u1,u2) are the two components of the velocity field
      implicit real*8 (a-h,o-z)

      integer error
      real *8 q(ns)
c     strength of Laplace DLP
      real *8 xs(ns),ys(ns)
c     location of the source/target points
      real *8 dir1(ns),dir2(ns)  
c     dirctional derivative direction.  Typically, it is the normal
c     direction
      real *8 pot(ns)
c     x and y components of the velocity field
      
      real *8, allocatable :: source(:,:)
c     location of sources/targets
      real *8, allocatable :: charge(:)
c     charge strength of single-layer potential term
      real *8, allocatable :: dipstr(:),dipvec(:,:)
c     charge strength and direction of the 
c     double-layer potential term
      real *8, allocatable :: grad(:,:)
c     room for two components that have to be summed
c     to form the velocity field

      allocate(source(2,ns),stat=error)
      if (stat .ne. 0) then
        print*,'ERROR WITH ALLOCATING ERROR'
      endif
      allocate(charge(ns),stat=error)
      if (stat .ne. 0) then
        print*,'ERROR WITH ALLOCATING ERROR'
      endif
      allocate(dipstr(ns),stat=error)
      if (stat .ne. 0) then
        print*,'ERROR WITH ALLOCATING ERROR'
      endif
      allocate(dipvec(2,ns),stat=error)
      if (stat .ne. 0) then
        print*,'ERROR WITH ALLOCATING ERROR'
      endif
      allocate(grad(2,ns),stat=error)
      if (stat .ne. 0) then
        print*,'ERROR WITH ALLOCATING ERROR'
      endif
c     allocate memory for temporary variables
      twopi = 8.d0*datan(1.d0)

      iprec = 4 ! enough for 12 digits

      ifpot = 1 ! don't need the single-layer potential
      ifgrad = 0 ! need the gradient
      ifhess = 0 ! don't need the Hessian
      do i=1,ns
        source(1,i) = xs(i)
        source(2,i) = ys(i)
      enddo
c     set charge locations


c     START OF FORMING FIRST COMPONENET OF VELCOTIY
      ifcharge = 0 ! need a charge component
      ifdipole = 1 ! need a dipole componenet

      do i=1,ns
        charge(i) = 0.d0
        dipstr(i) = q(i)
        dipvec(1,i) = dir1(i)
        dipvec(2,i) = dir2(i)
      enddo

      call rfmm2dpart(ierr,iprec,ns,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,ifgrad,grad,ifhess,hess)
c     compute the first component in the Laplace DLP

      do i=1,ns
        pot(i) = pot(i)/twopi
      enddo

          




      end

