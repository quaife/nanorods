      program testCall
      implicit real*8 (a-h,o-z)

      parameter (npts = 100000)


      dimension xs(npts),ys(npts)
      dimension q(npts)
      dimension dir1(npts),dir2(npts)
      dimension pot(npts),potExact(npts)
      dimension rfield(npts),cfield(npts)


      twopi = 8.d0*datan(1.d0)
      dx = twopi/dble(npts)

      do i=1,npts
        theta = dble(i-1)*dx
        xs(i) = dcos(theta)
        ys(i) = dsin(theta)
        dir1(i) = dcos(theta)
        dir2(i) = dsin(theta)
        q(i) = dcos(theta)/dble(npts)
      enddo

c      call cpu_time(t0) 
      print*,'FMM START'
      call laplaceDLP(q,xs,ys,npts,dir1,dir2,pot)
      print*,'FMM END'
c      call cpu_time(t1)

c      call cpu_time(t0) 
      print*,'DIRECT START'
      do i=1,npts
        potExact(i) = 0.d0
        do j=1,npts
          if (j .ne. i) then
            rdotn = (xs(i) - xs(j))*dir1(j) + 
     $              (ys(i) - ys(j))*dir2(j)
            rho2 = (xs(i) - xs(j))**2.d0 +
     $             (ys(i) - ys(j))**2.d0
            potExact(i) = potExact(i) - q(j)*
     $        rdotn/rho2
          endif
        enddo
        potExact(i) = potExact(i)/twopi
      enddo
      print*,'DIRECT END'
c      call cpu_time(t1)
c      print*,t1-t0


      error = 0.d0
      do i = 1,npts
        error = error + (pot(i) - potExact(i))**2.d0
      enddo
      print*,error

      end





