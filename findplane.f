c     Subroutine to find the galactic plane in hidrodynamic
c     simulations (HART)
c     We used a geometric technique from pages 9 to 11 of the book
c     "Selected exercises in galactic astronomy" of I.Atanasijevic
c     D.Reidel Publishing company dordrecht-holland (astrophysics and
c     space science library)

      subroutine findplane(x1,y1,z1,vx1,vy1,vz1,npart,age,tgas)

      include 'parameters.h'

      real*8 pl,pm,pn,plf,pmf,pnf,pl1,pm1,pn1,pl2,pm2,pn2,pl3,pm3,pn3
      dimension x1(npartm+nchmax),y1(npartm+nchmax),z1(npartm+nchmax)
      dimension vx1(npartm+nchmax),vy1(npartm+nchmax),vz1(npartm+nchmax)
      dimension vxn(npartm+nchmax),vyn(npartm+nchmax),vzn(npartm+nchmax)
      dimension npart(5),age(npartm),xn(npartm+nchmax),yn(npartm+nchmax)
      dimension zn(npartm+nchmax),tgas(nchmax),plf(maxiter+10)
      dimension pmf(maxiter+10),pnf(maxiter+10)
c  ,x2(2000000),y2(2000000),z2(2000000)
      integer npart,n

      xmean=0.0d0
      ymean=0.0d0
      zmean=0.0d0
      ntest=0


      do i=1,npart(2)
         dist=dsqrt(x1(i)**2.0d0+y1(i)**2.0d0+z1(i)**2.0d0)
     

         if((dist.lt.15.0d0).and.(dist.gt.5.0d0)) then
         xmean=xmean+x1(i)
         ymean=ymean+y1(i)
         zmean=zmean+z1(i)
         ntest=ntest+1
         endif

      enddo
       if(ntest.ne.0) then
         xmean=xmean/dble(ntest)
         ymean=ymean/dble(ntest)
         zmean=zmean/dble(ntest)
         else
         write(*,*)'No particles in 15kpc region'
         xmean=0.0d0
         ymean=0.0d0
         zmean=0.0d0
       endif


       write(*,100)xmean,ymean,zmean


      do i=1,npart(1)+npart(5)

         x1(i)=x1(i)-xmean
         y1(i)=y1(i)-ymean
         z1(i)=z1(i)-zmean
         xn(i)=x1(i)
         yn(i)=y1(i)
         zn(i)=z1(i)

      enddo

      do k=1,maxiter
         plf(k)=0.0d0
         pmf(k)=0.0d0
         pnf(k)=0.0d0
      enddo


c First we compute eigen values
      do k=1,maxiter
      a11=0.0d0
      a22=0.0d0
      a33=0.0d0
      a12=0.0d0
      a13=0.0d0
      a23=0.0d0
      

      do i=1,npart(2)

          x2=xn(i)
          y2=yn(i)
          z2=zn(i)


        if(k.gt.1) then
         distz=abs(z2)
         else
         distz=0.0d0

        endif

         if(distz.lt.2.5d0) then

          dist=dsqrt(x2**2.0d0+y2**2.0d0+z2**2.0d0)
          if ((dist.lt.15.0d0).and.(dist.gt.5.0d0)) then
  
            a11=a11+x2**2.0d0
            a22=a22+y2**2.0d0
            a33=a33+z2**2.0d0
            a12=a12+x2*y2
            a13=a13+x2*z2
            a23=a23+y2*z2

         endif
        endif
      enddo

      A=-(a11+a22+a33)
      B=-(a13**2.0d0+a23**2.0d0+a12**2.0d0-a11*a22-a11*a33-a22*a33)
      C=-(2.0d0*a12*a13*a23+a11*a22*a33-a12**2.0d0*a33-a13**2.0d0*a22-
     +a23**2.0d0*a11)

c     We solve the 3rd order equation to find the eigenvalues using
c     the analytical expression

cc      Q=(A**2.0d0-3.0d0*B)/9.0d0
cc      R=(2.0d0*A**3.0d0-9.0d0*A*B+27.0d0*C)/54.d0

cc      if(R**2.0d0.lt.Q**3.0d0) then
cc      roi=1
c      write(*,*)'1: r^2 < Q^3'
cc      theta=dacos(R/dsqrt(Q**3.0d0))
 
c There are three possible eigenvalues, we will select the one
c that minimizes the mean distance to the plane

cc      alambd1=-2.0d0*sqrt(Q)*dcos(theta/3.0d0)-A/3.0d0
cc      alambd2=-2.0d0*sqrt(Q)*dcos((theta+2.0d0*pi)/3.0d0)-A/3.0d0
cc      alambd3=-2.0d0*sqrt(Q)*dcos((theta-2.0d0*pi)/3.0d0)-A/3.0d0

cc      else
cc      roi=2
c      write(*,*)'2: r^2 > Q^3'
 
cc      alambd1=(-R+dsqrt(R**2.0d0-Q**3.0d0))**(0.33333d0)+
cc     +(-R-dsqrt(R**2.0d0-Q**3.0d0))**(0.33333d0)-A/3.0d0

cc      endif

c   We solve the 3rd order equation to find the eigenvalues using 
c  an iterative Newton-Rapson method

c We will use as a first guess the elements on the diagonal, a11,a22,a33
       alambd1=a11
       alambd2=a22
       alambd3=a33

      do kiter=1,maxiter
         alambd1_0=-(alambd1**3.0d0+A*alambd1**2.0d0+B*alambd1+C)/
     +(3.0d0*alambd1**2.0d0+2.0d0*A*alambd1+B)
         alambd2_0=-(alambd2**3.0d0+A*alambd2**2.0d0+B*alambd2+C)/
     +(3.0d0*alambd2**2.0d0+2.0d0*A*alambd2+B)
         alambd3_0=-(alambd3**3.0d0+A*alambd3**2.0d0+B*alambd3+C)/     
     +(3.0d0*alambd3**2.0d0+2.0d0*A*alambd3+B)
         alambd1=alambd1+alambd1_0
         alambd2=alambd2+alambd2_0
         alambd3=alambd3+alambd3_0
      enddo

c we compute the n,m,l values with the obtained values for the
c  eighenvalue

      pn1=(a12*a23-a13*(a22-alambd1))/dsqrt(((a22-alambd1)*(a33-alambd1)
     +-a23**2.0d0)**2.0d0+(a12*(a33-alambd1)-a13*a23)**2.0d0+((a12*a23-
     +a13*(a22-alambd1))**2.0d0))
      pm1=-pn1*(a12*(a33-alambd1)-a13*a23)/(a12*a23-a13*(a22-alambd1))
      pl1=pn1*((a22-alambd1)*(a33-alambd1)-a23**2.0d0)/(a12*a23-a13*(a22
     +-alambd1))

cc      If(roi.eq.1) then

      pn2=(a12*a23-a13*(a22-alambd2))/dsqrt(((a22-alambd2)*(a33-alambd2)
     +-a23**2.0d0)**2.0d0+(a12*(a33-alambd2)-a13*a23)**2.0d0+((a12*a23-
     +a13*(a22-alambd2))**2.0d0))
      pm2=-pn2*(a12*(a33-alambd2)-a13*a23)/(a12*a23-a13*(a22-alambd2))
      pl2=pn2*((a22-alambd2)*(a33-alambd2)-a23**2.0d0)/(a12*a23-a13*(a22
     +-alambd2))

      pn3=(a12*a23-a13*(a22-alambd3))/dsqrt(((a22-alambd3)*(a33-alambd3)
     +-a23**2.0d0)**2.0d0+(a12*(a33-alambd3)-a13*a23)**2.0d0+((a12*a23-
     +a13*(a22-alambd3))**2.0d0))
      pm3=-pn3*(a12*(a33-alambd3)-a13*a23)/(a12*a23-a13*(a22-alambd3))
      pl3=pn3*((a22-alambd3)*(a33-alambd3)-a23**2.0d0)/(a12*a23-a13*(a22
     +-alambd3))     

c  Select the solution that minimizes the distance
c we use only young stars and the ones in the inner region.
      dist1=0.d0
      dist2=0.d0
      dist3=0.d0

      do i=1,npart(2)
        x2=xn(i)
        y2=yn(i)
        z2=zn(i)

        dist=dsqrt(x2**2.0d0+y2**2.0d0+z2**2.0d0)
        if(k.gt.1) then
         distz=abs(z2)
         else
         distz=0.0d0
        endif
        if(distz.lt.2.5d0) then
        if ((dist.lt.15.0d0).and.(dist.gt.5.0d0)) then
c         if (age(i).lt.1.0) then
          dist1=dist1+(pl1*x2+pm1*y2+pn1*z2)**2.0d0
          dist2=dist2+(pl2*x2+pm2*y2+pn2*z2)**2.0d0
          dist3=dist3+(pl3*x2+pm3*y2+pn3*z2)**2.0d0
         endif
c        endif
       endif
      enddo

c      write(*,*)dist1,dist2,dist3

      If(dist1.lt.dist2) then
        If(dist1.lt.dist3) then
          pl=pl1
          pm=pm1
          pn=pn1
c         write(*,*)'dist1',dist1
          else
           pl=pl3
           pm=pm3
           pn=pn3
c         write(*,*)'dist3',dist3
        endif
        else
         If(dist2.lt.dist3) then
          pl=pl2
          pm=pm2
          pn=pn2
c         write(*,*)'dist2',dist2
          else
           pl=pl3
           pm=pm3
           pn=pn3
c         write(*,*)'dist3',dist3
         endif
      endif

       write(*,*)'l,m,n=',pl,pm,pn
       dplmn=dsqrt(pl**2.0d0+pm**2.0d0+pn**2.0d0)
       If((dplmn.lt.0.9).or.(dplmn.gt.1.1)) then
      write(*,*)'ERROR: l^2+m^2+n^2 = ',pl**2.0d0+pm**2.0d0+pn**2.0d0,',
     +has to be 1 and is not.'
       goto 21
       endif
c We define the new coordinates with respect to the new plane to
c continue with the next iteration


      do i=1,npart(2)
         x2=xn(i)
         y2=yn(i)
         z2=zn(i)     

      xn(i)=(pl*pn*x2-pm*y2)/dsqrt(pl**2.0d0+pm**2.0d0)-pl*z2
      yn(i)=(pm*pn*x2+pl*y2)/dsqrt(pl**2.0d0+pm**2.0d0)-pm*z2
      zn(i)=dsqrt(pl**2.0d0+pm**2.0d0)*x2+pn*z2
      enddo
    
c Convergence conditions (tolerance factor)
 
      plf(k)=pl
      pmf(k)=pm
      pnf(k)=pn

       
      if((dabs(dabs(pnf(k)*1.0d0)-1.0d0)).lt.tolerance) goto 12
      if(k.eq.maxiter) then
       write(*,*)'Max iter reached without reaching tolerance'
       goto 12
      endif
      enddo
       
12    continue
 
       write(*,*)'The finding iteration converged to l,m,n=',plf(k),
     +pmf(k),pnf(k),' while the ideal flat case is 0,0,1'

       do l=1,k
       do i=1,npart(1)+npart(5)

      xn(i)=(plf(l)*pnf(l)*x1(i)-pmf(l)*y1(i))/dsqrt(plf(l)**2.0d0+
     +pmf(l)**2.0d0)-plf(l)*z1(i)
      yn(i)=(pmf(l)*pnf(l)*x1(i)+plf(l)*y1(i))/dsqrt(plf(l)**2.0d0+
     +pmf(l)**2.0d0)-pmf(l)*z1(i)
      zn(i)=dsqrt(plf(l)**2.0d0+pmf(l)**2.0d0)*x1(i)+pnf(l)*z1(i)
      vxn(i)=(plf(l)*pnf(l)*vx1(i)-pmf(l)*vy1(i))/dsqrt(plf(l)**2.0d0+
     +pmf(l)**2.0d0)-plf(l)*vz1(i)
      vyn(i)=(pmf(l)*pnf(l)*vx1(i)+plf(l)*vy1(i))/dsqrt(plf(l)**2.0d0+
     +pmf(l)**2.0d0)-pmf(l)*vz1(i)
      vzn(i)=dsqrt(plf(l)**2.0d0+pmf(l)**2.0d0)*vx1(i)+pnf(l)*vz1(i)
         x1(i)=xn(i)
         y1(i)=yn(i)
         z1(i)=zn(i)
         vx1(i)=vxn(i)
         vy1(i)=vyn(i)
         vz1(i)=vzn(i)
       enddo
       enddo

  21   return
  10   format(6(f10.4,1x))
 100   format('Checking the input coordinates are centered at 0,0,0',/,
     +'(xmean,ymean,zmean)= ',3(f20.6,1x))
       end 
