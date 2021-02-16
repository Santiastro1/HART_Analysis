c -------------------------------------------------------------------
c       Computation of fourier modes m=1,10 as a function of radius
c       It computes the amplitude and angular position of the modes
c       in cylindrical shells equispaced (0.2 Kpc) with a deltaR growing
c       in radius
c       The output files (fouriermods0.XXXX.out) includes R,angle(rad),
c       angle(deg), amplitude (Am/A0), for m=1 to 10.
c       The only needs are the ascii files from simulations
c       (PMcrs0a0.XXXX.asc)
c --------------------------------------------------------------------
      SUBROUTINE Fouriermod(n,numstep,sntimer,atime,x,y,z,numtot,ndisk,
     +armangle)
C------------------------------------------

      INCLUDE 'parameters.h'
      Real*8 aget,age,angle,xcc,ycc,zcc,radarm,armangle,aa,a,b
      Character*6 sntime
      character*3 sntimeG
      character*5 sntimeH
      dimension x(npartm),y(npartm),z(npartm)
      dimension aget(0:numstep),radarm(0:narmbinsmax)
      dimension armangle(0:6,0:numstep,0:narmbinsmax)
      dimension aa(0:6,0:numstep,0:narmbinsmax)
      dimension a(0:6,0:numstep),b(0:6,0:numstep)

    
 
      IF(simclass.eq.'A') then
        call trans(sntime,sntimer)
        aget(n)=dble((age(dble(sntimer))-age(0.6d0)))
c        write(*,*)sntimer,age(dble(sntimer)),age(0.6d0)
      else
      if(simclass.eq.'G') then
      ntimei=int(sntimer)
      call trans2(sntimeG,ntimei)
        aget(n)=dble(atime*1e6)
       else
      call trans3(sntimeH,sntimer)
        aget(n)=dble((age(dble(sntimer))))       
      endif
      endif

         xcc=0.
         ycc=0.
         zcc=0.
         num1=0
         do i=1,ndisk
          if((abs(x(i)).lt.1.0).and.(abs(y(i)).lt.1.0).and.
     +(abs(z(i)).lt.1.0)) then
            xcc=xcc+x(i)
            ycc=ycc+y(i)
            zcc=zcc+z(i)
            num1=num1+1
           endif
         enddo
       xcc=xcc/dble(num1)
       ycc=ycc/dble(num1)
       zcc=zcc/dble(num1)

       write(*,*)xcc,ycc,zcc

       ncorrang=0
c------------------------------------------------------------------------------
c Arms
c------------------------------------------------------------------------------
     
      narmbins=nbins

      write(*,*)'Number of radial bins',narmbins

      Do i=0,Narmbins            ! set bins for arm structure
          radarm(i)=dble(i)*steparm
c      write(*,*)'Bin number ',i,' centered at ',radarm(i),' Kpc'
      EndDo


      do m=0,6
       Do i=1,Narmbins
         R3a  = (Radarm(i)+dble(i)*0.03d0)**2 ! variable width in R
         R1a  = (Radarm(i)-dble(i)*0.03d0)**2

         A(m,n)=0.
         B(m,n)=0.

         Do jp=1,Ndisk !loop to select particles in a shell

           dd = ((X(jp)-Xcc)**2+(Y(jp)-Ycc)**2)
           dz = dabs(Z(jp)-Zcc)
           If(dd.lt.R3a.and.dd.gt.R1a.and.dz.lt.dZ_max)Then

c now compute fourier components of second harmonic


            acosin=dacos((X(jp)-Xcc)/dsqrt((Y(jp)-Ycc)**2
     ++(X(jp)-Xcc)**2))
            asinus=dasin((Y(jp)-Ycc)/dsqrt((Y(jp)-Ycc)**2
     ++(X(jp)-Xcc)**2))
            If(asinus.lt.0.0) then
             acosin=2.0d0*pi-acosin
            endif
             A(m,n) = A(m,n)+dcos(dble(m)*acosin)
             B(m,n) = B(m,n)+dsin(dble(m)*acosin) 
           EndIf
         EndDo
         AA(m,n,i)=dsqrt(A(m,n)**2+B(m,n)**2) ! amplitud m harmonic
c         write(*,*)AA(0,i)
         if (m.gt.0) then
          AA(m,n,i)=AA(m,n,i)/AA(0,n,i)
c          write(*,*)'Amplitude of mode ',m,'at bin ',i,' =',aa(m,n,i)
         endif
         if(m.gt.0) armangle(m,n,i)=atan(A(m,n)/B(m,n)) ! Angle, as a maximization of fourier mode, like in Valenzuela et al. 2003
         if(m.eq.0) armangle(m,n,i)=0.
      control=acos(B(m,n)/dsqrt(B(m,n)**2.+A(m,n)**2.))
      if (armangle(m,n,i).gt.0.and.control.gt.pi/2.) then
             armangle(m,n,i)=armangle(m,n,i)-pi

       else
        if(armangle(m,n,i).lt.0.and.control.gt.pi/2.) then
             armangle(m,n,i)=armangle(m,n,i)+pi
        endif
      endif

c conditions to take into account all possible cuadrant changes
c we skip the first point (is not usefull when we use i-1), and the
c second one because some times it diverges (it is too in the center)

      if(i.gt.2) then

c You can put a threshold for the arms amplitude (AA2< or equal to that value),
c it means that the angle will be the one of the previous bin

        If (aa(m,n,i).le.0.000) then
           if(aa(m,n,i-1).gt.0.000) then
           armangle(m,n,i)=armangle(m,n,i-1)
           aa(m,n,i)=aa(m,n,i-1)
             else
              if(aa(m,n,i+1).gt.0.000) then
                armangle(m,n,i)=armangle(m,n,i+1)
                aa(m,n,i)=aa(m,n,i+1)
                else
                write(*,*)'No detecta brazos en el entorno del bin',i
              endif
           endif
           write(*,*) 'no detecta brazos, bin=',i,' snapshot a =',atime
     
        endif

c We start taking into account the case when the difference between previous
c and next angle is larger or equal to 360 deg. We also take into account
c if in the previous bin it happened and we applied the correction.

        if((abs(armangle(m,n,i-1)-armangle(m,n,i)).gt.0.8*2.*pi/dble(m)
     +).or.ncorrang.eq.1) then

c we indicate that the correction has been applied

          ncorrang=1

c If we have negative values, if the difference between bins is higher or equal 
c to 4 pi or 2 pi, we correct it substracting 2pi or 4 pi (it is usual to have 
c to correct 4 pi when in the previous bin we corrected for 2pi).

          if(armangle(m,n,i-1).gt.armangle(m,n,i)) then
             armangle(m,n,i)=armangle(m,n,i)+2.*pi/dble(m)
            else            
            armangle(m,n,i)=armangle(m,n,i)-2.*pi/dble(m)
           endif
        endif
       endif
      EndDo
      enddo 

      write(*,*)'Before opening Fourier output files'

      IF(simclass.eq.'A') then
      open(1,file=path//'fouriermods'//sntime//'.out',status='unknown')
      else
      If(simclass.eq.'G') then
      open(1,file=path//'fouriermods'//sntimeG//'.out',status='unknown')
      else
      open(1,file=path//'fouriermods'//sntimeH//'.out',status='unknown')
      endif
      endif
      Do i=2,Narmbins-1
         write(1,997)(i*steparm),(armangle(m,n,i),armangle(m,n,i)*180.
     +/3.141592,aa(m,n,i),m=0,6)
      EndDo
      close(1)
c      write(*,*)'n',n,numstep-1 
      if(n.eq.numstep-1) then
        write(*,*)'Starting fourier analyisis, function of time'

        open(20,file=path//'fourier_tm0.out',status='unknown')
        open(21,file=path//'fourier_tm1.out',status='unknown')
        open(22,file=path//'fourier_tm2.out',status='unknown')
        open(23,file=path//'fourier_tm3.out',status='unknown')
        open(24,file=path//'fourier_tm4.out',status='unknown')
        open(25,file=path//'fourier_tm5.out',status='unknown')
        open(26,file=path//'fourier_tm6.out',status='unknown')
c        open(27,file='fourier_tm7.out',status='unknown')
c        open(28,file='fourier_tm8.out',status='unknown')
c        open(29,file='fourier_tm9.out',status='unknown')
c        open(30,file='fourier_tm10.out',status='unknown')

      do l=0,numstep-1

        write(20,919)aget(l)/1e9,(aa(0,l,j),j=2,Narmbins-1)
        write(21,919)aget(l)/1e9,(aa(1,l,j),j=2,Narmbins-1)
        write(22,919)aget(l)/1e9,(aa(2,l,j),j=2,Narmbins-1)
        write(23,919)aget(l)/1e9,(aa(3,l,j),j=2,Narmbins-1)
        write(24,919)aget(l)/1e9,(aa(4,l,j),j=2,Narmbins-1)
        write(25,919)aget(l)/1e9,(aa(5,l,j),j=2,Narmbins-1)
        write(26,919)aget(l)/1e9,(aa(6,l,j),j=2,Narmbins-1)
c        write(27,919)aget(l)/1e9,(aa(7,l,j),j=2,Narmbins-1)
c        write(28,919)aget(l)/1e9,(aa(8,l,j),j=2,Narmbins-1)
c        write(29,919)aget(l)/1e9,(aa(9,l,j),j=2,Narmbins-1)
c        write(30,919)aget(l)/1e9,(aa(10,l,j),j=2,Narmbins-1)

      enddo
 
      close(20)
      close(21)
      close(22)
      close(23)
      close(24)
      close(25)
      close(26)
      close(27)
      close(28)
      close(29)
      close(30)
 
      
      endif
      write(*,*)'Back to the main program from fourier analyisis'

 919  format(300(f15.3,1x))
 997  format(34(f15.3,1x))
      return
      End
