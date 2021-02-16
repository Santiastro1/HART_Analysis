C------------------------------------------
C                Find length and position angle of bar 
      SUBROUTINE Barprop(n,numstep,sntimer,atime,x,y,z,numtot,ndisk)
C------------------------------------------

      INCLUDE 'parameters.h'
      Real*8 aget,age,angle,xcc,ycc,zcc,bar_ecc,barangle,eccentr
     +,angmedian,xx,yy,xy,bar_length
      character*5 sntimeH
      dimension x(npartm),y(npartm),z(npartm)
      dimension radbar(0:nbins),eccentr(0:nbins),barangle(0:nbins)
      dimension ecc_copy(0:nbins),bar_copy(0:nbins)
      dimension bar_angle(0:numstep),bar_angle1(0:numstep)
      dimension aget(0:numstep),bar_ecc(0:numstep)
      dimension nBarin(0:nbins),amCos(0:nbins),amSin(0:nbins)
      dimension velbart(0:numstep),A2_Int(0:numstep),bar_amp(0:numstep)
      dimension nmom(0:numstep),bar_length(0:numstep),ecc_max(0:numstep)
      integer num1,nbinsbar

      call trans3(sntimeH,sntimer)
        aget(n)=dble((age(dble(sntimer))))

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

      step =  steparm/2.0     
      
      Nbinsbar=int(5.0/step)

      Do i=0,Nbinsbar            ! set bins for bar structure
         RadBar(i) =dble(i)*step ! radii in Kpc
      EndDo

      ecc_max(n) =0. !Loop over different radii:find bar angle and length

      Do i=1,Nbinsbar
         
         R3  = RadBar(i)**2
         R1  = RadBar(i-1)**2
      
      eccentr(i) =0.0   ! initial guess
      angle = 0.

      Do iter =1,15
         Xcos = dcos(angle)**2
         Xsin = dsin(angle)**2
         Xcs  = dsin(2.*angle)
         Xx =0.                    ! find moments of inertia for angle
         Yy =0.
         Zz =0.
         XY =0.
         nmom(n) =0

         Do jp=1,ndisk
           dd = ((X(jp)-Xcc)**2*(1.-eccentr(i)*Xcos) +
     &          (Y(jp)-Ycc)**2*(1.-eccentr(i)*Xsin) -
     &          (X(jp)-Xcc)*(Y(jp)-Ycc)*Xcs*eccentr(i)) /
     &           (1.-eccentr(i))
c           write(*,*)'r1,r3,dd='r1,r3,dd
           dz = dabs(Z(jp)-Zcc)
           If(dd.lt.R3.and.dd.gt.R1.and.dz.lt.dZ_max)Then
              nmom(n) = nmom(n) + 1
             Xx =Xx + (X(jp)-Xcc)**2
             Yy =YY + (Y(jp)-Ycc)**2
             XY =XY + (X(jp)-Xcc)*(Y(jp)-Ycc)
           EndIf
         EndDo

         CALL  Inertia(Xx,Yy,XY,angle,eccentr(i))

         If(angle.lt.0.or.angle.gt.pi)Then
c            write (*,*)  ' angle=',angle
c            write (*,*)  ' inretia=',Xx,Yy,XY
c            write (*,*)  ' iter=',iter,' bin=',i
            stop 'error of angle'
         endIf
      EndDo   ! end iter
      
      barangle(i) = angle

               bar_par = nmom(n)

      EndDo
c            Make another iteration to find bar angle and length:
c   
      ecc_max(n) =0.

      Do i=3,Nbinsbar
         dfi =0.
         If(barangle(i).gt.barangle(i-1)+pi/2.)dfi=-pi
         If(barangle(i).lt.barangle(i-1)-pi/2.)dfi=pi
         barangle(i) = barangle(i) +dfi
         If(eccentr(i).gt.ecc_max(n))Then
            ecc_max(n) = eccentr(i)
            i_max_ecc = i
         EndIf
      EndDo
      ecc_copy(1) =(eccentr(1)+eccentr(2))/2.
      bar_copy(1) =(barangle(1)+barangle(2))/2.
      Do i=2,Nbinsbar-1
         ecc_copy(i) =(eccentr(i-1)+eccentr(i)+eccentr(i+1))/3.
         bar_copy(i) =(barangle(i-1)+barangle(i)+barangle(i+1))/3.
      EndDo
      Do i=2,Nbinsbar-1
         eccentr(i) =ecc_copy(i)
         barangle(i)=bar_copy(i)
      EndDo

      i_in = max(i_max_ecc-4,4)
      i_out= min(i_in+5,Nbinsbar)
c      write(*,*)'imaxecc=',i_max_ecc-4,4,i_in
c      write(*,*)'iin=',i_in+5,Nbins,i_out
      N_angle = i_out-i_in+1
      angle = AngMedian(barangle(i_in),N_angle) ! take median angle
      ecc_out = AngMedian(eccentr(Nbinsbar-5),4)
      ecc_limit = (ecc_max(n) + ecc_out)/2.
      iout = Nbinsbar
      Do i=i_max_ecc,Nbinsbar
         if(eccentr(i).lt.ecc_limit)Then
            iout = i   ! outer limit for bar length
            goto 20
         endIf
      EndDo
      iout =i_max_ecc   ! eccentricity is too low,
                        ! use radius of max eccentricity
 20   iout =max(iout,3)
      Do j=iout,i_max_ecc-4,-1
         dev_max = 0.    ! find max angle deviation inside radius
         i = j
c         Do i=3,j
            dev  = dabs(barangle(i)-angle)
            dev1 = dabs(barangle(i)-angle-pi)
            dev2 = dabs(barangle(i)-angle+pi)
            dev_max = max(min(dev,dev1,dev2),dev_max)
c         EndDo
         write (*,'(" BarRadius=",f8.3," BarAngle+deviations=",5f8.3)')
     &    RadBar(j),
     &     barangle(i)*deg,dev_max*deg,
     &     dev*deg,dev1*deg,dev2*deg
         If(RadBar(j).lt.bar_max   ! if inside allowed R and MaxAngle
     &     .and.RadBar(j).gt.bar_min
     & .and.dev_max.lt.angLimit*pi/180.)Then ! is small => bar length
               ecc_max(n) =eccentr(j)
               bar_rad = RadBar(j)
               write (*,*)  '     found:',j, iout,i_max_ecc
               write (*,*)  '     found:',bar_rad,dev_max*deg,nmom(n)
               goto 30
         EndIf
      EndDo
      bar_angle(n)= angle*1000.    ! did not find the bar:
      ecc_max(n)  = eccentr(iout) ! use data from max_eccentricity
      bar_rad  = RadBar(iout)
      write (*,*)  ' Did not find the bar. iout=',iout,i_max_ecc

 30   bar_length(n) = bar_rad
      If(angle.lt.0.)angle=angle+pi
      If(angle.gt.pi)angle=angle-pi
      bar_angle(n)= angle
      bar_ecc(n)    = ecc_max(n)
      nmom(n)       = bar_par


c La condicion que hay a continuación se impone segun si la velocidad de 
c la barra es alta o baja (teniendo en cuenta cuantos radianes recorre en un
c paso de tiempo) y segun la direccion de rotacion de la misma
c el caso programado es por una rotación antihoraria, y de un paso de 30 
c grados
       
          bar_angle1(0)=bar_angle(0)
          if(n.gt.0) then
            if(bar_angle(n-1).gt.bar_angle(n)) then
              bar_angle1(n)=bar_angle(n)+pi 
             do i=1,1000
              if(bar_angle1(n-1).gt.bar_angle1(n)) then
                 bar_angle1(n)=bar_angle1(n)+pi
              endif
             enddo
            else
            bar_angle1(n)=bar_angle(n)
             do i=1,1000
              if(bar_angle1(n-1).gt.bar_angle1(n)) then
                 bar_angle1(n)=bar_angle1(n)+pi
              endif
             enddo
           endif  
          endif




c          write (*,'(" Angle =",f9.2," axial ratio=",f8.4,
c     &            "  Nparticles=",i8," axes=",3f8.3)')
c     &            bar_angle1(n)*deg,dsqrt(1-bar_ecc),nmom,
c     &            bar_length,bar_length*dsqrt(1-bar_ecc)


      Do i=1,Nbinsbar
         barangle(i) = barangle(i) -bar_angle(n)
      EndDo

           
c--------------------------------------------------------------------  
C               accumulate for Bar statistics:  
C               Use statistics from BarPosition 
C------------------------------------------ 

      A2_Int_real = 0.0
      A2_Int_im = 0.0
      NA2 = 0


      Do i=0,Nbinsbar
         nBarin(i) =0
         amSin(i)  =0.
         amCos(i)  =0.
      EndDo
      Do jp=1,ndisk                ! get amplitudes of sin ans cos 
         Xp= (X(jp)- Xcc)   ! scale to kpc inits 
         Yp= (Y(jp)- Ycc)
         Zp= (Z(jp)- Zcc)
            rd = Xp**2 + Yp**2
            rs = sqrt(rd)
             ind = INT(rs/step)
             ind = min(ind,Nbinsbar)
             nBarin(ind)  = nBarin(ind)  + 1
             cs  = Xp/max(rs,1.e-5)
             sn  = Yp/max(rs,1.e-5)
             amSin(ind) =amSin(ind) +2.*sn*cs
             amCos(ind) =amCos(ind) +(cs-sn)*(cs+sn)
             if(rs.le.10.0)then
                A2_Int_real =  A2_Int_real + (cs-sn)*(cs+sn)
                A2_Int_im = A2_Int_im + 2.*sn*cs
                NA2 = NA2 + 1
             endif

      EndDo

      bar_amp(n) =0.
      Do i=1,Nbinsbar
       If(nBarin(i).gt.2.and.RadBar(i).gt.bar_min
     &    .and.RadBar(i).lt.bar_max)Then
         acs = amCos(i)/max(nBarin(i),1)
         asn = amSin(i)/max(nBarin(i),1)
         amp = sqrt(acs**2+asn**2)
         if(amp.gt.bar_amp(n))Then
            i_amp =i
            bar_amp(n) = amp

         EndIf
         phi = acos(acs/amp)
         phb = asin(asn/amp)
         ppp = phb
         if(acs.lt.0.)phb = -pi-phb
       EndIf

      EndDo

      A2_Int(n) = sqrt(A2_Int_real**2 +  A2_Int_im**2)/NA2


      bar_amp(n) = 0.
      Do i=i_amp-1,i_amp+1
         acs = amCos(i)/max(nBarin(i),1)
         asn = amSin(i)/max(nBarin(i),1)
         amp = sqrt(acs**2+asn**2)
         bar_amp(n) =bar_amp(n) + amp
      EndDo
      bar_amp(n) =bar_amp(n)/3.

      write(*,*)'Bar_amp =  ',bar_amp(n),'A2_Int =  ',A2_Int(n)
      open(5,file=path//'barproperties.out',status='unknown')
      write(5,*)'# Snapage(Gyr) BarA2 A2_Int Barlen1(kpc) Barlen2(kpc) '
      write(5,*)(aget(n))/1.e9,bar_amp(n),A2_Int(n),bar_length(n),
     +bar_length(n)*dsqrt(1-bar_ecc(n))
      close(5)
      If(n.eq.numstep-1) then
        write(*,*)'Starting bar analyisis, function of time'
      velbar=0.

      open(2,file=path//'barlength.out',status='unknown')
      open(3,file=path//'barvel.out',status='unknown')
      open(4,file=path//'baramp.out',status='unknown')

      
      do l=1,numstep
       n=l
       if(l.eq.numstep-2) goto 950
      if(numstep.eq.1) then
       n=0
       goto 949
      endif
      velbart(n)=((bar_angle1(n)-bar_angle1(n-1))/(aget(n)-aget(n-1))
     ++(bar_angle1(n+1)-bar_angle1(n))/(aget(n+1)-aget(n)))/2.0d0
c      velbar=velbar+velbart(l)
 949  continue
      write(2,999)bar_angle1(n),(aget(n))/1.e9,nmom(n),bar_length(n)
     &,bar_length(n)*dsqrt(1-bar_ecc(n))
      write(3,998)(aget(n))/1.e9,velbart(n)*(3.08518e16/
     +3.1536e7) !We transform years to Gyrs and rad/yr to Km/s/Kpc
      write(4,998)(aget(n))/1.e9,A2_Int(n),bar_amp(n)
     
      enddo
      
 950  close(2)
      close(3)
      close(4)
      endif
      write(*,*)'Back to main program from bar analysis'
 998  format(f15.3,1x,f15.3,1x,f15.3)
 999  format(f15.3,1x,f15.3,1x,I7,f15.6,1x,f15.6)
      return
      End


C------------------------------------------
C                                                      Read particles
      FUNCTION AngMedian(a,N)
C------------------------------------------
      implicit real*8 (a-h,o-z)
    
      DIMENSION a(N)
      
      AngMedian = a(1)   ! make this just in case
      Do i=1,N
         aa = a(i)
         less = 0
         Do j=1,N
           If(i.ne.j.and.a(j).le.aa)less =less+1
         EndDo
         If(less.eq.INT(N/2))AngMedian=aa
      EndDo

      Return
      End
c-------------------------------------------------------------------- 
       SUBROUTINE Inertia(Xx,Yy,XY,angle,eccentr)
c-------------------------------------------------------------------- 
c               use moments of inertia to find position angle,
c               amplitude of the second harmonic, and eccentricity                 
       implicit real*8(a-h,o-z)
       PARAMETER (pi =3.14159265) ! 
c       Real*8 Xx,Yy,XY,s1,s2,ss
c       write(*,*)'C',xx,yy,xy
       If(dabs(Xx).lt.1.d-5.and.dabs(Yy).lt.1.d-5)Then ! empty bin
          angle =0.
          eccentr =0.
          return
       EndIf
         If(dabs(Xx-Yy).lt.1.d-3*(Xx+Yy))Then
            If(XY.gt.0.)Then
               angle = pi/4.
            else
               angle =3.*pi/4.
            EndIf
         Else
            angle = dATAN(2.*XY/(Xx-Yy))/2.
            If(Xx.lt.Yy)then
               angle =angle +pi*0.5
            Else
               If(XY.lt.0.)angle =angle +pi
            EndIf
         EndIf
c         write(*,*)'C',angle
         ss = Xx+Yy
         amplitude = 4.*dsqrt( (Xy*(Xy-ss)+Xx*Yy) )/ss
         s1 = 0.5*(ss -dsqrt(ss**2 -4.*(Xx*Yy-XY**2)))
         s2 = 0.5*(ss +dsqrt(ss**2 -4.*(Xx*Yy-XY**2)))
         eccentr = 1.-s1/s2    ! = 1-b**2/a**2
       Return
       End

