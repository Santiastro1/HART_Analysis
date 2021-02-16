C_________Subroutine components_____________________________
C                 21-03-2014
C   S.Roca-Fàbrega Universitat de Barcelona
C__________________________________________________________________
c
c      Selects one or several stellar galactic components depending on the
c     value of variable 'compo' in parameters.h. 1=all, 2=disk, 3=bulge, 4=halo, 5=thin disk, 
c     6=thick disk, 7=disk+bulge, 8=dark matter, 0=gas
c
c
c     WARNING!!!!
c     For the moment only 0, 1, 2, 4(you need to select your own cuts) and 8 have been implemented
c
c Modifications
c 01-04-2014: File names
c 04-03-2015: Compute cooling and X-ray luminosity
 
      Subroutine components(sntimer,
     +npart,x2,y2,z2,vx2,vy2,vz2,mass2,zII,zIa,rhoh,thcell,tbirth,rshift
     +,cell_len,xc,yc,zc,massc,vxc,vyc,vzc,vcirc,ncomp,nbtot)

       include 'parameters.h'

      real Thcell(nchmax),alphaFec(nstarmax)
      real zII(npartm+nchmax),zIa(npartm+nchmax),fec(nstarmax)
      dimension x2(npartm+nchmax),y2(npartm+nchmax),z2(npartm+nchmax)
      dimension vx2(npartm+nchmax),vy2(npartm+nchmax),vz2(npartm+nchmax)
      dimension mass2(npartm+nchmax),npart(5),tbirth(npartm),
     +rshift(npartm),agec(nstarmax),vcirc(nbmax+1,6),vci(npartm+nchmax)
      dimension xc(npartm),yc(npartm),zc(npartm)
     +,vxc(npartm),vyc(npartm),vzc(npartm)
     +,massc(npartm),vc(npartm),rhoh(nchmax),cell_len(nchmax)
      real*8 mass2,massc,rhoh,lambda,ztotcool,rhocool,tcellcool,lumX
      real*8 ne_nh,masscool,masshidro_gr,mass1200,cell_len
      character*5 sntimeH

      call trans3(sntimeH,sntimer) 

c We add the variable Vcirc. Vcirc is the circular velocity of a star
c if it was only supported by rotation (vc=sqrt(G*M_totin/r))
c To compute Vcirc for each star we interpolate the values from
c the ones obtained for cylindrical shells in Vcir_comp.f


      ne_nh=1.0
      masshidro_gr=1.0d0/(6.022e23)


      do i=1,npart(1)+npart(5)
       vci(i)=0.0
      enddo

      do j=2,nbtot
         if(vcirc(j,2).ne.0.0d0) then
         nbmx=j
         else
         if(vcirc(j,1).eq.0.0d0) then
          write(*,*)'No enough particles in bin',j
         else
         write(*,*)'WARNING! Vcirc = 0.0 at step ',j
         goto 2
         endif
         endif
      enddo

 2    continue
      mass1200=0.0d0
      do i=1,npart(1)+npart(5)
         rad=dsqrt(x2(i)**2.0d0+y2(i)**2.0d0+z2(i)**2.0d0)  
         if(rad.le.1200.d0) mass1200=mass1200+mass2(i) 
       do j=1,nbmx
        if((vcirc(j,2).ne.0.0d0).and.(vcirc(j+1,2).ne.0.0d0)) then
         if((rad.gt.vcirc(j,1)).and.(rad.le.vcirc(j+1,1))) then
         vci(i)=((vcirc(j+1,1)-rad)*vcirc(j,2)+(rad-vcirc(j,1))*
     +vcirc(j+1,2))/(vcirc(j+1,1)-vcirc(j,1))
         endif
         endif
       enddo
      enddo

      write(99,*)'Rho1200=', mass1200/(4.0d0/3.0d0*3.141592653*1.2
     +**3.0d0)*0.7*0.7/2.775e11/0.3

      ncomp=0

      open(3,file=path//'all_stars_a'//sntimeH//'.out',status='unknown')
      do i=1,npart(2)
       if(zIa(i).le.0.0000001d0) then
      write(3,13)x2(i),y2(i),z2(i),vx2(i),vy2(i),vz2(i),mass2(i),
     +log10(zIa(i)/0.00178),log10(zII(i)/zIa(i))-log10(0.0161/0.00178)
     +,tbirth(i),rshift(i),zIa(i),999999.000d0,vci(i)
       else
      write(3,13)x2(i),y2(i),z2(i),vx2(i),vy2(i),vz2(i),mass2(i),
     +log10(zIa(i)/0.00178),log10(zII(i)/zIa(i))-log10(0.0161/0.00178)
     +,tbirth(i),rshift(i),zIa(i),zII(i)/zIa(i),vci(i)
       endif
      enddo
      close(3)
 
      if(compo.eq.8) then
 
      Open(5,file=path//'DM_component_a'//sntimeH//'.out',status=
     +'unknown')
        ncomp=0
c        write(*,*)npart(2),npart(1),npartm,x2(12519381),y2(12519381),
c     +z2(12519381),vx2(12519381),vy2(12519381),vz2(12519381),
c     +mass2(12519381),vci(12519381),x2(12519481)
        do i=npart(2)+1,npart(1)
         if(dsqrt(x2(i)**2.0d0+y2(i)**2.0d0+z2(i)**2.0d0).lt.Rvir*1.5)
     + then
              ncomp=ncomp+1
              xc(ncomp)=x2(i)
              yc(ncomp)=y2(i)
              zc(ncomp)=z2(i)
              vxc(ncomp)=vx2(i)
              vyc(ncomp)=vy2(i)
              vzc(ncomp)=vz2(i)
              massc(ncomp)=mass2(i)
              vc(ncomp)=vci(i)
c              write(*,*)i,ncomp
       write(5,12)xc(ncomp),yc(ncomp),zc(ncomp),vxc(ncomp),vyc(ncomp),
     +vzc(ncomp),massc(ncomp),vc(ncomp)
c       write(*,*)i,'3'
       endif
        enddo
       close(5)
       goto 10
      endif

      Open(1,file=path//'disk_component_a'//sntimeH//'.out',status=
     +'unknown')
      open(2,file=path//'halo_component_a'//sntimeH//'.out',status=
     +'unknown')
 
      if(compo.ne.0.and.compo.ne.8) then
c  Find rotation direction
      Vsum=0.
      do i=1,npart(2)
        if(abs(z2(i)).lt.0.5d0) then
         V=vx2(i)*y2(i)/dsqrt(x2(i)**2.0d0+y2(i)**2.0d0)-
     +vy2(i)*x2(i)/dsqrt(x2(i)**2.0d0+y2(i)**2.0d0)
        Vsum=Vsum+V
        endif
      enddo
      do i=1,npart(2)
        if(compo.eq.2) then

c  Conditions to be a disk particle
c  1- to be in the plane +/- 0.5Kpc
c  2- to be inside 15 Kpc
c  3- to have [Fe/H] higher than -0.3
c  4- to have a Jz/Jc between 0.5 and 1.25
c  5- to have a cos(alpha) > 0.7 where alpha is the angle
c between Jz of each star and Jtotal.
c  6- Age selection

          If(dsqrt(x2(i)**2.0d0+y2(i)**2.0d0).lt.25.0d0) then
           if(abs(z2(i)).lt.0.5d0) then
            V=vx2(i)*y2(i)/dsqrt(x2(i)**2.0d0+y2(i)**2.0d0)-
     +vy2(i)*x2(i)/dsqrt(x2(i)**2.0d0+y2(i)**2.0d0)
c           if(log10(zIa(i)/0.00178).gt.-0.3d0) then
             if(Vsum.lt.0) V=-V !Correct galaxy rotation
             angmomjzjc=V/vci(i)
            if((angmomjzjc.gt.0.5d0).and.(angmomjzjc.lt.1.5d0)) then
        cosalpha=V*dsqrt(x2(i)**2.0d0+y2(i)**2.0d0)/dsqrt((vz2(i)*y2(i)-
     +vy2(i)*z2(i))**2.0d0+(vx2(i)*z2(i)-vz2(i)*x2(i))**2.0d0+(vy2(i)*x2
     +(i)-vx2(i)*y2(i))**2.0d0)
             if (cosalpha.gt.0.7) then
c             if((tbirth(i).gt.8.1).and.(tbirth(i).lt.12.0)) then
              ncomp=ncomp+1
              xc(ncomp)=x2(i)
              yc(ncomp)=y2(i)
              zc(ncomp)=z2(i)
              vxc(ncomp)=vx2(i)
              vyc(ncomp)=vy2(i)
              vzc(ncomp)=vz2(i)
              massc(ncomp)=mass2(i)
              vc(ncomp)=vci(i)
       Fec(ncomp)=log10(zIa(i)/0.00178) ! From Asplund 2009 A&A Review
       alphaFec(ncomp)=log10((zII(i))/zIa(i))-log10(0.0161/0.00178)
              agec(ncomp)=tbirth(i)
c       write(*,*)'No tinc clar com calculo la alpha/Fe!!' ! Hi ha Fe i
c Ni també en les explosions de SNIa!!
      write(1,11)xc(ncomp),yc(ncomp),zc(ncomp),vxc(ncomp),vyc(ncomp),
     +vzc(ncomp),massc(ncomp),fec(ncomp),alphaFec(ncomp),agec(ncomp),
     +10.0d0**(fec(ncomp))*0.00178,10.0d0**(alphafec(ncomp)+fec(ncomp))
     +*0.0161,vc(ncomp)
               endif
c              endif
             endif
            endif
c           endif
          endif
        
        else
        if(compo.eq.3) then
          continue
        endif
        if(compo.eq.4)then
c  Conditions to be a halo particle
c  0- To be inside a sphere of R=100Kpc
c  1- to have low metallicity (<-0.1)
c  2- to have a high velocity dispersion in z
c  3- to be older than 12 Gyr
      if(dsqrt(x2(i)**2.0d0+y2(i)**2.0d0+z2(i)**2.0d0).lt.Rvir) then
cc         if(log10((zII(i)+zIa(i))/0.01989).lt.-1.0) then
cc           if(tbirth(i).gt.12.0d0) then
cc             if(vz2(i).gt.200.d0) then
            ncomp=ncomp+1
            xc(ncomp)=x2(i)
            yc(ncomp)=y2(i)
            zc(ncomp)=z2(i)
            vxc(ncomp)=vx2(i)
            vyc(ncomp)=vy2(i)
            vzc(ncomp)=vz2(i)
            massc(ncomp)=mass2(i)
            vc(ncomp)=vci(i)
       Fec(ncomp)=log10(zIa(i)/0.00178) ! From Asplund 2009 A&A Review
       alphaFec(ncomp)=log10(zII(i)/zIa(i))-log10(0.0161/0.00178)
            agec(ncomp)=tbirth(i)

      write(2,11)xc(ncomp),yc(ncomp),zc(ncomp),vxc(ncomp),vyc(ncomp),
     +vzc(ncomp),massc(ncomp),fec(ncomp),alphaFec(ncomp),agec(ncomp),
     +10.0d0**(fec(ncomp))*0.00178,10.0d0**(alphafec(ncomp)+fec(ncomp))
     +*0.0161,vc(ncomp)

cc             endif
cc           endif
cc         endif
          endif
        endif
        endif

      enddo
      

      else
      
       if(compo.eq.8) then

      Open(6,file=path//'gas_component_rvir_a'//sntimeH//'.out',status=
     +'unknown')


c Compute X-ray luminosity following Crain et al 2010
c Using CLOUDY for cooling rates

      redshiftnow=1.0/sntimer-1.0

      call Set_Cooling ()
      call Set_Cooling_Rate_rs ( rs )

      do i=1,npart(5)

       if(dsqrt(x2(i+npart(1))**2.0d0+y2(i+npart(1))**2.0d0+
     +z2(i+npart(1))**2.0d0).lt.Rvir) then 

      ztotcool=log10((zIa(i+npart(1))+zII(i+npart(1)))/
     +(0.0161+0.00178))
      tcellcool=thcell(i)/10e4
      rhocool=log10(rhoh(i)/(3.086e18)**3.0d0*1.989e33*6.022e23)
c      write(*,*)'rho=',rhocool,'gr/cm^3  ','T=',tcellcool,'10^4K'
c      write(*,*)'met=',ztotcool
      call cooling(rhocool,tcellcool,ztotcool,lambda)
c      write(*,*)'Cooling func.=',lambda
      masscool=mass2(i+npart(1))*1.989e33
      rhocool=10**rhocool/(6.022e23)

      lumX=ne_nh*masscool*rhocool*lambda/(masshidro_gr)**2.0d0/1e31

c     Units, x,y,z kpc, vx,vy,vz km/s, mass Msun, rho, Msun/kpc^3, temp K
      write(6,13)x2(i+npart(1)),y2(i+npart(1)),z2(i+npart(1))
     +,vx2(i+npart(1)),vy2(i+npart(1)),vz2(i+npart(1)),mass2(i+npart(1))
     +,rhoh(i)*1.e9,log10(zIa(i+npart(1))/0.00178),
     +log10(zII(i+npart(1))/zIa(i+npart(1)))-log10(0.0161/0.00178),
     +thcell(i),vc(i),cell_len(i)
       endif
      enddo
      close(6)
   
      endif

      endif

 10   continue

      Open(6,file=path//'gas_component_5Rvir_a'//sntimeH//'.out',status=
     +'unknown')
      do i=1,npart(5)

      if(dsqrt(x2(i+npart(1))**2.0d0+y2(i+npart(1))**2.0d0+
     +z2(i+npart(1))**2.0d0).lt.5.0*Rvir) then

      ztotcool=log10((zIa(i+npart(1))+zII(i+npart(1)))/
     +(0.0161+0.00178))
      tcellcool=thcell(i)/10e4
      rhocool=log10(rhoh(i)/(3.086e18)**3.0d0*1.989e33*6.022e23)
      call cooling(rhocool,tcellcool,ztotcool,lambda)
      masscool=mass2(i+npart(1))*1.989e33
      rhocool=10**rhocool/(6.022e23)
      lumX=ne_nh*masscool*rhocool*lambda/(masshidro_gr)**2.0d0/1e31
      write(6,13)x2(i+npart(1)),y2(i+npart(1)),z2(i+npart(1))
     +,vx2(i+npart(1)),vy2(i+npart(1)),vz2(i+npart(1)),mass2(i+npart(1))
     +,rhoh(i)*1e9,log10(zIa(i+npart(1))/0.00178),
     +log10(zII(i+npart(1))/zIa(i+npart(1)))-log10(0.0161/0.00178),
     +thcell(i),vci(i+npart(1)),cell_len(i)
      endif
      enddo
      close(6) 
      write(*,*)'Number of stars selected=', ncomp
      write(*,*)'End of components selection, back to the Main program'

   
      close(1)
      close(2)
   11 format(13(f19.8,1x))
   12 format(8(f19.8,1x))
   13 format(15(f19.8,1x))
      return
      END


