c ---------------------------------------------------------------------
c       Compute circular velocity and surface density in
c       equispaced cylindrical shells if shell variable set to 1
c       or in equispaced espherical shells if shell variable is set to 2
c ---------------------------------------------------------------------        
c --------------------------------------------------------------------
      SUBROUTINE Vradrd(n,sntimer,atime,x,y,z,vxx,vyy,massp,npart
     +,thcell,xc,yc,zc,massc,vxc,vyc,ncomp)
C------------------------------------------

      INCLUDE 'parameters.h'

      Real*8 r,massp,massc
      character*5 sntimeH
      dimension x(npartm+nchmax),y(npartm+nchmax),z(npartm+nchmax)
      dimension npart(5),dens(nbins,2),massp(npartm+nchmax)
      dimension xc(npartm),yc(npartm),zc(npartm)
      dimension vxc(npartm),vyc(npartm),massc(npartm)
      dimension vxx(npartm+nchmax),vyy(npartm+nchmax)
      dimension vrot(nbins,2),vtg(nchmax)
      dimension vtc(npartm),posc(npartm,5),thcell(nchmax)
      dimension rwout(nbins),densm(nbins)

      nbarbins=nbins
      numtot=npart(1)+npart(5)
      ndisk=npart(2)
      nb=npart(3)
      ngas=npart(5)
      call trans3(sntimeH,sntimer)

c     Selected component

      do i=1,ncomp
       posc(i,1)=xc(i)
       posc(i,2)=yc(i)
       posc(i,3)=zc(i)
       posc(i,4)=dsqrt(xc(i)**2.d0+yc(i)**2.d0)
       posc(i,5)=dsqrt(xc(i)**2.d0+yc(i)**2.d0+zc(i)**2.d0)
       vtc(i)=-vyc(i)*posc(i,1)/posc(i,4)+vxc(i)*posc(i,2)/
     +posc(i,4) 
      enddo

c     gas component

      do i=npart(1)+1,npart(1)+npart(5)
      vtg(i-npart(1))=(-vyy(i)*x(i)+vxx(i)*y(i))/
     +dsqrt(x(i)**2.0d0+y(i)**2.0d0)
      enddo

      do i=0,nbins
        densm(i)=0.0d0
        do j=1,2
         vrot(i,j)=0.0d0
         dens(i,j)=0.0d0
        enddo
      enddo


      If(shell.eq.1) then

      step = steparm      ! dR in kpc
      j=0
      do r=rmin,rmax-step,step
      j=j+1
      rwout(j)=r+step/2

c     Selected component, vrot

         numpdi=0
         superficie=(pi*(r+step)**2.0d0-pi*r**2.0d0)*1e6 !Surface in pc^2

         do i=1,ncomp
          if(posc(i,4).gt.r.and.posc(i,4).le.(r+step)) then
           if(abs(posc(i,3)).lt.0.5d0) then 
              numpdi=numpdi+1
              vrot(j,1)=vrot(j,1)+vtc(i)!*massc(i)
              dens(j,1)=dens(j,1)+massc(i)
           endif 
          endif
         if(numpdi.eq.0) then 
           dens(j,1)=1.0d0
           numpdi=1
         endif
         enddo
c          vrot(j,1)=vrot(j,1)/dens(j,1)
          vrot(j,1)=vrot(j,1)/dble(numpdi)
          dens(j,1)=dens(j,1)/superficie

c     Selected Gas component
         
         numpgas=0

         do i=1,npart(5)
          rgas=dsqrt(x(npart(1)+i)**2.0d0+y(npart(1)+i)**2.0d0)
          if(rgas.gt.r.and.rgas.le.(r+step)) then
           if(abs(z(npart(1)+i)).lt.0.5d0) then   
            if(Thcell(i).le.3.0d4) then 
              numpgas=numpgas+1
              vrot(j,2)=vrot(j,2)+vtg(i)*massp(i+npart(1))
              dens(j,2)=dens(j,2)+massp(i+npart(1))
            endif
           endif
          endif
         if(numpgas.eq.0) then
           dens(j,2)=1.0d0 
           numpgas=1
         endif
         enddo
         vrot(j,2)=vrot(j,2)/dens(j,2)
         dens(j,2)=dens(j,2)/superficie
         jmax=j
      enddo


        vrmeanc=0.
        vrmeang=0.
      do j=1,jmax
        vrmeanc=vrmeanc+vrot(j,1)
        vrmeang=vrmeans+vrot(j,2)
      enddo
        vrmeanc=vrmeanc/dble(jmax)
        vrmeang=vrmeang/dble(jmax)
        
      write(*,*)'Mean Vrot, comp, gas: ',vrmeanc,vrmeang

c Correct the galactic rotation direction
      if(vrmeanc.lt.0.0) then
        do j=1,jmax
          do k=1,2
             vrot(j,k)=-vrot(j,k)
          enddo
        enddo
      endif

      open(2,file=path//'vrotdens'//sntimeH//'.out',status='unknown')

      do i=1,nbins
      write(2,997)rwout(i),(vrot(i,k),k=1,2),(dens(i,k),k=1,2)
      enddo
      close(2)

      else
 
      If(shell.eq.2) then

      stp = steparm      ! dR in kpc
      
      r=rmin
      do k=1,nbins
      r=exp(log(r)+(log(rmax-stp)-log(rmin))/dble(nbins))
      step=exp(log(r)+(log(rmax-stp)-log(rmin))/dble(nbins))-r
      rwout(k)=r+step/2

c     Selected component, dens (usually suposed DM component, compo=8

         numpdm=0
         numpdmm=0
       volume=4./3.*(pi*(r+step)**3.0d0-pi*r**3.0d0) !Volume in Kpc^3
       volume1=4./3.*pi*(r+step/2.0d0)**3.0d0
         do i=1,ncomp
          if(posc(i,5).gt.r.and.posc(i,5).le.(r+step)) then
              numpdm=numpdm+1
              dens(k,1)=dens(k,1)+massc(i)
c              write(*,*)massc(i)
          endif
         if(numpdm.eq.0) then
           dens(k,1)=1.0d0
           numpdm=1
         endif
         enddo
          dens(k,1)=dens(k,1)/volume

c computing Rvir from DM component

         do i=npart(2)+1,npart(1)
         
          if(dsqrt(x(i)**2.0d0+y(i)**2.0d0+z(i)**2.0d0).le.r+step/2.0d0)
     + then
            densm(k)=densm(k)+massp(i)
            numpdmm=numpdmm+1
          endif
         if(numpdmm.eq.0) then
           densm(k)=1.0d0
         endif
         enddo
          densm(k)=densm(k)/volume1

          if(k.gt.1) then 
          if(densm(k).lt.denscv.and.densm(k-1).ge.denscv) then
            rv=rwout(k)
          endif
          if(densm(k).lt.densc200.and.densm(k-1).ge.densc200) then
            r200=rwout(k)
          endif 
          endif
c     Gas component (dens)

         numpgas=0

         do i=1,npart(5)
            rgas=dsqrt(x(npart(1)+i)**2.0d0+y(npart(1)+i)**2.0d0+
     +z(npart(1)+i)**2.0d0)
          if(rgas.gt.r.and.rgas.le.(r+step)) then
            if(Thcell(i).ge.3.0d5) then
c              write(*,*)'T=',Thcell(i),massp(i+npart(1))
              numpgas=numpgas+1
              dens(k,2)=dens(k,2)+massp(i+npart(1))
            endif
           endif
         if(numpgas.eq.0) then
           dens(k,2)=1.0d0
           numpgas=1
         endif
         enddo
         dens(k,2)=dens(k,2)/volume
      enddo
         
c     Total stellar component (Mvir,Nvir,M200,N200)

         numstarv=0
         numstar200=0
         virms=0.0d0
         sm200=0.0d0
         vmaxmstar=0.0
         vminmstar=1000000000.

         do i=1,npart(2)
          radsph=dsqrt(x(i)**2.0d0+y(i)**2.0d0+z(i)**2.0d0)
          if(radsph.le.rv) then
            virms=virms+massp(i)
            numstarv=numstarv+1
            vmaxmstar=dmax1(vmaxmstar,massp(i))
            vminmstar=dmin1(vminmstar,massp(i))
          endif

          if(radsph.le.r200) then
            sm200=sm200+massp(i)
            numstar200=numstar200+1
          endif
         enddo

         write(*,*)'Minimum stellar mass=',vminmstar,'Mo'
         write(*,*)'Maximum stellar mass=',vmaxmstar,'Mo'

c     Total DM component (Mvir,Nvir,M200,N200)

         numdmv=0
         numdm200=0
         virmdm=0.0d0
         dmm200=0.0d0

         do i=npart(2)+1,npart(1)
          radsph=dsqrt(x(i)**2.0d0+y(i)**2.0d0+z(i)**2.0d0)
          if(radsph.le.rv) then
            virmdm=virmdm+massp(i)
            numdmv=numdmv+1
          endif

          if(radsph.le.r200) then
            dmm200=dmm200+massp(i)
            numdm200=numdm200+1
          endif
         enddo

c     Gas component (Mvir,Nvir,M200,N200)

         virmg=0.0d0
         virmghot=0.0d0
         virmgcold=0.0d0
         gm200=0.0d0
         gmhot200=0.0d0
         gmcold200=0.0d0

         do i=1,npart(5)
          rgas=dsqrt(x(npart(1)+i)**2.0d0+y(npart(1)+i)**2.0d0+
     +z(npart(1)+i)**2.0d0)

          if(rgas.le.rv) then
            virmg=virmg+massp(i+npart(1))
           if(Thcell(i).ge.3.0d5) then
            virmghot=virmghot+massp(i+npart(1))
           endif
           if(Thcell(i).le.3.0d4) then
            virmgcold=virmgcold+massp(i+npart(1))
           endif
          endif
          if(rgas.le.r200) then
            gm200=gm200+massp(i+npart(1))
           if(Thcell(i).ge.3.0d5) then
            gmhot200=gmhot200+massp(i+npart(1))
           endif
           if(Thcell(i).le.3.0d4) then
            gmcold200=gmcold200+massp(i+npart(1))
           endif

          endif
         enddo

      endif
      open(2,file=path//'Dens'//sntimeH//'.out',status='unknown')

      do i=1,nbins
      write(2,998)rwout(i),(dens(i,k),k=1,2)
      enddo
      close(2)
      write(*,*)'------------------------------------------------------'
      write(*,*)'Rvir= ',rv,'; R200= ',r200
      write(*,*)'  '
      write(*,*)'DM: Mvir= ',virmdm,'; Nvir= ',numdmv
      write(*,*)'    M200= ',dmm200,'; N200= ',numdm200
      write(*,*)'Stars: Mvir= ',virms,'; Nvir= ',numstarv
      write(*,*)'       M200= ',sm200,'; N200= ',numstar200
      write(*,*)'All Gas:             Mvir=',virmg
      write(*,*)'                     M200=',gm200
      write(*,*)'Hot Gas (T>3·10^5):  Mvir=',virmghot
      write(*,*)'                     M200=',gmhot200
      write(*,*)'Cold Gas (T<3·10^4): Mvir=',virmgcold
      write(*,*)'                     M200=',gmcold200
      write(*,*)'   '
      write(*,*)'TOTAL: Mvir=',virmdm+virms+virmg,'; Nvir=',numdmv+numst
     +arv
      write(*,*)'      M200=',dmm200+sm200+gm200,'; N200=',numstar200+nu
     +mdm200
      write(*,*)'------------------------------------------------------'



      endif

      write(*,*)'Back to the main program from vrot, dens computation'

 997  format(5(f17.3,1x))
 998  format(3(f17.3,1x))
      Return
      End
