c            Subroutine Vcircomp
C                 21-03-2014
C   S.Roca-Fàbrega Universitat de Barcelona
C__________________________________________________________________
C This routine computes Vcir of stellar, dark, gas and total components
C The output is returned to the main program as variable vcirc 

      Subroutine Vcircomp(n,sntimer,x,y,z,massp,npart,vcirc,nbtot)

      Include 'parameters.h' !'dmstargas_par.h'
c      Include 'a_tree.h'

c      common / RUNPARAM / boxh, Om0, Oml0, Omb0, hubble,
c     &                    aexpn, ainit, gamma

      dimension x(npartm+nchmax),y(npartm+nchmax),z(npartm+nchmax)
      dimension massp(npartm+nchmax)
      dimension npart(5)
      dimension vcirc(nbmax+1,6)
      Character*6 sntime
      character*5 sntimeH
      character*3 sntimeG

      real*8 massp
      Real*8 xt(npartm+nchmax),yt(npartm+nchmax),zt(npartm+nchmax)
      real*8 sumst,sumcell
      Real mdm(0:nbmax),mst(0:nbmax),mb(0:nbmax),mg(0:nbmax)
      real rext(0:nbmax),rmean(0:nbmax)
      Integer ndm(0:nbmax),nst(0:nbmax),nb(0:nbmax),nga(0:nbmax),npart

      IF(simclass.eq.'A') then
      call trans(sntime,sntimer)

      Open(21,file=path//'rca'//sntime//'.out',status='unknown')

      xx =-(1.-Om0)*AEXPN**3/(Om0+(1.-Om0)*AEXPN**3)
      Deltavir =(178.+82.*xx-39.*xx**2)/(1.+xx)
      Write(*,*) 'Om0= ', Om0, 'AEXPN= ', AEXPN, 'Deltavir= ', Deltavir

c.... comoving size of the zeroth-level cell in /h Mpc
      r0= hubble/aexpn/1000.

      else
      If(simclass.eq.'G') then
      ntimei=int(sntimer)
      call trans2(sntimeG,ntimei)

      Open(21,file=path//'rca'//sntimeG//'.out',status='unknown')

      aexpn=1.0
      hubble=1.0
      r0=1.0/1000.
      else
 
      call trans3(sntimeH,sntimer)
      r0= hubble/aexpn/1000.0

      Open(21,file=path//'rca'//sntimeH//'.out',status='unknown')

      endif
      endif

c.... v0 - velocity units in km/s

      v0 = 1.0
      vscale= v0
      xmaxx1=-99999
      xminn1=99999

      do i=1,npart(1)+npart(5)
         xt(i)=r0*x(i)
         yt(i)=r0*y(i)
         zt(i)=r0*z(i)
        xmaxx1=max(xmaxx1,xt(i))
        xminn1=min(xminn1,xt(i))
      enddo
         nstars=npart(2)
         nbulge=npart(3)
         ndark=npart(4)
         ngas=npart(5) 
     
      write(*,*) 'Max x and Min x= ',xmaxx1,xminn1 
      Write(*,*) 'Number of dm particles= ', ndark
      wRITE(*,*) 'NUmber of stars=', nstars
      write(*,*) 'Number of bulge particles=',nbulge
      write(*,*) 'Number of gas cells=',ngas


      Do i=0, nbmax
         rext(i)= 0.0
         rmean(i)= 0.0
         ndm(i)= 0
         nst(i)= 0
         mst(i)= 0.0
         mdm(i)= 0.0
         mb(i)= 0.0
         nb(i)= 0
         mg(i)= 0.0
         nga(i)= 0  
      End do
         
      roff=101.0 - width*rinf
      rlw=10**rinf
      Do ib=0, nbmax
         If(ib .eq. 0)then
            rext(ib)= rlw
         Else
            rext(ib)= 10**(rinf + (1./width)*ib)
         End if
      End do
      

      Do ib=0, nbmax
         If(ib .eq. 0)then
            rmean(ib)= 0.5*rlw
         Else
            rmean(ib)= 0.5*(rext(ib) + rext(ib-1))
         End if
      End do

      ndmh= 0
      dx= 0.
      dy= 0.
      dz= 0.
      dvx= 0.
      dvy= 0.
      dvz= 0.
      Halo_mass=0.
      Do i=npart(2)+npart(3)+1, npart(1)
         dx= xt(i)
         dy= yt(i)
         dz= zt(i)

         dis2= dx*dx+dy*dy+dz*dz
         dis= 1.e+3*sqrt(dis2) ! kpc/h
         If(dis .le. Rvir)then
            ndmh= ndmh + 1
            Halo_Mass=Halo_Mass+massp(i)
            If(dis .lt. rlw)then
               ibin=0
            Else
               ibin=int(width*log10(dis) + roff)-100
            End if
            ndm(ibin)= ndm(ibin) + 1
            mdm(ibin)= mdm(ibin) + massp(i)
         Endif    
      End do
      Write(*,*) 'Number of dm particles inside halo=', ndmh
      Write(*,*) 'Halo mass(Msun)= ', Halo_Mass
      nstarh= 0
      nstarRgal= 0
      mstarRgal= 0.0
      Rgal= 0.1*Rvir ! comovil kpc/h: 07/03/2011
      Do i=1, npart(2)
         dx= xt(i)
         dy= yt(i)
         dz= zt(i)

         dis2= dx*dx+dy*dy+dz*dz
         dis= 1.e+3*sqrt(dis2) ! kpc/h
         If(dis .le. Rvir)then
            nstarh= nstarh + 1
            If(dis .lt. rlw)then
               ibin=0
            Else
               ibin=int(width*log10(dis) + roff)-100
            End if
            nst(ibin)= nst(ibin) + 1
            mst(ibin)= mst(ibin) +  massp(i)
            If(dis .le. Rgal)then
               nstarRgal= nstarRgal + 1
               mstarRgal= mstarRgal + massp(i)
            End if
         Endif    
      End do
      Write(*,*) 'Number of star particles inside halo=', nstarh
      Write(*,*) 'Number of star particles inside Rgal=', nstarRgal
      ngash= 0
      ngasRgal= 0
      mgasRgal= 0.0
      Rgal= 0.1*Rvir ! comovil kpc/h: 07/03/2011
      Do i=npart(1)+1,npart(1)+npart(5)
         dx= xt(i)
         dy= yt(i)
         dz= zt(i)

         dis2= dx*dx+dy*dy+dz*dz
         dis= 1.e+3*sqrt(dis2) ! kpc/h
         If(dis .le. Rvir)then
            ngash= ngash + 1
            If(dis .lt. rlw)then
               ibin=0
            Else
               ibin=int(width*log10(dis) + roff)-100
            End if
            nga(ibin)= nga(ibin) + 1
            mg(ibin)= mg(ibin) +  massp(i)
            If(dis .le. Rgal)then
               ngasRgal= ngasRgal + 1
               mgasRgal= mgasRgal + massp(i)
            End if
         Endif
      End do
      Write(*,*) 'Number of cells inside halo=', ngash
      Write(*,*) 'Number of cells inside Rgal=', ngasRgal      
      nbulgeh= 0
      nbulgeRgal= 0
      mbulgeRgal= 0.0
      Rgal= 0.1*Rvir ! comovil kpc/h: 07/03/2011
      Do i=npart(2)+1, npart(2)+npart(3)
         dx= xt(i)
         dy= yt(i)
         dz= zt(i)

         dis2= dx*dx+dy*dy+dz*dz
         dis= 1.e+3*sqrt(dis2) ! kpc/h
         If(dis .le. Rvir)then
            nbulgeh= nbulgeh + 1
            If(dis .lt. rlw)then
               ibin=0
            Else
               ibin=int(width*log10(dis) + roff)-100
            End if
            nb(ibin)= nb(ibin) + 1
            mb(ibin)= mb(ibin) +  massp(i)
            If(dis .le. Rgal)then
               nbulgeRgal= nbulgeRgal + 1
               mbulgeRgal= mbulgeRgal + massp(i)
            End if
         Endif
      End do

      Write(*,*) 'Number of bulge particles inside halo=', nbulgeh
      Write(*,*) 'Number of bulge particles inside Rgal=', nbulgeRgal

      totmdm=0.0
      totmstar= 0.0
      totmbulge= 0.0
      totmgas= 0.0
      Vcirdm= 0.0
      Vcirstar= 0.0
      Vcirbulge= 0.0
      Vcirgas= 0.0
      Vcir= 0.0
      Do i=0, nbmax
         totmdm= mdm(i) + totmdm
         totmstar= mst(i) + totmstar
         totmbulge= mb(i) + totmbulge
         totmgas= mg(i) + totmgas
         rphy= aexpn*(rext(i)/hubble)   ! kpc
         If(ndm(i).gt.0 .and. nst(i).gt.0 )then       
            vcirdm= 1.e-5*sqrt(6.446e+11*G*totmdm/rphy)
            vcirstar= 1.e-5*sqrt(6.446e+11*G*totmstar/rphy)
            vcirbulge= 1.e-5*sqrt(6.446e+11*G*totmbulge/rphy)
            vcirgas= 1.e-5*sqrt(6.446e+11*G*totmgas/rphy)
           vcir= sqrt(vcirdm**2 + vcirstar**2 + vcirbulge**2+vcirgas**2)
            Write(21,101) i,rphy,(rext(i)/Rgal),rext(i),vcirdm,
     &      vcirstar,vcirbulge,vcirgas,vcir
            vcirc(i+1,1)=rphy
            vcirc(i+1,2)=vcir
            vcirc(i+1,3)=vcirdm
            vcirc(i+1,4)=vcirstar
            vcirc(i+1,5)=vcirbulge
            vcirc(i+1,6)=vcirgas
            nbtot=i
         End if
      End do
      close(21)

      npart(1)=nstars+nbulge+ndark
      npart(2)=nstars
      npart(3)=nbulge
      npart(4)=ndark
      npart(5)=ngas


 101  Format(i4,8(f11.4,1X))
      return
      End

