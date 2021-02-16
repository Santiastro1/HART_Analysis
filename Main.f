C_________MAIN N-BODY ANALYSIS PROGRAM_____________________________
C                 00-00-2013
C__________________________________________________________________

c       implicit real*8(a-h,o-z)

       INCLUDE 'parameters.h'

      dimension massp(npartm+nchmax)
      dimension npart(5)
      dimension xt(nchmax+npartm),yt(nchmax+npartm),zt(nchmax+npartm)
      dimension xc(npartm),yc(npartm)
      dimension zc(npartm),massc(npartm)
      dimension vxc(npartm),vyc(npartm)
      dimension vzc(npartm) 
      dimension Thcell(nchmax),rhoh(nchmax),vxt(nchmax+npartm)
      dimension vyt(nchmax+npartm),vzt(nchmax+npartm)
      dimension zII(nchmax+npartm),zIa(nchmax+npartm)
      dimension vcirc(nbmax+10,5),kappaq(npartm),ome2q(npartm)
      dimension armanglet(0:6,0:500,0:narmbinsmax)
      dimension tbirth(npartm),zshift(npartm),sntimerarray(1000)
      dimension cell_len(nchmax)
      real*8 massp,massc,rhoh
      integer n

      numstep=0
      write(*,*)'An inputfiles.out file with snapshot times is needed'
      write(*,*)'inputfiles.out needs to be located at ',path1
      write(*,*)'The format of snapshots time needs to be X.XXX in '
      write(*,*)'scale factor units, and decreasing from highest a.'
      open(1,file=path1//'inputfiles.out',status='old')
      do i=1,1000
       read(1,*,end=3)sntimerarray(i)
       numstep=numstep+1
      enddo
  3   close(1)

      write(*,*)'inputfiles.out has been read correctly.'
      write(*,*)numstep,' snapshots will be analysed.'

       do i=1,numstep


          write(*,*)'Before reading file'

c           sntimer=aini+dble(i)*asteps
          sntimer=sntimerarray(i)

       CALL READ_HART_cosmo(i,sntimer,npart,xt,yt,zt,vxt,vyt,vzt,massp,
     +zII,zIa,rhoh,thcell,tbirth,zshift,cell_len)



          write(*,*)'File aexp=',sntimer,' readed'

      nstar=npart(2)
      ndark=npart(4)
      ntot=npart(1)
      nbulge=npart(3)
      ngas=npart(5)
      
      if(ntot.gt.npartm) write(*,*)'ERROR npart > npartmax'
c      Spectrograms

cc       write(*,*)'Start spectrograms analyisis'

cc       CALL Spectrograms(i,sntimer,atime,x,y,z,vxx,vyy,vzz,npart)

c      Rotation curve

       write(*,*)'Start Vcirc analysis'

       CALL Vcircomp(i,sntimer,xt,yt,zt,massp,npart,vcirc,nbtot) 

cc      write(*,*)'Mark and write tipsy marked'

cc      call marktip(npart,tbirth,thcell,massp,xt,yt,zt,vxt,vyt,vzt,zII
cc     +,zIa,vcirc,nbtot)

c     Separete the components (thin-thick disk, bulge, halo)

      write(*,*)'Start components identification'

      If (compo.eq.1) then
      CALL components(sntimer,npart,xt,yt,zt,vxt,vyt,vzt,massp,zII,zIa,
     +rhoh,thcell,tbirth,zshift,cell_len,xc,yc,zc,massc,vxc,vyc,vzc,
     +vcirc,ncomp,nbtot)
        do j=1,npart(2)
           xc(j)=xt(j)
           yc(j)=yt(j)
           zc(j)=zt(j)
           vxc(j)=vxt(j)
           vyc(j)=vyt(j)
           vzc(j)=vzt(j)
           massc(j)=massp(j)
        enddo
           ncomp=npart(2)
           goto 15

         else

      CALL components(sntimer,npart,xt,yt,zt,vxt,vyt,vzt,massp,zII,zIa,
     +rhoh,thcell,tbirth,zshift,cell_len,xc,yc,zc,massc,vxc,vyc,vzc,
     +vcirc,ncomp,nbtot)

      endif
  15  continue

  
c      Disk parameters

       write(*,*)'Start disk analysis'
       CALL dispQ(i,numstep,sntimer,atime,xc,yc,zc,vxc,vyc,vzc,massc,
     +ncomp,nbtot,vcirc,kappaq,ome2q,j)

c      Resonance computation

c       write(*,*)'Start resonance computation'
c       CALL reson(sntimer,kappaq,ome2q,j)

     

c      Compute density and velocity rotation curve
       
       write(*,*)'Start vr and dens computation'
       Call Vradrd(i,sntimer,atime,xt,yt,zt,vxt,vyt,massp,npart,
     +thcell,xc,yc,zc,massc,vxc,vyc,ncomp)

c      Compute vertex deviation and density 


cc       write(*,*)'Start cartesian lv analysis'
cc       Call lv_analysis(i,sntimer,atime,x,y,z,vxx,vyy,vzz,npart,vcirc
cc     +,nbtot)


cc       write(*,*)'Start polar lv analysis'
cc       Call lv_analysis_polar(i,sntimer,atime,x,y,z,vxx,vyy,vzz,npart
cc     +,vcirc,nbtot)


c      Fourier analysis
     
cc       write(*,*)ncomp

       write(*,*)'Start Fourier analyisis'
       Call Fouriermod(i,numstep,sntimer,atime,xc,yc,zc,ncomp,ncomp,
     +armanglet)

c     Finding the bar
      
       write(*,*)'Start Bar analyisis'
       call Barprop(i,numstep,sntimer,atime,xc,yc,zc,ncomp,ncomp)

c      Spiral arms velocity curve

cc       If(i.eq.numstep-1) then

cc       write(*,*)'Start arms velocity computation'
cc       Call armsvel(armanglet)

cc       endif

c      Spiral arms pitch angle

cc       write(*,*)'Start pitch angle computation'
cc       Call pitchangle(i,sntimer,armanglet)
 

       enddo
       stop
       end
