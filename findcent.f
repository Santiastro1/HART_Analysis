c      Subroutine to compute the galactic center in Hidrodinamical
c      simulations (HART)
c      Modified by S.Roca-Fabrega 28-04-2014
c               Option of selecting other system than the main one
c               as the center


       subroutine findcent (n,x1,y1,z1,npart,xc,yc,zc)

       include 'parameters.h'

       dimension x1(npartm),y1(npartm),z1(npartm),npart(5)
       dimension x11(npartm),y11(npartm),z11(npartm)
c       dimension den(-40:40,-40:40,-40:40)

      xc=0.
      yc=0.
      zc=0.
      xct=0.
      yct=0.
      zct=0.
      if(optcent.eq.1) then
      xcenti=xcentini*1000.
      ycenti=ycentini*1000.
      zcenti=zcentini*1000.
      else
      xcenti=0.0
      ycenti=0.0
      zcenti=0.0
      endif
      nxbord1=0
      nybord1=0
      nzbord1=0
      nxbord2=0
      nybord2=0
      nzbord2=0
      boxhe=boxh/aexpn*1000.

      if (n.eq.1) then

c We need to take into account that a resimulated halo can
c fall in the border limits of the simulated box
c We use the stars' position for that

c      do i=npart(2)+1,npart(2)+1+lspecies(2)
       do i=1,npart(2)
         If(x1(i).lt.boxhe/0.7d0/13.0d0) then
c           write(*,*)'low',x1(i)
           nxbord1=nxbord1+1
           else
           if(x1(i).gt.boxhe/0.7*(1.0d0-1.0d0/13.0d0)) then
c           write(*,*)'high',x1(i)
           nxbord2=nxbord2+1
           endif
         endif
         If(y1(i).lt.boxhe/0.7d0/13.0d0) then
           nybord1=nybord1+1
           else
           if(y1(i).gt.boxhe/0.7*(1.0d0-1.0d0/13.0d0)) then
           nybord2=nybord2+1
           endif
         endif
         If(z1(i).lt.boxhe/0.7d0/13.0d0) then
           nzbord1=nzbord1+1
           else
           if(z1(i).gt.boxhe/0.7*(1.0d0-1.0d0/13.0d0)) then
           nzbord2=nzbord2+1
           endif
         endif
      enddo
      write(*,*)nxbord1,nxbord2,nybord1,nybord2,nzbord1,nzbord2
      If(nxbord1.gt.nxbord2) then
        if(dble(nxbord1).gt.dble(npart(2))/6.0d0) boxcutx=1
        else
        if(dble(nxbord2).gt.dble(npart(2))/6.0d0) boxcutx=2
      endif
      If(nybord1.gt.nybord2) then
        if(dble(nybord1).gt.dble(npart(2))/6.0d0) boxcuty=1
        else
        if(dble(nybord2).gt.dble(npart(2))/6.0d0) boxcuty=2
      endif
      If(nzbord1.gt.nzbord2) then
        if(dble(nzbord1).gt.dble(npart(2))/6.0d0) boxcutz=1
        else
        if(dble(nzbord2).gt.dble(npart(2))/6.0d0) boxcutz=2
      endif

      do i=1,npart(1)
       if(nxbord1.ne.0.and.nxbord2.ne.0) then
        If(boxcutx.eq.1) then
        If(x1(i).gt.boxhe/0.7d0/2.0d0) x1(i)=x1(i)-boxhe
     +/0.7d0/2.0d0
        else
         if(boxcutx.eq.2) then
        If(x1(i).lt.boxhe/0.7d0/2.0d0) x1(i)=x1(i)+boxhe
     +/0.7d0/2.0d0
        endif
        endif
       endif
       if(nybord1.ne.0.and.nybord2.ne.0) then
        If(boxcuty.eq.1) then
        If(y1(i).gt.boxhe/0.7d0/2.0d0) y1(i)=y1(i)-boxhe
     +/0.7d0/2.0d0
        else
         if(boxcuty.eq.2) then
        If(y1(i).lt.boxhe/0.7d0/2.0d0) y1(i)=y1(i)+boxhe
     +/0.7d0/2.0d0
        endif
        endif
       endif
       if(nzbord1.ne.0.and.nzbord2.ne.0) then
        If(boxcutz.eq.1) then
        If(z1(i).gt.boxhe/0.7d0/2.0d0) z1(i)=z1(i)-boxhe
     +/0.7d0/2.0d0
        else
         if(boxcutz.eq.2) then
        If(z1(i).lt.boxhe/0.7d0/2.0d0) z1(i)=z1(i)+boxhe
     +/0.7d0/2.0d0
        endif
        endif
       endif
      enddo
      else
      xct=xcprev
      yct=ycprev
      zct=zcprev
      boxhe=100.0d0
      endif

      do l=1,40

       xct=xct+xc
       yct=yct+yc
       zct=zct+zc

       xcenti=xcenti-xc
       ycenti=ycenti-yc
       zcenti=zcenti-zc

       do i=1,npart(2)+1+lspecies(2)
         if(l.eq.1) then
          x11(i)=x1(i)
          y11(i)=y1(i)
          z11(i)=z1(i)
          else
         x11(i)=x11(i)-xc !
         y11(i)=y11(i)-yc !
         z11(i)=z11(i)-zc !
         endif
       enddo

       xc=0.!15800.
       yc=0.!12400.
       zc=0.!62950.
       xcs=xc
       ycs=yc
       zcs=zc
       num=0
       nums=0

      If(optcent.eq.0) then
c First aproach to compute galactic center using, first,dm 1sp meanposition
c      write(*,*)'Centering using dm 1sp'
      iini=npart(2)
      do i=iini,npart(2)+1+lspecies(2)
      if(l.eq.1) then
        rlim=2.0*boxhe/0.7d0
       else 
        rlim=boxhe/2.0d0/dble(l)
      endif
c      write(*,*)dsqrt(x11(i)**2.0d0+y11(i)**2.0d0+z11(i)**2.0d0)
      if(dsqrt(x11(i)**2.0d0+y11(i)**2.0d0+z11(i)**2.0d0).lt.rlim) 
     +  then
c        write(*,*)'in',l,rlim
        xc=xc+x11(i)
        yc=yc+y11(i)
        zc=zc+z11(i)
        num=num+1
      endif
      enddo
        xc=xc/dble(num)
        yc=yc/dble(num)
        zc=zc/dble(num)
        write(*,*)'xc,yc,zc',xc,yc,zc
      else
c If we are looking for other system than the one in the center of the
c high resolution area, we use input center coordinates from parameter.h
      if(optcent.eq.1) then
       xc=xcenti
       yc=ycenti
       zc=zcenti
      endif
      endif
c      write(*,*)'Center of DM 1st species: ',xc,yc,zc

c Second approach to compute galactic center using stars

      if(l.eq.1) then
       xcs=0.0
       ycs=0.0
       zcs=0.0
      else
      xcs=xc
      ycs=yc
      zcs=zc
      endif

c       do i=1,npart(2)+1+lspecies(2)
c         x11(i)=x11(i)-xc !
c         y11(i)=y11(i)-yc !
c         z11(i)=z11(i)-zc !
c       enddo
c      write(*,*)'Centering using stars'

      iini=1
      do i=iini,npart(2)
      if(l.eq.1) then
c        if (optcent.eq.1) then
c        rlim=boxhe/0.7d0/5.0d0
c        else
        rlim=2*boxhe/0.7d0
c        endif
       else
        rlim=boxhe/10.0d0/dble(l)
      endif
      if(dsqrt(x11(i)**2.0d0+y11(i)**2.0d0+z11(i)**2.0d0).lt.rlim)
     +  then
        xcs=xcs+x11(i)
        ycs=ycs+y11(i)
        zcs=zcs+z11(i)
        nums=nums+1
      endif
      enddo

        if(nums.ne.0) then
        xcs=xcs/dble(nums)
        ycs=ycs/dble(nums)
        zcs=zcs/dble(nums)
        else
         write(*,*)'ERROR: No particles in the centering area'
         write(*,*)'maybe you need to reduce max l value, now l=',l
         write(*,*)'Non principal systems need large l values >200'
         write(*,*)'Main systems usually need small l values, 10-50'
         stop
        endif

        if(optcent.eq.0) then
c        write(*,*)'Stellar center',xcs,ycs,zcs
        xc=(xc+xcs)/2.0d0
        yc=(yc+ycs)/2.0d0
        zc=(zc+zcs)/2.0d0
cc        xc=xcs
cc        yc=ycs
cc        zc=zcs
        else
         if(optcent.eq.1) then
          rlim=10000.0d0/dble(l)
          if(rlim.lt.200.0) then
           xc=xcs
           yc=ycs
           zc=zcs
          else
           xc=(xc+xcs)/2.0d0
           yc=(yc+ycs)/2.0d0
           zc=(zc+zcs)/2.0d0
          endif 
         endif
        endif

c        write(*,*)l,xct,yct,zct

      enddo
        xc=xct
        yc=yct
        zc=zct
        xcprev=xc
        ycprev=yc
        zcprev=zc
      write(*,*)'Center of DM 1st species + star /2: ',xc,yc,zc

c Second aproach: compute galactic center using stars max density

c Computing DM 1sp density aroud the first center found.

c       do i=-40,40
c        do j=-40,40
c         do k=-40,40
c            den(i,j,k)=0
c         enddo
c        enddo
c       enddo
c
c       do l=1,npart(2)
c
c       i=int((x1(l)-xc)/5.0)
c       j=int((y1(l)-yc)/5.0)
c       k=int((z1(l)-zc)/5.0)
c
c        if((abs(i).le.40).and.(abs(j).le.40).and.(abs(k).le.40)) then
c        
c         den(i,j,k)=den(i,j,k)+1
c
c        endif
c
c       enddo
c
c        den(0,0,0)=dble(den(0,0,0))/8.0d0
c       do i=-40,40
c         den(0,0,i)=dble(den(0,0,i))/4.0d0
c         den(i,0,0)=dble(den(i,0,0))/4.0d0
c         den(0,i,0)=dble(den(0,i,0))/4.0d0
c        do j=-40,40
c           den(0,i,j)=dble(den(0,i,j))/2.0d0
c           den(i,0,j)=dble(den(i,0,j))/2.0d0
c           den(i,j,0)=dble(den(i,j,0))/2.0d0
c        enddo
c       enddo
c   
c       denmax=0.0d0             
c 
c       do i=-40,40
c        do j=-40,40
c         do k=-40,40
c
c           If(den(i,j,k).gt.denmax) then
c             denmax=den(i,j,k)
c             imax=i
c             jmax=j
c             kmax=k 
c           endif
cc          write(*,*)i,j,k,den(i,j,k)
c         enddo
c        enddo
c       enddo
c
c       if(imax.ne.0) then
c         xc=xc+(dble(imax)+0.5)*5.0d0
c         else
c         xc=xc+(dble(imax))*5.0d0
c       endif
c       if(jmax.ne.0) then
c         yc=yc+(dble(jmax)+0.5)*5.0d0
c         else
c         yc=yc+(dble(jmax))*5.0d0
c       endif
c       if(zmax.ne.0) then
c         zc=zc+(dble(kmax)+0.5)*5.0d0
c         else
c         zc=zc+(dble(kmax))*5.0d0
c       endif
c

       return
       end
