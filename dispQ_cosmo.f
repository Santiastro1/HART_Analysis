c S. Roca-Fàbrega 2013
c Last update (buggs correction) August 2016

        Subroutine dispQ(n,numstep,sntimer,atime,x,y,z,vxx,vyy,vzz,massp
     +,ndisk,nbtot,vcirc,kappaq,ome2q,j)

c-------------------------------------------------------------------------
c        Compute sigmaR,sigmaPHI,sigmaZ,Kappa,Ome,Dens sup,Q,X
c        as a function of radius, in cylindrical shells
c        Needs the rotation curves 
c-------------------------------------------------------------------------        
      INCLUDE 'parameters.h'
      character*5 sntimeH
      integer i, n0, j, kk,ngas,nstar,ndim,numcap,numrotcurv
      integer nmean
      real*8 mass,radivi,curvrot,vlrs,sigma,ro,massp,
     +kappa,kappaq,kappast,kappasb,ome2q,ome2t,ome2b,ome2,
     +sigUrmean,sigVrmean,sigWrmean,kappaqrmean,ome2qrmean,
     +rormean,toomrermean,xirmean,aget,age,xmean
      dimension XRSR(nbmax+1),YRSR(nbmax+1)
      dimension ZRSR(nbmax+1),sigma(nbmax+1,3),xlv(nbmax+1)
      dimension xmean(nbmax+1,3)
      dimension umean(nbmax+1),vmean(nbmax+1),nmean(nbmax+1)
      dimension kappa(nbmax+1),ro(nbmax+1)
      dimension pos(npartm,3),vel(npartm,6),mass(npartm)
      dimension kappaq(nbmax+1),ome2(nbmax+1),dome2(nbmax+1),
     +toomre(nbmax+1),Xi(nbmax+1),ome2q(nbmax+1)
      dimension x(npartm),y(npartm),z(npartm),massp(npartm)
      dimension vxx(npartm),vyy(npartm),vzz(npartm)
      dimension vcirc(nbmax+1,5),curvrot(nbmax+1,5)
      dimension sigUrmean(0:numstep),sigVrmean(0:numstep),
     +sigWrmean(0:numstep),kappaqrmean(0:numstep),ome2qrmean(0:numstep),
     +rormean(0:numstep),toomrermean(0:numstep),xirmean(0:numstep),
     +aget(0:numstep)

        n0=0
        numrotcurv=0
        sigUrmean(n)=0.0
        sigVrmean(n)=0.0
        sigWrmean(n)=0.0
        kappaqrmean(n)=0.0
        ome2qrmean(n)=0.0
        rormean(n)=0.0
        toomrermean(n)=0.0
        xirmean(n)=0.0


      call trans3(sntimeH,sntimer)
        aget(n)=dble((age(dble(sntimer))))
          numrotcurv=nbtot+1         

        do i=1,numrotcurv
           curvrot(i,1)=vcirc(i,1)
           curvrot(i,2)=vcirc(i,2)
           ome2(i)=(curvrot(i,2)/curvrot(i,1))**2.d0
        enddo
           
           do i=4,numrotcurv-4

      dome2(i)=(ome2(i-4)/280.0d0-4.d0*ome2(i-3)/105.d0+ome2(i-2)/5.0d0-
     +4.0d0*ome2(i-1)/5.0d0+4.0d0*ome2(i+1)/5.0d0-ome2(i+2)/5.0d0+4.0d0*
     +ome2(i+3)/105.0d0-ome2(i+4)/280.0d0)/
     +dabs((curvrot(i,1)-curvrot(i+1,1)))

           enddo

      ome2(numrotcurv)=(curvrot(numrotcurv,2)/curvrot(numrotcurv-1,1))**
     $2.d0
      ome2(numrotcurv-1)=(curvrot(numrotcurv-1,2)/curvrot(numrotcurv-2,1
     &))**2.d0


      dome2(1)=(-49.d0/20.d0*ome2(1)+6.0d0*ome2(2)-15.d0/2.d0*ome2(3)+
     +20.d0/3.d0*ome2(4)-15.d0/4.d0*ome2(5)+6.d0/5.d0*ome2(6)-1.d0/6.d0*
     +ome2(7))/(curvrot(i,1)-curvrot(i+1,1))
      dome2(2)=(-49.d0/20.d0*ome2(2)+6.0d0*ome2(3)-15.d0/2.d0*ome2(4)+
     +20.d0/3.d0*ome2(5)-15.d0/4.d0*ome2(6)+6.d0/5.d0*ome2(7)-1.d0/6.d0*
     +ome2(8))/(curvrot(i,1)-curvrot(i+1,1))
      dome2(3)=(-49.d0/20.d0*ome2(3)+6.0d0*ome2(4)-15.d0/2.d0*ome2(5)+
     +20.d0/3.d0*ome2(6)-15.d0/4.d0*ome2(7)+6.d0/5.d0*ome2(8)-1.d0/6.d0*
     +ome2(9))/(curvrot(i,1)-curvrot(i+1,1))
      dome2(numrotcurv)=(49.d0/20.d0*ome2(numrotcurv)-6.0d0*
     +ome2(numrotcurv-1)+15.d0/2.d0*ome2(numrotcurv-2)-20.d0/3.d0*
     +ome2(numrotcurv-3)+15.d0/4.d0*ome2(numrotcurv-4)-6.d0/5.d0*
     +ome2(numrotcurv-5)+1.d0/6.d0*ome2(numrotcurv-6))/(curvrot(i,1)
     +-curvrot(i+1,1))
      dome2(numrotcurv-1)=(49.d0/20.d0*ome2(numrotcurv-1)-6.0d0*
     +ome2(numrotcurv-2)+15.d0/2.d0*ome2(numrotcurv-3)-20.d0/3.d0*
     +ome2(numrotcurv-4)+15.d0/4.d0*ome2(numrotcurv-5)-6.d0/5.d0*
     +ome2(numrotcurv-6)+1.d0/6.d0*ome2(numrotcurv-7))/(curvrot(i,1)
     +-curvrot(i+1,1))
      dome2(numrotcurv-2)=(49.d0/20.d0*ome2(numrotcurv-2)-6.0d0*
     +ome2(numrotcurv-3)+15.d0/2.d0*ome2(numrotcurv-4)-20.d0/3.d0*
     +ome2(numrotcurv-5)+15.d0/4.d0*ome2(numrotcurv-6)-6.d0/5.d0*
     +ome2(numrotcurv-7)+1.d0/6.d0*ome2(numrotcurv-8))/(curvrot(i,1)
     +-curvrot(i+1,1))
      dome2(numrotcurv-3)=(49.d0/20.d0*ome2(numrotcurv-3)-6.0d0*
     +ome2(numrotcurv-4)+15.d0/2.d0*ome2(numrotcurv-5)-20.d0/3.d0*
     +ome2(numrotcurv-6)+15.d0/4.d0*ome2(numrotcurv-7)-6.d0/5.d0*
     +ome2(numrotcurv-8)+1.d0/6.d0*ome2(numrotcurv-9))/(curvrot(i,1)
     +-curvrot(i+1,1))

c      open(123,file='kappaprova.out',
c     +status='unknown')

        do i=1,numrotcurv
           kappa(i)=dsqrt(dabs(curvrot(i,1)*dome2(i)+4.d0*ome2(i)))
c       write(123,*)curvrot(i,1),dsqrt(ome2(i))/61.009465d0,
c     +kappa(i)/61.009465d0
       enddo

       nstar=ndisk

         do i=1,ndisk
            pos(i,1)=x(i)
            pos(i,2)=y(i)
            pos(i,3)=z(i)
            vel(i,1)=vxx(i)
            vel(i,2)=vyy(i)
            vel(i,3)=vzz(i)
            mass(i)=massp(i)
         enddo


c Un cop llegides les dades per una de les captures, escollim la regió a
c analitzar, el radi fins on analitzarem, el radi mínim, el pas (en kpc), i la
c amplada de la regió (en cas que sigui una regió cuadrada, la aresta

        pas=steparm
        j=0
        do r=rmin,rmax-pas,pas
           j=j+1
        rdsota=-999.d0
        rdsobre=999.d0

           do l=1,numrotcurv
              rdif=curvrot(l,1)-r-pas/2.d0
              If(rdif.lt.0.) then
                if (rdif.gt.rdsota) then
                   rdsota=rdif
                   kappast=kappa(l)
                   ome2t=ome2(l)
                   rsota=curvrot(l,1)
                endif
              else
               if (rdif.gt.0.) then
                  if(rdif.lt.rdsobre) then
                    rdsobre=rdif
                    kappasb=kappa(l)
                    ome2b=ome2(l)
                    rsobre=curvrot(l,1)
                  endif
                else
                if(rdif.eq.0.) then
                   kappaq(j)=kappa(l)
                   ome2q(j)=ome2(l)
                   goto 11
                endif
               endif
             endif
            enddo
           kappaq(j)=kappasb*abs(rdsota/(rsobre-rsota))+kappast*abs
     &(rdsobre/(rsobre-rsota))
           ome2q(j)=ome2b*abs(rdsota/(rsobre-rsota))+ome2t*abs
     &(rdsobre/(rsobre-rsota))

 11      continue

c escollim analitzar només centrant-nos al pla

         h=0.0d0

c inicialitzem les variables que farem servir per calcular lv
        enddo
        do j=1,nbmax+1
           ro(j)=0
           nmean(j)=0
           do i=1,2
              xmean(j,i)=0.d0
           enddo
           do i=1,3
              sigma(j,i)=0.d0
           enddo
        enddo
        j=0

        do r=rmin,rmax-pas,pas
        j=j+1 
c comencem a calcular la lv i la densitat
c        write(*,*)r  
        nskew=0
 23        continue
        
        do kk=1,ndisk

         X2=pos(kk,1)
         Y2=pos(kk,2)
         Z2=pos(kk,3)

         VX2=vel(kk,1)
         VY2=vel(kk,2)
         VZ2=vel(kk,3)
         dist1=dsqrt(X2**2.d0+Y2**2.d0)

          if(dist1.le.(r+pas).and.dist1.gt.r)then
         
              UU=-1.d0/(dsqrt(x2**2.d0+y2**2.d0))*(vx2*x2+vy2*y2)
              VV=1.d0/(dsqrt(x2**2.d0+y2**2.d0))*(vy2*x2-vx2*y2)
            if(nskew.eq.0) then

              nmean(j)=nmean(j)+1
              xmean(j,1)=xmean(j,1)+UU
              xmean(j,2)=xmean(j,2)+VV
              xmean(j,3)=xmean(j,3)+vz2
              ro(j)=ro(j)+mass(kk)
                        
            else

              sigma(j,1)=sigma(j,1)+(UU-xmean(j,1))**2.d0
              sigma(j,2)=sigma(j,2)+(VV-xmean(j,2))**2.d0
              sigma(j,3)=sigma(j,3)+(vz2-xmean(j,3))**2.d0

       
             endif
          endif
         enddo
          if (nskew.eq.0) then
              if(nmean(j).le.10) then 
               nmean(j)=1
               xmean(j,1)=-9999
               xmean(j,2)=-9999
               xmean(j,3)=-9999
              endif
              xmean(j,1)=xmean(j,1)/dble(nmean(j))
              xmean(j,2)=xmean(j,2)/dble(nmean(j))
              xmean(j,3)=xmean(j,3)/dble(nmean(j))
              nskew=1
              goto 23
            else
              if(nmean(j).le.10) then
                nmean(j)=1
                sigma(j,1)=0.0d0
                sigma(j,2)=0.0d0
                sigma(j,3)=0.0d0
              endif
              sigma(j,1)=sigma(j,1)/dble(nmean(j))
              sigma(j,2)=sigma(j,2)/dble(nmean(j))
              sigma(j,3)=sigma(j,3)/dble(nmean(j))

          endif
c          write(*,*)xmean(j,1),xmean(j,2),xmean(j,3) 
          ro(j)=ro(j)/(pi*(r+pas)**2.0d0-pi*r**2.0d0)
          if(ro(j).eq.0.) ro(j)=1.0
      Toomre(j)=dsqrt(sigma(j,1))*kappaq(j)/(3.36*6.67e-11*ro(j)/1.e6)  
          toomre(j)=toomre(j)*3.08568015/1.9891*1e-11
      Xi(j)=kappaq(j)**2.d0*(r+pas/2.0d0)/(2.0d0*pi*6.67e-11*ro(j)/1.e6)
      mpar=2
      Xi(j)=Xi(j)*3.08568015/1.9891*1e-11/dble(mpar)
        enddo
        open(3,file=path//'diskpar'//sntimeH//'.out',status='unknown')
        write(3,153)
        l=0
        do i=1,j
       write(3,152)i*pas-pas/2.d0,dsqrt(sigma(i,1)),dsqrt(sigma(i,2)
     +),dsqrt(sigma(i,3)),kappaq(i)/61.009465,dsqrt(ome2q(i))/61.009465
     +,ro(i)/1.e6,toomre(i),Xi(i)
c           totaltot=totaltot+nmean(i)
c           write(*,*)totaltot

         If((i*pas-pas/2.0d0).le.trmax.and.(i*pas-pas/2.0d0).ge.trmin)
     + then

           sigUrmean(n)=sigUrmean(n)+dsqrt(sigma(i,1))
           sigVrmean(n)=sigVrmean(n)+dsqrt(sigma(i,2))
           sigWrmean(n)=sigWrmean(n)+dsqrt(sigma(i,3))
           kappaqrmean(n)=kappaqrmean(n)+kappaq(i)/61.009465
           ome2qrmean(n)=ome2qrmean(n)+dsqrt(ome2q(i))/61.009465
           rormean(n)=rormean(n)+ro(i)/1.e6
           toomrermean(n)=toomrermean(n)+toomre(i)
           xirmean(n)=xirmean(n)+xi(i)
           l=l+1

         endif
        enddo

           sigUrmean(n)=sigUrmean(n)/dble(l)
           sigVrmean(n)=sigVrmean(n)/dble(l)
           sigWrmean(n)=sigWrmean(n)/dble(l)
           kappaqrmean(n)=kappaqrmean(n)/dble(l)
           ome2qrmean(n)=ome2qrmean(n)/dble(l)
           rormean(n)=rormean(n)/dble(l)
           toomrermean(n)=toomrermean(n)/dble(l)
           xirmean(n)=xirmean(n)/dble(l)   
                
        close(2)
        close(1)    
        close(112)   


        If(n.eq.numstep-1) then

        write(*,*)'Starting disk analyisis, function of time'

        open(2,file=path//'diskpar_t.out',status='unknown')
        do l=0,numstep-1
         write(2,152)aget(l),sigUrmean(l),sigVrmean(l),sigWrmean(l),
     +kappaqrmean(l),ome2qrmean(l),rormean(l),toomrermean(l),xirmean(l)
        enddo
        close(2)
        endif

        write(*,*)'Back to main program from disk analysis'

151     FORMAT(8(F12.5,1x))
152     format(9(f18.5,1x))
153     format('# Rad (kpc) ','sigmaR (km/s)','sigmaPHI ',
     +'sigmaZ ','Kappa ','Ome ','Dsup (msun/pc^2) ','Q ','X')
 234    format(7(f15.5,1x))
    
        return 
        end
