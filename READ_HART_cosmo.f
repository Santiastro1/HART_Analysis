C  Subroutine Read_HART_cosmo.f
C                 21-03-2014
C   S.Roca-Fàbrega Universitat de Barcelona
C__________________________________________________________________
C This subroutine reads all particles and mesh properties (star,dark matter and gas)
C from the HART binary output files.
C The subroutine also finds the main galaxy (findcent.f), centers all coordinates and
C velocities (velcent.f) and finds the disk plane (findplane.f)
C Subroutines getrow.f, trans.f, trans1.f and trans2.f are needed for
C the computations 
C
C          Scales internal PM coordinates and velocities 
C          to  coordinates in Mpc/h and km/s
C          In PM the coordinates are in the range 1 - (NGRID+1)
C                     velocities are P = a_expansion*V_pec/(x_0H_0)
C                     where     x_0 = comoving cell_size=Box/Ngrid
C                                    H_0 = Hubble at z=0
C                   
C	     NROW = number of particles in 1D
C	     NGRID= number of cells        in 1D
C...........................................................................
C ---------------modifications---------------------------------------------
C ---------S. Roca-Fabrega 13-05-2014---------------------------------
C Add the possibility of not reading all DM particles, only 1st and 2nd
C species if the simulation file is too big for computer internal memory
c----------S. Roca-Fabrega 28-05-2014----------------------------------
C We commented line 109 to avoid redundance in parameter nbyteword
c----------S. Roca-Fabrega 1-09-2016-----------------------------------
c Add a counter for the number of analysed snapshots. It is usefull to
c use the center found in the previous snapshot as a first guess for the
c galactic center computation of the new one.

       Subroutine Read_HART_cosmo(n,sntimer,npart,x2,y2,z2,vx2,vy2,vz2,
     +mass2,zII,zIa,rhoh,thcell,tbirth,zshift,cell_len)

       include 'parameters.h'

       real*8 a2b,b2a,age1,f_b2a,f_a2b,aexpnpro,coord,vcoord,x,y,z,massa
       real*8 massp,mass2,rhoh,Cell_len
c     GAS parameters

      common / RUNPARAM / boxh1
      integer n
      real Thcell(nchmax),mhcell(nchmax)
      Real methcellsnII(nchmax),zII(npartm+nchmax),methcellsnIa(nchmax)
     &,zIa(npartm+nchmax)
      dimension coord(3,nchmax),vcoord(3,nchmax),cell_len(nchmax)
      dimension x(npartm),y(npartm),z(npartm),massp(npartm)
      dimension vxx(npartm),vyy(npartm),vzz(npartm),npart(5)
     &,lspecies1(10),x1(npartm),y1(npartm),z1(npartm),
     &vx1(npartm),vy1(npartm),vz1(npartm),mass2(npartm+nchmax)
      dimension massa(10),wpar(10),lsp(10),nsp(npartm,2),pw(npartm)
     &,pw0(npartm),tbirth(npartm),zstII(npartm),zstIa(npartm),
     +zshift(npartm),rhoh(nchmax)
      dimension x2(npartm+nchmax),y2(npartm+nchmax),z2(npartm+nchmax)
      dimension vx2(npartm+nchmax),vy2(npartm+nchmax),vz2(npartm+nchmax)

       Character  FileASCII*50
       character sntime*5  !para simulaciones con HART
       character tail*4
       character tail1*4
       character head1*6
       character head2*7
       character head4*7


      Call trans3(sntime,sntimer)

      filegas=path1//gashead//sntime//'.d'     
      tail='.DAT'
      tail1='.dat'
      head1='PMcrda'
      head2='PMcrs0a'
      head4='stars_a'



      call Read_Gas_File(filegas,coord,vcoord,Thcell,mhcell,rhoh,ncellh,
     &boxh1,methcellsnII,methcellsnIa,cell_len)

cc      do i=1,ncellh
cc       if(rhoh(i).eq.0.0) write(*,*)'ERROR! rho=0'
cc      enddo

      Box  =dble(boxh1) !(box size Mpc/h)                                                                                                               
      scrip=1 !cada cuantas particulas se van a escribir en el archivo de salida tipsy
      xmax =-1.e+9
      xmin = 1.e+9
      vmax =-1.e+9
      vmin = 1.e+9


C...................................................................
C			Read data and open files
      open(3,file=path1//head1//sntime//tail,form='unformatted',
     &status='unknown')
      read      (3,err=10,end=10) HEADER,
     &              AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW,
     &              TINTG,EKIN,EKIN1,EKIN2,AU0,AEU0,
     &              NROWC,NGRIDC,Nspecies,Nseed,Om0,Oml0,hubble,Wp5
     &              ,Ocurv,extras

 
c      ScaleV = BoxV/AEXPN/NGRID  ! scale factor for Velocities                                                                            \
c      ScaleV=Box*50.0d0/NGRID*sqrt(Om0)/AEXPN                       
c      ScaleC = Box*AEXPN/NGRID         ! scale factor for Coordinates

      ScaleT = 3.0856e19/(50.*hubble*sqrt(Om0)*3.15e7)!Scale factor fortime      
      ScaleC = Box/dble(NGRID)/dble(hubble)         ! scale factor for Coordinates
      ScaleV = Scalec/(ScaleT*3.15e7/3.0856e19)  ! scale factor for Velocities
      ScaleV=ScaleV/aexpn  !Correct by scale factor (velocity is independent of RS)

 
      close(3)
 10   continue
      nbyte  = nrecl * 4
c      nbytewordh=4!  defines length of direct-access record, 4 para HART, 1 para el resto
      nacces = nbyte / nbyteword
      xn=float(ngrid)+1.-1.E-7
      yn=float(ngrid)
c      nrowc=1024
      aexpn=dble(aexpn)
      open ( 1 , file =path1//head2//sntime//tail, access = 'direct',
     &           status = 'unknown', recl = nacces      )

         N_particles =lspecies(Nspecies+1)   ! Total number of particles                                                                    
         do i=1,nspecies
          write(*,*)'Particulas de la especie',i,' = ',lspecies(i+1),'we
     &eight= ',wspecies(i+1)
         enddo
         do i=1,nspecies
         massa(i)=wspecies(i+1)*(box/ngrid)**3*Om0/hubble/(3.64e-12)
         enddo
         write(*,*)'La especie',nspecies,' es la estelar',lspecies
     &(nspecies+1)-lspecies(nspecies)
         
         Npages = (N_particles -1)/NPAGE +1
         N_in_last=N_particles -NPAGE*(Npages-1)
      write (*,*) ' Nparticles=',N_particles

      if (optbigfile.eq.1) then
       write(*,*)'Big file option selected, only particles of DMsp 1 
     +will be readed.'
      endif

      DO  IROW=1, Npages         ! Loop over particle pages                                                                                
            In_page =NPAGE
            If(IROW.eq.Npages)In_page =N_in_last
            iL = NPAGE*(IROW-1)
         CALL GETROW(IROW,1) ! read in a page of particles                                                                                 
         DO  IN=1, In_page          ! Loop over particles                                                                                  
                ip =IN+iL                     ! current particle number                                                                    
                  WPAR =iWeight(ip)   ! particles weight                                                                                  
c check rounding errors (REAL*8 <=> REAL*4)                                                                                                


c If the files are too big we can read only part of DM particles
                if (optbigfile.eq.1) then
        if((ip.gt.lspecies(2)).and.(ip.le.lspecies(nspecies))) then
                   go to 99 
                  else
                  if (ip.gt.lspecies(nspecies)) then
                   ip=ip-lspecies(nspecies)+lspecies(2)
                  endif
                 endif
                else
                  continue
                endif
                x(ip) = dble(xpar(in))
                IF(x(ip).LT.1.) x(ip) = x(ip) + yn
                IF(x(ip).GT.xn) x(ip) = x(ip) - yn
                y(ip) = dble(ypar(in))
                IF(y(ip).LT.1.) y(ip) = y(ip) + yn
                IF(y(ip).GT.xn) y(ip) = y(ip) - yn
                z(ip) = dble(zpar(in))
                IF(z(ip).LT.1.) z(ip) = z(ip) + yn
                IF(z(ip).GT.xn) z(ip) = z(ip) - yn
             If(  x(ip).ge.xn)Then
                x(ip) = x(ip) -1.d-4
                write (*,*) x(ip),ip,1
             endif
             If(  y(ip).ge.xn)Then
                y(ip) = y(ip) -1.d-4
                write (*,*) y(ip),ip,2
             endif
             If(  z(ip).ge.xn)Then
                z(ip) = z(ip) -1.d-4
                write (*,*) z(ip),ip,3
             endif
                vxx(ip) = vx(in)
                vyy(ip) = vy(in)
                vzz(ip) = vz(in)
           X(ip)  =dble(ScaleC)* (dble(X(Ip))-1.d0)         
           Y(ip)  =dble(ScaleC)* (dble(Y(Ip))-1.d0)
           Z(ip)  =dble(ScaleC)* (dble(Z(Ip))-1.d0)          
           Vxx(ip)=dble(ScaleV)* dble(VXx(Ip))
           Vyy(ip)=dble(ScaleV)* dble(VYy(Ip))
           Vzz(ip)=dble(ScaleV)* dble(VZz(Ip)) 
           xmax =MAX(xmax,X(ip),Y(ip),Z(ip))
           xmin =MIN(xmin,X(ip),Y(ip),Z(ip))
           vmax =MAX(vmax,Vxx(ip),Vyy(ip),Vzz(ip))
           vmin =MIN(vmin,Vxx(ip),Vyy(ip),Vzz(ip))
 99      continue
         Enddo
       Enddo

        if (optbigfile.eq.1) then
         lspecies(nspecies+1)=lspecies(2)+lspecies(nspecies+1)-
     +lspecies(nspecies)
         lspecies(nspecies)=lspecies(2)
        endif

      write(*,*)'Xmax,xmin= ',xmax,xmin
c Read stars information

      write(*,*)'Reading stars'

      open ( 60 ,file =path1//head4//sntime//tail1,
     &       form = 'unformatted',status = 'old' )


      read(60) tdum, adum
      read(60) nstars
      
      if ( nstars .ne. (lspecies(nspecies+1)-lspecies(nspecies))) then
        write(*,*) '1  : In Read_Stellar_Data: something is wrong:'
        write(*,*) 'nstars =',nstars,' is iconsistent with'
        write(*,*) 'lspecies(nspecies+1)-lspecies(nspecies)',
     &               lspecies(nspecies+1)-lspecies(nspecies)
        write(*,*) 'stopping...'
        stop
      endif
      if ( nstars .eq. 0 ) goto 1234
      read (60) ws_old, ws_oldi
      read (60) (pw(ic1),ic1=1,nstars) ! weights
      read (60) (pw0(ic1),ic1=1,nstars)      ! initial masses
      read (60) (tbirth(ic1),ic1=1,nstars)   ! birth times
c#ifdef ENRICH
      read (60) (zstII(ic1),ic1=1,nstars)    ! metallicity 
c#endif
c#ifdef ENRICH_SNIa
      read (60) (zstIa(ic1),ic1=1,nstars)    ! metallicity 
c#endif
      close( 60 )
c
 1234 continue
      
c-------------------------------------------------------------------------------------------------
c-----------------------------------------------------------------------------------------------
      
      hfact=hubble
      t0=2.0*3.0856e19/(3.15e7*100*hubble*sqrt(Om0))
       aexpnpro=dble(aexpn)


      write(*,*)'Defining star variables values'

      do i=lspecies(nspecies)+1,lspecies(nspecies+1)
     &,scrip
        j=i-lspecies(nspecies)
        pw(j)=pw(j)*(box/ngrid)**3*2.776e11*Om0/hubble
        pw0(j)=pw0(j)*(box/ngrid)**3*(2.776e11)*Om0/hubble
        massp(j)=pw(j)
        x1(j)=(x(i))*1000.
        y1(j)=(y(i))*1000.
        z1(j)=(z(i))*1000.
        vx1(j)=vxx(i)
        vy1(j)=vyy(i)
        vz1(j)=vzz(i)
        zstII(j)=zstII(j)*pw(j)
        zstIa(j)=zstIa(j)*pw(j)
        zshift(j)=1.0d0/b2a(tbirth(j)*1.d0)-1.0d0
        tbirth(j)=age1(a2b(aexpnpro)*1.d0)-age1(tbirth(j)*1.d0)
      enddo
    

      npart(1)=lspecies(nspecies+1)
      npart(2)=lspecies(nspecies+1)-lspecies(nspecies)
      npart(3)=0
      npart(4)=lspecies(nspecies)
      npart(5)=ncellh
      write(*,*)'Total number of particles= ',npart(1)
      write(*,*)'Stellar particles= ',npart(2)
      write(*,*)'Bulge particles= ',npart(3)
      write(*,*)'DM particles= ',npart(4)
      write(*,*)'Gas cells= ',ncellh 



c      do l=1,ncellh
c       if(rhoh(l).eq.0.0) then
c        write(*,*)'ERROR! rho=0'
c        write(*,*)'Inside READ_HART, before DMdef'
c        stop
c       endif
c      enddo


      write(*,*)'Defining DM variables values. N1sp=',lspecies(2)

c Dark Matter
      
      xcenter=0.
      ycenter=0.
      zcenter=0.

      tbirth(npart(2)+1)=age1(a2b(aexpnpro)*1.d0)
      do i=npart(2)+1,npart(1)
        massp(i)=massa(1)
      if (optbigfile.ne.1) then
        do j=1,nspecies-2
        If((i-npart(2)).gt.lspecies(j+1).and.((i-npart(2)).
     +le.lspecies(j+2))) then
        massp(i)=massa(j+1)
        endif
        enddo
      endif
        x1(i)=(x(i-npart(2)))*1000.
        y1(i)=(y(i-npart(2)))*1000.
        z1(i)=(z(i-npart(2)))*1000.
        vx1(i)=vxx(i-npart(2))
        vy1(i)=vyy(i-npart(2))
        vz1(i)=vzz(i-npart(2))
        tbirth(i)=tbirth(npart(2)+1)
        zstII(i)=0.0
        zstIa(i)=0.0

      enddo


c      do l=1,ncellh
c       if(rhoh(l).eq.0.0) then
c        write(*,*)'ERROR! rho=0'
c        write(*,*)'Inside READ_HART, before finding galcent'
c        stop
c       endif
c      enddo


c Find galactic center
 
      write(*,*)'Finding the galactic center and centering'

      call findcent(n,x1,y1,z1,npart,xcenter,ycenter,zcenter)
      
c      xcenter=16166.629
c      ycenter=2819.809
c      zcenter=23347.849

      do i=1,npart(1)
         x2(i)=x1(i)-xcenter
         y2(i)=y1(i)-ycenter
         z2(i)=z1(i)-zcenter
         mass2(i)=massp(i)
         zII(i)=zstII(i)/mass2(i)
         zIa(i)=zstIa(i)/mass2(i)
c      if (i.lt.npart(2)) write(*,*)x1(i),x2(i),y1(i),y2(i),z1(i),z2(i)
      enddo

      write(*,*)'Finding the galaxy velocity with respect to the box'

      call findvelc(x2,y2,z2,vx1,vy1,vz1,npart,vxcent,vycent,vzcent)

      do i=1,npart(1)
         vx2(i)=vx1(i)-vxcent
         vy2(i)=vy1(i)-vycent
         vz2(i)=vz1(i)-vzcent
      enddo

c Gas

      write(*,*)'Defining Gas variables values, ncell=',npart(5)

      do i=1,ncellh
        Cell_len(i)=Cell_len(i)*1000.
        x2(i+npart(1))=coord(1,i)*1000.-xcenter
        y2(i+npart(1))=coord(2,i)*1000.-ycenter
        z2(i+npart(1))=coord(3,i)*1000.-zcenter
        vx2(i+npart(1))=vcoord(1,i)-vxcent
        vy2(i+npart(1))=vcoord(2,i)-vycent
        vz2(i+npart(1))=vcoord(3,i)-vzcent
        Thcell(i)=Thcell(i)
        rhoh(i)=rhoh(i)
        mass2(i+npart(1))=mhcell(i)
        zII(i+npart(1))=methcellsnII(i)
        zIa(i+npart(1))=methcellsnIa(i)
c        write(*,*)x2(i+npart(1)),y2(i+npart(1)),z2(i+npart(1))
      enddo


c      do l=1,ncellh
c       if(rhoh(l).eq.0.0) then
c        write(*,*)'ERROR! rho=0'
c        write(*,*)'Inside READ_HART, before finding galplane'
c        stop
c       endif
c      enddo


c Find galactic plane

      write(*,*)'Finding the galactic plane and defining as x,y plane'

      call findplane(x2,y2,z2,vx2,vy2,vz2,npart,tbirth,thcell)


c      do l=1,ncellh
c       if(rhoh(l).eq.0.0) then
c        write(*,*)'ERROR! rho=0'
c        write(*,*)'Inside READ_HART, before return'
c        stop
c       endif
c      enddo


      return
      END


c     ---------------------------
      real*8 function age1 ( td )
c     ---------------------------
c
c     returns age of the Universe at td (time in code units)
c
c     uses Oleg's formula for flat LCDM 
c
c

      real*8 td, b2a
      INCLUDE 'parameters.h'

      dOm0=dble(Om0)
      dOml0=dble(Oml0)
c      write(*,*)Om0,Oml0
      as = b2a ( td ) ! convert code time to expansion factor
      ff = dOm0/(1.-dOm0)/as**3
c
c.... calculate age of the universe in Gyrs 
c
      age1 = 9.779/hubble * 2./3./
     &      sqrt(1.-dOm0)*log((1.+sqrt(1.+ff))/sqrt(ff))

      return
      end
c     --------------------------
      real*8 function b2a ( bt )
c     --------------------------
c
      real*8 bt
      real*8 f_b2a, fp(1)
      real zbrent
      external zbrent, f_b2a
      include 'parameters.h'
c
c
      dOm0=dble(Om0)
      dOml0=dble(Oml0)
      IF ( (dOm0 .eq. 1.0) .and. (dOml0 .eq. 0.0) ) THEN
        b2a = (1.d0 / (1.d0 - bt))**2
      ELSE
        if ( (bt .lt. -190.d0) .or. (bt .gt. 1.d-1) ) then
          write(*,*) ' b2a : something is wrong : bt =',bt
          stop
        endif
        if ( (dOm0 .eq. 0.0) .and. ( dOml0 .eq. 0.0 ) ) then
          write(*,*) 'b2a : something is wrong : Om0 = Oml0 = 0.0'
          stop
        endif
        fp(1) = bt
        b2a = zbrent ( f_b2a , fp , 1 , 1.d-4 , 1.1d0 , 1.d-9 )
      ENDIF
c
      return
      end
c
c     --------------------------------------
      real*8 function f_b2a ( at , fp , np )
c     --------------------------------------
c
c     at - expansion factor
c     fp(1) = b
c
      integer np
      real*8 at , fp(np)
      real*8 a2b
      external a2b
c
      f_b2a = a2b ( at ) - fp(1)
c      
      return
      end
c
c     -------------------------------------
      FUNCTION zbrent(func,fp,np,x1,x2,tol)
c     -------------------------------------
c
      INTEGER ITMAX
      REAL zbrent
      integer np
      real*8 fp(np)
      real*8 tol,x1,x2,EPS
      real*8 func
      EXTERNAL func
      PARAMETER (ITMAX=100,EPS=3.e-8)
      INTEGER iter
      REAL*8 a,b
      REAL c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
      a=x1
      b=x2
      fa=func(a,fp,np)
      fb=func(b,fp,np)
      if((fa.gt.0..and.fb.gt.0.).or.(fa.lt.0..and.fb.lt.0.)) then
        write(*,*) 'root must be bracketed for zbrent'
        write(*,*) 'fa =',fa,'  fb =', fb
        pause
      endif
      c=b
      fc=fb
      do 11 iter=1,ITMAX
        if((fb.gt.0..and.fc.gt.0.).or.(fb.lt.0..and.fc.lt.0.))then
          c=a
          fc=fa
          d=b-a
          e=d
        endif
        if(abs(fc).lt.abs(fb)) then
          a=b
          b=c
          c=a
          fa=fb
          fb=fc
          fc=fa
        endif
        tol1=2.*EPS*abs(b)+0.5*tol
        xm=.5*(c-b)
        if(abs(xm).le.tol1 .or. fb.eq.0.)then
          zbrent=b
          return
        endif
        if(abs(e).ge.tol1 .and. abs(fa).gt.abs(fb)) then
          s=fb/fa
          if(a.eq.c) then
            p=2.*xm*s
            q=1.-s
          else
            q=fa/fc
            r=fb/fc
            p=s*(2.*xm*q*(q-r)-(b-a)*(r-1.))
            q=(q-1.)*(r-1.)*(s-1.)
          endif
          if(p.gt.0.) q=-q
          p=abs(p)
          if(2.*p .lt. min(3.*xm*q-abs(tol1*q),abs(e*q))) then
            e=d
            d=p/q
          else
            d=xm
            e=d
          endif
        else
          d=xm
          e=d
        endif
        a=b
        fa=fb
        if(abs(d) .gt. tol1) then
          b=b+d
        else
          b=b+sign(tol1,xm)
        endif
        fb=func(b,fp,np)
11    continue
      pause 'zbrent exceeding maximum iterations'
      zbrent=b
      return
      END

c     --------------------------
      real*8 function a2b ( at )
c     --------------------------
c
c     translates expansion factor at into hydro time variable b
c
      real*8 at
      real*8 fp(2), atst
      real*8 INTEGRATE , f_a2b
      external INTEGRATE , f_a2b
      include 'parameters.h' 
c      
c
      dOm0=dble(Om0)
      dOml0=dble(Oml0)
      if ( (at .lt. 0.0) ) then
        write(*,*) 'a2b : something is wrong : at = ', at
        stop
      endif
      IF ( (dOm0 .eq. 1.0) .and. (dOml0 .eq. 1.0) ) THEN
        a2b = (1.d0 - 1.d0/sqrt(at))
      ELSE
        if ( (dOm0 .eq. 0.0) .and. ( dOml0 .eq. 0.0 ) ) then
          write(*,*) 'a2b : something is wrong : Om0 = Oml0 = 0.0'
          stop
        endif
        fp(1) = dOm0
        fp(2) = dOml0
        atst  = 1.d-1 * (1.d0 - at)
        a2b = INTEGRATE( f_a2b , fp , 2 , 1.d0 , at , atst , 1.d-9 )
      ENDIF
c
      return
      end




c     --------------------------------------------------------------
      double precision function INTEGRATE(FUNC,fp,np,a,b,dxinit,eps)
c     --------------------------------------------------------------
c
c     Quadrature using fifth order Runge-Kutta with adaptive step size.
c     Based on Press et al, Numerical Recipes in C, 2nd ed, pp 719-722.
c
c     Runge-Kutta driver with adaptive stepsize control.  Integrate starting
c     value y from a to b with accuracy eps, storing intermediate results in
c     global variables.  dxinit should be set as a guessed first stepsize.
c
c     Pass a vector of parameters of length np to FUNC in fp(np).
c
c     Copyright (c) 1997 Michael A. K. Gross.  You may use this program for
c     personal, educational or research purposes.  Commercial purposes require
c     special arrangements. If you publish a paper that depends upon this code,
c     please cite it appropriately.
c
c     Questions and/or comments may be sent to gross@fozzie.gsfc.nasa.gov.
c
c     slight modification by A.Kravtsov to allow to pass a vector 
c     of parameters to FUNC

      implicit none
      integer np
      double precision a, b, eps, dxinit, FUNC, fp(np)
      external FUNC

      integer maxsteps
      parameter(maxsteps=100000000)

      double precision x, dx, dxnext, y, dydx, yscale
      integer  Nstep

      x     = a
      dx    = dxinit
      y     = 0.d0
      Nstep = 0

      do while ((x-b)*(b-a).lt.0.d0.and.Nstep.lt.maxsteps)
        Nstep = Nstep + 1
        dydx = FUNC(x,fp,np)
c
c       yscale is the scaling used to monitor accuracy.  This general-purpose
c       choice can be modified if need be.
c
        yscale = max(abs(y) + abs(dx*dydx), 1.d-20)
        if ((x+dx-b)*(x+dx-a).gt.0.d0)  ! If stepsize overshoots, decrease it.
     1    dx = b - x

        call RUNGE5VAR(y,dydx,x,dx,eps,yscale,dxnext,FUNC,fp,np)

        dx = dxnext
      end do

      if (Nstep.ge.maxsteps)
     1  write (*,*) 'WARNING: failed to converge in INTEGRATE.'

      INTEGRATE = y

      return
      end
c
c     -------------------------------------------------------------
      SUBROUTINE RUNGE5VAR(y,dydx,x,htry,eps,yscale,hnext,DERIVS,
     1                     fp,np)
c     -------------------------------------------------------------
c
c     Fifth-order Runge-Kutta step with monitoring of local truncation error
c     to ensure accuracy and adjust stepsize.  Input are the dependent
c     variable y and its derivative dydx at the starting value of the
c     independent variable x.  Also input are the stepsize to be attempted
c     htry, the required accuracy eps, and the value yscale, against which the
c     error is scaled.  On output, y and x are replaced by their new values.
c     hdid is the stepsize that was actually accomplished, and hnext is the
c     estimated next stepsize.  DERIVS is the user-supplied routine that
c     computes right-hand-side derivatives.  The argument fp is a vector 
c     of parameters (np parameters) to be passed to DERIVS 
c     (NOT integrated over).
c
c
c     Copyright (c) 1997 Michael A. K. Gross.  You may use this program for
c     personal, educational or research purposes.  Commercial purposes require
c     special arrangements. If you publish a paper that depends upon this code,
c     please cite it appropriately.
c
c     Questions and/or comments may be sent to gross@fozzie.gsfc.nasa.gov.
c
c     slight modification by A.Kravtsov to allow to pass a vector 
c     of parameters to FUNC
c
      implicit none
      integer np
      double precision eps,hnext,htry,x,dydx,y,yscale,DERIVS,fp(np)
      external DERIVS

      double precision errmax,h,hold,htemp,xnew,yerr,ytemp

      double precision safety,pgrow,pshrink,errcon
      parameter (safety  =  0.9d0)
      parameter (pgrow   = -0.2d0)
      parameter (pshrink = -0.25d0)
      parameter (errcon  =  1.89d-4)

      h = htry                         ! Set stepsize to initial accuracy.
      errmax = 10.d0
      do while (errmax.gt.1.d0)
        call RUNGE(y,dydx,x,h,ytemp,yerr,DERIVS,fp,np)

        errmax = abs(yerr/yscale)/eps   ! Scale relative to required accuracy.
        if (errmax.gt.1.d0) then        ! Truncation error too large; reduce h
          htemp = safety*h*(errmax**pshrink)
          hold = h
          h = sign(max(abs(htemp),0.1d0*abs(h)),h)  ! No more than factor of 10
          xnew = x + h
          if (xnew.eq.x) then
            write (*,*) 'WARNING: ',
     1                  'Stepsize underflow in RUNGE5VAR().'
            h = hold
            errmax = 0.d0
          end if
        end if
      end do
c
c     Step succeeded.  Compute estimated size of next step.
c
      if (errmax.gt.errcon) then
        hnext = safety*h*(errmax**pgrow)
      else
        hnext = 5.d0 * h                ! No more than factor of 5 increase.
      end if
      x = x + h

      y = ytemp

      return
      end


c     ---------------------------------------------------
      SUBROUTINE RUNGE(y,dydx,x,h,yout,yerr,DERIVS,fp,np)
c     ---------------------------------------------------
c
c     Given values for a variable y and its derivative dydx known at x, use
c     the fifth-order Cash-Karp Runge-Kutta method to advance the solution
c     over an interval h and return the incremented variables as yout.  Also
c     return an estimate of the local truncation error in yout using the
c     embedded fourth order method.  The user supplies the routine
c     DERIVS(x,y,dydx), which returns derivatives dydx at x.
c
c     Copyright (c) 1997 Michael A. K. Gross.  You may use this program for
c     personal, educational or research purposes.  Commercial purposes require
c     special arrangements. If you publish a paper that depends upon this code,
c     please cite it appropriately.
c
c     Questions and/or comments may be sent to gross@fozzie.gsfc.nasa.gov.
c
c     slight modification by A.Kravtsov to allow to pass a vector 
c     of parameters to FUNC
c
      implicit none

      integer np
      double precision h,x,dydx,y,yerr,yout,DERIVS,fp(np)

      external DERIVS

      double precision ak3, ak4, ak5 ,ak6

      double precision a2,a3,a4,a5,a6
      double precision c1,c3,c4,c6,dc1,dc3,dc4,dc5,dc6
      parameter(a2  =    0.2d0)
      parameter(a3  =    0.3d0)
      parameter(a4  =    0.6d0)
      parameter(a5  =    1.d0)
      parameter(a6  =    0.875d0)
      parameter(c1  =   37.d0/378.d0)
      parameter(c3  =  250.d0/621.d0)
      parameter(c4  =  125.d0/594.d0)
      parameter(c6  =  512.d0/1771.d0)
      parameter(dc1 = c1 -  2825.d0/27648.d0)
      parameter(dc3 = c3 - 18575.d0/48384.d0)
      parameter(dc4 = c4 - 13525.d0/55296.d0)
      parameter(dc5 = -277.d0/14336.d0)
      parameter(dc6 = c6 -     0.25d0)
c      write(*,*)fp(1),fp(2)
      ak3 = DERIVS(x+a3*h,fp,np)
      ak4 = DERIVS(x+a4*h,fp,np)
      ak5 = DERIVS(x+a5*h,fp,np)
      ak6 = DERIVS(x+a6*h,fp,np)
c
c     Estimate the fifth order value.
c
      yout = y + h*(c1*dydx + c3*ak3 + c4*ak4  + c6*ak6)
c
c     Estimate error as difference between fourth and fifth order
c
      yerr = h*(dc1*dydx + dc3*ak3 + dc4*ak4 + dc5*ak5 + dc6*ak6)

      return
      end

c     -------------------------------------
      real*8 function f_a2b ( x , fp , np )
c     -------------------------------------
c
c     input : x     - is expansion factor variable
c             fp(1) - Om0  = present-day matter density 
c             fp(2) - Oml0 = present-day vacuum contribution
c
      integer np
      real*8 x , fp(np), d1
c
      d1 = x**3
      f_a2b = .5d0 * sqrt(fp(1)) / d1
     &             / sqrt(fp(1)/d1 + fp(2) + (1.d0-fp(1)-fp(2))/x**2)
c
      return
      end

c    -----------------------------------------------------------------
      subroutine Read_Gas_File(filegas,coord,vcoord,Thcell,mhcell,rhoh,
     &numcellh,boxh1,methcellsnII,methcellsnIa,Cell_len)
c    ------------------------------------------------------------------

      include 'parameters.h'
      real*8 dt,coord,vcoord,rhoh,cell_len
      Real methcellsnIa(nchmax),methcellsnII(nchmax)
      Real Thcell(nchmax),mhcell(nchmax)
      dimension coord(3,nchmax),vcoord(3,nchmax),rhoh(nchmax)
      dimension cell_len(nchmax)

c      Write(6,'(''Enter (x,v) halo coordinates and Rvir(kpc/h), fact'',$
c     &)')
c      Read(5,*) xc,yc,zc,vxc,vyc,vzc,Rvir,fact

      Do i=1, nchmax  ! Initialize arrays
         Do j=1, 3
            coord(j,i) = 0.0d0
            vcoord(j,i)= 0.0d0
         End do
         cell_len(i)=0.0
         Thcell(i)= 0.0
         mhcell(i)= 0.0
         methcellsnIa(i)= 0.0
         methcellsnII(i)= 0.0
         rhoh(i)= 0.0
      End do

      open  ( 19 , file = filegas, form = 'unformatted' )
      read ( 19 ) jname

      read ( 19 ) istep , t , dt, aexpn, ainit
 20   format ('istep =',i4,1x,' t =',g15.8,1x,
     &        'dt =',g15.8,1x,' a =',g15.8,' ainit =',g15.8)

      do Level = MinLevel , MaxLevel
        aexp(Level) = aexpn
      enddo

c     if having trouble on linux, try inserting a dummy first in the read
c     list below

      read ( 19 ) boxh1, Om0, Oml0, Omb0, hubble
c
c.... rho0 = 3H0^2 * Om0 / (8*pi*G) - unit of density in Msun/Mpc^3
c
      rho0 = 2.776e11 * hubble**2 * Om0
      rhoscale= (rho0/aexpn**3)*1.e-18       ! [Msun/pc^3]
      aM0= rho0 * (boxh1/hubble)**3 / ncell0  ! [Msun]   
c
c.... T_0 = unit of temperature in K and in keV)

c.... v0 - velocity units in km/s
c
      r0= dble(boxh1)/dble(nrowg)
      v0 = 50. * r0 * sqrt(Om0)
      vscale= v0/aexpn

      T_0 = 3.03e5 * r0**2 * wmu * Om0 ! [K]
      Tscale= T_0/aexpn**2

      read ( 19 ) nextras
      if ( nextras .gt. nextra ) then
        write(*,*) '* error: Read_Gas_Binary : nextras > nextra :',
     &             nextras, nextra
        write(*,*) '* (check nextra parameter in a_control.h) '
        stop
      endif
      read ( 19 ) (extra(i), i=1,nextras)
      read ( 19 ) (lextra(i), i=1,nextras)

      read ( 19 ) MinLev , MaxLevelNow
      if ( MaxLevelNow .gt. MaxLevel ) then
         write(*,*)
     &        '* error: Read_Gas_Binary : MaxLevelNow > MaxLevel :',
     &        MaxLevelNow, MaxLevel
         stop
      endif

      read ( 19 ) ( tl(i)     , i = MinLevel , MaxLevelNow )
      read ( 19 ) ( dtl(i)    , i = MinLevel , MaxLevelNow )
      read ( 19 ) ( tlold(i)  , i = MinLevel , MaxLevelNow )
      read ( 19 ) ( dtlold(i) , i = MinLevel , MaxLevelNow )
      read ( 19 ) ( iSO(i)    , i = MinLevel , MaxLevelNow )
      read ( 19 ) ncell
      if ( ncell .ne. ncell0 ) then
         write(*,*)
     &        '* error: Read_Gas_Binary : ncell not equal ncell0 :',
     &        ncell, ncell0
         stop
      endif
      read ( 19 ) (iOctCh(i), i=1,ncell)
      read ( 19 ) ((hvar(ivar,i),ivar=1,nhvar), i=1,ncell)
      nleave= 0
      totmass= 0.
      Tcellmin= 1.e+10
      Rvir1= 1.e-3*Rvir1 ! Mpc/h
      Do i=1, ncell
         If(iOctCh(i) .eq. nil)then
            nleave= nleave + 1
            totmass= totmass + hvar(1,i)
            Tcell= hvar(6,i)/hvar(1,i)
            Tcellmin= min(Tcell,Tcellmin)
         End if
      End do
      Write(*,*) 'Number of leafs at zero level', nleave, 'totmass',
     & totmass

      read ( 19 ) ((var(ivar,i),ivar=2,3), i=1,ncell)
      ntot = 0
      if ( MaxLevelNow .gt. MinLevel ) then
         read ( 19 ) iOctFree , nOct
         write(*,*) 'iOctFree =',iOctFree,' nOct =', nOct
         if ( ncell+nchild*nOct .gt. mcell ) then
            write(*,*)
     &   '* error: Read_Gas_Binary : ncell+nchild*nOct > mcell :',
     &           n0ell+nchild*nOct, mcell
            write(*,*)
     &           '* (size of the input exceeds max. of cells (mcell))'
            write(*,*) '* increase mcell in a_setup.h and try again...'
            stop
         endif

         Tcellmax= 0.0  
         Tcellhmax= 0.0
         vcellhmax= 0.0
         numcellh= 0
         pmax= 0.0
         rhohax= 0.0
         rhomax= -1.0
         rhomin= 1.e+19
         cellhalomass= 0.0
         cellmin= 1.e+19
         cellmax= 0.0
         do Lev = MinLevel+1 , MaxLevelNow
            read ( 19 ) Level , iNOLL(Level) , iHOLL(Level)
            iOct   = iHOLL(Level)
            nLevel = iNOLL(Level)
            nLevCells = nLevel*nchild
            if ( nLevCells .gt. nclmax ) then
               write(*,*) 'error : L =',Lev,' nLevCell =',nLevCell,
     &              ' > nclmax =',nclmax,' set in a_setup.h'
               write(*,*) '=> increase nclmax and rerun.'
               close ( 19 )
               stop
            endif
            ntot = ntot + nLevel
            write(*,*) 'reading tree (oct) data for level ', Lev, '  
     &nlevel=', nlevel
            do ic1 = 1 , nLevel
               read(19) (iOctPs(i,iOct),i=1,3),(iOctNb(i,iOct),i=1,6),
     &              iOctPr(iOct), iOctLv(iOct), iOctLL1(iOct),
     &              iOctLL2(iOct)
               iOct = iOctLL1(iOct)
            enddo
            write(*,*) 'reading cell data: ncells =', nLevel*nchild

            do ic1 = 1 , nLevel*nchild
               CellVolume= 1.0 * 2.0**(-3.0*Level)
               Cellsize= 1.0 * 2.0**(-Level)
               read ( 19 ) idc, iOctCh(idc), (hvar(i,idc),i=1,nhvar),
     &              (var(i,idc), i=2,3)
               if ( iOctCh(idc) .eq. nil ) then
                  call Ps(idc, Posx,Posy,Posz ) ! finds position

                  xcell= dble(r0)*(dble(Posx) - 1.0d0)
                  ycell= dble(r0)*(dble(Posy) - 1.0d0)
                  zcell= dble(r0)*(dble(Posz) - 1.0d0)
                  Cellsizeg=Cellsize*dble(r0)
c                 If(idc .eq. 9858863)then
c                    Write(*,*) xcell, ycell, zcell
c                 End if

                  rx= dble(xcell)-dble(boxh1)/2.d0  ! Define x & v coordinates of halo
                  ry= dble(ycell)-dble(boxh1)/2.d0  ! with respect to cm of halo
                  rz= dble(zcell)-dble(boxh1)/2.d0
                  rhomax= max(hvar(1,idc),rhomax)

                  If(abs(rx) .gt. boxh1/2.)then  ! periodic conditions
                     If(rx .gt. 0.0)then
                        xcell=dble(xcell) - dble(boxh1)
                     Else
                        xcell=dble(xcell) + dble(boxh1)
                     End if
                  End if
                  If(abs(ry) .gt. boxh1/2.)then
                     If(ry .gt. 0.0)then
                        ycell=dble(ycell) - dble(boxh1)
                     Else
                        ycell=dble(ycell) + dble(boxh1)
                     End if
                  End if
                  If(abs(rz) .gt. boxh1/2.)then
                     If(rz .gt. 0.0)then
                        zcell=dble(zcell) - dble(boxh1)
                     Else
                        zcell=dble(zcell) + dble(boxh1)
                     End if
                  End if
c                  dis2= (xcell-xc)**2 + (ycell-yc)**2 + (zcell-zc)**2
c                  dis= sqrt(dis2)
                  cellmass = hvar(1,idc)*CellVolume
                  rhocell= hvar(1,idc)
                  Tcell = hvar(6,idc)/hvar(1,idc)
                  vxcell= hvar(3,idc)/hvar(1,idc)
                  vycell= hvar(4,idc)/hvar(1,idc)
                  vzcell= hvar(5,idc)/hvar(1,idc)
                  vcell2= vxcell**2 + vycell**2 + vzcell**2
                  vcell= sqrt(vcell2)
                  rhomin= min(hvar(1,idc),rhomin)
                  Tcellmax = max(Tcell,Tcellmax)
                  Tcellmin= min(Tcell,Tcellmin)
                  totmass = totmass + cellmass
                  pressure= hvar(6,idc)

c                 If(dis .le. fact*Rvir .and. Tcell*Tscale .lt. TGasmin)then
c                  If(dis .le. fact*Rvir)then

                     Cell_len(numcellh)=dble(Cellsizeg)/dble(hubble)
                     numcellh= numcellh + 1
                     coord(1,numcellh)= dble(xcell)/dble(hubble)
                     coord(2,numcellh)= dble(ycell)/dble(hubble)
                     coord(3,numcellh)= dble(zcell)/dble(hubble)

                     vcoord(1,numcellh)= vxcell*vscale
                     vcoord(2,numcellh)= vycell*vscale
                     vcoord(3,numcellh)= vzcell*vscale
           methcellsnIa(numcellh)=hvar(10,idc)/hvar(1,idc)
           methcellsnII(numcellh)=hvar(9,idc)/hvar(1,idc)
                     Thcell(numcellh)= Tcell*Tscale
                     mhcell(numcellh)= hvar(1,idc)*CellVolume*aM0
                     rhoh(numcellh)= hvar(1,idc)*rhoscale
                     rhohmax= max(rhocell,rhohmax)
                     rhohmin= min(rhocell,rhohmin)
                     cellhalomass= hvar(1,idc)*CellVolume + cellhalomass
                     cellmin= min(cellmass,cellmin)
                     cellmax= max(cellmass,cellmax)
                     Tcellhmax= max(Tcell,Tcellhmax)
                     vcellhmax= max(vcell,vcellhmax)
                     pmax= max(pressure,pmax)
c                     End if
                  Endif
            Enddo
         Enddo
      Endif
      close (19)
      write(*,*) 'done reading',ncell0+nchild*ntot,' cells...'
      tiempo= seconds()
      Write(*,*) 'tiempo(seg)=  ', tiempo

      write(*,100)'min & max density code units=', rhohmin, rhohmax
      Write(*,*) 'max pressure inside halo=', pmax
      Write(*,*) 'rhomax,min=', rhomax*rhoscale, rhomin*rhoscale
      Write(*,*) 'rho0=', rho0, 'totmass=', totmass, 'wmu=', wmu
      write(*,100)'totmass=', aM0*totmass
      write(*,100)'min & max cell mass=', aM0*cellmin, aM0*cellmax
      Write(*,101)'Tmax & Tmin (K)=',  Tcellmax*Tscale, Tcellmin*Tscale
      Write(*,*) 'ncell in halo=', numcellh
      Write(*,100)'Halo gas mass (Msun/h)=', aM0*cellhalomass*hubble
      Write(*,102)'T halo max [K]=', Tcellhmax*Tscale
      Write(*,102)'V halo max [km/s]=', vcellhmax*vscale
      Write(*,*) 'rho halo max [Msun/pc^3]', rhohmax*rhoscale
      Write(*,*) 'rho halo max [cm^-3]', rhohmax*rhoscale*30.37

 100  Format(1x,a24,2g12.3)
 101  Format(1x,a17,2f12.2)
 102  Format(1x,a17,f12.2)

      Return
      End

c     ------------------------
      real function seconds ()
c     ------------------------
c
c     purpose: returns elapsed time in seconds
c     uses   : xlf90 utility function timef()
c              (in this case - to be compiled with xlf90 or with -lxlf90)
c              or dtime (Exemplar) or etime (Power Challenge)
c      real*8 timef
      real*8 dummy
      real tarray(2)

      dummy = 0.0

c.... for IBM SP
      dummy   = timef ()
      seconds = dummy / 1000

c.... for exemplar
c      dummy = dtime(tarray) 
c      seconds = dtime(tarray)

c.... for hitachi
c      dummy = dtime(tarray)
c       seconds = second()

c.... for power challenge
      dummy = etime(tarray)
      seconds = etime(tarray)

      return
      end
c     -------------------------------------
      subroutine Ps ( iC , Posx,Posy,Posz )      ! finds position
c     -------------------------------------
c     purpose : returns coordinates of cell center
c     Input   : iC     - pointer to a cell 
c     Output  : Posxyz    - x,y,z positions
c
      include 'parameters.h'
      integer iC
      integer idelta(8,3)
      data      idelta / -1,  1, -1,  1, -1,  1, -1,  1,
     &                   -1, -1,  1,  1, -1, -1,  1,  1,
     &                   -1, -1, -1, -1,  1,  1,  1,  1 /
      if ( iC .le. ncell0 ) then
        iC_ = iC - 1
        i = iC_ / ng2               ! ng2 = ng**2 
        j = ( iC_ - i*ng2 ) / ng
        k =   iC_ - i*ng2 - j*ng
        Posx = dble(i) + 1.5d0
        Posy = dble(j) + 1.5d0
        Posz = dble(k) + 1.5d0
      else
        iC_ = iC + nbshift
        iO = ishft ( iC_ , - ndim )
        id = ishft ( 1, MaxLevel - iOctLv(iO) )
        j  = iC_ + 1 - ishft( iO , ndim )
        Posx=dble(d_x) * (dble(iOctPs(1,iO))+dble(sign(id,idelta(j,1))))
        Posy=dble(d_x) * (dble(iOctPs(2,iO))+dble(sign(id,idelta(j,2))))
        Posz=dble(d_x) * (dble(iOctPs(3,iO))+dble(sign(id,idelta(j,3))))
      endif
      return
      end


c     -----------------------------
       Function Age(aexp)
c     -----------------------------
c                                  age for LCDM model 
       implicit real*8(a-h,o-z)
       PARAMETER (Om0 =0.3)
       PARAMETER (tconst =9.31e+9)

       x = ((1.-Om0)/Om0)**0.333333*aexp
       Age = tconst /sqrt(1.-Om0)*(log(x**1.5+sqrt(1.+x**3)))

       Return
       End

