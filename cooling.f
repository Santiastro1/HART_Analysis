
c     ---------------------------------------------------------
      subroutine Cooling ( rhogl, T_g , Z_met , coolrate )
      
c     ---------------------------------------------------------
c     input:  
c         rhogl - log10(n_H), where n_H is hydrogen number density in cm^-3
c         T_g  - gas temperature in units of 10^4 K
c         Z_met - log10([rho_metals/rho_gas/solar)
c         if ENRICH is not defined Z_met will be ignored 
c         if SD93_COOLING is defined, rhogl will be ignored 
c
c     output: coolrate - net output cooling/heating rate in code units
c        de/dt = e - coolrate * rho_g**2 / aexpn
c        where e is internal energy per unit volume and rho_g is 
c        gas density (e and rho_g are in code units)
c
c#     include "a_def.h"
c      include 'a_tree.h'
c      include 'a_control.h'
       include 'parameters.h'
      real*8 rhogl, rho_g, T_g, Z_met, coolrate
      real*8 Tlog,rs
      real*8 ac, bc, ah, bh
c      write(*,*)'T,Z,rho=',T_g,Z_met,rhogl     
c      hubble=0.7
c      Om0=0.3
c      rhogl=10.0
c      T_G=100.0
c      Z_met=0.02
c      rs=0.00

c      call Set_Cooling ()
c      call Set_Cooling_Rate_rs ( rs )

c
c.... use net cooling/heating rate tabulated using Cloudy
c

      Tlog = log10(T_g) + 4.d0 
      it1 = int((Tlog - tlmin)*dlti) + 1
      it2 = it1 + 1
      id1 = int((rhogl - dlmin)*dldi) + 1
      id2 = id1 + 1
      iz1 = int((Z_met - zlmin)*dlzi) + 1
      iz2 = iz1 + 1
      it1 = max(it1,1) 
      it1 = min(it1,nlt)
      it2 = max(it2,1) 
      it2 = min(it2,nlt)
      id1 = max(id1,1) 
      id1 = min(id1,nld)
      id2 = max(id2,1) 
      id2 = min(id2,nld)
      iz1 = max(iz1,1) 
      iz1 = min(iz1,nlz)
      iz2 = max(iz2,1) 
      iz2 = min(iz2,nlz)

      td = tlmin + dlt * (it1 - 1) + dlt
      dd = dlmin + dld * (id1 - 1) + dld
      zd = zlmin + dlz * (iz1 - 1) + dlz
      t1 = (td - Tlog) * dlti
      d1 = 1.d0 - t1
      t2 = (dd - rhogl) * dldi
      d2 = 1.d0 - t2
      t3 = (zd - Z_met) * dlzi 
      d3 = 1.d0 - t3
c
c.... trilinear (CIC) interpolation 
c
      coolrate = t1*t2*t3 * ccl_rs(it1,id1,iz1) +
     &           d1*t2*t3 * ccl_rs(it2,id1,iz1) + 
     &           t1*d2*t3 * ccl_rs(it1,id2,iz1) + 
     &           d1*d2*t3 * ccl_rs(it2,id2,iz1) +
     &           t1*t2*d3 * ccl_rs(it1,id1,iz2) + 
     &           d1*t2*d3 * ccl_rs(it2,id1,iz2) + 
     &           t1*d2*d3 * ccl_rs(it1,id2,iz2) + 
     &           d1*d2*d3 * ccl_rs(it2,id2,iz2)

c      write(*,*)'t1,t2,t3',t1,t2,t3
c      write(*,*)'d1,d2,d3',d1,d2,d3
c      write(*,*)ccl_rs(it1,id1,iz1),ccl_rs(it2,id1,iz1),
c     +ccl_rs(it1,id2,iz1),ccl_rs(it2,id2,iz1),ccl_rs(it1,id1,iz2),
c     +ccl_rs(it2,id1,iz2),ccl_rs(it1,id2,iz2),ccl_rs(it2,id2,iz2)
c      write(*,*)'coolrate',coolrate
c      stop
      return
      end


c     -------------------------------------
      subroutine Set_Cooling_Rate_rs ( rs ) 
c     -------------------------------------
c
c     prepare cooling rate table for a given redshift rs
c     from the Cloudy cooling rate table 
c
c     real*8 rs - redshift 
c     
c     this routine is only used if CLOUDY_COOLING is defined in a_def.h
c     it should be called only in the beginning of every step 
c
c
      include 'parameters.h' 
      real*8 rs, rsd
      real*8 ac, bc, ah, bh

      irs = int((rs - rsmin)*drsi) + 1
      irs1 = max(irs,1)
      irs1 = min(irs1,nrs)
      irs2 = min(irs+1,nrs)
      irs2 = max(irs2,1)

      if ( irs1 .eq. irs2 ) then 
        do ilt = 1 , nlt
        do ild = 1 , nld
        do ilz = 1 , nlz 
          ccl_rs(ilt,ild,ilz) = coolcl(ilt,ild,ilz,irs1)
        enddo
        enddo
        enddo
      else
        rs1 = rsmin + drs*(irs1-1)
        rs2 = rsmin + drs*(irs2-1)

        do ilt = 1 , nlt
        do ild = 1 , nld
        do ilz = 1 , nlz 
          ac = (coolcl(ilt,ild,ilz,irs2) - coolcl(ilt,ild,ilz,irs1)) /
     &         (rs2 - rs1)
          bc = coolcl(ilt,ild,ilz,irs1) - ac * rs1 
          ccl_rs(ilt,ild,ilz) = ac * rs + bc 
        enddo
        enddo
        enddo
      endif

      return
      end


c     -------------------------
      subroutine Set_Cooling ()
c     -------------------------
c
c     tabulate cooling curve as a function of T in units of 10^4 K
c     table entry i corresponds to i = 100*T^1/4
c     AL_0 must be set in Set_Units prior to call to this routine

      include 'parameters.h'
      real*8 cdum, hdum, ct_crit
c.... use net cooling/heating rate tabulated using Cloudy
c
      open ( 40 , file = 'clcool.dat' )
      read(40,*) tlmin, tlmax, dlt, nlt
      read(40,*) dlmin, dlmax, dld, nld
      read(40,*) Zlmin, Zlmax, dlZ, nlz
      read(40,*) rsmin, rsmax, drs, nrs

c      write(*,*) tlmin, tlmax, dlt, nlt
c      write(*,*) dlmin, dlmax, dld, nld
c      write(*,*) Zlmin, Zlmax, dlZ, nlz
c      write(*,*) rsmin, rsmax, drs, nrs

      dlti = 1.d0 / dlt
      dldi = 1.d0 / dld
      dlzi = 1.d0 / dlz
      drsi = 1.d0 / drs 

      do irs = 1, nrs
        do ilz = 1 , nlz
          do ild = 1 , nld 

            ct_crit = 0.0
            do ilt = 1 , nlt
              read(40,*) d1, d2, d3, d4, d5, d6, d7, d8, d9, cdum, hdum
              if ( d1 .ge. 3.2 .and. ct_crit .eq. 0.0 ) then
                ct_crit = cdum  
              endif
              cdum = max ( cdum , smallrate )
              cdum = max ( cdum , ct_crit   )  ! fix for a trough in
c                                                equilibrium H2 cooling curve

cc       AL_SD = 1.6625d0 * (1.d0 - Y_p)**2 * hubble /
cc     &         sqrt(Om0) / (boxh/ng)**2 ! in ergs/s/cm^3

              hdum = max ( hdum , smallrate ) 
              coolcl(ilt,ild,ilz,irs) = (cdum - hdum) !* 1.d23 * AL_SD
c              write(*,*)coolcl(ilt,ild,ilz,irs)
              f_ion(ilt,ild,ilz,irs) = d6 / 10.d0**(d2)
            enddo
          enddo
        enddo
      enddo
      close ( 40 ) 

      return
      end

