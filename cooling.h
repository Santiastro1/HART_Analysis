      integer nlt, nld, nlz, nrs
      real*8 tlmin, tlmax, dlt
      real*8 dlmin, dlmax, dld
      real*8 Zlmin, Zlmax, dlz
      real*8 rsmin, rsmax, drs
      real*8 dlti, dldi, dlzi, drsi
      common / CLCOOL01 / tlmin, tlmax, dlt, nlt
      common / CLCOOL02 / dlmin, dlmax, dld, nld                
      common / CLCOOL03 / Zlmin, Zlmax, dlz, nlz
      common / CLCOOL04 / rsmin, rsmax, drs, nrs
      common / CLCOOL05 / dlti, dldi, dlzi, drsi

      real*8 smallrate
      parameter ( smallrate = 1.d-30 ) 
      parameter ( nltmax = 71 ) 
      parameter ( nldmax = 11 ) 
      parameter ( nlzmax = 9  ) 
      parameter ( nrsmax = 21 ) 
      real*8 coolcl(nltmax,nldmax,nlzmax,nrsmax)
      real*8 ccl_rs(nltmax,nldmax,nlzmax)
      real*8 f_ion(nltmax,nldmax,nlzmax,nrsmax)
      common / CLCOOL05 / coolcl
      common / CLCOOL06 / ccl_rs
      common / CLCOOL07 / f_ion
