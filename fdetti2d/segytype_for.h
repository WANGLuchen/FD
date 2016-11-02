      type :: segy
      integer*4 tracl,tracr,fldr,tracf,ep,cdp,cdpt
      integer*2 trid,nvs,nhs,duse
      integer*4 offset,gelev,selev,sdepth,gdel,sdel,swdep,gwdep
      integer*2 scalel,scalco
      integer*4 sx,sy,gx,gy
      integer*2 counit,wevel,swevel,sut,gut,sstat,gstat,
     1            tstat,laga,lagb,delrt,muts,mute,ns,dt
      integer*2 gain,igc,igi,corr,sfs,sfe,slen,styp,stas,
     1            stae,tatyp,afilf,afils,nofilf,nofils,lcf,hcf,
     1            lcs,hcs,year,day,hour,minute,sec,timbas,
     1            trwf,grnors,grnofr,grnlof,gaps,otrav
      real*4 d1,f1,d2,f2,ungpow,unscale
      integer*2 mark,unass(17)
      end type segy
