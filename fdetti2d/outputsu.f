c
      subroutine outputsu(fn_su,
     +                    n_sample,n_trace,nsi,dx,dz,dt,
     +                    ishot,newxs,newxr,newzr,sei_ps)
      implicit none
      include 'segy_simulation.h'

      integer   i,k,gdx,gdz
      integer   n_sample,n_trace,nsi,ishot,newxs
      integer:: newxr(1:n_trace),newzr(1:n_trace)
      real      dx,dz,dt
      real ::   sei_ps(1:n_trace,1:n_sample)
      character(*) fn_su

      include 'setsegyhead.h'

      open(30,file=trim(fn_su),form='unformatted',access='direct',
     +        recl=(60+n_sample)*4,status='replace')   
  
         do i=1,n_trace  
         ns=n_sample
         dt_head=nsi*dt*1000*1000
         tracl=i
         tracr=i
         fldr=ishot
         sx=(newxs-1)*dx
         gx=(newxr(i)-1)*dx
         offset=gx-sx
         if(i .eq. 1)then
             gdx=newxr(2)-newxr(1)
             gdz=newzr(2)-newzr(1)
         else
             gdx=newxr(i+1)-newxr(i)
             gdz=newzr(i+1)-newzr(i)
         endif
         d2=((gdx*dx)**2+(gdz*dz)**2)**0.5
         
        write(30,rec=i)
     1   tracl,tracr,fldr,tracf,ep,cdp,cdpt,
     1   trid,nvs,nhs,duse,offset,gelev,selev,sdepth,
     1   gdel,sdel,swdep,gwdep,scalel,scalco,
     1   sx,sy,gx,gy,counit,wevel,swevel,sut,gut,sstat,gstat,
     1   tstat,laga,lagb,delrt,muts,mute,ns,dt_head,
     1   gain,igc,igi,corr,sfs,sfe,slen,styp,stas,
     1   stae,tatyp,afilf,afils,nofilf,nofils,lcf,hcf,
     1   lcs,hcs,year,day,hour,minute,sec,timbas,
     1   trwf,grnors,grnofr,grnlof,gaps,otrav,
     1   d1,f1,d2,f2,ungpow,unscale,mark,(unass(k),k=1,17),
     1   (sei_ps(i,k),k=1,n_sample)
         enddo
       close(30)

       return
       end subroutine outputsu
