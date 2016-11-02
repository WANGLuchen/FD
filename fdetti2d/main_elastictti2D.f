c***********************************************************************
c           2D tti pseudo pure P wave equation forward simulation
c
c           modified by Wang Luchen
c           2013.10.27
c***********************************************************************

c***********************************************************************
c*************    program simulation start    **************************
c***********************************************************************
      program simulation_forward_elastic2D
      implicit none
      include 'segy_simulation.h'    !the different with segy.h is the dt change to dt_head 
      include 'mpif.h'  

      real pi
      parameter(pi=3.14159265)
c----------  indata  ---------------------------------
      integer source_flag,geometry_flag
      integer mx,mz,natte,nt
      integer shotbegx,shotbegz,sdx,sdz,n_shot
      integer gathbegx,gathbegz,gdx,gdz,n_trace
      real    dx,dz,dt,fp,pr,alpha,ppml
      integer mxl,mxr,dsg
      integer nn_shot_beg,n_sample,nsi
      integer mrlx,mrrx
      real    epsilons,delta,theta
      character *80 forvpfile,invvpfile,fepsilon,fdelta,ftheta
c-----------------------------------------------------
      integer nx,nz,ix,iz
      integer i,k,n,nts,nts1,sample_count,MM
      integer xl,xr
      
      integer newxs,newzs
      integer,allocatable::shotn(:),shotx(:),shotz(:)
      integer,allocatable::sta_x(:),sta_z(:)
      integer,allocatable::newxr(:),newzr(:)
c------------------------------------------------------
      integer nshot,nshotcom,nshotres,isb,ise,ishot
      integer p,me,info,status(mpi_status_size)
c------------------------------------------------------

      real b,t2
      real cpu,tel(2),dtime
      real dens
c--------- moment tensor and source wavelet function ------
      real m_tensor(3,3)
      real,allocatable::wtime(:)
c----------------------------------------------------------
      integer xlexpand,xrexpand
c      
      real,allocatable::sei_p(:,:),sei_r(:,:),sei_sum(:,:)
      real,allocatable::vp(:,:),vs(:,:),Aepsilon(:,:),Adelta(:,:),
     +                  Atheta(:,:)
c      real,allocatable::ps(:,:),rs(:,:)    !for staggered grid
c      real,allocatable::vpx(:,:),vpz(:,:),vrx(:,:),vrz(:,:)
      real,allocatable::pxx(:,:),pzz(:,:),pxz(:,:),vx(:,:),vz(:,:)
      real,allocatable::psfilter(:,:),rsfilter(:,:)
c
      character *80 fn_sour,fn_sta              ! source file & station file
      character *20 shotfile2
      character *15 filehead2
      character *13 shotfile1
      character *8  filehead1
      character *1  no

      include 'setsegyhead.h'

      no='0'
!      p=1
!      me=0

      filehead1='data/ccc'
      fn_sour='map_source'
      fn_sta='map_station'
!      filehead2='shotr_gather/rp'
c-----------------------------------------------------
c initialize mpi
      call mpi_init(info)
c get number of processes and processes rank.
      call mpi_comm_size(mpi_comm_world,p,info)
      call mpi_comm_rank(mpi_comm_world,me,info)
c-----------------------------------------------------
      cpu=dtime(tel)
c-----------------------------------------------------

c--------------------- open files ------------------------  
      open(1,file='indat',status='old')
      rewind(1)
c--------------------- read indat--------------------------
      read(1,*) source_flag,geometry_flag
      read(1,*) mx,mz,natte,nt,n_shot,n_trace
      read(1,*) shotbegx,shotbegz,sdx,sdz
      read(1,*) gathbegx,gathbegz,gdx,gdz
      read(1,*) dx,dz,dt,fp,pr,alpha,ppml
      read(1,*) mxl,mxr,dsg
      read(1,*) nn_shot_beg,n_sample,nsi      
      read(1,*) mrlx,mrrx
      read(1,*) forvpfile
      read(1,*) fepsilon
      read(1,*) fdelta
      read(1,*) ftheta
      read(1,*) invvpfile
      close(1)

      write(*,*) 'read in parameter OK'
      write(*,*) 'nt = ',nt
c-----------------------------------------------------
      nshot=n_shot
      nx=mx+natte+natte+mxl+mxr
      nz=mz+natte+natte

      allocate (  wtime(1:nt) )
      allocate (  shotn(1:n_shot),
     +            shotx(1:n_shot),
     +            shotz(1:n_shot) )
      allocate (  sta_x(1:n_trace),
     +            sta_z(1:n_trace) )

      allocate (  newxr(1:n_trace),
     +            newzr(1:n_trace) )

c   source_flag---> 0=calculated source, 1=source file as singleshot, 2=source file as multishot
      if (source_flag .ne. 0) then
          open(27,file=fn_sour,status='unknown')
          do i=1,n_shot
             read(27,*) shotn(i),shotx(i),shotz(i)
          enddo
          close(27)
          write(*,*) 'read in source file, OK'
          if (source_flag .eq. 2)then
              nshot=1
          endif
      endif
c   geometry_flag---> 1=station file,  2=calculated station
      if (geometry_flag .eq. 1) then
          open(28,file=fn_sta,status='unknown')
          do i=1,n_trace
             read(28,*) sta_x(i),sta_z(i)
          enddo
          close(28)
          write(*,*) 'read in station file, OK'
      endif

      nts=nint(1./fp/dt)
      nts1=2*nts+1   
   
c--------------- compute elastic parameter------------------
      dens=2500.0
c*************** create source wavelet *********************
c--------------- Ricker wavelet-----------------------------
      b=(pi*fp)**2    

      do n=1,nts1
         t2=((n-1-nts)*dt)**2
         wtime(n)=(1.-2.*b*t2)*exp(-b*t2)
      enddo
      do n=nts1+1,nt
         wtime(n)=0.0
      enddo
c--------------------------------------------------------------
c
c mpi distribution
       nshotcom=int(nshot/p)
       nshotres=nshot-(int(nshot/p))*p
       if (nshotres.eq.0) then
          isb=me*nshotcom+1
          ise=(me+1)*nshotcom
       endif
       if (nshotres.gt.0) then
          if ((me.ge.0).and.(me.le.(nshotres-1))) then
             isb=me*(nshotcom+1)+1
             ise=(me+1)*(nshotcom+1)
          endif
          if (me.gt.(nshotres-1)) then 
             isb=me*nshotcom+nshotres+1
             ise=(me+1)*nshotcom+nshotres
          endif
       endif
c--------------------------------------------------------------
c***************************mpi loop**************************
      do 999 ishot=isb+nn_shot_beg,ise+nn_shot_beg,1
c--------------set geometry ----------------------------------
        if(source_flag.eq.0) then
           newxs=(ishot-1)*sdx+shotbegx
           newzs=(ishot-1)*sdz+shotbegz+natte
c        elseif(source_flag.eq.1) then
        else
           newxs=shotx(ishot)
           newzs=shotz(ishot)+natte
        endif
        if(geometry_flag.eq.1) then
           do i=1,n_trace
              newxr(i)=sta_x(i)
              newzr(i)=sta_z(i)+natte
           enddo 
        elseif(geometry_flag.eq.2) then
           do i=1,n_trace
c              newzr(i)=401+gni*(i-1+dsg)+natte                  !receiver direction: single
c              newzr(i)=newxs-(n_trace-1)*gni/2+gni*(i-1+dsg)  !receiver direction: double
              newxr(i)=(i-1+dsg)*gdx+gathbegx
              newzr(i)=(i-1+dsg)*gdz+gathbegz+natte
           enddo 
        else
            write(*,*) 'the geometry_flag must is 1 or 2'
            stop
        endif
      
      xl=newxs-mrlx
      xr=newxs+mrrx

      allocate ( sei_p(1:n_trace,1:n_sample),
     +           sei_r(1:n_trace,1:n_sample),
     +         sei_sum(1:n_trace,1:n_sample) )
c
      allocate ( vp(xl-natte:xr+natte,nz),
     +           vs(xl-natte:xr+natte,nz),
     +     Aepsilon(xl-natte:xr+natte,nz),
     +       Atheta(xl-natte:xr+natte,nz),
     +       Adelta(xl-natte:xr+natte,nz) )
cc   for staggered grid
c      allocate ( ps(xl-natte:xr+natte,nz),
c     +           rs(xl-natte:xr+natte,nz), 
c     +          vpx(xl-natte:xr+natte,nz),
c     +          vpz(xl-natte:xr+natte,nz), 
c     +          vrx(xl-natte:xr+natte,nz),
c     +          vrz(xl-natte:xr+natte,nz) ) 
      allocate ( pxx(xl-natte:xr+natte,nz),
     +           pzz(xl-natte:xr+natte,nz), 
     +     psfilter(xl-natte:xr+natte,nz), 
     +     rsfilter(xl-natte:xr+natte,nz), 
     +          pxz(xl-natte:xr+natte,nz),
     +          vx(xl-natte:xr+natte,nz), 
     +          vz(xl-natte:xr+natte,nz) ) 
c
c****************  model velocity parameter  ********************************
c
c------ P wave velocity --------------------


      xlexpand=max(1-xl,0)
      xrexpand=max(xr-mx,0)
    
      do i=xl+xlexpand,xr-xrexpand
         do k=1+natte,nz-natte
            vp(i,k)=3200.0
            vs(i,k)=1847.57
            Aepsilon(i,k)=0.24
            Adelta(i,k)=0.1
            Atheta(i,k)=45.0
         enddo
      enddo

c*************************************************************************
c------ set values around absorbing area--------------
      call  setboundaries2d
     +          (xl,xr,1+natte,mz+natte,mx,mz,natte,vp)
      call  setboundaries2d
     +          (xl,xr,1+natte,mz+natte,mx,mz,natte,vs)
      call  setboundaries2d
     +          (xl,xr,1+natte,mz+natte,mx,mz,natte,Aepsilon)
      call  setboundaries2d
     +          (xl,xr,1+natte,mz+natte,mx,mz,natte,Adelta)
      call  setboundaries2d
     +          (xl,xr,1+natte,mz+natte,mx,mz,natte,Atheta)

c************************ forward *****************************
c--------------- initial variable values-------------------
cc   for staggered grid
c       ps=0.0
c       rs=0.0
c      vpx=0.0
c      vpz=0.0
c      vrx=0.0
c      vrz=0.0
       psfilter=0.0
       rsfilter=0.0
      pxx=0.0
      pzz=0.0
      pxz=0.0
      vx=0.0
      vz=0.0

      sei_p=0.0
      sei_r=0.0
      sei_sum=0.0
c
c---------------------------
      sample_count=0
c---------------forward computing---------------------------

c---------------loop over time------------------------------
      do 100 n=1,nt
c*************** add source(Ricker wavelet) force **************************
c
      if(source_flag.eq.2) then ! multishot
          do i=1,n_shot
          pxx(shotx(i),shotz(i)+natte)= pxx(shotx(i),shotz(i)+natte)
     +                                -wtime(n)
          pzz(shotx(i),shotz(i)+natte)= pzz(shotx(i),shotz(i)+natte)
     +                                -wtime(n)
          enddo
      else! singleshot
          ! let source contribute within limited distance
          do ix=max(1,newxs-5),min(mx,newxs+5)
          do iz=max(1+natte,newzs-5),min(mz+natte,newzs+5)
         pzz(ix,iz)= pzz(ix,iz)
     +             +wtime(n)*exp(-1.0*(ix-newxs)**2-1.0*(iz-newzs)**2)
         pxx(ix,iz)= pxx(ix,iz)
     +             +wtime(n)*exp(-1.0*(ix-newxs)**2-1.0*(iz-newzs)**2)
          enddo
          enddo
      endif
c***********************************************************
c
c###################  compute the field  ####################################
c
c------------  use staggered grid to compute the velocity and stress
c
      call forwarding_tti_elastic_2d_staggered_10order
     +    (xl,xr,1+natte,mz+natte,natte,vp,vs,Aepsilon,Adelta,Atheta,
     +     vx,vz,pxx,pzz,pxz,dx,dz,dt,4)

      call spongeatten_2d_2(xl,xr,1+natte,mz+natte,natte,alpha,vx)
      call spongeatten_2d_2(xl,xr,1+natte,mz+natte,natte,alpha,vz)

c----------------------------------------------------------------------
c
      if (mod((n-1),nsi) .eq. 0 )then
      sample_count=sample_count+1
      do i=1,n_trace
       sei_p(i,sample_count)  = pxx(newxr(i),newzr(i))
      enddo 
      endif 
c
c---------------  make snapshot  ------------------------------------
      if ((n.gt.0).and.(n.le.nt))then
         if ( mod((n-1),50) .eq. 0 )then
c      call filter_pseudo_pureP_tti_2d
c     +         (xl,xr,1+natte,mz+natte,natte,vp,Aepsilon,Adelta,Atheta,
c     +          ps,rs,psfilter,rsfilter,dx,dz,dt,4)
      call snapshots2D(ishot,n,xl,xr,1+natte,mz+natte,
     +                 mx,mz,natte,pxx,'fps')
      call snapshots2D(ishot,n,xl,xr,1+natte,mz+natte,
     +                 mx,mz,natte,pxx+pzz,'fpxz')
c      call snapshots2D(ishot,n,xl,xr,1+natte,mz+natte,
c     +                 mx,mz,natte,psfilter+rsfilter,'fss')
c      call snapshots2D(ishot,n,xl,xr,1+natte,mz+natte,
c     +                 mx,mz,natte,rsfilter,'frs')
         endif
      endif
100   continue
c---------------------output the shotgather -------------------------
      call read_shot_filename(filehead1,ishot,1,shotfile1)     
      
      write(*,*)'shotfile1=',shotfile1
      call outputsu(trim(shotfile1)//'_ps.su',
     +              n_sample,n_trace,nsi,dx,dz,dt,
     +              ishot,newxs,newxr,newzr,sei_p) 

c---------------deallocate all the arrays -----------------------
      deallocate ( sei_p,sei_r,sei_sum )
      deallocate ( vp,vs,Aepsilon,Adelta )
cc   for staggered grid
cc      deallocate ( vpx,vpz,vrx,vrz )
c   for regular grid
      deallocate ( pxx,pzz,pxz,vx,vz )
      deallocate ( psfilter,rsfilter )

999   continue

      deallocate ( shotn,shotx,shotz )
      deallocate ( sta_x,sta_z )
      deallocate ( newxr,newzr )

c------------------------------------
      cpu=dtime(tel)
c------------------------------------
      write(*,*) 'computing finished!'

      write(*,*)'node:',me+1,'no. of tables:',ise-isb+1

      write(*,*)'total cpu time used=',cpu,'(sec)'

      call mpi_barrier(mpi_comm_world,info)
      call mpi_finalize(info)

      stop
      end
c***************************************************************
c***********                                       *************
c***********   main program sm_for_elastic2D end   *************
c***********                                       *************
c***************************************************************
