c***********************************************************************
c           3D tti elastic wave equation forward simulation
c
c           modified by Wang Luchen
c           2015.07.06
c***********************************************************************

c***********************************************************************
c*************    program simulation start    **************************
c***********************************************************************
      program simulation_forward_tti_elastic_3D
      implicit none
      include 'segy_simulation.h'    !the different with segy.h is the dt change to dt_head 
      include 'mpif.h'  

      real pi
      parameter(pi=3.14159265)
c----------  indata  ---------------------------------
      integer source_flag,geometry_flag,model_flag
      integer mx,my,mz,natte,nt
      integer shotbegx,shotbegy,shotbegz,sdx,sdy,sdz,n_shot
      integer gathbegx,gathbegy,gathbegz,gdx,gdy,gdz,n_trace
      real    dx,dy,dz,dt,fp,pr,alpha,ppml
      integer mxl,mxr,myl,myr,dsg
      integer nn_shot_beg,n_sample,nsi
      integer mrlx,mrrx,mrly,mrry
      real    epsilons,delta,gammas,theta,phi
      character *80 forvpfile,invvpfile
c-----------------------------------------------------
      integer nx,ny,nz,ix,iy,iz
      integer xbeg,xend,ybeg,yend,zbeg,zend
      integer i,j,k,n,nts,nts1,sample_count,MM
      
      integer newxs,newys,newzs
      integer,allocatable::shotn(:),shotx(:),shoty(:),shotz(:)
      integer,allocatable::sta_n(:),sta_x(:),sta_y(:),sta_z(:)
      integer,allocatable::newxr(:),newyr(:),newzr(:)
c------------------------------------------------------
      integer nshot,nshotcom,nshotres,isb,ise,ishot
      integer p,me,info,status(mpi_status_size),access,filestatus
c------------------------------------------------------

      real b,t2
      real cpu,tel(2),dtime
      real dens
c--------- moment tensor and source wavelet function ------
      real m_tensor(3,3)
      real,allocatable::wtime(:)
c----------------------------------------------------------
      integer xlexpand,xrexpand,ylexpand,yrexpand
c      
      real,allocatable::
     +       sei_vx(:,:),sei_vy(:,:),sei_vz(:,:)
      real,allocatable::
     +       vp(:,:,:),vs(:,:,:),
     +       Aepsilon(:,:,:),Adelta(:,:,:),Agamma(:,:,:),
     +       Atheta(:,:,:),Aphi(:,:,:)
c      real,allocatable::xzvp(:,:),xzAepsilon(:,:),xzAdelta(:,:)
c      !for staggered grid
      real,allocatable::
     +        vx(:,:,:), vy(:,:,:), vz(:,:,:),
     +       pxx(:,:,:),pyy(:,:,:),pzz(:,:,:),
     +       pxy(:,:,:),pxz(:,:,:),pyz(:,:,:)
c  save the summation of wave-field snapshot
!      real,allocatable::
!     +       ss_vx(:,:,:), ss_vy(:,:,:), ss_vz(:,:,:)
c  stiffness matrix coefficient
      real,allocatable::
     +       c11(:,:,:),c12(:,:,:),c13(:,:,:),c14(:,:,:),c15(:,:,:),
     +       c16(:,:,:),c22(:,:,:),c23(:,:,:),c24(:,:,:),c25(:,:,:),
     +       c26(:,:,:),c33(:,:,:),c34(:,:,:),c35(:,:,:),c36(:,:,:),
     +       c44(:,:,:),c45(:,:,:),c46(:,:,:),c55(:,:,:),c56(:,:,:),
     +       c66(:,:,:)
c
      character *80 fn_sour,fn_sta              ! source file & station file
      character *80 shotfile1
      character *80 filehead1

      include 'setsegyhead.h'

      filehead1='data/ccc'
      fn_sour='map_source'
      fn_sta='map_station'
      m_tensor=0.0
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
      read(1,*) source_flag,geometry_flag,model_flag
      read(1,*) mx,my,mz,natte,nt,n_shot,n_trace
      read(1,*) shotbegx,shotbegy,shotbegz,sdx,sdy,sdz
      read(1,*) gathbegx,gathbegy,gathbegz,gdx,gdy,gdz
      read(1,*) dx,dy,dz,dt,fp,pr,alpha,ppml
      read(1,*) mxl,mxr,myl,myr,dsg
      read(1,*) nn_shot_beg,n_sample,nsi      
      read(1,*) mrlx,mrrx,mrly,mrry
      read(1,*) theta,phi
      read(1,*) m_tensor(1,1),m_tensor(2,2),m_tensor(3,3),
     +          m_tensor(2,3),m_tensor(1,3),m_tensor(1,2)
c      read(1,*) forvpfile
c      read(1,*) fepsilon
c      read(1,*) fdelta
c      read(1,*) invvpfile
      close(1)

      write(*,*) 'read in parameter OK'
      write(*,*) 'nt=',nt
c-----------------------------------------------------
      nshot=n_shot
      nx=mx+natte+natte+mxl+mxr
      ny=my+natte+natte+myl+myr
      nz=mz+natte+natte

      allocate (  wtime(1:nt) )
      allocate (  shotn(1:nn_shot_beg+n_shot),
     +            shotx(1:nn_shot_beg+n_shot),
     +            shoty(1:nn_shot_beg+n_shot),
     +            shotz(1:nn_shot_beg+n_shot) )
      allocate (  sta_n(1:n_trace),
     +            sta_x(1:n_trace),
     +            sta_y(1:n_trace),
     +            sta_z(1:n_trace) )

      allocate (  newxr(1:n_trace),
     +            newyr(1:n_trace),
     +            newzr(1:n_trace) )

c   source_flag---> 0=calculated source, 1=source file as singleshot, 2=source file as multishot
      if (source_flag .ne. 0) then
          open(27,file=fn_sour,status='unknown')
          do i=1,nn_shot_beg+n_shot
             read(27,*) shotn(i),shotx(i),shoty(i),shotz(i)
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
             read(28,*) sta_n(i),sta_x(i),sta_y(i),sta_z(i)
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
       write(*,*)'loop start!'
c***************************mpi loop**************************
      do 999 ishot=isb+nn_shot_beg,ise+nn_shot_beg,1
        write(*,*)'isbeg=',isb,'isend=',ise,'ishot=',ishot
c----------- check output file name ----------------------
      call read_shot_filename(filehead1,ishot,1,shotfile1)
      filestatus=access(trim(shotfile1)//'_vx.su',' ')
      if(filestatus.eq.0)then
           write(*,*)'file:',trim(shotfile1),'already exists'
           goto 999
      endif
c--------------set geometry ----------------------------------
        if(source_flag.eq.0) then
           newxs=(ishot-1)*sdx+shotbegx
           newys=(ishot-1)*sdy+shotbegy
           newzs=(ishot-1)*sdz+shotbegz
        else
           newxs=shotx(ishot)
           newys=shoty(ishot)
           newzs=shotz(ishot)
        endif
        if(geometry_flag.eq.1) then
           do i=1,n_trace
              newxr(i)=sta_x(i)
              newyr(i)=sta_y(i)
              newzr(i)=sta_z(i)
           enddo 
        elseif(geometry_flag.eq.2) then
           do i=1,n_trace
c              newzr(i)=401+gni*(i-1+dsg)                  !receiver direction: single
c              newzr(i)=newxs-(n_trace-1)*gni/2+gni*(i-1+dsg)  !receiver direction: double
              newxr(i)=(i-1+dsg)*gdx+gathbegx
              newyr(i)=(i-1+dsg)*gdy+gathbegy
              newzr(i)=(i-1+dsg)*gdz+gathbegz
           enddo 
        else
            write(*,*) 'the geometry_flag must is 1 or 2'
            stop
        endif
      
        write(*,*)'shotx=',newxs,'shoty=',newys,'shotz=',newzs

      xbeg=newxs-mrlx
      xend=newxs+mrrx
      ybeg=newys-mrly
      yend=newys+mrry
      zbeg=1
      zend=mz

      allocate ( sei_vx(1:n_trace,1:n_sample),
     +           sei_vy(1:n_trace,1:n_sample),
     +           sei_vz(1:n_trace,1:n_sample) )
c
      allocate ( 
     +          vp(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +             zbeg-natte:zend+natte),
     +          vs(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +             zbeg-natte:zend+natte),
     +    Aepsilon(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +             zbeg-natte:zend+natte),
     +      Adelta(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +             zbeg-natte:zend+natte),
     +      Agamma(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +             zbeg-natte:zend+natte),
     +      Atheta(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +             zbeg-natte:zend+natte),
     +        Aphi(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +             zbeg-natte:zend+natte) )
!      allocate ( xzvp(xbeg-natte:xend+natte,zbeg-natte:zend+natte),
!     +     xzAepsilon(xbeg-natte:xend+natte,zbeg-natte:zend+natte),
!     +       xzAdelta(xbeg-natte:xend+natte,zbeg-natte:zend+natte) )
cc   for staggered grid
      allocate ( 
     +          vx(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +             zbeg-natte:zend+natte),
     +          vy(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +             zbeg-natte:zend+natte), 
     +          vz(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +             zbeg-natte:zend+natte), 
     +         pxx(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +             zbeg-natte:zend+natte),
     +         pyy(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +             zbeg-natte:zend+natte), 
     +         pzz(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +             zbeg-natte:zend+natte),
     +         pxy(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +             zbeg-natte:zend+natte), 
     +         pxz(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +             zbeg-natte:zend+natte),
     +         pyz(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +             zbeg-natte:zend+natte) ) 

      allocate ( c11(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +               zbeg-natte:zend+natte),
     +           c12(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +               zbeg-natte:zend+natte), 
     +           c13(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +               zbeg-natte:zend+natte), 
     +           c14(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +               zbeg-natte:zend+natte),
     +           c15(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +               zbeg-natte:zend+natte), 
     +           c16(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +               zbeg-natte:zend+natte),
     +           c22(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +               zbeg-natte:zend+natte), 
     +           c23(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +               zbeg-natte:zend+natte), 
     +           c24(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +               zbeg-natte:zend+natte), 
     +           c25(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +               zbeg-natte:zend+natte), 
     +           c26(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +               zbeg-natte:zend+natte), 
     +           c33(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +               zbeg-natte:zend+natte), 
     +           c34(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +               zbeg-natte:zend+natte), 
     +           c35(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +               zbeg-natte:zend+natte), 
     +           c36(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +               zbeg-natte:zend+natte), 
     +           c44(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +               zbeg-natte:zend+natte), 
     +           c45(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +               zbeg-natte:zend+natte), 
     +           c46(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +               zbeg-natte:zend+natte), 
     +           c55(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +               zbeg-natte:zend+natte), 
     +           c56(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +               zbeg-natte:zend+natte),
     +           c66(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +               zbeg-natte:zend+natte) ) 
!      allocate ( ss_vx(xbeg-natte:xend+natte,ybeg-natte:yend+natte,zbeg-natte:zend+natte),
!     +           ss_vy(xbeg-natte:xend+natte,ybeg-natte:yend+natte,zbeg-natte:zend+natte),
!     +           ss_vz(xbeg-natte:xend+natte,ybeg-natte:yend+natte,zbeg-natte:zend+natte) )
c
c****************  model velocity parameter  ********************************
      ! This is a two layers model, an isotropic media and a tti media
      do k=1,40
          vp(:,:,k)=3200.0
          vs(:,:,k)=1800.0
          Aepsilon(:,:,k)=0.0
          Adelta(:,:,k)=0.0
      enddo
      do k=41,85
          vp(:,:,k)=4000.0
          vs(:,:,k)=2300.0
          Aepsilon(:,:,k)=0.24
          Adelta(:,:,k)=0.1
      enddo
      do k=86,zend
          vp(:,:,k)=4200.0
          vs(:,:,k)=2500.0
          Aepsilon(:,:,k)=0.05
          Adelta(:,:,k)=0.05
      enddo
      Agamma=0.0
      Atheta=theta
      Aphi=phi

ccc  this part script is correct, wanglc 2015.07.06  **********
c
c------ P wave velocity --------------------
c      write(*,*) 'for_vp=',forvpfile
c      write(*,*) 'for_epsilon=',fepsilon
c      write(*,*) 'for_delta=',fdelta
c      xlexpand=max(1-xl,0)
c      xrexpand=max(xr-mx,0)
c      ylexpand=max(1-yl,0)
c      yrexpand=max(yr-my,0)
c    
c      if(model_flag.eq.3)then
cc   read in 3D model file
c        call input_vel3D_zxy(trim(forvpfile),xl,xr,yl,yr,
c     +                     1+natte,mz+natte,mx,my,mz,natte,vp)
c        call input_vel3D_zxy(trim(fepsilon),xl,xr,yl,yr,
c     +                     1+natte,mz+natte,mx,my,mz,natte,Aepsilon)
c        call input_vel3D_zxy(trim(fdelta),xl,xr,yl,yr,
c     +                     1+natte,mz+natte,mx,my,mz,natte,Adelta)
c      elseif(model_flag.eq.2)then
cc   read in 2D model file and extrapolate to 3D
c        open(2,file=trim(forvpfile),
c     1    form='unformatted',access='direct',
c     2    recl=mz*4,status='old')
c        open(3,file=trim(fepsilon),
c     1    form='unformatted',access='direct',
c     2    recl=mz*4,status='old')
c        open(4,file=trim(fdelta),
c     1    form='unformatted',access='direct',
c     2    recl=mz*4,status='old')
c
c         do i=xl+xlexpand,xr-xrexpand
c            read(2,rec=i) (xzvp(i,k),k=1+natte,mz+natte)
c            read(3,rec=i) (xzAepsilon(i,k),k=1+natte,mz+natte)
c            read(4,rec=i) (xzAdelta(i,k),k=1+natte,mz+natte)
c
c            do j=yl+ylexpand,yr-yrexpand
c            do k=1+natte,nz-natte
c                  vp(i,j,k)=xzvp(i,k)
c            Aepsilon(i,j,K)=xzAepsilon(i,k)
c              Adelta(i,j,k)=xzAdelta(i,k)
c            enddo
c            enddo
c
c         enddo
c         close(2)
c         close(3)
c         close(4)
c      endif

c*************************************************************************
c------ set values around absorbing area--------------
      call  setboundaries3d
     +          (xbeg,xend,ybeg,yend,zbeg,zend,mx,my,mz,natte,vp)
      call  setboundaries3d
     +          (xbeg,xend,ybeg,yend,zbeg,zend,mx,my,mz,natte,vs)
      call  setboundaries3d
     +          (xbeg,xend,ybeg,yend,zbeg,zend,mx,my,mz,natte,Aepsilon)
      call  setboundaries3d
     +          (xbeg,xend,ybeg,yend,zbeg,zend,mx,my,mz,natte,Adelta)
      call  setboundaries3d
     +          (xbeg,xend,ybeg,yend,zbeg,zend,mx,my,mz,natte,Agamma)
      call  setboundaries3d
     +          (xbeg,xend,ybeg,yend,zbeg,zend,mx,my,mz,natte,Atheta)
      call  setboundaries3d
     +          (xbeg,xend,ybeg,yend,zbeg,zend,mx,my,mz,natte,Aphi)

c************************ forward *****************************
c--------------- initial variable values-------------------
cc   for staggered grid
      vx=0.0
      vy=0.0
      vz=0.0
      pxx=0.0
      pyy=0.0
      pzz=0.0
      pxy=0.0
      pxz=0.0
      pyz=0.0

      sei_vx=0.0
      sei_vy=0.0
      sei_vz=0.0

c
c---------------------------
      sample_count=0
c********************************************************************
c*************  forward computing  **********************************
c********************************************************************

      call stiffness_matrix_coefficient
     +   ( xbeg,xend,ybeg,yend,zbeg,zend,natte,
     +     vp,vs,Aepsilon,Adelta,Agamma,Atheta,Aphi,
     +     c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,
     +     c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)

c---------------loop over time------------------------------
      do 100 n=1,nt
c*************** add source(Ricker wavelet) force **************************
c
!      if(source_flag.eq.2) then ! multishot
!          do i=1,n_shot
!             ps(shotx(i),shoty(i),shotz(i))= 
!     +      +ps(shotx(i),shoty(i),shotz(i))-wtime(n)/dt
!          enddo
!      else! singleshot
         call source_loading(newxs,newys,newzs,wtime,
     +        xbeg,xend,ybeg,yend,zbeg,zend,natte,nt,dx,dy,dz,dt,
     +        m_tensor,pxx,pyy,pzz,pyz,pxz,pxy,vx,vy,vz,n,1)
!      endif

c***********************************************************
c
c###################  compute the field  ####################################
c
c------------  use staggered grid to compute the velocity and stress
c
       call forwarding_tti_elastic_3d_staggered_10order_2
     +    ( xbeg,xend,ybeg,yend,zbeg,zend,natte,
     +      c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,
     +      c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66,
     +      vx,vy,vz,pxx,pyy,pzz,pxy,pxz,pyz,dx,dy,dz,dt,4 )


      call spongeatten_3d_2(xbeg,xend,ybeg,yend,
     +                      zbeg,zend,natte,alpha,vx)
      call spongeatten_3d_2(xbeg,xend,ybeg,yend,
     +                      zbeg,zend,natte,alpha,vy)
      call spongeatten_3d_2(xbeg,xend,ybeg,yend,
     +                      zbeg,zend,natte,alpha,vz)
c----------------------------------------------------------------------
c
      if (mod(n,nsi) .eq. 0 )then
      sample_count=sample_count+1
      do i=1,n_trace
       sei_vx(i,sample_count)  = vx(newxr(i),newyr(i),newzr(i))
       sei_vy(i,sample_count)  = vy(newxr(i),newyr(i),newzr(i))
       sei_vz(i,sample_count)  = vz(newxr(i),newyr(i),newzr(i))
      enddo 
      endif 
c
c---------------  make snapshot  ------------------------------------
      if ( mod((n-1),100) .eq. 0 )then
          write(*,*)'snap, n=',n
      call snapshots3D(ishot,n,xbeg,xend,ybeg,yend,zbeg,zend,
     +                 mx,my,mz,natte,newys,newzs,vx,'fvx')
      call snapshots3D(ishot,n,xbeg,xend,ybeg,yend,zbeg,zend,
     +                 mx,my,mz,natte,newys,newzs,vy,'fvy')
      call snapshots3D(ishot,n,xbeg,xend,ybeg,yend,zbeg,zend,
     +                 mx,my,mz,natte,newys,newzs,vz,'fvz')
      endif
c-----------  output the 3D wavefield  ---------------------
      if ( mod(n,500) .eq. 0 )then
      call outputbin3D(n,'vx',xbeg,xend,ybeg,yend,zbeg,zend,
     +               mx,my,mz,natte,vx)
      call outputbin3D(n,'vy',xbeg,xend,ybeg,yend,zbeg,zend,
     +               mx,my,mz,natte,vy)
      call outputbin3D(n,'vz',xbeg,xend,ybeg,yend,zbeg,zend,
     +               mx,my,mz,natte,vz)
      endif
c--------------------------------------------------------------------
100   continue
c---------------------output the shotgather -------------------------
      call read_shot_filename(filehead1,ishot,1,shotfile1)

      write(*,*)'shotfile1=',shotfile1
      call outputsu3(trim(shotfile1)//'_vx.su',
     +               n_sample,n_trace,nsi,dx,dy,dt,
     +               ishot,newxs,newys,newxr,newyr,sei_vx)
      call outputsu3(trim(shotfile1)//'_vy.su',
     +               n_sample,n_trace,nsi,dx,dy,dt,
     +               ishot,newxs,newys,newxr,newyr,sei_vy)
      call outputsu3(trim(shotfile1)//'_vz.su',
     +               n_sample,n_trace,nsi,dx,dy,dt,
     +               ishot,newxs,newys,newxr,newyr,sei_vz)

c---------------deallocate all the arrays -----------------------
      deallocate ( sei_vx,sei_vy,sei_vz )
!      deallocate ( ss_vx,ss_vy,ss_vz )
      deallocate ( vp,vs,Aepsilon,Adelta,Agamma,Atheta,Aphi )
      deallocate ( c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,
     +             c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66 )
!      deallocate ( xzvp,xzAepsilon,xzAdelta )
cc   for staggered grid
      deallocate ( vx,vy,vz,pxx,pyy,pzz,pxy,pxz,pyz )

999   continue

      deallocate ( shotn,shotx,shoty,shotz )
      deallocate ( sta_n,sta_x,sta_y,sta_z )
      deallocate ( newxr,newyr,newzr )

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
c*********                                           ***********
c*********   main program main_forttipseudop3d end   ***********
c*********                                           ***********
c***************************************************************
