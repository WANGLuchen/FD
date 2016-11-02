c******************************************************************************
c******************   TTI elastic 2D forwarding staggered grid   ***********
c
      subroutine forwarding_tti_elastic_2d_staggered_10order
     1  ( xbeg,xend,zbeg,zend,natte,vp,vs,Aepsilon,Adelta,Atheta,
     2      vx,vz,pxx,pzz,pxz,dx,dz,dt,inf ) 
      implicit none

      integer xbeg,xend,zbeg,zend,natte
      real theta,dx,dz,dt,dens
      integer inf
      real:: vp(xbeg-natte:xend+natte,zbeg-natte:zend+natte),
     1       vs(xbeg-natte:xend+natte,zbeg-natte:zend+natte),
     1 Aepsilon(xbeg-natte:xend+natte,zbeg-natte:zend+natte),
     1   Adelta(xbeg-natte:xend+natte,zbeg-natte:zend+natte),
     1   Atheta(xbeg-natte:xend+natte,zbeg-natte:zend+natte),
     1       vx(xbeg-natte:xend+natte,zbeg-natte:zend+natte),
     1       vz(xbeg-natte:xend+natte,zbeg-natte:zend+natte),
     1      pxx(xbeg-natte:xend+natte,zbeg-natte:zend+natte),
     1      pzz(xbeg-natte:xend+natte,zbeg-natte:zend+natte),
     1      pxz(xbeg-natte:xend+natte,zbeg-natte:zend+natte)
      real,parameter:: PI=3.14159265357257
c   staggered grid differential coefficients
      real,parameter:: a5_1=1.2112427,
     1                 a5_2=-0.089721680,
     1                 a5_3=0.013842773,
     1                 a5_4=-0.0017656599,
     1                 a5_5=0.00011867947
      real vp2,vs2,vpx2,vpn2,ep,de,cos2,sin2,sina,cosa,cossin,coef
      real c011,c013,c033,c055
      real c11,c13,c15,c33,c35,c55
      real dpx_x,dpz_z,dpxz_z,dpxz_x,dvx_x,dvz_z,dvx_z,dvz_x
      real dtx,dtz
      integer i,k,k_zbeg,k_zend,i_xbeg,i_xend

      dtx=dt/dx
      dtz=dt/dz

c----------------- define the computing area ----------------------------------
      if( inf == 4 )then                  !all area 
            k_zbeg=zbeg-natte+5
            k_zend=zend+natte-5
            i_xbeg=xbeg-natte+5
            i_xend=xend+natte-5
      elseif( inf == 3 )then                  !the up is free boundary 
            k_zbeg=zbeg
            k_zend=zend+natte-5
            i_xbeg=xbeg-natte+5
            i_xend=xend+natte-5
      elseif( inf == 0 )then                  !the around is free boundary 
            k_zbeg=zbeg
            k_zend=zend
            i_xbeg=xbeg
            i_xend=xend
      else 
            write(*,*)" The simulation subrountine(inf) must be 4,3,0 "
            stop
      end if

c----------------- compute velocity--------------------------------------------
      do k=k_zbeg,k_zend
      do i=i_xbeg,i_xend
      dpx_x= a5_1*(pxx(i+1,k)-pxx(i,k))   + a5_2*(pxx(i+2,k)-pxx(i-1,k))
     1      +a5_3*(pxx(i+3,k)-pxx(i-2,k)) + a5_4*(pxx(i+4,k)-pxx(i-3,k))
     2      +a5_5*(pxx(i+5,k)-pxx(i-4,k))
      dpxz_z=a5_1*(pxz(i,k)-pxz(i,k-1))   + a5_2*(pxz(i,k+1)-pxz(i,k-2))
     1      +a5_3*(pxz(i,k+2)-pxz(i,k-3)) + a5_4*(pxz(i,k+3)-pxz(i,k-4))
     2      +a5_5*(pxz(i,k+4)-pxz(i,k-5))

      vx(i,k)=vx(i,k)+(dtx*dpx_x+dtz*dpxz_z)

      dpz_z= a5_1*(pzz(i,k+1)-pzz(i,k))   + a5_2*(pzz(i,k+2)-pzz(i,k-1))
     1      +a5_3*(pzz(i,k+3)-pzz(i,k-2)) + a5_4*(pzz(i,k+4)-pzz(i,k-3))
     2      +a5_5*(pzz(i,k+5)-pzz(i,k-4))
      dpxz_x=a5_1*(pxz(i,k)-pxz(i-1,k))   + a5_2*(pxz(i+1,k)-pxz(i-2,k))
     1      +a5_3*(pxz(i+2,k)-pxz(i-3,k)) + a5_4*(pxz(i+3,k)-pxz(i-4,k))
     2      +a5_5*(pxz(i+4,k)-pxz(i-5,k))

      vz(i,k)=vz(i,k)+(dtz*dpz_z+dtx*dpxz_x)

      enddo
      enddo
                                        
c------------------- compute stress--------------------------------------------          
      do k=k_zbeg,k_zend
      do i=i_xbeg,i_xend
      dvx_x=a5_1*(vx(i,k)-vx(i-1,k))    +a5_2*(vx(i+1,k)-vx(i-2,k))
     1     +a5_3*(vx(i+2,k)-vx(i-3,k))  +a5_4*(vx(i+3,k)-vx(i-4,k))
     2     +a5_5*(vx(i+4,k)-vx(i-5,k))
      dvx_z=a5_1*(vx(i,k+1)-vx(i,k))    +a5_2*(vx(i,k+2)-vx(i,k-1))
     1     +a5_3*(vx(i,k+3)-vx(i,k-2))  +a5_4*(vx(i,k+4)-vx(i,k-3))
     2     +a5_5*(vx(i,k+5)-vx(i,k-4))
      dvz_z=a5_1*(vz(i,k)-vz(i,k-1))    +a5_2*(vz(i,k+1)-vz(i,k-2))
     1     +a5_3*(vz(i,k+2)-vz(i,k-3))  +a5_4*(vz(i,k+3)-vz(i,k-4))
     2     +a5_5*(vz(i,k+4)-vz(i,k-5))
      dvz_x=a5_1*(vz(i+1,k)-vz(i,k))    +a5_2*(vz(i+2,k)-vz(i-1,k))
     1     +a5_3*(vz(i+3,k)-vz(i-2,k))  +a5_4*(vz(i+4,k)-vz(i-3,k))
     2     +a5_5*(vz(i+5,k)-vz(i-4,k))      

         vp2=vp(i,k)**2
         vs2=vs(i,k)**2
         ep=1+2*Aepsilon(i,k)
         de=1+2*Adelta(i,k)
         vpx2=vp2*ep
         vpn2=vp2*de
         theta=Atheta(i,k)
         cos2=cos(theta*PI/180.0)*cos(theta*PI/180.0)
         sin2=sin(theta*PI/180.0)*sin(theta*PI/180.0)
         sina=sin(theta*PI/180.0)*cos(theta*PI/180.0)*2
         cosa=cos2-sin2
         c011=vpx2
         c033=vp2
         c055=vs2
         c013=sqrt((vp2-vs2)*(vpn2-vs2))-vs2
         c11= (cos2*c011+sin2*c013)*cos2+(cos2*c013+sin2*c033)*sin2
     +       +sina**2*c055
         c13= (cos2*c011+sin2*c013)*sin2+(cos2*c013+sin2*c033)*cos2
     +       -sina**2*c055
         c15=-0.5*(cos2*c011+sin2*c013)*sina
     +       +0.5*(cos2*c013+sin2*c033)*sina+sina*cosa*c055
         c33= (sin2*c011+cos2*c013)*sin2+(sin2*c013+cos2*c033)*cos2
     +       +sina**2*c055
         c35=-0.5*(sin2*c011+cos2*c013)*sina
     +       +0.5*(sin2*c013+cos2*c033)*sina-sina*cosa*c055
         c55= 0.25*sina**2*(c011-2*c013+c033)+cosa**2*c055

      pxx(i,k)= pxx(i,k)
     1         +dtx*dvx_x*c11+(dtx*dvz_x+dtz*dvx_z)*c15+dtz*dvz_z*c13
      pzz(i,k)= pzz(i,k)
     1         +dtx*dvx_x*c13+(dtx*dvz_x+dtz*dvx_z)*c35+dtz*dvz_z*c33
      pxz(i,k)=pxz(i,k)
     1         +dtx*dvx_x*c15+(dtx*dvz_x+dtz*dvx_z)*c55+dtz*dvz_z*c35

      enddo
      enddo 
         
      return
      end subroutine forwarding_tti_elastic_2d_staggered_10order
c******************   end of TTI elastic staggered grid   ******************
c
c
c************************************************************************
c******************attenuate boundaries(sponge absorbing)****************
c************************************************************************
      subroutine spongeatten_2d( xbeg,xend,zbeg,zend,natte,alpha,p_str )
c**natte should be set to 50 ,alpha =0.0001~0.00006       
      implicit none    
      integer xbeg,xend,zbeg,zend,natte
      real alpha,dis2,aa
      integer i,k     
      real::p_str(-natte+xbeg:xend+natte,-natte+zbeg:zend+natte) 
c--------- initial variable values--------------------------
      aa=0.0
c----------- damp left and right boundaries--------------------------
      do k=zbeg-natte,zend+natte

           do i=xbeg-1,xbeg-natte,-1
                dis2=(xbeg-i)**2
                aa=exp(-alpha*dis2)
                p_str(i,k)=p_str(i,k)*aa
           enddo
 
           do i=xend+1,xend+natte,1
                dis2=(i-xend)**2
                aa=exp(-alpha*dis2)
                p_str(i,k)=p_str(i,k)*aa
           enddo

      enddo
c------------------- damp middle boundaries-------------------
      do i=xbeg-natte,xend+natte

           do k=zbeg-1,zbeg-natte,-1
                dis2=(zbeg-k)**2
                aa=exp(-alpha*dis2)
                p_str(i,k)=p_str(i,k)*aa
           enddo

           do k=zend+1,zend+natte,1
                dis2=(k-zend)**2
                aa=exp(-alpha*dis2)
                p_str(i,k)=p_str(i,k)*aa
           enddo

      enddo

      return
      end subroutine spongeatten_2d
c**************************************************************************

c************************************************************************
c******************attenuate boundaries(sponge absorbing)****************
c************************************************************************
      subroutine spongeatten_2d_2
     1         ( xbeg,xend,zbeg,zend,natte,alpha,p_str )
c    used distance/natte to match the alpha        
c    natte=50 ,alpha =0.25 ~0.15
c    natte=100, alpha=0.1   is ok 
c    natte=200, alpha=0.05 is prefect
      implicit none    
      integer xbeg,xend,zbeg,zend,natte
      real alpha,dis2,aa
      integer i,k     
      real::p_str(-natte+xbeg:xend+natte,-natte+zbeg:zend+natte) 
c--------- initial variable values--------------------------
      aa=0
c----------- damp left and right boundaries--------------------------
      do k=zbeg-natte,zend+natte

           do i=xbeg-1,xbeg-natte,-1
                dis2=(xbeg-i)**2
            dis2=dis2/natte**2
                aa=exp(-alpha*dis2)
                p_str(i,k)=p_str(i,k)*aa
           enddo
 
           do i=xend+1,xend+natte,1
                dis2=(i-xend)**2
            dis2=dis2/natte**2
                aa=exp(-alpha*dis2)
                p_str(i,k)=p_str(i,k)*aa
           enddo

      enddo
c------------------- damp middle boundaries-------------------
      do i=xbeg-natte,xend+natte

           do k=zbeg-1,zbeg-natte,-1
                dis2=(zbeg-k)**2
            dis2=dis2/natte**2
                aa=exp(-alpha*dis2)
                p_str(i,k)=p_str(i,k)*aa
           enddo

           do k=zend+1,zend+natte,1
                dis2=(k-zend)**2
            dis2=dis2/natte**2
                aa=exp(-alpha*dis2)
                p_str(i,k)=p_str(i,k)*aa
           enddo

      enddo

      return
      end subroutine spongeatten_2d_2
c***************************************************************************
c***************************************************************************
c********************** Ricker wavelet *************************************
c***************************************************************************
        subroutine Rickerwavelet(fp,nlenth,dt,Ricker)
        
        implicit none
        real,parameter::pi=3.14159265
        integer fp,nlenth,nts,nts1,n
        real dt,b,t2
        real Ricker(nlenth)
        
        b=(pi*fp)**2    
        nts=nint(1./fp/dt)
        nts1=2*nts+1 
        do n=1,nts1
        t2=((n-1-nts)*dt)**2
        Ricker(n)=(1.-2.*b*t2)*exp(-b*t2)
        enddo 
        
        return
        end subroutine Rickerwavelet
c***************************************************************************

c***************************************************************************
c******************** Read in velocity *************************************
c***************************************************************************
        subroutine readinvelocity2d
     1     (vfile,xbeg,xend,zbeg,zend,mx,mz,natte,vp)
!       the xbeg ,xend should in <1,nx_model>   
        implicit none
        integer nc_file
        character(len=*) vfile
        integer xbeg,xend,zbeg,zend,mx,mz,natte
        integer i,k
        real vp(-natte+xbeg:xend+natte,-natte+zbeg:zend+natte)


        open(1,file=trim(vfile),form='unformatted',access='direct'
     1 ,recl=4*mz,status='old')

        do i=max(xbeg,1),min(xend,mx)
          read(1,rec=i,err=100) (vp(i,k),k=1,mz)
        enddo
        goto 200

100   write(*,*) "read velocity file fault: ",vfile  
       
200        close(1)
        
        return
        end subroutine readinvelocity2d
c*****************************************************************************

c**************************************************************************        
c****************** set boundaries (2d) ***********************************
c**************************************************************************
        subroutine setboundaries2d
     1     (xbeg,xend,zbeg,zend,mx,mz,natte,vp)

        implicit none
        integer xbeg,xend,zbeg,zend,mx,mz,natte
        real::vp(-natte+xbeg:xend+natte,-natte+zbeg:zend+natte) 
        integer k,i
c-------------------- set values around absorbing area--------------     
        do k=zbeg,zend
                do i=max(xbeg,1)-1,xbeg-natte,-1
                        vp(i,k)=vp( max(xbeg,1) , k)
                enddo

                do i=min(xend,mx)+1,xend+natte,1
                        vp(i,k)=vp( min(xend,mx) , k)
                enddo
        enddo

        do i=xbeg-natte,xend+natte
                do k=max(zbeg,1)-1,zbeg-natte,-1
                        vp(i,k)=vp( i , max(zbeg,1) )
                enddo
        
                do k=min(zend,mz)+1,zend+natte,1
                        vp(i,k)=vp( i , min(zend,mz) )
                enddo
        enddo

!        write(*,*) "set values around absorbing area finished"
        
        return 
        end subroutine setboundaries2d
c***************************************************************************

c*****************************************************************************
c******************* Assigning the tasks by mpi environment ******************
c*****************************************************************************
      subroutine assign_task_mpi
     1   (ntask,task_beg,p,me,itask_beg,itask_end)

      implicit none

      integer ntask,task_beg,p,me,itask_beg,itask_end
      integer ntask_com,ntask_res

      ntask_com=int(ntask/p)
      ntask_res=ntask-(ntask_com*p)

      if ( ntask_res == 0 )then
      itask_beg=me*ntask_com+task_beg+1
      itask_end=itask_beg+ntask_com-1
      endif 

      if ( ntask_res >= 0 )then 
       if ( me>=0 .and. me<=(ntask_res-1) )then
       itask_beg=me*(ntask_com+1)+task_beg+1
       itask_end=itask_beg+ntask_com
       endif
       if( me > (ntask_res-1) )then
       itask_beg=me*ntask_com+ntask_res+task_beg+1
       itask_end=itask_beg+ntask_com-1
       endif
      endif

      return
      end subroutine assign_task_mpi
c*****************************************************************************
c************************* initialize the tracehead ***************************
c******************************************************************************
      subroutine initialize_segy_tracehead( n_trace, tracehead )
      implicit none

      include 'segytype_for.h' 
      integer n_trace
      type(segy) :: tracehead(n_trace)
      
      tracehead=segy( 0,0,0,0,0,0,0,0,0,0,
     1                        0,0,0,0,0,0,0,0,0,0,
     2                        0,0,0,0,0,0,0,0,0,0,
     3                        0,0,0,0,0,0,0,0,0,0,
     4                        0,0,0,0,0,0,0,0,0,0,
     5                        0,0,0,0,0,0,0,0,0,0,
     6                        0,0,0,0,0,0,0,0,0,0,
     7                        0,0,0,0,0,0,0,0,0 )      
c                  the tracehead%unass(17) just be counted as 1 element          

!     tracehead%tracl=0      !      tracehead%tracr=0      !      tracehead%fldr=0      !      tracehead%tracf=0
!      tracehead%ep=0            !      tracehead%cdp=0            !      tracehead%cdpt=0      !      tracehead%trid=0
!      tracehead%nvs=0            !      tracehead%nhs=0            !      tracehead%duse=0      !      tracehead%offset=0
!      tracehead%gelev=0      !      tracehead%selev=0      !      tracehead%sdepth=0      !     tracehead%gdel=0
!      tracehead%sdel=0      !      tracehead%swdep=0      !      tracehead%gwdep=0      !      tracehead%scalel=0
!      tracehead%scalco=0      !      tracehead%sx=0            !      tracehead%sy=0            !      tracehead%gx=0
!      tracehead%gy=0            !      tracehead%counit=0      !      tracehead%wevel=0      !      tracehead%swevel=0
!      tracehead%sut=0            !      tracehead%gut=0            !      tracehead%sstat=0      !      tracehead%gstat=0
!      tracehead%tstat=0      !      tracehead%laga=0      !      tracehead%lagb=0      !      tracehead%delrt=0
!      tracehead%muts=0      !      tracehead%mute=0      !      tracehead%ns=0            !      tracehead%dt=0
!      tracehead%gain=0      !      tracehead%igc=0            !      tracehead%igi=0            !      tracehead%corr=0
!      tracehead%sfs=0            !      tracehead%sfe=0            !      tracehead%slen=0      !      tracehead%styp=0
!      tracehead%stas=0      !      tracehead%stae=0      !      tracehead%tatyp=0      !      tracehead%afilf=0
!      tracehead%afils=0      !      tracehead%nofilf=0      !      tracehead%nofils=0      !      tracehead%lcf=0
!      tracehead%hcf=0            !      tracehead%lcs=0            !      tracehead%hcs=0            !      tracehead%year=0
!      tracehead%day=0            !      tracehead%hour=0      !      tracehead%minute=0      !      tracehead%sec=0
!      tracehead%timbas=0      !      tracehead%trwf=0      !      tracehead%grnors=0      !      tracehead%grnofr=0
!      tracehead%grnlof=0      !      tracehead%gaps=0      !      tracehead%otrav=0      !      tracehead%d1=0
!      tracehead%f1=0            !      tracehead%d2=0            !      tracehead%f2=0            !      tracehead%ungpow=0
!      tracehead%unscale=0      !      tracehead%mark=0      !      tracehead%unass(1)=0!      tracehead%unass(2)=0
!      tracehead%unass(3)=0!      tracehead%unass(4)=0!      tracehead%unass(5)=0!      tracehead%unass(6)=0
!      tracehead%unass(7)=0!      tracehead%unass(8)=0!      tracehead%unass(9)=0!      tracehead%unass(10)=0
!      tracehead%unass(11)=0!      tracehead%unass(12)=0!      tracehead%unass(13)=0!      tracehead%unass(14)=0
!      tracehead%unass(15)=0!      tracehead%unass(16)=0!      tracehead%unass(17)=0

      return
      end subroutine initialize_segy_tracehead

c******************************************************************************

c******************************************************************************
c********************* elastic wave forwarding with staggered gird ************
c******************************* 10 order *************************************
      subroutine forwarding_elastic_2d_staggered_10order
     1  ( xbeg,xend,zbeg,zend,natte,vp,vs,
     2      vx,vz,pxx,pzz,pxz,dx,dz,dt,inf ) 
      implicit none

      integer xbeg,xend,zbeg,zend,natte
      real dx,dz,dt,dens
      integer inf
      real:: vp(xbeg-natte:xend+natte,zbeg-natte:zend+natte),
     1       vs(xbeg-natte:xend+natte,zbeg-natte:zend+natte),
     1       vx(xbeg-natte:xend+natte,zbeg-natte:zend+natte),
     1       vz(xbeg-natte:xend+natte,zbeg-natte:zend+natte),
     1      pxx(xbeg-natte:xend+natte,zbeg-natte:zend+natte),
     1      pzz(xbeg-natte:xend+natte,zbeg-natte:zend+natte),
     1      pxz(xbeg-natte:xend+natte,zbeg-natte:zend+natte)
      real,parameter:: a5_1=1.2112427,
     1                 a5_2=-0.089721680,
     1                 a5_3=0.013842773,
     1                 a5_4=-0.0017656599,
     1                 a5_5=0.00011867947
      real dpx_x,dpz_z,dpxz_z,dpxz_x,dvx_x,dvz_z,dvx_z,dvz_x
      real dtx,dtz
      integer i,k,k_zbeg,k_zend,i_xbeg,i_xend

      dtx=dt/dx
      dtz=dt/dz

c----------------- define the computing area ----------------------------------
      if( inf == 4 )then                  !all area 
            k_zbeg=zbeg-natte+5
            k_zend=zend+natte-5
            i_xbeg=xbeg-natte+5
            i_xend=xend+natte-5
      elseif( inf == 3 )then                  !the up is free boundary 
            k_zbeg=zbeg
            k_zend=zend+natte-5
            i_xbeg=xbeg-natte+5
            i_xend=xend+natte-5
      elseif( inf == 0 )then                  !the around is free boundary 
            k_zbeg=zbeg
            k_zend=zend
            i_xbeg=xbeg
            i_xend=xend
      else 
            write(*,*)" The simulation subrountine(inf) must be 4,3,0 "
            stop
      end if

c----------------- compute velocity--------------------------------------------
      do k=k_zbeg,k_zend
      do i=i_xbeg,i_xend
      dpx_x= a5_1*(pxx(i+1,k)-pxx(i,k))   + a5_2*(pxx(i+2,k)-pxx(i-1,k))
     1      +a5_3*(pxx(i+3,k)-pxx(i-2,k)) + a5_4*(pxx(i+4,k)-pxx(i-3,k))
     2      +a5_5*(pxx(i+5,k)-pxx(i-4,k))
      dpxz_z=a5_1*(pxz(i,k)-pxz(i,k-1))   + a5_2*(pxz(i,k+1)-pxz(i,k-2))
     1      +a5_3*(pxz(i,k+2)-pxz(i,k-3)) + a5_4*(pxz(i,k+3)-pxz(i,k-4))
     2      +a5_5*(pxz(i,k+4)-pxz(i,k-5))

      vx(i,k)=vx(i,k)+(dtx*dpx_x+dtz*dpxz_z)

      dpz_z= a5_1*(pzz(i,k+1)-pzz(i,k))   + a5_2*(pzz(i,k+2)-pzz(i,k-1))
     1      +a5_3*(pzz(i,k+3)-pzz(i,k-2)) + a5_4*(pzz(i,k+4)-pzz(i,k-3))
     2      +a5_5*(pzz(i,k+5)-pzz(i,k-4))
      dpxz_x=a5_1*(pxz(i,k)-pxz(i-1,k))   + a5_2*(pxz(i+1,k)-pxz(i-2,k))
     1      +a5_3*(pxz(i+2,k)-pxz(i-3,k)) + a5_4*(pxz(i+3,k)-pxz(i-4,k))
     2      +a5_5*(pxz(i+4,k)-pxz(i-5,k))

      vz(i,k)=vz(i,k)+(dtz*dpz_z+dtx*dpxz_x)

      enddo
      enddo
                                        
c------------------- compute stress--------------------------------------------          
!        do k=zbeg,zend+natte-5
      do k=k_zbeg,k_zend
      do i=i_xbeg,i_xend
      dvx_x=a5_1*(vx(i,k)-vx(i-1,k))    +a5_2*(vx(i+1,k)-vx(i-2,k))
     1     +a5_3*(vx(i+2,k)-vx(i-3,k))  +a5_4*(vx(i+3,k)-vx(i-4,k))
     2     +a5_5*(vx(i+4,k)-vx(i-5,k))
      dvx_z=a5_1*(vx(i,k+1)-vx(i,k))    +a5_2*(vx(i,k+2)-vx(i,k-1))
     1     +a5_3*(vx(i,k+3)-vx(i,k-2))  +a5_4*(vx(i,k+4)-vx(i,k-3))
     2     +a5_5*(vx(i,k+5)-vx(i,k-4))
      dvz_z=a5_1*(vz(i,k)-vz(i,k-1))    +a5_2*(vz(i,k+1)-vz(i,k-2))
     1     +a5_3*(vz(i,k+2)-vz(i,k-3))  +a5_4*(vz(i,k+3)-vz(i,k-4))
     2     +a5_5*(vz(i,k+4)-vz(i,k-5))
      dvz_x=a5_1*(vz(i+1,k)-vz(i,k))    +a5_2*(vz(i+2,k)-vz(i-1,k))
     1     +a5_3*(vz(i+3,k)-vz(i-2,k))  +a5_4*(vz(i+4,k)-vz(i-3,k))
     2     +a5_5*(vz(i+5,k)-vz(i-4,k))      

      pxx(i,k)=pxx(i,k)+(dtx*dvx_x+dtz*dvz_z)*vp(i,k)*vp(i,k)-
     1             (dtz*dvz_z)*2*vs(i,k)*vs(i,k)
      pzz(i,k)=pzz(i,k)+(dtx*dvx_x+dtz*dvz_z)*vp(i,k)*vp(i,k)-
     1             (dtx*dvx_x)*2*vs(i,k)*vs(i,k)
      pxz(i,k)=pxz(i,k)+(dtz*dvx_z+dtx*dvz_x)*vs(i,k)*vs(i,k)

      enddo
      enddo 
         
      return
      end subroutine forwarding_elastic_2d_staggered_10order
c******************************************************************************

c  save the 2D wave field snapshot
c
      subroutine snapshots2D(ishot,n,xbeg,xend,zbeg,zend,
     +                       mx,mz,natte,P,Pname)
      implicit none
      integer ishot,n,xbeg,xend,zbeg,zend,mx,my,mz,natte
      real:: P(xbeg-natte:xend+natte,zbeg-natte:zend+natte)
      character(*) Pname

      integer i,j,k,countkey
      character *5  tmp1,tmp2
      character *80 fn_su,snapname_xz
c  make a new dir to save the snapshot of eace shot (not ready!!)
c
      call read_shot_filename('snap/'//Pname,ishot,1,fn_su)
c      call system('mkdir $fn_su')
c
      if ((n.ge.1).and.(n.le.9)) then
          tmp1='0000'
          write(tmp2,'(I1)') n
          snapname_xz=trim(fn_su)//'_xz'//trim(tmp1)//tmp2
      else if((n.ge.10).and.(n.le.99)) then
          tmp1='000'
          write(tmp2,'(I2)') n
          snapname_xz=trim(fn_su)//'_xz'//trim(tmp1)//tmp2
      else if((n.ge.100).and.(n.le.999)) then
          tmp1='00'
          write(tmp2,'(I3)') n
          snapname_xz=trim(fn_su)//'_xz'//trim(tmp1)//tmp2
      else if((n.ge.1000).and.(n.le.9999)) then 
          tmp1='0'
          write(tmp2,'(I4)') n
          snapname_xz=trim(fn_su)//'_xz'//trim(tmp1)//tmp2
      else
          tmp1=''
          write(tmp2,'(I5)') n
          snapname_xz=trim(fn_su)//'_xz'//trim(tmp1)//tmp2
      endif

      open(77,file=snapname_xz,form='unformatted',access='direct',
     +     recl=mz*4, status='replace')
      countkey=1
      do i=1,mx
      if(i.lt.xbeg .or. i.gt.xend)then
         write(77,rec=countkey) (P(i,k)*0.0,k=zbeg,zend)
      else
         write(77,rec=countkey) (P(i,k),k=zbeg,zend)
      endif
         countkey=countkey+1
      enddo
      close(77)

      return
      end subroutine snapshots2D
c*************************************************************************

c*******************  read shot filename  ***********************
c
      subroutine read_shot_filename(ccc,number,key,outnam)

      character *80 outnam
      character (*) ccc
      character *80 outn
      character *80 outa
      character *1 outno
      js=number
      outno='0'

      if(key.eq.1) outa=trim(ccc)

      if(js.eq.0) then
        outn=trim(outa)//outno//outno//outno//outno//outno
      endif

      if((js.ge.1).and.(js.le.9)) then
        outn=trim(outa)//outno//outno//outno//outno//char(js+48)
      endif

      if((js.ge.10).and.(js.le.99)) then
        ia=int(js/10)
        ib=mod(js,10)
        outn=trim(outa)//outno//outno//outno//char(ia+48)//char(ib+48)
      endif

      if((js.ge.100).and.(js.le.999)) then
        ia=int(js/100)
        ib=mod(int(js/10),10)
        ic=mod(js,10)
        outn=trim(outa)//outno//outno//char(ia+48)//
     +                   char(ib+48)//char(ic+48)
      endif

      if((js.ge.1000).and.(js.le.9999)) then
        ia=int(js/1000)
        ib=mod(int(js/100),10)
        ic=mod(int(js/10),10)
        id=mod(js,10)
        outn=trim(outa)//outno//char(ia+48)//char(ib+48)//char(ic+48)
     1 //char(id+48)
      endif 

      outnam=trim(outn)

      return
      end
c****************************************************************************
c******************************************************************************
      subroutine unify_filename
     1           (filehead,nc_head,number,nc_number,outname)
      implicit none

      integer nc_head,number,nc_number
      character(len=nc_head)::filehead
      character(len=nc_number)::char_number
      character(len=nc_head+nc_number)::outname
      integer i
     
!      write(char_number,"(I5.5)") number
      
      do i=1,nc_number 
         char_number(i:i)=
     1       char( mod( int( number/(10**(nc_number-i)) ),10 )+48 )
      enddo
 
      outname=filehead//char_number

      return 
      end subroutine unify_filename
c******************************************************************************
