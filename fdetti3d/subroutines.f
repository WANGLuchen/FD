cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c***********   TTI elastic 3D forwarding staggered grid   ******
c
      !  this is a simple version
      subroutine forwarding_tti_elastic_3d_staggered_10order_2
     1  ( xbeg,xend,ybeg,yend,zbeg,zend,natte,
     1    c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,
     1    c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66,
     1    vx,vy,vz,pxx,pyy,pzz,pxy,pxz,pyz,dx,dy,dz,dt,inf )
      implicit none

      integer xbeg,xend,ybeg,yend,zbeg,zend,natte
      real theta,phi,dx,dy,dz,dt,dens
      integer inf
      real:: 
     1       vx(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1       vy(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1       vz(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1      pxx(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1      pyy(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1      pzz(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1      pxy(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1      pxz(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1      pyz(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte)
      real::
     1      c11(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1      c12(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1      c13(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1      c14(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1      c15(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1      c16(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1      c22(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1      c23(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1      c24(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1      c25(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1      c26(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1      c33(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1      c34(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1      c35(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1      c36(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1      c44(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1      c45(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1      c46(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1      c55(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1      c56(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1      c66(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte)

      real,parameter:: PI=3.14159265357257
c   staggered grid differential coefficients
      real,parameter:: a5_1=1.2112427,
     1                 a5_2=-0.089721680,
     1                 a5_3=0.013842773,
     1                 a5_4=-0.0017656599,
     1                 a5_5=0.00011867947

      real dpxx_x,dpxy_y,dpxz_z,
     1     dpxy_x,dpyy_y,dpyz_z,
     1     dpxz_x,dpyz_y,dpzz_z,
     2     dvx_x,dvy_y,dvz_z,
     2     dvy_z,dvz_y,dvx_z,dvz_x,dvx_y,dvy_x
      real dtx,dty,dtz
      integer i,j,k,k_zbeg,k_zend,j_ybeg,j_yend,i_xbeg,i_xend

      dtx=dt/dx
      dty=dt/dy
      dtz=dt/dz

c----------------- define the computing area ----------------------------------
      if( inf == 4 )then                  !all area 
            k_zbeg=zbeg-natte+5
            k_zend=zend+natte-5
            j_ybeg=ybeg-natte+5
            j_yend=yend+natte-5
            i_xbeg=xbeg-natte+5
            i_xend=xend+natte-5
      elseif( inf == 3 )then                  !the up is free boundary 
            k_zbeg=zbeg
            k_zend=zend+natte-5
            j_ybeg=ybeg-natte+5
            j_yend=yend+natte-5
            i_xbeg=xbeg-natte+5
            i_xend=xend+natte-5
      elseif( inf == 0 )then                  !the around is free boundary 
            k_zbeg=zbeg
            k_zend=zend
            j_ybeg=ybeg
            j_yend=yend
            i_xbeg=xbeg
            i_xend=xend
      else 
            write(*,*)" The simulation subrountine(inf) must be 4,3,0 "
            stop
      end if

c----------------- compute velocity--------------------------------------------
      do k=k_zbeg,k_zend
      do j=j_ybeg,j_yend
      do i=i_xbeg,i_xend
      dpxx_x= 
     1 a5_1*(pxx(i  ,j,k)-pxx(i-1,j,k))+a5_2*(pxx(i+1,j,k)-pxx(i-2,j,k))
     2+a5_3*(pxx(i+2,j,k)-pxx(i-3,j,k))+a5_4*(pxx(i+3,j,k)-pxx(i-4,j,k))
     3+a5_5*(pxx(i+4,j,k)-pxx(i-5,j,k))
      dpxy_y= 
     1 a5_1*(pxy(i,j+1,k)-pxy(i,j  ,k))+a5_2*(pxy(i,j+2,k)-pxy(i,j-1,k))
     2+a5_3*(pxy(i,j+3,k)-pxy(i,j-2,k))+a5_4*(pxy(i,j+4,k)-pxy(i,j-3,k))
     3+a5_5*(pxy(i,j+5,k)-pxy(i,j-4,k))
      dpxz_z= 
     1 a5_1*(pxz(i,j,k  )-pxz(i,j,k-1))+a5_2*(pxz(i,j,k+1)-pxz(i,j,k-2))
     2+a5_3*(pxz(i,j,k+2)-pxz(i,j,k-3))+a5_4*(pxz(i,j,k+3)-pxz(i,j,k-4))
     3+a5_5*(pxz(i,j,k+4)-pxz(i,j,k-5))

      vx(i,j,k)=vx(i,j,k)+(dtx*dpxx_x+dty*dpxy_y+dtz*dpxz_z)

      dpxy_x= 
     1 a5_1*(pxy(i+1,j,k)-pxy(i  ,j,k))+a5_2*(pxy(i+2,j,k)-pxy(i-1,j,k))
     2+a5_3*(pxy(i+3,j,k)-pxy(i-2,j,k))+a5_4*(pxy(i+4,j,k)-pxy(i-3,j,k))
     3+a5_5*(pxy(i+5,j,k)-pxy(i-4,j,k))
      dpyy_y= 
     1 a5_1*(pyy(i,j  ,k)-pyy(i,j-1,k))+a5_2*(pyy(i,j+1,k)-pyy(i,j-2,k))
     2+a5_3*(pyy(i,j+2,k)-pyy(i,j-3,k))+a5_4*(pyy(i,j+3,k)-pyy(i,j-4,k))
     3+a5_5*(pyy(i,j+4,k)-pyy(i,j-5,k))
      dpyz_z= 
     1 a5_1*(pyz(i,j,k  )-pyz(i,j,k-1))+a5_2*(pyz(i,j,k+1)-pyz(i,j,k-2))
     2+a5_3*(pyz(i,j,k+2)-pyz(i,j,k-3))+a5_4*(pyz(i,j,k+3)-pyz(i,j,k-4))
     3+a5_5*(pyz(i,j,k+4)-pyz(i,j,k-5))

      vy(i,j,k)=vy(i,j,k)+(dtx*dpxy_x+dty*dpyy_y+dtz*dpyz_z)

      dpxz_x= 
     1 a5_1*(pxz(i+1,j,k)-pxz(i  ,j,k))+a5_2*(pxz(i+2,j,k)-pxz(i-1,j,k))
     2+a5_3*(pxz(i+3,j,k)-pxz(i-2,j,k))+a5_4*(pxz(i+4,j,k)-pxz(i-3,j,k))
     3+a5_5*(pxz(i+5,j,k)-pxz(i-4,j,k))
      dpyz_y= 
     1 a5_1*(pyz(i,j+1,k)-pyz(i,j  ,k))+a5_2*(pyz(i,j+2,k)-pyz(i,j-1,k))
     2+a5_3*(pyz(i,j+3,k)-pyz(i,j-2,k))+a5_4*(pyz(i,j+4,k)-pyz(i,j-3,k))
     3+a5_5*(pyz(i,j+5,k)-pyz(i,j-4,k))
      dpzz_z= 
     1 a5_1*(pzz(i,j,k+1)-pzz(i,j,k  ))+a5_2*(pzz(i,j,k+2)-pzz(i,j,k-1))
     2+a5_3*(pzz(i,j,k+3)-pzz(i,j,k-2))+a5_4*(pzz(i,j,k+4)-pzz(i,j,k-3))
     3+a5_5*(pzz(i,j,k+5)-pzz(i,j,k-4))

      vz(i,j,k)=vz(i,j,k)+(dtx*dpxz_x+dty*dpyz_y+dtz*dpzz_z)

      enddo
      enddo
      enddo
                                        
c------------------- compute stress--------------------------------------------          
      do k=k_zbeg,k_zend
      do j=j_ybeg,j_yend
      do i=i_xbeg,i_xend
      dvx_x=
     1 a5_1*(vx(i+1,j,k)-vx(i  ,j,k))+a5_2*(vx(i+2,j,k)-vx(i-1,j,k))
     2+a5_3*(vx(i+3,j,k)-vx(i-2,j,k))+a5_4*(vx(i+4,j,k)-vx(i-3,j,k))
     3+a5_5*(vx(i+5,j,k)-vx(i-4,j,k)) 
      dvy_y=
     1 a5_1*(vy(i,j+1,k)-vy(i,j  ,k))+a5_2*(vy(i,j+2,k)-vy(i,j-1,k))
     2+a5_3*(vy(i,j+3,k)-vy(i,j-2,k))+a5_4*(vy(i,j+4,k)-vy(i,j-3,k))
     3+a5_5*(vy(i,j+5,k)-vy(i,j-4,k))
      dvz_z=
     1 a5_1*(vz(i,j,k  )-vz(i,j,k-1))+a5_2*(vz(i,j,k+1)-vz(i,j,k-2))
     2+a5_3*(vz(i,j,k+2)-vz(i,j,k-3))+a5_4*(vz(i,j,k+3)-vz(i,j,k-4))
     3+a5_5*(vz(i,j,k+4)-vz(i,j,k-5))

      dvy_z=
     1 a5_1*(vy(i,j,k+1)-vy(i,j,k  ))+a5_2*(vy(i,j,k+2)-vy(i,j,k-1))
     2+a5_3*(vy(i,j,k+3)-vy(i,j,k-2))+a5_4*(vy(i,j,k+4)-vy(i,j,k-3))
     3+a5_5*(vy(i,j,k+5)-vy(i,j,k-4))
      dvz_y=
     1 a5_1*(vz(i,j  ,k)-vz(i,j-1,k))+a5_2*(vz(i,j+1,k)-vz(i,j-2,k))
     2+a5_3*(vz(i,j+2,k)-vz(i,j-3,k))+a5_4*(vz(i,j+3,k)-vz(i,j-4,k))
     3+a5_5*(vz(i,j+4,k)-vz(i,j-5,k))

      dvx_z=
     1 a5_1*(vx(i,j,k+1)-vx(i,j,k  ))+a5_2*(vx(i,j,k+2)-vx(i,j,k-1))
     2+a5_3*(vx(i,j,k+3)-vx(i,j,k-2))+a5_4*(vx(i,j,k+4)-vx(i,j,k-3))
     3+a5_5*(vx(i,j,k+5)-vx(i,j,k-4))
      dvz_x=
     1 a5_1*(vz(i  ,j,k)-vz(i-1,j,k))+a5_2*(vz(i+1,j,k)-vz(i-2,j,k))
     2+a5_3*(vz(i+2,j,k)-vz(i-3,j,k))+a5_4*(vz(i+3,j,k)-vz(i-4,j,k))
     3+a5_5*(vz(i+4,j,k)-vz(i-5,j,k))

      dvx_y=
     1 a5_1*(vx(i,j  ,k)-vx(i,j-1,k))+a5_2*(vx(i,j+1,k)-vx(i,j-2,k))
     2+a5_3*(vx(i,j+2,k)-vx(i,j-3,k))+a5_4*(vx(i,j+3,k)-vx(i,j-4,k))
     3+a5_5*(vx(i,j+4,k)-vx(i,j-5,k))
      dvy_x=
     1 a5_1*(vy(i  ,j,k)-vy(i-1,j,k))+a5_2*(vy(i+1,j,k)-vy(i-2,j,k))
     2+a5_3*(vy(i+2,j,k)-vy(i-3,j,k))+a5_4*(vy(i+3,j,k)-vy(i-4,j,k))
     3+a5_5*(vy(i+4,j,k)-vy(i-5,j,k))
     
c      c11=26390000.0
c      c12=1271000.0
c      c13=6110000.0
c      c14=   0.0000000    
c      c15=  -0.0000000    
c      c16=   0.0000000    
c      c22=   26390000.    
c      c23=   6110000.0    
c      c24=  -0.0000000    
c      c25=   0.0000000    
c      c26=   0.0000000    
c      c33=   15600000.    
c      c34=   0.0000000    
c      c35=   0.0000000    
c      c36=   0.0000000    
c      c44=   4380000.3    
c      c45=   0.0000000    
c      c46=   0.0000000    
c      c55=   4380000.3    
c      c56=   0.0000000    
c      c66=   6840000.0

c      c11=15155200.0
c      c12=4914600.0
c      c13=4369863.0
c      c14=   0.0000000    
c      c15=  -0.0000000    
c      c16=   0.0000000    
c      c22=   15155200.    
c      c23=   4369863.0    
c      c24=  -0.0000000    
c      c25=   0.0000000    
c      c26=   0.0000000    
c      c33=   10240000.    
c      c34=   0.0000000    
c      c35=   0.0000000    
c      c36=   0.0000000    
c      c44=   3413533.3    
c      c45=   0.0000000    
c      c46=   0.0000000    
c      c55=   3413533.3    
c      c56=   0.0000000    
c      c66=   5120300.0

      pxx(i,j,k)= pxx(i,j,k)
     1   +dtx*dvx_x*c11(i,j,k) + (dtz*dvy_z+dty*dvz_y)*c14(i,j,k)
     1   +dty*dvy_y*c12(i,j,k) + (dtx*dvz_x+dtz*dvx_z)*c15(i,j,k)
     1   +dtz*dvz_z*c13(i,j,k) + (dtx*dvy_x+dty*dvx_y)*c16(i,j,k)
      pyy(i,j,k)= pyy(i,j,k)
     1   +dtx*dvx_x*c12(i,j,k) + (dtz*dvy_z+dty*dvz_y)*c24(i,j,k)
     1   +dty*dvy_y*c22(i,j,k) + (dtx*dvz_x+dtz*dvx_z)*c25(i,j,k)
     1   +dtz*dvz_z*c23(i,j,k) + (dtx*dvy_x+dty*dvx_y)*c26(i,j,k)
      pzz(i,j,k)= pzz(i,j,k)
     1   +dtx*dvx_x*c13(i,j,k) + (dtz*dvy_z+dty*dvz_y)*c34(i,j,k)
     1   +dty*dvy_y*c23(i,j,k) + (dtx*dvz_x+dtz*dvx_z)*c35(i,j,k)
     1   +dtz*dvz_z*c33(i,j,k) + (dtx*dvy_x+dty*dvx_y)*c36(i,j,k)
      pyz(i,j,k)= pyz(i,j,k)
     1   +dtx*dvx_x*c14(i,j,k) + (dtz*dvy_z+dty*dvz_y)*c44(i,j,k)
     1   +dty*dvy_y*c24(i,j,k) + (dtx*dvz_x+dtz*dvx_z)*c45(i,j,k)
     1   +dtz*dvz_z*c34(i,j,k) + (dtx*dvy_x+dty*dvx_y)*c46(i,j,k)
      pxz(i,j,k)= pxz(i,j,k)
     1   +dtx*dvx_x*c15(i,j,k) + (dtz*dvy_z+dty*dvz_y)*c45(i,j,k)
     1   +dty*dvy_y*c25(i,j,k) + (dtx*dvz_x+dtz*dvx_z)*c55(i,j,k)
     1   +dtz*dvz_z*c35(i,j,k) + (dtx*dvy_x+dty*dvx_y)*c56(i,j,k)
      pxy(i,j,k)= pxy(i,j,k)
     1   +dtx*dvx_x*c16(i,j,k) + (dtz*dvy_z+dty*dvz_y)*c46(i,j,k)
     1   +dty*dvy_y*c26(i,j,k) + (dtx*dvz_x+dtz*dvx_z)*c56(i,j,k)
     1   +dtz*dvz_z*c36(i,j,k) + (dtx*dvy_x+dty*dvx_y)*c66(i,j,k)

      enddo
      enddo 
      enddo 
         
      return
      end subroutine forwarding_tti_elastic_3d_staggered_10order_2

      subroutine forwarding_tti_elastic_3d_staggered_10order
     1  ( xbeg,xend,ybeg,yend,zbeg,zend,natte,
     1    vp,vs,Aepsilon,Adelta,Agamma,Atheta,Aphi,
     1    vx,vy,vz,pxx,pyy,pzz,pxy,pxz,pyz,dx,dy,dz,dt,inf )
      implicit none

      integer xbeg,xend,ybeg,yend,zbeg,zend,natte
      real theta,phi,dx,dy,dz,dt,dens
      integer inf
      real:: 
     1       vp(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1       vs(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1 Aepsilon(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1   Adelta(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1   Agamma(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1   Atheta(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1     Aphi(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1       vx(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1       vy(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1       vz(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1      pxx(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1      pyy(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1      pzz(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1      pxy(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1      pxz(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1      pyz(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte)
      real,parameter:: PI=3.14159265357257
c   staggered grid differential coefficients
      real,parameter:: a5_1=1.2112427,
     1                 a5_2=-0.089721680,
     1                 a5_3=0.013842773,
     1                 a5_4=-0.0017656599,
     1                 a5_5=0.00011867947
      real vp2,vs2,vpx2,vpn2,ep,de,ga,
     +     costhe, sinthe, cosphi, sinphi,
     +     sin2the,cos2the,sin2phi,cos2phi

      real c011,c012,c013,c022,c023,c033,c044,c055,c066
      real c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,
     +     c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66
      real dpxx_x,dpxy_y,dpxz_z,
     1     dpxy_x,dpyy_y,dpyz_z,
     1     dpxz_x,dpyz_y,dpzz_z,
     2     dvx_x,dvy_y,dvz_z,
     2     dvy_z,dvz_y,dvx_z,dvz_x,dvx_y,dvy_x
      real dtx,dty,dtz
      integer i,j,k,k_zbeg,k_zend,j_ybeg,j_yend,i_xbeg,i_xend

      dtx=dt/dx
      dty=dt/dy
      dtz=dt/dz

c----------------- define the computing area ----------------------------------
      if( inf == 4 )then                  !all area 
            k_zbeg=zbeg-natte+5
            k_zend=zend+natte-5
            j_ybeg=ybeg-natte+5
            j_yend=yend+natte-5
            i_xbeg=xbeg-natte+5
            i_xend=xend+natte-5
      elseif( inf == 3 )then                  !the up is free boundary 
            k_zbeg=zbeg
            k_zend=zend+natte-5
            j_ybeg=ybeg-natte+5
            j_yend=yend+natte-5
            i_xbeg=xbeg-natte+5
            i_xend=xend+natte-5
      elseif( inf == 0 )then                  !the around is free boundary 
            k_zbeg=zbeg
            k_zend=zend
            j_ybeg=ybeg
            j_yend=yend
            i_xbeg=xbeg
            i_xend=xend
      else 
            write(*,*)" The simulation subrountine(inf) must be 4,3,0 "
            stop
      end if

c----------------- compute velocity--------------------------------------------
      do k=k_zbeg,k_zend
      do j=j_ybeg,j_yend
      do i=i_xbeg,i_xend
      dpxx_x= 
     1 a5_1*(pxx(i  ,j,k)-pxx(i-1,j,k))+a5_2*(pxx(i+1,j,k)-pxx(i-2,j,k))
     2+a5_3*(pxx(i+2,j,k)-pxx(i-3,j,k))+a5_4*(pxx(i+3,j,k)-pxx(i-4,j,k))
     3+a5_5*(pxx(i+4,j,k)-pxx(i-5,j,k))
      dpxy_y= 
     1 a5_1*(pxy(i,j+1,k)-pxy(i,j  ,k))+a5_2*(pxy(i,j+2,k)-pxy(i,j-1,k))
     2+a5_3*(pxy(i,j+3,k)-pxy(i,j-2,k))+a5_4*(pxy(i,j+4,k)-pxy(i,j-3,k))
     3+a5_5*(pxy(i,j+5,k)-pxy(i,j-4,k))
      dpxz_z= 
     1 a5_1*(pxz(i,j,k  )-pxz(i,j,k-1))+a5_2*(pxz(i,j,k+1)-pxz(i,j,k-2))
     2+a5_3*(pxz(i,j,k+2)-pxz(i,j,k-3))+a5_4*(pxz(i,j,k+3)-pxz(i,j,k-4))
     3+a5_5*(pxz(i,j,k+4)-pxz(i,j,k-5))

      vx(i,j,k)=vx(i,j,k)+(dtx*dpxx_x+dty*dpxy_y+dtz*dpxz_z)

      dpxy_x= 
     1 a5_1*(pxy(i+1,j,k)-pxy(i  ,j,k))+a5_2*(pxy(i+2,j,k)-pxy(i-1,j,k))
     2+a5_3*(pxy(i+3,j,k)-pxy(i-2,j,k))+a5_4*(pxy(i+4,j,k)-pxy(i-3,j,k))
     3+a5_5*(pxy(i+5,j,k)-pxy(i-4,j,k))
      dpyy_y= 
     1 a5_1*(pyy(i,j  ,k)-pyy(i,j-1,k))+a5_2*(pyy(i,j+1,k)-pyy(i,j-2,k))
     2+a5_3*(pyy(i,j+2,k)-pyy(i,j-3,k))+a5_4*(pyy(i,j+3,k)-pyy(i,j-4,k))
     3+a5_5*(pyy(i,j+4,k)-pyy(i,j-5,k))
      dpyz_z= 
     1 a5_1*(pyz(i,j,k  )-pyz(i,j,k-1))+a5_2*(pyz(i,j,k+1)-pyz(i,j,k-2))
     2+a5_3*(pyz(i,j,k+2)-pyz(i,j,k-3))+a5_4*(pyz(i,j,k+3)-pyz(i,j,k-4))
     3+a5_5*(pyz(i,j,k+4)-pyz(i,j,k-5))

      vy(i,j,k)=vy(i,j,k)+(dtx*dpxy_x+dty*dpyy_y+dtz*dpyz_z)

      dpxz_x= 
     1 a5_1*(pxz(i+1,j,k)-pxz(i  ,j,k))+a5_2*(pxz(i+2,j,k)-pxz(i-1,j,k))
     2+a5_3*(pxz(i+3,j,k)-pxz(i-2,j,k))+a5_4*(pxz(i+4,j,k)-pxz(i-3,j,k))
     3+a5_5*(pxz(i+5,j,k)-pxz(i-4,j,k))
      dpyz_y= 
     1 a5_1*(pyz(i,j+1,k)-pyz(i,j  ,k))+a5_2*(pyz(i,j+2,k)-pyz(i,j-1,k))
     2+a5_3*(pyz(i,j+3,k)-pyz(i,j-2,k))+a5_4*(pyz(i,j+4,k)-pyz(i,j-3,k))
     3+a5_5*(pyz(i,j+5,k)-pyz(i,j-4,k))
      dpzz_z= 
     1 a5_1*(pzz(i,j,k+1)-pzz(i,j,k  ))+a5_2*(pzz(i,j,k+2)-pzz(i,j,k-1))
     2+a5_3*(pzz(i,j,k+3)-pzz(i,j,k-2))+a5_4*(pzz(i,j,k+4)-pzz(i,j,k-3))
     3+a5_5*(pzz(i,j,k+5)-pzz(i,j,k-4))

      vz(i,j,k)=vz(i,j,k)+(dtx*dpxz_x+dty*dpyz_y+dtz*dpzz_z)

      enddo
      enddo
      enddo
                                        
c------------------- compute stress--------------------------------------------          
      do k=k_zbeg,k_zend
      do j=j_ybeg,j_yend
      do i=i_xbeg,i_xend
      dvx_x=
     1 a5_1*(vx(i+1,j,k)-vx(i  ,j,k))+a5_2*(vx(i+2,j,k)-vx(i-1,j,k))
     2+a5_3*(vx(i+3,j,k)-vx(i-2,j,k))+a5_4*(vx(i+4,j,k)-vx(i-3,j,k))
     3+a5_5*(vx(i+5,j,k)-vx(i-4,j,k)) 
      dvy_y=
     1 a5_1*(vy(i,j+1,k)-vy(i,j  ,k))+a5_2*(vy(i,j+2,k)-vy(i,j-1,k))
     2+a5_3*(vy(i,j+3,k)-vy(i,j-2,k))+a5_4*(vy(i,j+4,k)-vy(i,j-3,k))
     3+a5_5*(vy(i,j+5,k)-vy(i,j-4,k))
      dvz_z=
     1 a5_1*(vz(i,j,k  )-vz(i,j,k-1))+a5_2*(vz(i,j,k+1)-vz(i,j,k-2))
     2+a5_3*(vz(i,j,k+2)-vz(i,j,k-3))+a5_4*(vz(i,j,k+3)-vz(i,j,k-4))
     3+a5_5*(vz(i,j,k+4)-vz(i,j,k-5))

      dvy_z=
     1 a5_1*(vy(i,j,k+1)-vy(i,j,k  ))+a5_2*(vy(i,j,k+2)-vy(i,j,k-1))
     2+a5_3*(vy(i,j,k+3)-vy(i,j,k-2))+a5_4*(vy(i,j,k+4)-vy(i,j,k-3))
     3+a5_5*(vy(i,j,k+5)-vy(i,j,k-4))
      dvz_y=
     1 a5_1*(vz(i,j  ,k)-vz(i,j-1,k))+a5_2*(vz(i,j+1,k)-vz(i,j-2,k))
     2+a5_3*(vz(i,j+2,k)-vz(i,j-3,k))+a5_4*(vz(i,j+3,k)-vz(i,j-4,k))
     3+a5_5*(vz(i,j+4,k)-vz(i,j-5,k))

      dvx_z=
     1 a5_1*(vx(i,j,k+1)-vx(i,j,k  ))+a5_2*(vx(i,j,k+2)-vx(i,j,k-1))
     2+a5_3*(vx(i,j,k+3)-vx(i,j,k-2))+a5_4*(vx(i,j,k+4)-vx(i,j,k-3))
     3+a5_5*(vx(i,j,k+5)-vx(i,j,k-4))
      dvz_x=
     1 a5_1*(vz(i  ,j,k)-vz(i-1,j,k))+a5_2*(vz(i+1,j,k)-vz(i-2,j,k))
     2+a5_3*(vz(i+2,j,k)-vz(i-3,j,k))+a5_4*(vz(i+3,j,k)-vz(i-4,j,k))
     3+a5_5*(vz(i+4,j,k)-vz(i-5,j,k))

      dvx_y=
     1 a5_1*(vx(i,j  ,k)-vx(i,j-1,k))+a5_2*(vx(i,j+1,k)-vx(i,j-2,k))
     2+a5_3*(vx(i,j+2,k)-vx(i,j-3,k))+a5_4*(vx(i,j+3,k)-vx(i,j-4,k))
     3+a5_5*(vx(i,j+4,k)-vx(i,j-5,k))
      dvy_x=
     1 a5_1*(vy(i  ,j,k)-vy(i-1,j,k))+a5_2*(vy(i+1,j,k)-vy(i-2,j,k))
     2+a5_3*(vy(i+2,j,k)-vy(i-3,j,k))+a5_4*(vy(i+3,j,k)-vy(i-4,j,k))
     3+a5_5*(vy(i+4,j,k)-vy(i-5,j,k))

         vp2=vp(i,j,k)**2
         vs2=vs(i,j,k)**2
         ep=1+2*Aepsilon(i,j,k)
         de=1+2*Adelta(i,j,k)
         ga=1+2*Agamma(i,j,k)
         vpx2=vp2*ep
         vpn2=vp2*de

           theta=Atheta(i,j,k)
           phi=Aphi(i,j,k)
         costhe=cos(theta*PI/180.0)
         sinthe=sin(theta*PI/180.0)
         cosphi=cos(phi*PI/180.0)
         sinphi=sin(phi*PI/180.0)
         sin2the=2.0*sinthe*costhe
         cos2the=costhe**2-sinthe**2
         sin2phi=2.0*sinphi*cosphi
         cos2phi=cosphi**2-sinphi**2

         c011=vpx2
         c022=c011
         c033=vp2
         c055=vs2
         c044=c055
         c066=vs2*ga
         c012=c011-2*c066
         c013=sqrt((vp2-vs2)*(vpn2-vs2))-vs2
         c023=c013

      c11= sin2phi*(c066*sin2phi*costhe**2 + c044*sin2phi*sinthe**2) 
     1    +cosphi**2*(costhe**2*(c011*cosphi**2*costhe**2 
     1    +c013*cosphi**2*sinthe**2 + c012*sinphi**2) 
     1    +sinthe**2*(c013*cosphi**2*costhe**2
     1    +c033*cosphi**2*sinthe**2 + c023*sinphi**2) 
     1    +c055*cosphi**2*sin2the**2) 
     1    +sinphi**2*(c012*cosphi**2*costhe**2 
     2    +c023*cosphi**2*sinthe**2 + c022*sinphi**2)

      c12= sinphi**2*(costhe**2*(c011*cosphi**2*costhe**2 
     2    +c013*cosphi**2*sinthe**2 + c012*sinphi**2) 
     2    +sinthe**2*(c013*cosphi**2*costhe**2 
     2    +c033*cosphi**2*sinthe**2 + c023*sinphi**2) 
     2    +c055*cosphi**2*sin2the**2) -sin2phi*(c066*sin2phi*costhe**2
     2    +c044*sin2phi*sinthe**2)+cosphi**2*(c012*cosphi**2*costhe**2
     2    +c023*cosphi**2*sinthe**2 + c022*sinphi**2)

      c13= costhe**2*(c013*cosphi**2*costhe**2 
     2    +c033*cosphi**2*sinthe**2 + c023*sinphi**2)
     2    +sinthe**2*(c011*cosphi**2*costhe**2 
     2    +c013*cosphi**2*sinthe**2 + c012*sinphi**2) 
     2    -c055*cosphi**2*sin2the**2

      c14= cosphi*(c044*costhe*sinthe*sin2phi 
     4    -c066*costhe*sinthe*sin2phi) - sinphi*(0.5*sin2the
     4    *(c013*cosphi**2*costhe**2 + c033*cosphi**2*sinthe**2 
     4    +c023*sinphi**2) - 0.5*sin2the*(c011*cosphi**2*costhe**2 
     4    +c013*cosphi**2*sinthe**2 + c012*sinphi**2) 
     4    +c055*cosphi**2*cos2the*sin2the)

      c15=-sinphi*(c044*costhe*sinthe*sin2phi 
     5    -c066*costhe*sinthe*sin2phi) - cosphi*(0.5*sin2the
     5    *(c013*cosphi**2*costhe**2 + c033*cosphi**2*sinthe**2 
     5    +c023*sinphi**2) - 0.5*sin2the*(c011*cosphi**2*costhe**2 
     5    +c013*cosphi**2*sinthe**2 + c012*sinphi**2) 
     5    +c055*cosphi**2*cos2the*sin2the)

      c16= 0.5*sin2phi*(costhe**2*(c011*cosphi**2*costhe**2 
     1    +c013*cosphi**2*sinthe**2 + c012*sinphi**2) + sinthe**2
     1    *(c013*cosphi**2*costhe**2 + c033*cosphi**2*sinthe**2 
     1    +c023*sinphi**2) + c055*cosphi**2*sin2the**2) - cos2phi
     1    *(c066*sin2phi*costhe**2 + c044*sin2phi*sinthe**2) 
     1    -0.5*sin2phi*(c012*cosphi**2*costhe**2 
     1    +c023*cosphi**2*sinthe**2 + c022*sinphi**2)

      c22= sin2phi*(c066*sin2phi*costhe**2 + c044*sin2phi*sinthe**2)
     1    +sinphi**2*(costhe**2*(c012*cosphi**2 
     1    +c011*costhe**2*sinphi**2 + c013*sinphi**2*sinthe**2) 
     1    +sinthe**2*(c023*cosphi**2 + c013*costhe**2*sinphi**2 
     1    +c033*sinphi**2*sinthe**2) + c055*sinphi**2*sin2the**2) 
     1    +cosphi**2*(c022*cosphi**2 + c012*costhe**2*sinphi**2
     1    +c023*sinphi**2*sinthe**2)

      c23=  costhe**2*(c023*cosphi**2 + c013*costhe**2*sinphi**2 
     1     +c033*sinphi**2*sinthe**2) + sinthe**2*(c012*cosphi**2 
     1     +c011*costhe**2*sinphi**2 + c013*sinphi**2*sinthe**2) 
     1     -c055*sinphi**2*sin2the**2

      c24=-cosphi*(c044*costhe*sinthe*sin2phi 
     1    -c066*costhe*sinthe*sin2phi) - sinphi*(0.5*sin2the
     1    *(c023*cosphi**2 + c013*costhe**2*sinphi**2 
     1    +c033*sinphi**2*sinthe**2) - 0.5*sin2the*(c012*cosphi**2
     1    +c011*costhe**2*sinphi**2 + c013*sinphi**2*sinthe**2) 
     1    +c055*cos2the*sinphi**2*sin2the)

      c25= sinphi*(c044*costhe*sinthe*sin2phi 
     1    -c066*costhe*sinthe*sin2phi) - cosphi*(0.5*sin2the
     1    *(c023*cosphi**2 + c013*costhe**2*sinphi**2 
     1    +c033*sinphi**2*sinthe**2) - 0.5*sin2the*(c012*cosphi**2 
     1    +c011*costhe**2*sinphi**2 + c013*sinphi**2*sinthe**2) 
     1    +c055*cos2the*sinphi**2*sin2the)

      c26= cos2phi*(c066*sin2phi*costhe**2 + c044*sin2phi*sinthe**2) 
     1    +0.5*sin2phi*(costhe**2*(c012*cosphi**2 
     1    +c011*costhe**2*sinphi**2 + c013*sinphi**2*sinthe**2) 
     1    +sinthe**2*(c023*cosphi**2 + c013*costhe**2*sinphi**2 
     1    +c033*sinphi**2*sinthe**2) + c055*sinphi**2*sin2the**2) 
     1    -0.5*sin2phi*(c022*cosphi**2 + c012*costhe**2*sinphi**2 
     1    +c023*sinphi**2*sinthe**2)

      c33= c055*sin2the**2 +costhe**2*(c033*costhe**2+c013*sinthe**2) 
     1    +sinthe**2*(c013*costhe**2 + c011*sinthe**2)

      c34= sinphi*(0.5*sin2the*(c013*costhe**2 + c011*sinthe**2) 
     1    -0.5*sin2the*(c033*costhe**2 + c013*sinthe**2) 
     1    +c055*cos2the*sin2the)

      c35= cosphi*(0.5*sin2the*(c013*costhe**2 + c011*sinthe**2) 
     1    -0.5*sin2the*(c033*costhe**2 + c013*sinthe**2) 
     1    +c055*cos2the*sin2the)

      c36= 0.5*sin2phi*(costhe**2*(c013*costhe**2 + c011*sinthe**2) 
     1    -c055*sin2the**2 +sinthe**2*(c033*costhe**2+c013*sinthe**2))
     1    -0.5*sin2phi*(c023*costhe**2 + c012*sinthe**2)

      c44= cosphi*(c044*cosphi*costhe**2 + c066*cosphi*sinthe**2) 
     1    +sinphi*(c055*sinphi*cos2the**2 
     1    +0.5*sin2the*(0.5*c011*sinphi*sin2the 
     1    -0.5*c013*sinphi*sin2the) - 0.5*sin2the
     1    *(0.5*c013*sinphi*sin2the - 0.5*c033*sinphi*sin2the))

      c45= cosphi*(c055*sinphi*cos2the**2 + 0.5*sin2the
     1    *(0.5*c011*sinphi*sin2the - 0.5*c013*sinphi*sin2the) 
     1    -0.5*sin2the*(0.5*c013*sinphi*sin2the
     1    -0.5*c033*sinphi*sin2the)) - sinphi*(c044*cosphi*costhe**2 
     1    +c066*cosphi*sinthe**2)

      c46= 0.5*sin2phi*((0.5*c011*sinphi*sin2the 
     2    -0.5*c013*sinphi*sin2the)*costhe**2+(0.5*c013*sinphi*sin2the
     2    -0.5*c033*sinphi*sin2the)*sinthe**2 
     2    -c055*cos2the*sinphi*sin2the) - 0.5*sin2phi
     2    *(0.5*c012*sinphi*sin2the - 0.5*c023*sinphi*sin2the)-cos2phi
     2    *(c044*cosphi*costhe*sinthe-c066*cosphi*costhe*sinthe)

      c55= sinphi*(c044*sinphi*costhe**2 + c066*sinphi*sinthe**2) 
     3    +cosphi*(c055*cosphi*cos2the**2 + 0.5*sin2the
     3    *(0.5*c011*cosphi*sin2the - 0.5*c013*cosphi*sin2the)
     3    -0.5*sin2the*(0.5*c013*cosphi*sin2the 
     3    -0.5*c033*cosphi*sin2the))

      c56= cos2phi*(c044*costhe*sinphi*sinthe 
     5    -c066*costhe*sinphi*sinthe) - 0.5*sin2phi
     5    *(0.5*c012*cosphi*sin2the - 0.5*c023*cosphi*sin2the) 
     5    +0.5*sin2phi*(0.5*(c011*cosphi*sin2the-c013*cosphi*sin2the)
     5    *costhe**2 +(0.5*c013*cosphi*sin2the 
     5    -0.5*c033*cosphi*sin2the)*sinthe**2 
     5    -c055*cosphi*cos2the*sin2the)

      c66= cos2phi*(c066*cos2phi*costhe**2 + c044*cos2phi*sinthe**2)
     3    +0.5*sin2phi*(sinthe**2*(0.5*c013*sin2phi*costhe**2 
     3    +0.5*c033*sin2phi*sinthe**2 - 0.5*c023*sin2phi) 
     3    +costhe**2*(0.5*c011*sin2phi*costhe**2 
     3    +0.5*c013*sin2phi*sinthe**2 - 0.5*c012*sin2phi) 
     3    +0.5*c055*sin2phi*sin2the**2) - 0.5*sin2phi
     3    *(0.5*c012*sin2phi*costhe**2 + 0.5*c023*sin2phi*sinthe**2
     3    -0.5*c022*sin2phi)

      write(*,*)'c11=',c11
      write(*,*)'c12=',c12
      write(*,*)'c13=',c13
      write(*,*)'c14=',c14
      write(*,*)'c15=',c15
      write(*,*)'c16=',c16
      write(*,*)'c22=',c22
      write(*,*)'c23=',c23
      write(*,*)'c24=',c24
      write(*,*)'c25=',c25
      write(*,*)'c26=',c26
      write(*,*)'c33=',c33
      write(*,*)'c34=',c34
      write(*,*)'c35=',c35
      write(*,*)'c36=',c36
      write(*,*)'c44=',c44
      write(*,*)'c45=',c45
      write(*,*)'c46=',c46
      write(*,*)'c55=',c55
      write(*,*)'c56=',c56
      write(*,*)'c66=',c66
      stop

      pxx(i,j,k)= pxx(i,j,k)
     1           +dtx*dvx_x*c11 + (dtz*dvy_z+dty*dvz_y)*c14
     1           +dty*dvy_y*c12 + (dtx*dvz_x+dtz*dvx_z)*c15
     1           +dtz*dvz_z*c13 + (dtx*dvy_x+dty*dvx_y)*c16
      pyy(i,j,k)= pyy(i,j,k)
     1           +dtx*dvx_x*c12 + (dtz*dvy_z+dty*dvz_y)*c24
     1           +dty*dvy_y*c22 + (dtx*dvz_x+dtz*dvx_z)*c25
     1           +dtz*dvz_z*c23 + (dtx*dvy_x+dty*dvx_y)*c26
      pzz(i,j,k)= pzz(i,j,k)
     1           +dtx*dvx_x*c13 + (dtz*dvy_z+dty*dvz_y)*c34
     1           +dty*dvy_y*c23 + (dtx*dvz_x+dtz*dvx_z)*c35
     1           +dtz*dvz_z*c33 + (dtx*dvy_x+dty*dvx_y)*c36
      pyz(i,j,k)= pyz(i,j,k)
     1           +dtx*dvx_x*c14 + (dtz*dvy_z+dty*dvz_y)*c44
     1           +dty*dvy_y*c24 + (dtx*dvz_x+dtz*dvx_z)*c45
     1           +dtz*dvz_z*c34 + (dtx*dvy_x+dty*dvx_y)*c46
      pxz(i,j,k)= pxz(i,j,k)
     1           +dtx*dvx_x*c15 + (dtz*dvy_z+dty*dvz_y)*c45
     1           +dty*dvy_y*c25 + (dtx*dvz_x+dtz*dvx_z)*c55
     1           +dtz*dvz_z*c35 + (dtx*dvy_x+dty*dvx_y)*c56
      pxy(i,j,k)= pxy(i,j,k)
     1           +dtx*dvx_x*c16 + (dtz*dvy_z+dty*dvz_y)*c46
     1           +dty*dvy_y*c26 + (dtx*dvz_x+dtz*dvx_z)*c56
     1           +dtz*dvz_z*c36 + (dtx*dvy_x+dty*dvx_y)*c66

      enddo
      enddo 
      enddo 
         
      return
      end subroutine forwarding_tti_elastic_3d_staggered_10order
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       the sequence of velocity bin data must be x->y->z   
        subroutine input_vel3D_xyz(vfile,xbeg,xend,
     +                 ybeg,yend,zbeg,zend,mx,my,mz,natte,vel)
        implicit none

        character(len=*) vfile
        integer xbeg,xend,ybeg,yend,zbeg,zend,mx,my,mz,natte
        integer xlexpand,xrexpand,ylexpand,yrexpand
        integer i,j,k
        real::
     1      vel(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     1         zbeg-natte:zend+natte)
        real,allocatable::tmpvel(:,:)

        allocate( tmpvel(mx,my) )

        xlexpand=max(1-xbeg,0)
        xrexpand=max(xend-mx,0)
        ylexpand=max(1-ybeg,0)
        yrexpand=max(yend-my,0)
        write(*,*) " velocity model file: ",vfile

        open(20,file=trim(vfile),form='unformatted',access='direct',
     1  recl=4*mx*my,status='old')

       do k=1,mz
          tmpvel=0.0
          read(20,rec=k,err=100) ((tmpvel(i,j),j=1,my),i=1,mx)
            do j=ybeg+ylexpand,yend-yrexpand
            do i=xbeg+xlexpand,xend-xrexpand
               vel(i,j,zbeg+k-1)=tmpvel(i,j)
            enddo
            enddo
       enddo

        goto 200
100   write(*,*) "read velocity file fault: ",vfile  
200   close(20)
      deallocate( tmpvel )

      return
      end subroutine input_vel3D_xyz
c
c       the sequence of velocity bin data must be z->x->y   
        subroutine input_vel3D_zxy(vfile,xbeg,xend,
     +                 ybeg,yend,zbeg,zend,mx,my,mz,natte,vel)
        implicit none

        character(len=*) vfile
        integer xbeg,xend,ybeg,yend,zbeg,zend,mx,my,mz,natte
        integer i,j,k
        integer nc_file
        real::
     1      vel(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     1         zbeg-natte:zend+natte)

        write(*,*) " velocity model file: ",vfile

        open(20,file=trim(vfile),form='unformatted',access='direct'
     1 ,recl=4*mz,status='old')

        do j=max(ybeg,1),min(yend,my)
          do i=max(xbeg,1),min(xend,mx)
                nc_file=(j-1)*mx+i
             read(20,rec=nc_file,err=100) (vel(i,j,k),k=zbeg,zend)
          enddo
        enddo
        goto 200

100   write(*,*) "read velocity file fault: ",vfile  
       
200   close(20)
        
      return
      end subroutine input_vel3D_zxy
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    output the 3d array data to a binary file
c
c     required-----"wfield/" subdirectories
c
      subroutine outputbin3D(n,vfile,xbeg,xend,
     +              ybeg,yend,zbeg,zend,mx,my,mz,natte,ps)
      implicit none

      character(len=*) vfile
        integer n,xbeg,xend,ybeg,yend,zbeg,zend,mx,my,mz,natte
        real::
     1      ps(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     1         zbeg-natte:zend+natte)

        integer i,j,k
        integer nc_file
        character *4  tmp1,tmp2
        character *7  header
        character *80 outname

      header='wfield/'

      if ((n.ge.1).and.(n.le.9)) then
          tmp1='0000'
          write(tmp2,'(I1)') n
          outname=trim(header)//trim(vfile)//trim(tmp1)//tmp2
      else if((n.ge.10).and.(n.le.99)) then
          tmp1='000'
          write(tmp2,'(I2)') n
          outname=trim(header)//trim(vfile)//trim(tmp1)//tmp2
      else if((n.ge.100).and.(n.le.999)) then
          tmp1='00'
          write(tmp2,'(I3)') n
          outname=trim(header)//trim(vfile)//trim(tmp1)//tmp2
      else 
          tmp1='0'
          write(tmp2,'(I4)') n
          outname=trim(header)//trim(vfile)//trim(tmp1)//tmp2
      endif

        open(20,file=trim(outname)//'.bin',form='unformatted',
     1          access='direct',recl=4*mz,status='replace')

        do j=max(ybeg,1),min(yend,my)
          do i=max(xbeg,1),min(xend,mx)
                nc_file=(j-1)*mx+i
             write(20,rec=nc_file) (ps(i,j,k),k=zbeg,zend)
          enddo
        enddo

        close(20)
        return
      end subroutine outputbin3D
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  the outputsu3 is used to output 3D data
      subroutine outputsu3(fn_su,
     +                    n_sample,n_trace,nsi,dx,dy,dt,
     +                    ishot,newxs,newys,newxr,newyr,sei_ps)
      implicit none
      include 'segy_simulation.h'

      integer   i,k,gdx,gdy
      integer   n_sample,n_trace,nsi,ishot,newxs,newys
      integer:: newxr(1:n_trace),newyr(1:n_trace)
      real      dx,dy,dt
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
         sy=(newys-1)*dy
         gx=(newxr(i)-1)*dx
         gy=(newyr(i)-1)*dy
         offset=((gx-sx)**2+(gy-sy)**2)**0.5
         if(i .eq. 1)then
             gdx=newxr(2)-newxr(1)
             gdy=newyr(2)-newyr(1)
         else
             gdx=newxr(i+1)-newxr(i)
             gdy=newyr(i+1)-newyr(i)
         endif
         d2=((gdx*dx)**2+(gdy*dy)**2)**0.5
         
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
       end subroutine outputsu3
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  shot filename
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
      end subroutine read_shot_filename
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Creating a traditional Ricker source wavelet
        subroutine Rickerwavelet2(fp,nlenth,dt,Ricker)
        
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
        end subroutine Rickerwavelet2
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  set boundaries (3d)
      subroutine setboundaries3d
     +           (xbeg,xend,ybeg,yend,zbeg,zend,mx,my,mz,natte,vp)
        implicit none
        integer xbeg,xend,ybeg,yend,zbeg,zend,mx,my,mz,natte
        real::
     1      vp(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     1         zbeg-natte:zend+natte)
        integer k,j,i
c-------------------- set values around absorbing area--------------     
        do k=max(zbeg,1)-1,zbeg-natte,-1                ! the upper
        do j=ybeg,yend
        do i=xbeg,xend
           vp(i,j,k)=vp(i,j,max(zbeg,1))
        enddo
        enddo
        enddo

        do k=min(zend,mz)+1,zend+natte,1                ! the bottom
        do j=ybeg,yend
        do i=xbeg,xend
           vp(i,j,k)=vp(i,j,min(zend,mz))
        enddo
        enddo
        enddo

        do k=zbeg-natte,zend+natte                      ! the left
        do j=ybeg,yend 
        do i=max(xbeg,1)-1,xbeg-natte,-1                    
           vp(i,j,k)=vp(max(xbeg,1),j,k)
        enddo
        enddo
        enddo
      
        do k=zbeg-natte,zend+natte                      ! the right  
        do j=ybeg,yend
        do i=min(xend,mx)+1,xend+natte,1                               
           vp(i,j,k)=vp(min(xend,mx),j,k)
        enddo
        enddo
        enddo

        do k=zbeg-natte,zend+natte                      ! the front
        do j=max(ybeg,1)-1,ybeg-natte,-1                        
        do i=xbeg-natte,xend+natte
           vp(i,j,k)=vp(i,max(ybeg,1),k) 
        enddo
        enddo
        enddo

        do k=zbeg-natte,zend+natte                      ! the back
        do j=min(yend,my)+1,yend+natte                      
        do i=xbeg-natte,xend+natte
           vp(i,j,k)=vp(i,min(yend,my),k) 
        enddo
        enddo
        enddo

!        write(*,*) "set values around absorbing area finished"
        
      return 
      end subroutine setboundaries3d
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  save the 3D wave field snapshot
c
c     required-----"snap/" subdirectories
c
      subroutine snapshots3D(ishot,n,xbeg,xend,ybeg,yend,zbeg,zend,
     +                       mx,my,mz,natte,newys,newzs,P,Pname)
      implicit none
      integer ishot,n,xbeg,xend,ybeg,yend,zbeg,zend,mx,my,mz,natte
      integer newys,newzs
      real:: P(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +         zbeg-natte:zend+natte)
      character(*) Pname

      integer i,j,k,countkey
      character *4  tmp1,tmp2
      character *80 fn_su,snapname_xz,snapname_xy
c  make a new dir to save the snapshot of eace shot (not ready!!)
c
      call read_shot_filename('snap/'//Pname,ishot,1,fn_su)
c      call system('mkdir $fn_su')
c
      if ((n.ge.1).and.(n.le.9)) then
          tmp1='0000'
          write(tmp2,'(I1)') n
          snapname_xz=trim(fn_su)//'_xz'//trim(tmp1)//tmp2
          snapname_xy=trim(fn_su)//'_xy'//trim(tmp1)//tmp2
      else if((n.ge.10).and.(n.le.99)) then
          tmp1='000'
          write(tmp2,'(I2)') n
          snapname_xz=trim(fn_su)//'_xz'//trim(tmp1)//tmp2
          snapname_xy=trim(fn_su)//'_xy'//trim(tmp1)//tmp2
      else if((n.ge.100).and.(n.le.999)) then
          tmp1='00'
          write(tmp2,'(I3)') n
          snapname_xz=trim(fn_su)//'_xz'//trim(tmp1)//tmp2
          snapname_xy=trim(fn_su)//'_xy'//trim(tmp1)//tmp2
      else 
          tmp1='0'
          write(tmp2,'(I4)') n
          snapname_xz=trim(fn_su)//'_xz'//trim(tmp1)//tmp2
          snapname_xy=trim(fn_su)//'_xy'//trim(tmp1)//tmp2
      endif

      open(77,file=snapname_xz,form='unformatted',access='direct',
     +     recl=(zend-zbeg+1)*4, status='replace')
      countkey=1
      do i=xbeg,xend
         write(77,rec=countkey) (P(i,newys,k),k=zbeg,zend)
         countkey=countkey+1
      enddo
      close(77)

      open(33,file=snapname_xy,form='unformatted',access='direct',
     +     recl=(yend-ybeg+1)*4, status='replace')
      countkey=1
      do i=xbeg,xend
         write(33,rec=countkey) (P(i,j,newzs+natte),j=ybeg,yend)
         countkey=countkey+1
      enddo
      close(33)


      return
      end subroutine snapshots3D
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  add source with focal mechanism
      subroutine source_loading(newxs,newys,newzs,wtime,
     +           xbeg,xend,ybeg,yend,zbeg,zend,natte,nt,dx,dy,dz,dt,
     +           mt,pxx,pyy,pzz,pyz,pxz,pxy,vx,vy,vz,n,style )
    !!!  if style == 1, load source on pxx,pyy,...
    !!!  if style == 2, load source on vx,vy,vz
      implicit none

      integer newxs,newys,newzs,n
      integer xbeg,xend,ybeg,yend,zbeg,zend,natte,nt
      integer style
      real mt(3,3),dx,dy,dz,dt
      real ::
     1       wtime(1:nt),
     1       vx(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1       vy(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1       vz(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1      pxx(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1      pyy(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1      pzz(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1      pxy(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1      pxz(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1      pyz(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte)

      !  style==1, source on stress
         pxx(newxs,newys,newzs)= pxx(newxs,newys,newzs)
     +                          -wtime(n)*mt(1,1)
         pyy(newxs,newys,newzs)= pyy(newxs,newys,newzs)
     +                          -wtime(n)*mt(2,2)
         pzz(newxs,newys,newzs)= pzz(newxs,newys,newzs)
     +                          -wtime(n)*mt(3,3)
         pyz(newxs,newys,newzs)= pyz(newxs,newys,newzs)
     +                          -wtime(n)*mt(2,3)
         pxy(newxs,newys,newzs)= pxy(newxs,newys,newzs)
     +                          -wtime(n)*mt(1,2)
         pxz(newxs,newys,newzs)= pxz(newxs,newys,newzs)
     +                          -wtime(n)*mt(1,3)


      !  style==2, source on velocity
cc  ISO
c      vx(newxs+1,newys,newzs)= vx(newxs+1,newys,newzs)
c     +                        +( -wtime(n)*mt(1,1) )
c      vx(newxs,newys,newzs)= vx(newxs,newys,newzs)
c     +                        -( -wtime(n)*mt(1,1) )
c      vy(newxs,newys+1,newzs)= vy(newxs,newys+1,newzs)
c     +                        +( -wtime(n)*mt(2,2) )
c      vy(newxs,newys,newzs)= vy(newxs,newys,newzs)
c     +                        -( -wtime(n)*mt(2,2) )
c      vz(newxs,newys,newzs+1)= vz(newxs,newys,newzs+1)
c     +                        +( -wtime(n)*mt(3,3) )
c      vz(newxs,newys,newzs)= vz(newxs,newys,newzs)
c     +                        -( -wtime(n)*mt(3,3) )
cc  DC
c      vx(newxs+1,newys,newzs)= vx(newxs+1,newys,newzs)
c     +                        -( -wtime(n)*mt(1,3) )
c      vx(newxs+1,newys,newzs-1)= vx(newxs+1,newys,newzs-1)
c     +                        +( -wtime(n)*mt(1,3) )
c      vz(newxs,newys,newzs-1)= vz(newxs+1,newys,newzs-1)
c     +                        +( -wtime(n)*mt(1,3) )
c      vz(newxs+1,newys,newzs-1)= vz(newxs+1,newys,newzs-1)
c     +                        -( -wtime(n)*mt(1,3) )
cc  CLVD
c      vx(newxs+1,newys,newzs)= vx(newxs+1,newys,newzs)
c     +                        +( -wtime(n)*mt(1,1) )
c      vx(newxs,newys,newzs)= vx(newxs,newys,newzs)
c     +                        -( -wtime(n)*mt(1,1) )
c      vy(newxs,newys+1,newzs)= vy(newxs,newys+1,newzs)
c     +                        +( -wtime(n)*mt(2,2)*(-2) )
c      vy(newxs,newys,newzs)= vy(newxs,newys,newzs)
c     +                        -( -wtime(n)*mt(2,2)*(-2) )
c      vz(newxs,newys,newzs+1)= vz(newxs,newys,newzs+1)
c     +                        +( -wtime(n)*mt(3,3) )
c      vz(newxs,newys,newzs)= vz(newxs,newys,newzs)
c     +                        -( -wtime(n)*mt(3,3) )

      !  style==3, source on stress + velocity
c         pxx(newxs,newys,newzs)= !pxx(newxs,newys,newzs)
c     +                          -wtime(n)*mt(1,1)/dt*dx
c         pyy(newxs,newys,newzs)= !pyy(newxs,newys,newzs)
c     +                          -wtime(n)*mt(2,2)/dt*dy
c         pzz(newxs,newys,newzs)= !pzz(newxs,newys,newzs)
c     +                          -wtime(n)*mt(3,3)/dt*dz
cc  DC
c      vx(newxs,newys,newzs+1)= vx(newxs,newys,newzs+1)
c     +                        -( -wtime(n)*mt(1,3) )
c      vx(newxs,newys,newzs)  = vx(newxs,newys,newzs)
c     +                        +( -wtime(n)*mt(1,3) )
c      vz(newxs+1,newys,newzs)= vz(newxs+1,newys,newzs)
c     +                        +( -wtime(n)*mt(1,3) )
c      vz(newxs,newys,newzs)= vz(newxs,newys,newzs)
c     +                        -( -wtime(n)*mt(1,3) )

      end subroutine source_loading
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   3d attenuate boundaries(sponge absorbing)
      subroutine spongeatten_3d(xbeg,xend,ybeg,yend,zbeg,zend,
     1                          natte,alpha,p_str)    
      real alpha,dis,aa
      integer xbeg,xend,ybeg,yend,zbeg,zend,natte,i,j,k     
      real::p_str(xbeg-natte:xend+natte,
     1            ybeg-natte:yend+natte,zbeg-natte:zend+natte)    
c--------- initial variable values--------------------------
      aa=0
      dis=0.0
c---------attenuate left and right boundery ----------------------------
        
        do k=zbeg-natte,zend+natte
          do j=ybeg-natte,yend+natte

            do i=xbeg-1,xbeg-natte,-1  !---------attenuate left boundery 
             dis=xbeg-i
             aa=exp(-alpha*dis*dis)
             p_str(i,j,k)=p_str(i,j,k)*aa
            end do
                                   
            do i=xend+1,xend+natte,1   !---------attenuate right boundery 
             dis=i-xend
             aa=exp(-alpha*dis*dis)
             p_str(i,j,k)=p_str(i,j,k)*aa
            end do

          end do
        end do

c---------attenuate front and bebind boundery ----------------------------
        
        do k=zbeg-natte,zend+natte

          do j=ybeg-1,ybeg-natte,-1 !---------attenuate front boundery 
             dis=ybeg-j
            do i=xbeg-natte,xend+natte 
             aa=exp(-alpha*dis*dis)
             p_str(i,j,k)=p_str(i,j,k)*aa
            end do
          enddo
                                   
          do j=yend+1,yend+natte    !---------attenuate bebind boundery 
             dis=j-yend
            do i=xbeg-natte,xend+natte 
             aa=exp(-alpha*dis*dis)
             p_str(i,j,k)=p_str(i,j,k)*aa
            end do
          end do

        end do

c---------attenuate up and down boundery ----------------------------
        
        do k=zbeg-1,zbeg-natte,-1 !---------attenuate up boundery 
             dis=zbeg-k
          do j=ybeg-natte,yend+natte
            do i=xbeg-natte,xend+natte 
             aa=exp(-alpha*dis*dis)
             p_str(i,j,k)=p_str(i,j,k)*aa
            end do
          end do
         end do
                                   
        do k=zend+1,zend+natte,1 !---------attenuate up boundery 
             dis=k-zend
          do j=ybeg-natte,yend+natte
            do i=xbeg-natte,xend+natte 
             aa=exp(-alpha*dis*dis)
             p_str(i,j,k)=p_str(i,j,k)*aa
            end do
          end do
         end do


      return
      end subroutine spongeatten_3d

      subroutine spongeatten_3d_2(xbeg,xend,ybeg,yend,zbeg,zend,
     1                          natte,alpha,p_str)
c    used distance/natte to match the alpha        
c    natte=50 ,alpha =0.25 ~0.15
c    natte=100, alpha=0.1   is ok 
c    natte=200, alpha=0.05 is prefect
      real alpha,dis2,aa
      integer xbeg,xend,ybeg,yend,zbeg,zend,natte,i,j,k     
      real::p_str(xbeg-natte:xend+natte,
     1            ybeg-natte:yend+natte,zbeg-natte:zend+natte)    
c--------- initial variable values--------------------------
      aa=0
      dis2=0.0
c---------attenuate left and right boundery ----------------------------
        
        do k=zbeg-natte,zend+natte
          do j=ybeg-natte,yend+natte

            do i=xbeg-1,xbeg-natte,-1  !---------attenuate left boundery 
             dis2=(xbeg-i)**2
             dis2=dis2/natte**2
             aa=exp(-alpha*dis2)
             p_str(i,j,k)=p_str(i,j,k)*aa
            end do
                                   
            do i=xend+1,xend+natte,1   !---------attenuate right boundery 
             dis2=(i-xend)**2
             dis2=dis2/natte**2
             aa=exp(-alpha*dis2)
             p_str(i,j,k)=p_str(i,j,k)*aa
            end do

          end do
        end do

c---------attenuate front and bebind boundery ----------------------------
        
        do k=zbeg-natte,zend+natte

          do j=ybeg-1,ybeg-natte,-1 !---------attenuate front boundery 
             dis2=(ybeg-j)**2
             dis2=dis2/natte**2
            do i=xbeg-natte,xend+natte 
             aa=exp(-alpha*dis2)
             p_str(i,j,k)=p_str(i,j,k)*aa
            end do
          enddo
                                   
          do j=yend+1,yend+natte    !---------attenuate bebind boundery 
             dis2=(j-yend)**2
             dis2=dis2/natte**2
            do i=xbeg-natte,xend+natte 
             aa=exp(-alpha*dis2)
             p_str(i,j,k)=p_str(i,j,k)*aa
            end do
          end do

        end do

c---------attenuate up and down boundery ----------------------------
        
        do k=zbeg-1,zbeg-natte,-1 !---------attenuate up boundery 
             dis2=(zbeg-k)**2
             dis2=dis2/natte**2
          do j=ybeg-natte,yend+natte
            do i=xbeg-natte,xend+natte 
             aa=exp(-alpha*dis2)
             p_str(i,j,k)=p_str(i,j,k)*aa
            end do
          end do
         end do
                                   
        do k=zend+1,zend+natte,1 !---------attenuate up boundery 
             dis2=(k-zend)**2
             dis2=dis2/natte**2
          do j=ybeg-natte,yend+natte
            do i=xbeg-natte,xend+natte 
             aa=exp(-alpha*dis2)
             p_str(i,j,k)=p_str(i,j,k)*aa
            end do
          end do
         end do


      return
      end subroutine spongeatten_3d_2
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   calculating the coefficients of stiffness matrix
      subroutine stiffness_matrix_coefficient
     1  ( xbeg,xend,ybeg,yend,zbeg,zend,natte,
     1    vp,vs,Aepsilon,Adelta,Agamma,Atheta,Aphi,
     1    c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,
     1    c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66 )
      implicit none

      integer xbeg,xend,ybeg,yend,zbeg,zend,natte
      real theta,phi,dx,dy,dz,dt,dens
      integer inf
      real:: 
     1       vp(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1       vs(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1 Aepsilon(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1   Adelta(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1   Agamma(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1   Atheta(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1     Aphi(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte)
      real::
     1      c11(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1      c12(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1      c13(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1      c14(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1      c15(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1      c16(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1      c22(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1      c23(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1      c24(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1      c25(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1      c26(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1      c33(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1      c34(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1      c35(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1      c36(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1      c44(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1      c45(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1      c46(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1      c55(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1      c56(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte),
     1      c66(xbeg-natte:xend+natte,ybeg-natte:yend+natte,
     +          zbeg-natte:zend+natte)

      real vp2,vs2,vpx2,vpn2,ep,de,ga,
     +     costhe, sinthe, cosphi, sinphi,
     +     sin2the,cos2the,sin2phi,cos2phi

      real c011,c012,c013,c022,c023,c033,c044,c055,c066
      integer i,j,k,k_zbeg,k_zend,j_ybeg,j_yend,i_xbeg,i_xend

      real,parameter:: PI=3.14159265357257

        k_zbeg=zbeg-natte
        k_zend=zend+natte
        j_ybeg=ybeg-natte
        j_yend=yend+natte
        i_xbeg=xbeg-natte
        i_xend=xend+natte

      do k=k_zbeg,k_zend
      do j=j_ybeg,j_yend
      do i=i_xbeg,i_xend

         vp2=vp(i,j,k)**2
         vs2=vs(i,j,k)**2
         ep=1+2*Aepsilon(i,j,k)
         de=1+2*Adelta(i,j,k)
         ga=1+2*Agamma(i,j,k)
         vpx2=vp2*ep
         vpn2=vp2*de

           theta=Atheta(i,j,k)
           phi=Aphi(i,j,k)
         costhe=cos(theta*PI/180.0)
         sinthe=sin(theta*PI/180.0)
         cosphi=cos(phi*PI/180.0)
         sinphi=sin(phi*PI/180.0)
         sin2the=2.0*sinthe*costhe
         cos2the=costhe**2-sinthe**2
         sin2phi=2.0*sinphi*cosphi
         cos2phi=cosphi**2-sinphi**2

         c011=vpx2
         c022=c011
         c033=vp2
         c055=vs2
         c044=c055
         c066=vs2*ga
         c012=c011-2*c066
         c013=sqrt((vp2-vs2)*(vpn2-vs2))-vs2
         c023=c013

      c11(i,j,k)= 
     1     sin2phi*(c066*sin2phi*costhe**2 + c044*sin2phi*sinthe**2) 
     1    +cosphi**2*(costhe**2*(c011*cosphi**2*costhe**2 
     1    +c013*cosphi**2*sinthe**2 + c012*sinphi**2) 
     1    +sinthe**2*(c013*cosphi**2*costhe**2
     1    +c033*cosphi**2*sinthe**2 + c023*sinphi**2) 
     1    +c055*cosphi**2*sin2the**2) 
     1    +sinphi**2*(c012*cosphi**2*costhe**2 
     2    +c023*cosphi**2*sinthe**2 + c022*sinphi**2)

      c12(i,j,k)= 
     2     sinphi**2*(costhe**2*(c011*cosphi**2*costhe**2 
     2    +c013*cosphi**2*sinthe**2 + c012*sinphi**2) 
     2    +sinthe**2*(c013*cosphi**2*costhe**2 
     2    +c033*cosphi**2*sinthe**2 + c023*sinphi**2) 
     2    +c055*cosphi**2*sin2the**2) -sin2phi*(c066*sin2phi*costhe**2
     2    +c044*sin2phi*sinthe**2)+cosphi**2*(c012*cosphi**2*costhe**2
     2    +c023*cosphi**2*sinthe**2 + c022*sinphi**2)

      c13(i,j,k)= 
     3     costhe**2*(c013*cosphi**2*costhe**2 
     3    +c033*cosphi**2*sinthe**2 + c023*sinphi**2)
     3    +sinthe**2*(c011*cosphi**2*costhe**2 
     3    +c013*cosphi**2*sinthe**2 + c012*sinphi**2) 
     3    -c055*cosphi**2*sin2the**2

      c14(i,j,k)= 
     4     cosphi*(c044*costhe*sinthe*sin2phi 
     4    -c066*costhe*sinthe*sin2phi) - sinphi*(0.5*sin2the
     4    *(c013*cosphi**2*costhe**2 + c033*cosphi**2*sinthe**2 
     4    +c023*sinphi**2) - 0.5*sin2the*(c011*cosphi**2*costhe**2 
     4    +c013*cosphi**2*sinthe**2 + c012*sinphi**2) 
     4    +c055*cosphi**2*cos2the*sin2the)

      c15(i,j,k)=
     5    -sinphi*(c044*costhe*sinthe*sin2phi 
     5    -c066*costhe*sinthe*sin2phi) - cosphi*(0.5*sin2the
     5    *(c013*cosphi**2*costhe**2 + c033*cosphi**2*sinthe**2 
     5    +c023*sinphi**2) - 0.5*sin2the*(c011*cosphi**2*costhe**2 
     5    +c013*cosphi**2*sinthe**2 + c012*sinphi**2) 
     5    +c055*cosphi**2*cos2the*sin2the)

      c16(i,j,k)= 
     1     0.5*sin2phi*(costhe**2*(c011*cosphi**2*costhe**2 
     1    +c013*cosphi**2*sinthe**2 + c012*sinphi**2) + sinthe**2
     1    *(c013*cosphi**2*costhe**2 + c033*cosphi**2*sinthe**2 
     1    +c023*sinphi**2) + c055*cosphi**2*sin2the**2) - cos2phi
     1    *(c066*sin2phi*costhe**2 + c044*sin2phi*sinthe**2) 
     1    -0.5*sin2phi*(c012*cosphi**2*costhe**2 
     1    +c023*cosphi**2*sinthe**2 + c022*sinphi**2)

      c22(i,j,k)= 
     1     sin2phi*(c066*sin2phi*costhe**2 + c044*sin2phi*sinthe**2)
     1    +sinphi**2*(costhe**2*(c012*cosphi**2 
     1    +c011*costhe**2*sinphi**2 + c013*sinphi**2*sinthe**2) 
     1    +sinthe**2*(c023*cosphi**2 + c013*costhe**2*sinphi**2 
     1    +c033*sinphi**2*sinthe**2) + c055*sinphi**2*sin2the**2) 
     1    +cosphi**2*(c022*cosphi**2 + c012*costhe**2*sinphi**2
     1    +c023*sinphi**2*sinthe**2)

      c23(i,j,k)=  
     1      costhe**2*(c023*cosphi**2 + c013*costhe**2*sinphi**2 
     1     +c033*sinphi**2*sinthe**2) + sinthe**2*(c012*cosphi**2 
     1     +c011*costhe**2*sinphi**2 + c013*sinphi**2*sinthe**2) 
     1     -c055*sinphi**2*sin2the**2

      c24(i,j,k)=
     1    -cosphi*(c044*costhe*sinthe*sin2phi 
     1    -c066*costhe*sinthe*sin2phi) - sinphi*(0.5*sin2the
     1    *(c023*cosphi**2 + c013*costhe**2*sinphi**2 
     1    +c033*sinphi**2*sinthe**2) - 0.5*sin2the*(c012*cosphi**2
     1    +c011*costhe**2*sinphi**2 + c013*sinphi**2*sinthe**2) 
     1    +c055*cos2the*sinphi**2*sin2the)

      c25(i,j,k)= 
     1     sinphi*(c044*costhe*sinthe*sin2phi 
     1    -c066*costhe*sinthe*sin2phi) - cosphi*(0.5*sin2the
     1    *(c023*cosphi**2 + c013*costhe**2*sinphi**2 
     1    +c033*sinphi**2*sinthe**2) - 0.5*sin2the*(c012*cosphi**2 
     1    +c011*costhe**2*sinphi**2 + c013*sinphi**2*sinthe**2) 
     1    +c055*cos2the*sinphi**2*sin2the)

      c26(i,j,k)= 
     1     cos2phi*(c066*sin2phi*costhe**2 + c044*sin2phi*sinthe**2) 
     1    +0.5*sin2phi*(costhe**2*(c012*cosphi**2 
     1    +c011*costhe**2*sinphi**2 + c013*sinphi**2*sinthe**2) 
     1    +sinthe**2*(c023*cosphi**2 + c013*costhe**2*sinphi**2 
     1    +c033*sinphi**2*sinthe**2) + c055*sinphi**2*sin2the**2) 
     1    -0.5*sin2phi*(c022*cosphi**2 + c012*costhe**2*sinphi**2 
     1    +c023*sinphi**2*sinthe**2)

      c33(i,j,k)= 
     1     c055*sin2the**2 +costhe**2*(c033*costhe**2+c013*sinthe**2) 
     1    +sinthe**2*(c013*costhe**2 + c011*sinthe**2)

      c34(i,j,k)= 
     1     sinphi*(0.5*sin2the*(c013*costhe**2 + c011*sinthe**2) 
     1    -0.5*sin2the*(c033*costhe**2 + c013*sinthe**2) 
     1    +c055*cos2the*sin2the)

      c35(i,j,k)= 
     1     cosphi*(0.5*sin2the*(c013*costhe**2 + c011*sinthe**2) 
     1    -0.5*sin2the*(c033*costhe**2 + c013*sinthe**2) 
     1    +c055*cos2the*sin2the)

      c36(i,j,k)= 
     1     0.5*sin2phi*(costhe**2*(c013*costhe**2 + c011*sinthe**2) 
     1    -c055*sin2the**2 +sinthe**2*(c033*costhe**2+c013*sinthe**2))
     1    -0.5*sin2phi*(c023*costhe**2 + c012*sinthe**2)

      c44(i,j,k)= 
     1     cosphi*(c044*cosphi*costhe**2 + c066*cosphi*sinthe**2) 
     1    +sinphi*(c055*sinphi*cos2the**2 
     1    +0.5*sin2the*(0.5*c011*sinphi*sin2the 
     1    -0.5*c013*sinphi*sin2the) - 0.5*sin2the
     1    *(0.5*c013*sinphi*sin2the - 0.5*c033*sinphi*sin2the))

      c45(i,j,k)= 
     1     cosphi*(c055*sinphi*cos2the**2 + 0.5*sin2the
     1    *(0.5*c011*sinphi*sin2the - 0.5*c013*sinphi*sin2the) 
     1    -0.5*sin2the*(0.5*c013*sinphi*sin2the
     1    -0.5*c033*sinphi*sin2the)) - sinphi*(c044*cosphi*costhe**2 
     1    +c066*cosphi*sinthe**2)

      c46(i,j,k)= 
     2     0.5*sin2phi*((0.5*c011*sinphi*sin2the 
     2    -0.5*c013*sinphi*sin2the)*costhe**2+(0.5*c013*sinphi*sin2the
     2    -0.5*c033*sinphi*sin2the)*sinthe**2 
     2    -c055*cos2the*sinphi*sin2the) - 0.5*sin2phi
     2    *(0.5*c012*sinphi*sin2the - 0.5*c023*sinphi*sin2the)-cos2phi
     2    *(c044*cosphi*costhe*sinthe-c066*cosphi*costhe*sinthe)

      c55(i,j,k)= 
     3     sinphi*(c044*sinphi*costhe**2 + c066*sinphi*sinthe**2) 
     3    +cosphi*(c055*cosphi*cos2the**2 + 0.5*sin2the
     3    *(0.5*c011*cosphi*sin2the - 0.5*c013*cosphi*sin2the)
     3    -0.5*sin2the*(0.5*c013*cosphi*sin2the 
     3    -0.5*c033*cosphi*sin2the))

      c56(i,j,k)= 
     3     cos2phi*(c044*costhe*sinphi*sinthe 
     5    -c066*costhe*sinphi*sinthe) - 0.5*sin2phi
     5    *(0.5*c012*cosphi*sin2the - 0.5*c023*cosphi*sin2the) 
     5    +0.5*sin2phi*(0.5*(c011*cosphi*sin2the-c013*cosphi*sin2the)
     5    *costhe**2 +(0.5*c013*cosphi*sin2the 
     5    -0.5*c033*cosphi*sin2the)*sinthe**2 
     5    -c055*cosphi*cos2the*sin2the)

      c66(i,j,k)= 
     3     cos2phi*(c066*cos2phi*costhe**2 + c044*cos2phi*sinthe**2)
     3    +0.5*sin2phi*(sinthe**2*(0.5*c013*sin2phi*costhe**2 
     3    +0.5*c033*sin2phi*sinthe**2 - 0.5*c023*sin2phi) 
     3    +costhe**2*(0.5*c011*sin2phi*costhe**2 
     3    +0.5*c013*sin2phi*sinthe**2 - 0.5*c012*sin2phi) 
     3    +0.5*c055*sin2phi*sin2the**2) - 0.5*sin2phi
     3    *(0.5*c012*sin2phi*costhe**2 + 0.5*c023*sin2phi*sinthe**2
     3    -0.5*c022*sin2phi)

      enddo
      enddo
      enddo

c      write(*,*)'c11=',c11(101,101,201)
c      write(*,*)'c12=',c12(101,101,201)
c      write(*,*)'c13=',c13(101,101,201)
c      write(*,*)'c14=',c14(101,101,201)
c      write(*,*)'c15=',c15(101,101,201)
c      write(*,*)'c16=',c16(101,101,201)
c      write(*,*)'c22=',c22(101,101,201)
c      write(*,*)'c23=',c23(101,101,201)
c      write(*,*)'c24=',c24(101,101,201)
c      write(*,*)'c25=',c25(101,101,201)
c      write(*,*)'c26=',c26(101,101,201)
c      write(*,*)'c33=',c33(101,101,201)
c      write(*,*)'c34=',c34(101,101,201)
c      write(*,*)'c35=',c35(101,101,201)
c      write(*,*)'c36=',c36(101,101,201)
c      write(*,*)'c44=',c44(101,101,201)
c      write(*,*)'c45=',c45(101,101,201)
c      write(*,*)'c46=',c46(101,101,201)
c      write(*,*)'c55=',c55(101,101,201)
c      write(*,*)'c56=',c56(101,101,201)
c      write(*,*)'c66=',c66(101,101,201)
c      stop

      return
      end subroutine stiffness_matrix_coefficient
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*********** unify the filename by adding number ***************
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
