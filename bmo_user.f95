!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!               provided by the user
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       initialisation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
subroutine init(ndim,xinit,xmin,xmax)
implicit real*8 (a-h,o-z)
!
integer, parameter :: ndata=20000,nparam=40,ntype=5
real*8,  dimension (12,ndata) :: vdata,uncdata,uncdata_1
real*8, dimension (nparam) ::x0,xinit,xmin,xmax
integer*4, dimension(nparam) ::iflag
integer*4,dimension (ntype) :: npd
character(len=50), dimension(ntype) ::filename 



! Common blocks
common /fmodel/imodel
common /param/x0
common /param_flag/iflag
common /data_obs/vdata,uncdata_1
common /nt_data/ndt
common /ni_data/npd
common /iter/n_iter

write(*,*)'Which forward model? 1-Mogi, 2-McTigue, 3-Okada'
read(*,*)imodel

! Files
if (imodel .eq. 1) then
  filename(1)="data_gps_mogi.dat"
  filename(2)="data_insar_mogi.dat"
  filename(3)="data_gps_mogi.dat"
  filename(4)="data_level_mogi.dat"
  filename(5)="data_gravi_mogi.dat"
elseif (imodel .eq. 2) then
  filename(1)="data_gps_mctigue.dat"
  filename(2)="data_insar_mctigue.dat"
  filename(3)="data_tilt_mctigue.dat"
  filename(4)="data_level_mctigue.dat"
  filename(5)="data_gravi_mctigue.dat"
elseif (imodel .eq. 3) then
  filename(1)="data_gps_okada.dat"
  filename(2)="data_insar_okada.dat"
  filename(3)="data_tilt_okada.dat"
  filename(4)="data_level_okada.dat"
  filename(5)="data_gravi_okada.dat"
endif  
 
! Read data and source parameters
call read_init(imodel, x0, iflag, xinit, xmin, xmax, ndim)
call read_data(filename,vdata,uncdata,ndt,npd)
call inverse_m(vdata,uncdata,uncdata_1,ndt)
n_iter=0
return
end
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! functional definition
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
subroutine func(ndim,x,f)
implicit double precision (a-h,o-z)
double precision f
integer, parameter :: ndata=20000,nparam=40,ntype=5
real*8,  dimension (12,ndata) :: vdata,vmodel,uncdata_1
real*8, dimension (nparam) ::x0,x
real*8, dimension (ntype) ::F_i
integer*4, dimension(nparam) ::iflag
integer*4,dimension (ntype) :: npd

common /fmodel/imodel         
common/nbeval/nbeval,nbevalp
common /param/x0
common /param_flag/iflag
common /data_obs/vdata,uncdata_1
common /data_mod/vmodel
common /nt_data/ndt
common /ni_data/npd
common /iter/n_iter
common /f_type/F_i

pi=3.14159265
bigg = 6.67e-11
rad=pi/180
	  
n_iter=n_iter+1
	  
if(nbeval.ge.0) then
  nbeval=nbeval+1
else
  nbevalp=nbevalp+1
endif

icallext=0   ! =0 si f definie ici, =1 sinon

if(icallext.eq.1) then
  j=0
  do i=1,nparam
    if(iflag(i).eq.1) then
      j=j+1
      x0(i)=x(j)
    endif
  enddo
  ! Calculate the foward model
  !		
  ! InSAR unit vector		  
  !		  
  uxsar = -cos(x0(1) *rad)*sin(x0(2)*rad)  
  uysar = sin(x0(1) *rad)*sin(x0(2)*rad)  
  uzsar = cos(x0(2)*rad)  
  vmodel(:,:)=0.0

  do i=1,ndt
    vmodel(1,i)=vdata(1,i)	!data type
    vmodel(2,i)=1.0			!time
    vmodel(3,i)=vdata(3,i)	!data point x
    vmodel(4,i)=vdata(4,i)	!data point y
	vmodel(5,i)=vdata(5,i)	!data point z
    ! 
    !  Estimation of the variables on the data point grid
    !  
    do j=1,int(x0(5))
      if (imodel .eq. 1) then
	    jj=(j-1)*8
        call mogi(x0(3),x0(4),x0(jj+6),x0(jj+7),x0(jj+8),x0(jj+9),&
                       x0(jj+10),x0(jj+11),x0(jj+12),x0(jj+13),&
                       vdata(3,i),vdata(4,i),ux,uy,uz,uzx,uzy,dg)  
        vmodel(6,i) = vmodel(6,i) + ux
        vmodel(7,i) = vmodel(7,i) + uy
        vmodel(8,i) = vmodel(8,i) + uz
        vmodel(9,i) = vmodel(9,i) + ux*uxsar+uy*uysar+uz*uzsar
        vmodel(10,i) = vmodel(10,i) + uzx
        vmodel(11,i) = vmodel(11,i) + uzy
        !vmodel(12,i) = vmodel(12,i) + dg
	  elseif (imodel .eq. 2) then
        jj=(j-1)*8
        call mctigue87(x0(3),x0(4),x0(jj+6),x0(jj+7),x0(jj+8),x0(jj+9),&
                       x0(jj+10),x0(jj+11),x0(jj+12),x0(jj+13),&
                       vdata(3,i),vdata(4,i),ux,uy,uz,uzx,uzy,dg)   	
        vmodel(6,i) = vmodel(6,i) + ux
        vmodel(7,i) = vmodel(7,i) + uy
        vmodel(8,i) = vmodel(8,i) + uz
        vmodel(9,i) = vmodel(9,i) + ux*uxsar+uy*uysar+uz*uzsar
        vmodel(10,i) = vmodel(10,i) + uzx
        vmodel(11,i) = vmodel(11,i) + uzy
        !vmodel(12,i) = vmodel(12,i) + dg
      elseif (imodel .eq. 3) then
	    jj=(j-1)*10
        ALPHA=(x0(3)+x0(4))/(x0(3)+2*x0(4))
	    ! Determination of the observation point in the Okada system
	    !For the rotation matrix
	    phi=rad*(x0(jj+8))
	    cosa=cos(phi)
	    sina=sin(phi)
	    !For the translation matrix
	    xt=vdata(3,i)-x0(jj+6)
	    yt=vdata(4,i)-x0(jj+7)
	    !Coordinates for Okada Model
	    Xok=xt*cosa-yt*sina
	    Yok=xt*sina+yt*cosa
	    Zok=vdata(5,i)
	    !Dimensions of the fracture
	    AL1=-x0(jj+11)/2
	    AL2=x0(jj+11)/2
	    AW1=-x0(jj+12)/2
	    AW2=x0(jj+12)/2
        call DC3D(ALPHA,Xok,Yok,Zok,x0(jj+9),&
	              x0(jj+10),AL1,AL2,AW1,AW2,x0(jj+13),&
	     		    x0(jj+14),x0(jj+15),ux,uy,uz,uxx,uyx,uzx,uxy,&
	    		    uyy,uzy,uxz,uyz,uzz,RET)
        vmodel(6,i) = vmodel(6,i) + cosa*ux-sina*uy
        vmodel(7,i) = vmodel(7,i) - sina*ux+cosa*uy
        vmodel(8,i) = vmodel(8,i) + uz
        vmodel(9,i) = vmodel(9,i) + ux*uxsar+uy*uysar+uz*uzsar
        vmodel(10,i) = vmodel(10,i) + cosa*uzx+sina*uzy
        vmodel(11,i) = vmodel(11,i) - sina*uzx+cosa*uzy
        !vmodel(12,i) = vmodel(12,i) + dg
      endif
    enddo
  enddo
  ! Calculate the functional
  call cost_function(vdata,uncdata_1,ndt,npd,vmodel,f,F_i,n_iter)
  	  
else if(icallext.eq.0) then
  f=0  !ndim
  do ii=1,ndim
!  f=f+cos(2*x(ii))+exp(-(x(ii)/20)**2)-(2/x(ii))*sin(2*x(ii))+2
!  f=f+x(ii)**2-cos(18*x(ii))
!  f=f+(x(ii)-10)**2
  enddo
endif

end
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Linearization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
subroutine funcpex(ndim,x,fpx)
implicit double precision (a-h,o-z)
double precision fpx(*),x(*)
common/nbeval/nbeval,nbevalp
nbevalp=nbevalp+1
!
!  user gradient (adjoint)
!
epsfd=1.e-5
do ii=1,ndim       
  fpx(ii)=0
enddo
do ii=1,ndim       
  x00=x(ii)
  x(ii)=x00+epsfd
  call func(ndim,x,fp)
  x(ii)=x00-2*epsfd
  call func(ndim,x,fm)
  fpx(ii)=(fp-fm)/(2*epsfd)
  x(ii)=x00
enddo
!
! linearization mogi

end
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Results
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
subroutine result(ndim,xoptg,vinit)
implicit double precision (a-h,o-z)      
double precision vinit(*),xoptg(*)
integer, parameter :: ndata=20000,nparam=40
real*8,  dimension (12,ndata) :: vdata,vmodel,uncdata_1
real*8, dimension (nparam) ::x0
integer*4, dimension(nparam) ::iflag

common /fmodel/imodel   
common /param_flag/iflag
common /param/x0
common /data_mod/vmodel
common /data_obs/vdata,uncdata_1
common /nt_data/ndt
!
ii=int(ndim/x0(5))
open(1,file='contsave.dat',status='unknown')
do i=1,int(x0(5))
  j=(i-1)*ii
  write(1,*) xoptg(j+1),xoptg(j+2),xoptg(j+3),xoptg(j+4),vinit(j+1),vinit(j+2),vinit(j+3),vinit(j+4)
enddo
close(1)    
j=0
do i=1,nparam
  if(iflag(i).eq.1) then
    j=j+1
    x0(i)=xoptg(j)
  endif
enddo 

vmodel(:,:)=0

do i=1,ndt
  vmodel(1,i)=vdata(1,i)	!data type
  vmodel(2,i)=1.0			!time
  vmodel(3,i)=vdata(3,i)	!data point x
  vmodel(4,i)=vdata(4,i)	!data point y
  vmodel(5,i)=vdata(5,i)	!data point z
  ! 
  !  Estimation of the variables on the data point grid
  !  
  do j=1,int(x0(5))
    if (imodel .eq. 1) then
    jj=(j-1)*8
      call mogi(x0(3),x0(4),x0(jj+6),x0(jj+7),x0(jj+8),x0(jj+9),&
                     x0(jj+10),x0(jj+11),x0(jj+12),x0(jj+13),&
                     vdata(3,i),vdata(4,i),ux,uy,uz,uzx,uzy,dg)  
      vmodel(6,i) = vmodel(6,i) + ux
      vmodel(7,i) = vmodel(7,i) + uy
      vmodel(8,i) = vmodel(8,i) + uz
      vmodel(9,i) = vmodel(9,i) + ux*uxsar+uy*uysar+uz*uzsar
      vmodel(10,i) = vmodel(10,i) + uzx
      vmodel(11,i) = vmodel(11,i) + uzy
      !vmodel(12,i) = vmodel(12,i) + dg
    elseif (imodel .eq. 2) then
      jj=(j-1)*8
      call mctigue87(x0(3),x0(4),x0(jj+6),x0(jj+7),x0(jj+8),x0(jj+9),&
                     x0(jj+10),x0(jj+11),x0(jj+12),x0(jj+13),&
                     vdata(3,i),vdata(4,i),ux,uy,uz,uzx,uzy,dg)   	
      vmodel(6,i) = vmodel(6,i) + ux
      vmodel(7,i) = vmodel(7,i) + uy
      vmodel(8,i) = vmodel(8,i) + uz
      vmodel(9,i) = vmodel(9,i) + ux*uxsar+uy*uysar+uz*uzsar
      vmodel(10,i) = vmodel(10,i) + uzx
      vmodel(11,i) = vmodel(11,i) + uzy
      !vmodel(12,i) = vmodel(12,i) + dg
    elseif (imodel .eq. 3) then
    jj=(j-1)*10
      ALPHA=(x0(3)+x0(4))/(x0(3)+2*x0(4))
    ! Determination of the observation point in the Okada system
    !For the rotation matrix
    phi=rad*(x0(jj+8))
    cosa=cos(phi)
    sina=sin(phi)
    !For the translation matrix
    xt=vdata(3,i)-x0(jj+6)
    yt=vdata(4,i)-x0(jj+7)
    !Coordinates for Okada Model
    Xok=xt*cosa-yt*sina
    Yok=xt*sina+yt*cosa
    Zok=vdata(5,i)
    !Dimensions of the fracture
    AL1=-x0(jj+11)/2
    AL2=x0(jj+11)/2
    AW1=-x0(jj+12)/2
    AW2=x0(jj+12)/2
      call DC3D(ALPHA,Xok,Yok,Zok,x0(jj+9),&
              x0(jj+10),AL1,AL2,AW1,AW2,x0(jj+13),&
     		    x0(jj+14),x0(jj+15),ux,uy,uz,uxx,uyx,uzx,uxy,&
    		    uyy,uzy,uxz,uyz,uzz,RET)
      vmodel(6,i) = vmodel(6,i) + cosa*ux-sina*uy
      vmodel(7,i) = vmodel(7,i) - sina*ux+cosa*uy
      vmodel(8,i) = vmodel(8,i) + uz
      vmodel(9,i) = vmodel(9,i) + ux*uxsar+uy*uysar+uz*uzsar
      vmodel(10,i) = vmodel(10,i) + cosa*uzx+sina*uzy
      vmodel(11,i) = vmodel(11,i) - sina*uzx+cosa*uzy
      !vmodel(12,i) = vmodel(12,i) + dg
    endif
  enddo
enddo

open(2,file='data_model.dat',status='unknown')
do i=1,ndt
 write(2,*)vmodel(1,i),vmodel(3,i),vmodel(4,i),vmodel(6,i),vmodel(7,i),vmodel(8,i),vmodel(10,i),vmodel(11,i)
enddo
close(2)
!
end
