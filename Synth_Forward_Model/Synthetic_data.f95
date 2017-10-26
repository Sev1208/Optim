program forward_model
!
! Creation of synthetic data
! Options: 1- Mogi (1957), 2- McTigue (1987), 3- Okada (1992)
!
implicit real*8 (a-h,o-z)
!
integer, parameter :: ndata=20000,nparam=40,ntype=4
real*8,  dimension (12,ndata) :: vdata,synth
real*8, dimension (nparam) ::x0,xinit,xmin,xmax
real*8 ux,uy,uz,uxx,uyx,uzx,uxy,uyy,uzy,uxz,uyz,uzz
integer*4, dimension(nparam) ::iflag
integer*4,dimension (ntype) :: npd
character(len=50), dimension(ntype) ::filename 

! Common blocks

common /np_data/npd

! Files
filename(1)="point_gps.dat"
filename(2)="point_insar.dat"
filename(3)="point_tilt.dat"
filename(4)="point_level.dat"
!filename(5)="point_gravi.dat"

write(*,*)'Which forward model? 1-Mogi, 2-McTigue, 3-Okada'
read(*,*)imodel
 
! Read data and source parameters
call read_init(imodel, x0, iflag, xinit, xmin, xmax, ndim)
call read_data(filename,vdata,ndt,npd)
synth(:,:)=0.0
pi=3.14159265
rad=pi/180.
!		
! InSAR unit vector		
!		
uxsar = -cos(x0(1)*rad)*sin(x0(2)*rad)		
uysar = sin(x0(1)*rad)*sin(x0(2)*rad)		
uzsar = cos(x0(2)*rad)

do i=1,ndt
  synth(1,i)=vdata(1,i)	!data type
  synth(2,i)=1.0		!time
  synth(3,i)=vdata(3,i)	!data point x
  synth(4,i)=vdata(4,i)	!data point y
  synth(5,i)=vdata(5,i)	!data point z
  ! 
  !  Estimation of the variables on the data point grid
  !
  do j=1,int(x0(5))
    if (imodel .eq. 1) then
	  jj=(j-1)*8
      call mogi(x0(3),x0(4),x0(jj+6),x0(jj+7),x0(jj+8),x0(jj+9),&
                     x0(jj+10),x0(jj+11),x0(jj+12),x0(jj+13),&
                     vdata(3,i),vdata(4,i),ux,uy,uz,uzx,uzy,dg)  
      synth(6,i) = synth(6,i) + ux
      synth(7,i) = synth(7,i) + uy
      synth(8,i) = synth(8,i) + uz
      synth(9,i) = synth(9,i) + ux*uxsar+uy*uysar+uz*uzsar
      synth(10,i) = synth(10,i) + uzx
      synth(11,i) = synth(11,i) + uzy
      !synth(12,i) = synth(12,i) + dg
	elseif (imodel .eq. 2) then
      jj=(j-1)*8
      call mctigue87(x0(3),x0(4),x0(jj+6),x0(jj+7),x0(jj+8),x0(jj+9),&
                     x0(jj+10),x0(jj+11),x0(jj+12),x0(jj+13),&
                     vdata(3,i),vdata(4,i),ux,uy,uz,uzx,uzy,dg)   	
      synth(6,i) = synth(6,i) + ux
      synth(7,i) = synth(7,i) + uy
      synth(8,i) = synth(8,i) + uz
      synth(9,i) = synth(9,i) + ux*uxsar+uy*uysar+uz*uzsar
      synth(10,i) = synth(10,i) + uzx
      synth(11,i) = synth(11,i) + uzy
      !synth(12,i) = synth(12,i) + dg
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
	  X=xt*cosa-yt*sina
	  Y=xt*sina+yt*cosa
	  Z=vdata(5,i)
	  !Dimensions of the fracture
	  AL1=-x0(jj+11)/2
	  AL2=x0(jj+11)/2
	  AW1=-x0(jj+12)/2
	  AW2=x0(jj+12)/2
      call DC3D(ALPHA,X,Y,Z,x0(jj+9),&
	            x0(jj+10),AL1,AL2,AW1,AW2,x0(jj+13),&
	   		    x0(jj+14),x0(jj+15),ux,uy,uz,uxx,uyx,uzx,uxy,&
	  		    uyy,uzy,uxz,uyz,uzz,RET)
    synth(6,i) = synth(6,i) + cosa*ux-sina*uy
    synth(7,i) = synth(7,i) - sina*ux+cosa*uy
    synth(8,i) = synth(8,i) + uz
    synth(9,i) = synth(9,i) + ux*uxsar+uy*uysar+uz*uzsar
    synth(10,i) = synth(10,i) + cosa*uzx+sina*uzy
    synth(11,i) = synth(11,i) - sina*uzx+cosa*uzy
    !synth(12,i) = synth(12,i) + dg
    endif
  enddo
enddo 
if (imodel .eq. 3) then
  open(2,file='sources.txt')
  write(2,*)x0(jj+6)-AL1*cos(phi),x0(jj+7)-AL1*sin(phi)
  write(2,*)x0(jj+6)+AL1*cos(phi),x0(jj+7)+AL1*sin(phi)
  close(2)
endif
!
call write_data(imodel,synth)
end


