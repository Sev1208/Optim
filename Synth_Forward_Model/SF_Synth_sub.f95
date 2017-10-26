!======================================================================
! Read the file containing the source parameters.
! Each parameter must have its own xmin/xmax file.
! The source file contains the following information:
! shearm,poiss: shear modulus and Poisson's ratio
! nsource: nb of sources
! lsource: informs which parameters are available to calculate mogi model
! xs,ys,zs: sources coordinates
! rho: density of the surrounding medium
! dvol,dp,rs: parameters of the point source
!
! Severine Furst (11-May-2017)
!
subroutine read_init(imodel, x0, iflag, xinit, xmin, xmax, ndim)
!
implicit real*8 (a-h,o-z)
integer, parameter :: nparam=40
real*8, dimension (nparam) ::x0,x0min,x0max,xinit,xmin,xmax
integer*4, dimension(nparam) ::iflag

x0(:)=0.
xinit(:)=0.
xmin(:)=0.
xmax(:)=0.
iflag(:)=0

if ((imodel .eq. 1) .OR. (imodel .eq. 2)) then
  open(10,file='ini_param_sphere')
  read(10,*) x0(1),x0min(1),x0max(1),iflag(1) !ahorizsar
  read(10,*) x0(2),x0min(2),x0max(2),iflag(2) !avertsar
  read(10,*) x0(3),x0min(3),x0max(3),iflag(3) !shearm
  read(10,*) x0(4),x0min(4),x0max(4),iflag(4) !poiss
  read(10,*) x0(5),x0min(5),x0max(5),iflag(5) !nsource
  do i=1,int(x0(5))
    ii=(i-1)*8
    read(10,*) x0(ii+6),x0min(ii+6),x0max(ii+6),iflag(ii+6) !lsource
    read(10,*) x0(ii+7),x0min(ii+7),x0max(ii+7),iflag(ii+7) !xs(i)
    read(10,*) x0(ii+8),x0min(ii+8),x0max(ii+8),iflag(ii+8) !ys(i)
    read(10,*) x0(ii+9),x0min(ii+9),x0max(ii+9),iflag(ii+9) !zs(i)
    read(10,*) x0(ii+10),x0min(ii+10),x0max(ii+10),iflag(ii+10) !rho(i)
    read(10,*) x0(ii+11),x0min(ii+11),x0max(ii+11),iflag(ii+11) !dvol(i)
    read(10,*) x0(ii+12),x0min(ii+12),x0max(ii+12),iflag(ii+12) !dp(i)
    read(10,*) x0(ii+13),x0min(ii+13),x0max(ii+13),iflag(ii+13) !rs(i)
  enddo
  close(10)
elseif (imodel .eq. 3) then
  open(10,file='ini_param_disl')
  read(10,*) x0(1),x0min(1),x0max(1),iflag(1) !ahorizsar
  read(10,*) x0(2),x0min(2),x0max(2),iflag(2) !avertsar
  read(10,*) x0(3),x0min(3),x0max(3),iflag(3) !shearm
  read(10,*) x0(4),x0min(4),x0max(4),iflag(4) !lambda
  read(10,*) x0(5),x0min(5),x0max(5),iflag(5) !nsource
  do i=1,int(x0(5))
    ii=(i-1)*7
	read(10,*) x0(ii+6),x0min(ii+6),x0max(ii+6),iflag(ii+6) !xs(i) source position at the surface
	read(10,*) x0(ii+7),x0min(ii+7),x0max(ii+7),iflag(ii+7) !ys(i) source position at the surface
	read(10,*) x0(ii+8),x0min(ii+8),x0max(ii+8),iflag(ii+8) !az(i) fracture azimuth
    read(10,*) x0(ii+9),x0min(ii+9),x0max(ii+9),iflag(ii+9) !zs(i)
    read(10,*) x0(ii+10),x0min(ii+10),x0max(ii+10),iflag(ii+10) !dip(i)
    read(10,*) x0(ii+11),x0min(ii+11),x0max(ii+11),iflag(ii+11) !length(i)
    read(10,*) x0(ii+12),x0min(ii+12),x0max(ii+12),iflag(ii+12) !width(i)
    read(10,*) x0(ii+13),x0min(ii+13),x0max(ii+13),iflag(ii+13) !DISL1(i)
    read(10,*) x0(ii+14),x0min(ii+14),x0max(ii+14),iflag(ii+14) !DISL2(i)
    read(10,*) x0(ii+15),x0min(ii+15),x0max(ii+15),iflag(ii+15) !DISL3(i)
  enddo
  close(10)
endif
j=0
do i=1,nparam
  if (iflag(i) .eq. 1) then
    j=j+1
    xinit(j)=x0(i)
    xmin(j)=x0min(i)
    xmax(j)=x0max(i)
  endif 
  ndim=j
enddo
return
end
!

!======================================================================
! Read the file 'filename' containing the data points.
! data type = 1(GPS),2(InSAR),3(Tilt),4(Levelling),5(Gravimetric)
! npd : number of data points
! vdata : data vectors (type,time,x,y)
!
! Severine Furst (17-July-2017)
!
subroutine read_data(filename,vdata,ndt,npd)
!
implicit real*8 (a-h,o-z)
!
integer, parameter :: ndata=20000,ntype=4
real*8,  dimension (12,ndata) :: vdata
integer*4,dimension (ntype) :: npd
character(len=50), dimension(ntype) ::filename
logical :: there
!

vdata(:,:)=0.
n=1
do i=1,ntype
  ! Test if the given file exists
  inquire(FILE=filename(i), EXIST=there)
  if (there)  then
    ! Number of line in the file
    open(10,file=filename(i),iostat=ierr,status="old")
    nline = 0
	do while (ierr .eq. 0)
  	  read(10,*,iostat=ierr)
  	  if (ierr .eq. 0) then
  		  nline = nline + 1
  	  endif
    enddo
    close(10)
    ! Test if number of data is greater than the allocated dimension
    if (nline .gt. ndata) then 
      write (*,*) 'ERROR, too many data in file ',filename(i)
	  stop
    endif
    open(10,file=filename(i),iostat=ierr,status="old")
	read(10,*)it
	npd(i)=nline-5
	do j=1,4
	  read(10,*)
	enddo
	do j=1,npd(i)
	  vdata(1,n)=it
	  read(10,*) vdata(2,n),vdata(3,n),vdata(4,n),vdata(5,n)
	  n=n+1
	enddo
  endif
enddo
ndt=n-1
return
end
!

!======================================================================
! surface displacement induced by a source of volumetric strain at xs,ys,
! source depth zs, volume change dvol, pore pressure change dp, source radius rs
! displacement at x and y : ux,uy,uz 
! tilt tx,ty
! gravity change dg
! medium : Poisson ratio poiss
! lsource = 1 : dvol,rs used as input
! lsource = 2 : dp,rs   used as input
! lsource = 3 : dp,rs used as input
!
! Segall (2010), pp. 207 et 291
!
subroutine mogi (shearm,poiss,tsource,xs,ys,zs,&
		              rho,dvol,dp,rs,&
		              x, y, ux, uy, uz, tx, ty, dg)
!
implicit real*8 (a-h,o-z)
pi   = 3.141592
bigg = 6.67e-11
!
! input parameters choice
!
if (tsource .eq. 1.) then
   rs3 = rs**3
   dp = dvol*shearm / (pi*rs3)
elseif (tsource .eq. 2) then
   rs3 = rs**3
   dvol = pi*dp*rs3 / shearm
elseif (tsource .eq. 3) then
   rs3 = dvol*shearm / (pi*dp)
   rs = rs3**(0.333333333333)
endif
!
coef = (1. - poiss) * dvol / pi
!
dx = x - xs
dy = y - ys
dz = -zs
dist  = sqrt(dx*dx+dy*dy+dz*dz)
denom = dist**3
!
! deplacements
!
ux = coef * (dx / denom)
uy = coef * (dy / denom)
uz = coef * (dz / denom)
!
! tilts
!
tx = -uz * 3.0 * dx / dist**2
ty = -uz * 3.0 * dy / dist**2
!
! gravité
!
coefg = bigg * rho * dvol
dg = coefg * (dist / denom)
!
return
end
!
!======================================================================
! McTigue (1987) solution for the point source
! Analytical solution with higher order terms that account 
! for the finite shape of a spherical body.
!
! surface displacement induced by a source of volumetric strain at xs,ys,
! source depth zs, volume change dvol, pore pressure change dp, source radius rs
! displacement at x and y : ux,uy,uz 
! tilt tx,ty
! medium : Poisson ratio poiss
! lsource = 1 : dvol,rs used as input
! lsource = 2 : dp,rs   used as input
! lsource = 3 : dp,rs used as input
!
! McTigue (1987) and USGS report Battaglia et al (2013)
! 
! Severine Furst 29-Aug-2017
!
subroutine mctigue87 (shearm,poiss,tsource,xs,ys,zs,&
		              rho,dvol,dp,rs,&
		              x, y, ux, uy, uz, tx, ty, dg)
!
implicit real*8 (a-h,o-z)
pi   = 3.141592
A=shearm/pi
!
! input parameters choice
!
if (tsource .eq. 1.) then
   rs3 = rs*rs*rs
   B=(1+(rs/zs)*(rs/zs)*(rs/zs)*(rs/zs))
   dp = A*dvol/(rs3*B)
elseif (tsource .eq. 2.) then
   rs3 = rs*rs*rs
   B=(1+(rs/zs)*(rs/zs)*(rs/zs)*(rs/zs))
   dvol = rs3*dp*B/A
elseif (tsource .eq. 3.) then
   write(*,*)'Impossible to estimate source radius, try giving option 1 or 2'
   STOP
endif
!
! First factor of McTigue equation
!
coef = (1.-poiss)*dvol/(pi*B)
!
! Denominators for displacements and tilts 
!
dx=x-xs
dy=y-ys
dz=-zs
dist  = sqrt(dx*dx+dy*dy+dz*dz)
denom1 = dist*dist*dist !displacement equation
denom2 = denom1*dist*dist !tilt equation
!
! Last terms of McTigue equation
!
C1=(1./(7.-5.*poiss))*(rs/zs)*(rs/zs)*(rs/zs)
C2=(1.+poiss)
C3=(2.-poiss)*((zs/dist)*(zs/dist))
!
C=1.-C1*(0.5*C2-3.75*C3)
D=3.-C1*(3.5*C2-18.75*C3)
!
! Displacements
!
ux = coef * C * dx / denom1
uy = coef * C * dy / denom1
uz = coef * C * dz / denom1
!
! Tilts
!
tx = -coef * D * dx * dz / denom2
ty = -coef * D * dy * dz / denom2
!
return
end

!======================================================================
! Write synthetic data previously estimated according to the chosen forward model
! 1- Mogi model
! 2- Okada model
!
! Severine Furst (17-July-2017)
!
subroutine write_data (imodel,vdata)
implicit real*8 (a-h,o-z)

integer, parameter :: ndata=10000,ntype=4
real*8,  dimension (12,ndata) :: vdata
integer*4,dimension (ntype) :: npd

common /np_data/npd

if (imodel .eq. 1) then
  open(1,file='data_gps_mogi.dat')
  open(2,file='data_insar_mogi.dat')
  open(3,file='data_tilt_mogi.dat')
  open(4,file='data_level_mogi.dat')
  open(5,file='data_gravi_mogi.dat')
elseif (imodel .eq. 2) then
  open(1,file='data_gps_mctigue.dat')
  open(2,file='data_insar_mctigue.dat')
  open(3,file='data_tilt_mctigue.dat')
  open(4,file='data_level_mctigue.dat')
  open(5,file='data_gravi_mctigue.dat')
elseif (imodel .eq. 3) then
  open(1,file='data_gps_okada.dat')
  open(2,file='data_insar_okada.dat')
  open(3,file='data_tilt_okada.dat')
  open(4,file='data_level_okada.dat')
  open(5,file='data_gravi_okada.dat')
endif
j=0

do i=1,ntype
  j=j+npd(i)
  write(i,*)int(vdata(1,j)), 'data type'
  write(i,*)'#t'
  write(i,*)'#x	(UTM)'
  write(i,*)'#y	(UTM)'
  write(i,*)'#z	(m)'
  if(vdata(1,j) .eq. 1) then
    write(i,*)' #ux	(m)'
    write(i,*)' #uy	(m)'
    write(i,*)' #uz	(m)'
	do k=1,j
      write(i,*)vdata(2:8,k)
	enddo
  else if(vdata(1,j) .eq. 2) then
    write(i,*)' #usar	(m)'
	do k=j-npd(i)+1,j
      write(i,*)vdata(2:5,k),vdata(9,k)
	enddo
  else if(vdata(1,j) .eq. 3) then
    write(i,*)' #tx	(rad)'
    write(i,*)' #ty	(rad)'
	do k=j-npd(i)+1,j
      write(i,*)vdata(2:5,k),vdata(10:11,k)
	enddo
  else if(vdata(1,j) .eq. 4) then
    write(i,*)' #uz	(m)'
	do k=j-npd(i)+1,j
      write(i,*)vdata(2:5,k),vdata(8,k)
	enddo
  else if (vdata(1,j) .eq. 5) then
    write(i,*)' #tg	(gal)'
	do k=j-npd(i)+1,j
      write(i,*)vdata(2:5,k),vdata(12,k)
	enddo
  endif
enddo
return
end

