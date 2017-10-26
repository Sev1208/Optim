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
    ii=(i-1)*10
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
! Read the file 'filename' containing the geodetic data.
! data type = 1(GPS),2(InSAR),3(Tilt),4(Levelling),5(Gravimetric)
! npd : number of data points
! vdata : data vectors (type,time,x,y,ux,uy,uz,usar,tx,ty,dg)
! vdata is build to be compared to the vector of modelled data from the mogi program
! (ux,uy,uz,usar,tx,ty,dg)
!
! Severine Furst (11-May-2017)
!
subroutine read_data(filename,vdata,uncdata,ndt,npd)
!
implicit real*8 (a-h,o-z)
!
integer, parameter :: ndata=20000,nparam=40,ntype=5
real*8,  dimension (12,ndata) :: vdata,uncdata
integer*4,dimension (ntype) :: npd
character(len=50), dimension(ntype) ::filename
logical :: there
!
unit_dis=1.e0 !convert the mm to meters
unit_ang=1.e0 !convert the angles to radians
unit_gal=1.e0 !convert the microgals to gals 

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
	! Read  data according to its type
	read(10,*)it
	if (it .eq. 1) then
	  npd(i)=nline-17
	  ! Skip the header
	  do j=1,16
	    read(10,*)
	  enddo
	  do j=1,npd(i)
	    vdata(1,n)=it
	    read(10,*) vdata(2,n),vdata(3,n),vdata(4,n),vdata(5,n),&
		           vdata(6,n),vdata(7,n),vdata(8,n),&
				   uncdata(1,n),uncdata(2,n),uncdata(3,n),&
				   uncdata(4,n),uncdata(5,n),uncdata(6,n),&
				   uncdata(7,n),uncdata(8,n),uncdata(9,n)
		n=n+1
	  enddo
    elseif (it .eq. 2) then
      npd(i)=nline-6
      do j=1,5
        read(10,*)
      enddo
      do j=1,npd(i)
        vdata(1,n) = it
	    read(10,*) vdata(2,n),vdata(3,n),vdata(4,n),vdata(5,n),vdata(9,n)
		n=n+1
	  enddo
    elseif (it .eq. 3) then
      npd(i)=nline-9
      do j=1,8
        read(10,*)
      enddo
      do j=1,npd(i)
	    vdata(1,n) = it
	    read(10,*) vdata(2,n),vdata(3,n),vdata(4,n),vdata(5,n),&
		           vdata(10,n),vdata(11,n),uncdata(1,n),uncdata(4,n)
		vdata(10:11,n)=unit_ang*vdata(10:11,n) !data in microrad
		uncdata(1:4,n)=unit_ang*uncdata(1:4,n) !data in microrad
		n=n+1
	  enddo
    elseif (it .eq. 4) then
      npd(i)=nline-6
      do j=1,5
        read(10,*)
      enddo
      do j=1,npd(i)
	    vdata(1,n) = it
	    read(10,*) vdata(2,n),vdata(3,n),vdata(4,n),vdata(5,n),vdata(8,n)
		n=n+1
	  enddo
    elseif (it .eq. 5) then
      npd(i)=nline-6
      do j=1,5
        read(10,*)
      enddo
      do j=1,npd(i)
	    vdata(1,n)=it
	    read(10,*) vdata(2,n),vdata(3,n),vdata(4,n),vdata(5,n),vdata(12,n)
		n=n+1
	  enddo
    endif
    close(10)
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
!
!======================================================================
! Subroutine to calculate the functional of each data type.
! ntype = 5 number of data types
! data type = 1(GPS),2(InSAR),3(Tilt),4(Levelling),5(Gravimetric)
! vdata : data vector (type,time,x,y,ux,uy,uz,usar,tx,ty,dg)
! vmodel : modelled data vector (type,time,x,y,ux,uy,uz,usar,tx,ty,dg)
! Levelling measures are corrected from relative values
!
! Severine Furst (23-May-2017)
!
subroutine cost_function (vdata,uncdata_1,ndt,npd,vmodel,f,F_i,n_iter)
!
implicit real*8 (a-h,o-z)
integer, parameter :: ndata=20000,ntype=5
real*8,  dimension (12,ndata) :: vdata,vmodel,uncdata_1
real*8, dimension (ntype)::dat_m, cal_m, fi, alphai, F_i
integer*4,dimension (ntype) :: npd
double precision f
!
f=0
fi(:)=0.
F_i(:)=0.
dat_m(:)=0.
cal_m(:)=0.
!alphai(:)=1./(ntype)
alphai(1)=0.5
alphai(3)=0.5
!
! Mean of relative data
!
do i=1,ndt
  if (vdata(1,i) .eq. 4) then
    dat_m(4)=dat_m(4)+vdata(8,i)
	cal_m(4)=cal_m(4)+vmodel(8,i)
  endif
enddo
if (npd(4) .ne. 0) then
  dat_m(4)=dat_m(4)/npd(4)
  cal_m(4)=cal_m(4)/npd(4)
else 
  dat_m(4)=0
  cal_m(4)=0
endif
!
! Functional
!
do i=1,ndt
  if (vdata(1,i) .eq. 1) then 
    rx=vmodel(6,i)-vdata(6,i)
	ry=vmodel(7,i)-vdata(7,i)
	rz=vmodel(8,i)-vdata(8,i)
    res=rx*(rx*uncdata_1(1,i)+ry*uncdata_1(2,i)+rz*uncdata_1(3,i))+&
	   ry*(rx*uncdata_1(4,i)+ry*uncdata_1(5,i)+rz*uncdata_1(6,i))+&
	   rz*(rx*uncdata_1(7,i)+ry*uncdata_1(8,i)+rz*uncdata_1(9,i))
    fi(1)=fi(1)+res
  else if (vdata(1,i) .eq. 2) then
    !Modify for uncertainties
    res=(vdata(9,i)-vmodel(9,i))**2
    fi(2)=fi(2)+res
  else if (vdata(1,i) .eq. 3) then
    rx=vmodel(10,i)-vdata(10,i)
	ry=vmodel(11,i)-vdata(11,i)
	res=rx*(rx*uncdata_1(1,i)+ry*uncdata_1(2,i))+&
	    ry*(rx*uncdata_1(3,i)+ry*uncdata_1(4,i))
	fi(3)=fi(3)+res
  else if (vdata(1,i) .eq. 4) then
    !Modify for uncertainties
    res=(vdata(8,i)-dat_m(4)-vmodel(8,i)+cal_m(4))**2
    fi(4)=fi(4)+res
  else if (vdata(1,i) .eq. 5) then
    !Modify for uncertainties
    res=(vdata(12,i)-vmodel(12,i))**2
    fi(5)=fi(5)+res
  endif
enddo
!
do i=1,ntype
  if (npd(i) .ne. 0) then
    F_i(i)=alphai(i)*fi(i)/npd(i)
    f=f+F_i(i)
  endif
enddo
return
end
!
!======================================================================
! Calculate the inverse covariance matrix of 1x1, 2x2 or 3x3 matrices
!
! Severine Furst (21-Sept-2017)
!
subroutine inverse_m (vdata,uncdata,uncdata_1,ndt,npd)
!
implicit real*8 (a-h,o-z)
integer, parameter :: ndata=20000
real*8,  dimension (12,ndata) :: vdata,uncdata,uncdata_1
real*8,  dimension (3,3) :: A,C
real*8,  dimension (2,2) :: B,D
logical :: FLAG

do i=1,ndt
  l=1
  if (vdata(1,i) .eq. 1) then
    do j=1,3
      do k=1,3
        A(j,k)=uncdata(l,i)
	    l=l+1
	  enddo 
    enddo
	call M33INV (A, C, FLAG)
    l=1
    do j=1,3
      do k=1,3
	    uncdata_1(l,i)=(C(j,k)+C(k,j))/2
	    l=l+1
	  enddo 
    enddo
  else if (vdata(1,i) .eq. 3) then
     do j=1,2
      do k=1,2
        B(j,k)=uncdata(l,i)
	    l=l+1
	  enddo 
    enddo
	call M22INV (B, D, FLAG)
    l=1
    do j=1,2
      do k=1,2
	    uncdata_1(l,i)=(D(j,k)+D(k,j))/2
	    l=l+1
	  enddo 
    enddo
  else 
    uncdata_1(1,i)=1./uncdata(1,i)
  endif
enddo

return
end



!
!***********************************************************************************************************************************
!  M33INV  -  Compute the inverse of a 3x3 matrix.
!
!  A       = input 3x3 matrix to be inverted
!  AINV    = output 3x3 inverse of matrix A
!  OK_FLAG = (output) .TRUE. if the input matrix could be inverted, and .FALSE. if the input matrix is singular.
!
!  Programmer:   David G. Simpson
!                NASA Goddard Space Flight Center
!                Greenbelt, Maryland  20771
!***********************************************************************************************************************************

      SUBROUTINE M33INV (A, AINV, OK_FLAG)

      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN)  :: A
      DOUBLE PRECISION, DIMENSION(3,3), INTENT(OUT) :: AINV
      LOGICAL, INTENT(OUT) :: OK_FLAG

      DOUBLE PRECISION, PARAMETER :: EPS = 1.0D-20
      DOUBLE PRECISION :: DET
      DOUBLE PRECISION, DIMENSION(3,3) :: COFACTOR

      DET =   A(1,1)*A(2,2)*A(3,3)  &
            - A(1,1)*A(2,3)*A(3,2)  &
            - A(1,2)*A(2,1)*A(3,3)  &
            + A(1,2)*A(2,3)*A(3,1)  &
            + A(1,3)*A(2,1)*A(3,2)  &
            - A(1,3)*A(2,2)*A(3,1)

      IF (ABS(DET) .LE. EPS) THEN
         AINV = 0.0D0
         OK_FLAG = .FALSE.
         RETURN
      END IF

      COFACTOR(1,1) = +(A(2,2)*A(3,3)-A(2,3)*A(3,2))
      COFACTOR(1,2) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))
      COFACTOR(1,3) = +(A(2,1)*A(3,2)-A(2,2)*A(3,1))
      COFACTOR(2,1) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))
      COFACTOR(2,2) = +(A(1,1)*A(3,3)-A(1,3)*A(3,1))
      COFACTOR(2,3) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))
      COFACTOR(3,1) = +(A(1,2)*A(2,3)-A(1,3)*A(2,2))
      COFACTOR(3,2) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))
      COFACTOR(3,3) = +(A(1,1)*A(2,2)-A(1,2)*A(2,1))

      AINV = TRANSPOSE(COFACTOR) / DET

      OK_FLAG = .TRUE.

      RETURN

      END SUBROUTINE M33INV

!***********************************************************************************************************************************
!  M22INV  -  Compute the inverse of a 2x2 matrix.
!
!  A       = input 2x2 matrix to be inverted
!  AINV    = output 2x2 inverse of matrix A
!  OK_FLAG = (output) .TRUE. if the input matrix could be inverted, and .FALSE. if the input matrix is singular.
!***********************************************************************************************************************************

      SUBROUTINE M22INV (A, AINV, OK_FLAG)

      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(2,2), INTENT(IN)  :: A
      DOUBLE PRECISION, DIMENSION(2,2), INTENT(OUT) :: AINV
      LOGICAL, INTENT(OUT) :: OK_FLAG

      DOUBLE PRECISION, PARAMETER :: EPS = 1.0D-20
      DOUBLE PRECISION :: DET
      DOUBLE PRECISION, DIMENSION(2,2) :: COFACTOR


      DET =   A(1,1)*A(2,2) - A(1,2)*A(2,1)

      IF (ABS(DET) .LE. EPS) THEN
         AINV = 0.0D0
         OK_FLAG = .FALSE.
         RETURN
      END IF

      COFACTOR(1,1) = +A(2,2)
      COFACTOR(1,2) = -A(2,1)
      COFACTOR(2,1) = -A(1,2)
      COFACTOR(2,2) = +A(1,1)

      AINV = TRANSPOSE(COFACTOR) / DET

      OK_FLAG = .TRUE.

      RETURN

      END SUBROUTINE M22INV