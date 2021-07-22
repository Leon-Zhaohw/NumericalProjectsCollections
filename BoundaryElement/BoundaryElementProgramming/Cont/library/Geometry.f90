!     Last change:  CD   25 Mar 2000    8:10 pm
MODULE Geometry_lib
USE Utility_lib
CONTAINS
SUBROUTINE Normal_Jac(v3,Jac,xsi,eta,ldim,nodes,inci,coords)
!--------------------------------------------------------
!   Computes normal vector and Jacobian
!   at point with local coordinates xsi,eta
!--------------------------------------------------------
IMPLICIT NONE
REAL, INTENT(IN):: xsi,eta            ! intrinsic co-ordinates of point
INTEGER,INTENT(IN):: ldim             ! element dimension
INTEGER,INTENT(IN):: nodes            ! number of nodes
INTEGER,INTENT(IN):: inci(:)          ! element incidences
REAL, INTENT(IN)  :: coords(:,:)      ! node coordinates
REAL,INTENT(OUT):: v3(:)              ! Vector normal to point
REAL,INTENT(OUT):: Jac                ! Jacobian
REAL,ALLOCATABLE  :: DNi(:,:)         ! Derivatives of shape function
REAL,ALLOCATABLE  :: v1(:),v2(:)      ! Vectors in xsi,eta directions
INTEGER :: Cdim ,i                      ! Cartesian dimension
!    Cartesian dimension:
 Cdim= ldim+1
!    Allocate temporary arrays
ALLOCATE (DNi(nodes,ldim),V1(Cdim),V2(Cdim))
!    Compute derivatives of shape function
Call Serendip_deriv(DNi,xsi,eta,ldim,nodes,inci)
!    Compute vectors in xsi (eta) direction(s)
DO I=1,Cdim
	V1(I)= DOT_PRODUCT(DNi(:,1),COORDS(I,:))
	IF(ldim == 2) THEN
		V2(I)= DOT_PRODUCT(DNi(:,2),COORDS(I,:))
	END IF
END DO
!    Compute normal vector
IF(ldim == 1) THEN
	v3(1)= V1(2)
	v3(2)= -V1(1)
	ELSE
  V3= Vector_ex(v1,v2)
END IF
!    Normalise
CAll Vector_norm(V3,Jac)
DEALLOCATE (DNi,V1,V2)
RETURN
END SUBROUTINE Normal_Jac

SUBROUTINE Tangent(v1,v2,xsi,eta,ldim,nodes,inci,coords)
!--------------------------------------------------------
!   Computes vectors tangent to BE
!   at point with local coordinates xsi,eta
!--------------------------------------------------------
IMPLICIT NONE
INTEGER,INTENT(IN):: ldim             ! element dimension
REAL, INTENT(IN)  :: xsi,eta            ! intrinsic co-ordinates of point
INTEGER,INTENT(IN):: nodes            ! number of nodes
INTEGER,INTENT(IN):: inci(:)          ! element incidences
REAL, INTENT(IN)  :: coords(:,:)      ! node coordinates
REAL,INTENT(OUT)  :: v1(ldim+1),v2(ldim+1)              ! Vector normal to point
REAL,ALLOCATABLE  :: DNi(:,:)         ! Derivatives of shape function
!      REAL,ALLOCATABLE  :: v1(:),v2(:)      ! Vectors in xsi,eta directions
INTEGER :: Cdim ,i                      ! Cartesian dimension
!    Cartesian dimension:
Cdim= ldim+1
!    Allocate temporary arrays
ALLOCATE (DNi(nodes,ldim))
!    Compute derivatives of shape function
Call Serendip_deriv(DNi,xsi,eta,ldim,nodes,inci)
!    Compute vectors in xsi (eta) direction(s)
DO i=1,Cdim
	v1(i)= DOT_PRODUCT(DNi(:,1),COORDS(i,:))
	IF(ldim == 2) THEN
		v2(i)= DOT_PRODUCT(DNi(:,2),COORDS(i,:))
	END IF
END DO
DEALLOCATE (DNi)
RETURN
END SUBROUTINE Tangent

Subroutine Serendip_func(Ni,xsi,eta,ldim,nodes,inci)
!---------------------------------
!   Computes Serendipity shape functions  Ni(xsi,eta)
!   for one and two-dimensional (linear/parabolic) finite boundary elements
!---------------------------------
REAL,INTENT(OUT):: Ni(:)      ! Array with shape function values at xsi,eta
REAL, INTENT(IN):: xsi,eta    ! intrinsic co-ordinates
INTEGER,INTENT(IN):: ldim     ! element dimension
INTEGER,INTENT(IN):: nodes    ! number of nodes
INTEGER,INTENT(IN):: inci(:)  ! element incidences
REAL              :: mxs,pxs,met,pet   !  temporary variables
SELECT CASE (ldim)
	CASE (1)   !   one-dimensional element
		Ni(1)= 0.5*(1.0 - xsi) ;  Ni(2)= 0.5*(1.0 + xsi)
		IF(nodes == 2) RETURN  !  linear element finshed
		Ni(3)=  1.0 - xsi*xsi
		Ni(1)= Ni(1) - 0.5*Ni(3) ; Ni(2)= Ni(2) - 0.5*Ni(3)
	CASE(2)    !    two-dimensional element
		mxs= 1.0-xsi ; pxs= 1.0+xsi ; met= 1.0-eta ; pet= 1.0+eta
		Ni(1)= 0.25*mxs*met ; Ni(2)= 0.25*pxs*met
		Ni(3)= 0.25*pxs*pet ; Ni(4)= 0.25*mxs*pet
    IF(nodes == 4) RETURN   !  linear element finshed
    IF(Inci(5) > 0) THEN        !  zero node number means node is missing
			Ni(5)= 0.5*(1.0 -xsi*xsi)*met
			Ni(1)= Ni(1) - 0.5*Ni(5) ; Ni(2)= Ni(2) - 0.5*Ni(5)
    END IF
    IF(Inci(6) > 0) THEN
			Ni(6)= 0.5*(1.0 -eta*eta)*pxs
			Ni(2)= Ni(2) - 0.5*Ni(6) ;  Ni(3)= Ni(3) - 0.5*Ni(6)
    END IF
    IF(Inci(7) > 0) THEN
			Ni(7)= 0.5*(1.0 -xsi*xsi)*pet
			Ni(3)= Ni(3) - 0.5*Ni(7) ; Ni(4)= Ni(4) - 0.5*Ni(7)
    END IF
    IF(Inci(8) > 0) THEN
			Ni(8)= 0.5*(1.0 -eta*eta)*mxs
			Ni(4)= Ni(4) - 0.5*Ni(8) ; Ni(1)= Ni(1) - 0.5*Ni(8)
    END IF
CASE DEFAULT     !   error message
CALL Error_message('Element dimension not 1 or 2' )
END SELECT
RETURN
END SUBROUTINE Serendip_func

SUBROUTINE Serendip_deriv(DNi,xsi,eta,ldim,nodes,inci)
!---------------------------------
!   Computes Derivatives ofSerendipity shape functions  Ni(xsi,eta)
!   for one and two-dimensional (linear/parabolic) finite boundary elements
!---------------------------------
IMPLICIT NONE
REAL,INTENT(OUT):: DNi(:,:)    ! Array with shape function derivatives at xsi,eta
REAL, INTENT(IN):: xsi,eta     ! intrinsic co-ordinates
INTEGER,INTENT(IN):: ldim     ! element dimension
INTEGER,INTENT(IN):: nodes    ! number of nodes
INTEGER,INTENT(IN):: inci(:)  ! element incidences
REAL            :: mxs,pxs,met,pet   !  temporary variables
SELECT CASE (ldim)
	CASE (1)   !   one-dimensional element
		DNi(1,1)= -0.5
		DNi(2,1)= 0.5
		IF(nodes == 2) RETURN  !  linear element finshed
		DNi(3,1)=  -2.0*xsi
		DNi(1,1)= DNi(1,1) - 0.5*DNi(3,1)
		DNi(2,1)= DNi(2,1) - 0.5*DNi(3,1)
	CASE (2)   !    two-dimensional element
		mxs= 1.0-xsi
		pxs= 1.0+xsi
		met= 1.0-eta
		pet= 1.0+eta
		DNi(1,1)= -0.25*met
		DNi(1,2)= -0.25*mxs
		DNi(2,1)=  0.25*met
		DNi(2,2)= -0.25*pxs
		DNi(3,1)=  0.25*pet
		DNi(3,2)=  0.25*pxs
		DNi(4,1)= -0.25*pet
		DNi(4,2)=  0.25*mxs
		IF(nodes == 4) RETURN  !  linear element finshed
		IF(Inci(5) > 0) THEN   !  zero node number means node is missing
			DNi(5,1)= -xsi*met
			DNi(5,2)= -0.5*(1.0 -xsi*xsi)
			DNi(1,1)= DNi(1,1) - 0.5*DNi(5,1)
			DNi(1,2)= DNi(1,2) - 0.5*DNi(5,2)
			DNi(2,1)= DNi(2,1) - 0.5*DNi(5,1)
			DNi(2,2)= DNi(2,2) - 0.5*DNi(5,2)
		END IF
		IF(Inci(6) > 0) THEN
			DNi(6,1)= 0.5*(1.0 -eta*eta)
			DNi(6,2)= -eta*pxs
			DNi(2,1)= DNi(2,1) - 0.5*DNi(6,1)
			DNi(2,2)= DNi(2,2) - 0.5*DNi(6,2)
			DNi(3,1)= DNi(3,1) - 0.5*DNi(6,1)
			DNi(3,2)= DNi(3,2) - 0.5*DNi(6,2)
		END IF
		IF(Inci(7) > 0) THEN
			DNi(7,1)= -xsi*pet
			DNi(7,2)= 0.5*(1.0 -xsi*xsi)
			DNi(3,1)= DNi(3,1) - 0.5*DNi(7,1)
			DNi(3,2)= DNi(3,2) - 0.5*DNi(7,2)
			DNi(4,1)= DNi(4,1) - 0.5*DNi(7,1)
			DNi(4,2)= DNi(4,2) - 0.5*DNi(7,2)
		END IF
		IF(Inci(8) > 0) THEN
			DNi(8,1)= -0.5*(1.0 -eta*eta)
			DNi(8,2)= -eta*mxs
			DNi(4,1)= DNi(4,1) - 0.5*DNi(8,1)
			DNi(4,2)= DNi(4,2) - 0.5*DNi(8,2)
			DNi(1,1)= DNi(1,1) - 0.5*DNi(8,1)
			DNi(1,2)= DNi(1,2) - 0.5*DNi(8,2)
		END IF
	CASE DEFAULT     !   error message
CALL Error_message('Element dimension not 1 or 2' )
END SELECT
RETURN
END SUBROUTINE Serendip_deriv
 
SUBROUTINE Cartesian(Ccor,Ni,ldim,elcor)
!--------------------------------------------------------
!   Computes Cartesian coordinates
!   at point with local coordinates xsi,eta
!--------------------------------------------------------
IMPLICIT NONE
REAL,INTENT(OUT)  :: Ccor(:)            ! Cart. coords of point xsi,eta
REAL,INTENT(IN)   :: Ni(:)              ! Shape functions at xsi,eta
REAL, INTENT(IN)  :: elcor(:,:)       ! element coordinates
INTEGER           :: ldim, Cdim,I             ! Cartesian dimension
!    Cartesian dimension:
Cdim= ldim+1
!    Compute vectors in xsi (eta) direction(s)
DO I=1,Cdim
	Ccor(I)= DOT_PRODUCT(Ni(:),Elcor(I,:))
END DO
RETURN
END SUBROUTINE Cartesian

SUBROUTINE Vector_norm(v3,Vlen)
!----------------------------------------
!   Normalise vector
!----------------------------------------
IMPLICIT NONE
REAL, INTENT(INOUT)  :: V3(:)               !     Vector to be normalised
REAL, INTENT(OUT)    :: Vlen                !     length of vector
Vlen= SQRT( SUM(v3*v3))
IF(Vlen == 0.) RETURN
V3= V3/Vlen
RETURN
END SUBROUTINE Vector_norm

FUNCTION Vector_ex(v1,v2)
!----------------------------------------
!   Returns vector x-product v1xv2
!   Where v1 and v2 are dimension 3
!----------------------------------------
IMPLICIT NONE
REAL, INTENT(IN) :: V1(3),V2(3)               !     Input
REAL             :: Vector_ex(3)              !     Result
Vector_ex(1)=V1(2)*V2(3)-V2(2)*V1(3)
Vector_ex(2)=V1(3)*V2(1)-V1(1)*V2(3)
Vector_ex(3)=V1(1)*V2(2)-V1(2)*V2(1)
RETURN
END FUNCTION Vector_ex

REAL FUNCTION Dist(xa,xe,Cdim)
!----------------------------------------
!    Computes the distance between two points
!    with coordinates (xa,ya) and (xe,ye)
!------------------------------------------
IMPLICIT NONE
REAL, INTENT(IN)    :: xa(:)      !  coords of point 1
REAL, INTENT(IN)    :: xe(:)      !  coords of point 2
INTEGER, INTENT(IN) :: Cdim       !  Cartesian dimension
INTEGER             :: N
REAL 								:: SUMS
SUMS= 0.0
DO N=1,Cdim
	SUMS= SUMS + (xa(n)-xe(n))**2
END DO
Dist= SQRT(SUMS)
RETURN
END FUNCTION Dist

REAL FUNCTION Direc(xA,xE)
!--------------------------------------------------------
!  Computes the Direction-angle from point xA to point xE
!--------------------------------------------------------
IMPLICIT NONE
REAL, INTENT(IN)     :: xA(2)
REAL, INTENT(IN)     :: xE(2)
REAL                 :: pi=3.1415926536
Direc=ATAN2((xE(2)-xA(2)),(xE(1)-xA(1)))
IF (Direc < 0.00000000) THEN
	Direc= Direc + 2*pi
END IF
RETURN
END FUNCTION Direc

SUBROUTINE Ortho(v3,v1,v2)
!----------------------------------------------------------
!     DETERMINES ORTHOGONAL VECTORS
!     V1, V2, V3
!     USING
!     V2 = V3 X VX
!     OR
!     V2 = V3 X VY     IF V3=VX
!     AND
!     V1 = V2 X V3
!--------------------------------------------------------
REAL, INTENT (IN) ::  v3(:)         !  Normal vector
REAL              ::  V1(:), V2(:)  !  Orthogonal vectors
REAL              ::  vx(3),vy(3)   !  Vectors in coordinate directions
vx= (/1.0,0.0,0.0/) ; vy= (/0.0,1.0,0.0/)
IF(ABS(V3(1)) + 0.005 .GE. 1.0) THEN
	v2= Vector_ex(v3,vy)
ELSE
	v2= Vector_ex(v3,vx)
END IF
V1= Vector_ex(v2,v3)
RETURN
END SUBROUTINE Ortho

REAL FUNCTION Min_dist(Elcor,xPi,Nodel,ldim,inci)
IMPLICIT NONE
REAL,INTENT(IN)				::	Elcor(:,:)							!	Coordinates of Element
REAL,INTENT(IN)				::	xPi(:)									!	Coordinates of Collocation point
REAL									::	DET,A,B,C,D	
REAL									::	xsi,eta,Dxsi,Deta
REAL									::	DistPS,DistPS_N
REAL									::	L,ELengx,ELenge
REAL,ALLOCATABLE			::	Ni(:)										! Shape function
REAL,ALLOCATABLE			::	DNi(:,:)								! Derivatives of shape function
REAL,ALLOCATABLE			::	r(:),xS(:),Dxs(:,:)
INTEGER,INTENT(IN)		::	Nodel,ldim
INTEGER,INTENT(IN)		::	inci(:)
INTEGER								::	n,Cdim

Cdim= ldim + 1
ALLOCATE(Ni(Nodel),DNi(Nodel,ldim),r(Cdim),xS(Cdim),DxS(Cdim,ldim))
SELECT CASE(Cdim)
	CASE(2)
		xsi=0.0
		CALL Serendip_func(Ni,xsi,eta,ldim,Nodel,inci)
		xS(1)= DOT_PRODUCT(Ni(:),Elcor(1,:))
		xS(2)= DOT_PRODUCT(Ni(:),Elcor(2,:))
		r= xS-xPi
		CALL Elength(L,Elcor,Nodel,ldim)
		DistPS=Dist(xPi(:),xS(:),Cdim)
		IF(((DistPS-L/2)/L) > 4.)THEN
			Min_dist=DistPS
			RETURN
		END IF
		DO n=1,20
			IF(n > 1)DistPS= DistPS_N
			CALL Serendip_deriv(DNi,xsi,eta,ldim,nodel,inci)
			DxS(1,1)= DOT_PRODUCT(DNi(:,1),Elcor(1,:))
			DxS(2,1)= DOT_PRODUCT(DNi(:,1),Elcor(2,:))
			DET= DxS(1,1)**2+DxS(2,1)**2
			Dxsi= -1/DET * DOT_PRODUCT(DxS(:,1),r(:))
			xsi= xsi+ Dxsi
			CALL Serendip_func(Ni,xsi,eta,ldim,Nodel,inci)
			xS(1)= DOT_PRODUCT(Ni(:),Elcor(1,:))
			xS(2)= DOT_PRODUCT(Ni(:),Elcor(2,:))
			r= xS- xPi
			DistPS_N=Dist(xPi(:),xS(:),Cdim)
			IF(ABS((DistPS- DistPS_N)/DistPS_N) < 0.05)EXIT
		END DO
		IF(xsi >= -1.0 .and. xsi <= 1.0)THEN
			Min_dist=DistPS_N
			RETURN
		ELSE IF(xsi < -1.0)THEN
			Min_dist=Dist(Elcor(:,1),xPi(:),Cdim)
			RETURN
		ELSE IF(xsi > -1.0)THEN
			Min_dist=Dist(Elcor(:,2),xPi(:),Cdim)
			RETURN
		END IF
	CASE(3)
		xsi=0.0
		eta=0.0
		CALL Serendip_func(Ni,xsi,eta,ldim,Nodel,inci)
		xS(1)= DOT_PRODUCT(Ni(:),Elcor(1,:))
		xS(2)= DOT_PRODUCT(Ni(:),Elcor(2,:))
		xS(3)= DOT_PRODUCT(Ni(:),Elcor(3,:))
		r= xS-xPi
		ELengx= Dist((Elcor(:,3)+Elcor(:,2))/2.,(Elcor(:,4)+Elcor(:,1))/2.,Cdim)  ! Lxsi
		ELenge= Dist((Elcor(:,2)+Elcor(:,1))/2.,(Elcor(:,3)+Elcor(:,4))/2.,Cdim)  ! Leta
		DistPS=Dist(xPi(:),xS(:),Cdim)
		IF(((DistPS-ELengx/2)/ELengx) > 4. .and.((DistPS-ELenge/2)/ELenge)> 4.)THEN
			Min_dist=DistPS
			RETURN
		END IF
		DO n=1,40
			IF(n > 1)DistPS= DistPS_N
			CALL Serendip_deriv(DNi,xsi,eta,ldim,nodel,inci)
			DxS(1,1)= DOT_PRODUCT(DNi(:,1),Elcor(1,:))
			DxS(2,1)= DOT_PRODUCT(DNi(:,1),Elcor(2,:))
			DxS(3,1)= DOT_PRODUCT(DNi(:,1),Elcor(3,:))
			DxS(1,2)= DOT_PRODUCT(DNi(:,2),Elcor(1,:))
			DxS(2,2)= DOT_PRODUCT(DNi(:,2),Elcor(2,:))
			DxS(3,2)= DOT_PRODUCT(DNi(:,2),Elcor(3,:))
			A= DxS(1,1)**2+DxS(2,1)**2+DxS(3,1)**2
			B= DxS(1,1)*DxS(1,2)+DxS(2,1)*DxS(2,2)+DxS(3,1)*DxS(3,2)
			C=B
			D= DxS(1,2)**2+DxS(2,2)**2+DxS(3,2)**2
			DET= A*D - C*B
			Dxsi = -1/DET * DOT_PRODUCT(DxS(:,1),r(:))
			Deta = -1/DET * DOT_PRODUCT(DxS(:,2),r(:))
			xsi= xsi+ Dxsi
			eta= eta+ Deta
			CALL Serendip_func(Ni,xsi,eta,ldim,Nodel,inci)
			xS(1)= DOT_PRODUCT(Ni(:),Elcor(1,:))
			xS(2)= DOT_PRODUCT(Ni(:),Elcor(2,:))
			xS(3)= DOT_PRODUCT(Ni(:),Elcor(3,:))
			r= xS- xPi
			DistPS_N=Dist(xPi(:),xS(:),Cdim)
			IF(ABS((DistPS- DistPS_N)/DistPS_N) < 0.01)EXIT
		END DO
		IF(xsi >= -1.0 .and. xsi <= 1.0 .and. eta >= -1.0 .and. eta <= 1.0)THEN
			Min_dist=DistPS_N
			RETURN
		ELSE IF(xsi < -1.0 .and. eta < -1.0)THEN
			Min_dist=Dist(Elcor(:,1),xPi(:),Cdim)
			RETURN
		ELSE IF(xsi > 1.0 .and. eta > 1.0)THEN
			Min_dist=Dist(Elcor(:,3),xPi(:),Cdim)
			RETURN
		ELSE IF(xsi > 1.0 .and. eta < -1.0)THEN
			Min_dist=Dist(Elcor(:,2),xPi(:),Cdim)
			RETURN
		ELSE IF(xsi < -1.0 .and. eta > 1.0)THEN
			Min_dist=Dist(Elcor(:,4),xPi(:),Cdim)
			RETURN
		ELSE IF(xsi >= -1.0 .and. xsi <= 1.0 .and. eta < -1.0)THEN
			eta=-1.0
			CALL Serendip_func(Ni,xsi,eta,ldim,Nodel,inci)
			xS(1)= DOT_PRODUCT(Ni(:),Elcor(1,:))
			xS(2)= DOT_PRODUCT(Ni(:),Elcor(2,:))
			xS(3)= DOT_PRODUCT(Ni(:),Elcor(3,:))
			Min_dist=Dist(xS(:),xPi(:),Cdim)
			RETURN
		ELSE IF(xsi >= -1.0 .and. xsi <= 1.0 .and. eta > 1.0)THEN
			eta=1.0
			CALL Serendip_func(Ni,xsi,eta,ldim,Nodel,inci)
			xS(1)= DOT_PRODUCT(Ni(:),Elcor(1,:))
			xS(2)= DOT_PRODUCT(Ni(:),Elcor(2,:))
			xS(3)= DOT_PRODUCT(Ni(:),Elcor(3,:))
			Min_dist=Dist(xS(:),xPi(:),Cdim)
			RETURN
		ELSE IF(eta >= -1.0 .and. eta <= 1.0 .and. xsi > 1.0)THEN
			xsi=1.0
			CALL Serendip_func(Ni,xsi,eta,ldim,Nodel,inci)
			xS(1)= DOT_PRODUCT(Ni(:),Elcor(1,:))
			xS(2)= DOT_PRODUCT(Ni(:),Elcor(2,:))
			xS(3)= DOT_PRODUCT(Ni(:),Elcor(3,:))
			Min_dist=Dist(xS(:),xPi(:),Cdim)
			RETURN
		ELSE IF(eta >= -1.0 .and. eta <= 1.0 .and. xsi < -1.0)THEN
			xsi=-1.0
			CALL Serendip_func(Ni,xsi,eta,ldim,Nodel,inci)
			xS(1)= DOT_PRODUCT(Ni(:),Elcor(1,:))
			xS(2)= DOT_PRODUCT(Ni(:),Elcor(2,:))
			xS(3)= DOT_PRODUCT(Ni(:),Elcor(3,:))
			Min_dist=Dist(xS(:),xPi(:),Cdim)
			RETURN
		END IF
END SELECT

END FUNCTION Min_dist

REAL FUNCTION Min_dist1(Elcor,xPi,Nodel,inci,ELengx,Elenge,ldim)
IMPLICIT NONE
REAL,INTENT(IN)		::	Elcor(:,:)			!	Coordinates of Element
REAL,INTENT(IN)		::	xPi(:)				!	Coordinates of Collocation point
REAL,INTENT(IN)		::	ELengx,ELenge       !   Elementlength xsi and eta
REAL							::	DET,A,B,C,D,F1,F2	
REAL							::	xsi,eta,Dxsi,Deta
REAL							::	DistPS,DistPS_N
REAL							::	L
REAL							::	ERR
REAL,ALLOCATABLE  ::	Ni(:)				! Shape function
REAL,ALLOCATABLE	::	DNi(:,:)			! Derivatives of shape function
REAL,ALLOCATABLE	::	r(:),xS(:),Dxs(:,:)
INTEGER,INTENT(IN)::	Nodel
INTEGER,INTENT(IN)::	inci(:)
INTEGER,INTENT(IN)::	ldim
INTEGER							::	n,cdim
cdim= ldim + 1
ALLOCATE(Ni(Nodel),DNi(Nodel,ldim),r(Cdim),xS(Cdim),DxS(Cdim,ldim))
SELECT CASE(Cdim)
  CASE(2)
	  xsi=0.0
		CALL Serendip_func(Ni,xsi,eta,ldim,Nodel,inci)  
    CALL Cartesian(xS,Ni,ldim,Elcor)
		r= xS-xPi
		DistPS=Dist(xPi(:),xS(:),Cdim)
		IF(((DistPS-ELengx/2)/Elengx) > 4.)THEN
			Min_dist1=DistPS-Elengx/2
			RETURN
		END IF
    DO n=1,40
      IF(n > 1)DistPS= DistPS_N
      CALL Serendip_deriv(DNi,xsi,eta,ldim,nodel,inci)
      DxS(1,1)= DOT_PRODUCT(DNi(:,1),Elcor(1,:))
      DxS(2,1)= DOT_PRODUCT(DNi(:,1),Elcor(2,:))
      DET= DxS(1,1)**2+DxS(2,1)**2
      Dxsi= -1/DET * DOT_PRODUCT(DxS(:,1),r(:))
      xsi= xsi+ Dxsi
      IF(ABS(xsi) > 1.) xsi= xsi/ABS(xsi)
      CALL Serendip_func(Ni,xsi,eta,ldim,Nodel,inci)
      CALL Cartesian(xS,Ni,ldim,Elcor)
      r= xS- xPi
      DistPS_N= Dist(xPi(:),xS(:),Cdim)
      IF(DistPS_N > DistPS)THEN
!                xsi=xsi- Dxsi
        Min_dist1= DistPS
        RETURN
      END IF
      ERR= (DistPS- DistPS_N)/DistPS_N
      IF(ERR < 0.05)THEN 
        Min_dist1= DistPS_N
!	      WRITE(2,*)'n=',n
        RETURN
      END IF
		END DO
	CASE(3)
		xsi=0.0
		eta=0.0
		CALL Serendip_func(Ni,xsi,eta,ldim,Nodel,inci)
    CALL Cartesian(xS,Ni,ldim,Elcor)
		r= xS-xPi
		DistPS=Dist(xPi(:),xS(:),Cdim)
		IF(((DistPS-ELengx/2)/ELengx) > 4. .and.((DistPS-ELenge/2)/ELenge)> 4.)THEN
			Min_dist1=DistPS
			RETURN
		END IF
    DO n=1,40
			IF(n > 1)DistPS= DistPS_N
			CALL Serendip_deriv(DNi,xsi,eta,ldim,nodel,inci)
			DxS(1,1)= DOT_PRODUCT(DNi(:,1),Elcor(1,:))
			DxS(2,1)= DOT_PRODUCT(DNi(:,1),Elcor(2,:))
			DxS(3,1)= DOT_PRODUCT(DNi(:,1),Elcor(3,:))
			DxS(1,2)= DOT_PRODUCT(DNi(:,2),Elcor(1,:))
			DxS(2,2)= DOT_PRODUCT(DNi(:,2),Elcor(2,:))
			DxS(3,2)= DOT_PRODUCT(DNi(:,2),Elcor(3,:))
			A= DxS(1,1)**2+DxS(2,1)**2+DxS(3,1)**2
			B= DxS(1,1)*DxS(1,2)+DxS(2,1)*DxS(2,2)+DxS(3,1)*DxS(3,2)
			C=B
			D= DxS(1,2)**2+DxS(2,2)**2+DxS(3,2)**2
			DET= A*D - C*B
			F1= DOT_PRODUCT(DxS(:,1),r(:))
      F2= DOT_PRODUCT(DxS(:,2),r(:))
      Dxsi = -1/DET * (F1*D - F2*B)
			Deta = -1/DET * (F2*A - F1*C)
			xsi= xsi+ Dxsi
			eta= eta+ Deta
      IF(ABS(xsi) > 1.) xsi= xsi/ABS(xsi)
      IF(ABS(eta) > 1.) eta= eta/ABS(eta)
			CALL Serendip_func(Ni,xsi,eta,ldim,Nodel,inci)
      CALL Cartesian(xS,Ni,ldim,Elcor)
			r= xS- xPi
			DistPS_N=Dist(xPi(:),xS(:),Cdim)
      IF(DistPS_N > DistPS)THEN
      xsi=xsi- Dxsi
      eta=eta- Deta
      Min_dist1= DistPS
!		  WRITE(2,*)'n=',n
!     WRITE(2,*)'XSI=',xsi
!     WRITE(2,*)'ETA=',eta
      RETURN
    END IF
    ERR= (DistPS- DistPS_N)/DistPS_N
    IF(ERR < 0.05)THEN
      Min_dist1= DistPS_N
  !   WRITE(2,*)'n=',n
  !		WRITE(2,*)'XSI=',xsi
  !		WRITE(2,*)'ETA=',eta
      RETURN
    END IF
	END DO
END SELECT
DEALLOCATE(Ni,DNi,r,xS,DxS)
END FUNCTION Min_dist1


SUBROUTINE Elength(L,Elcor,nodes,ldim)
!------------------------------------------------
!   Computes the length of a boundary element
!----------------------------------------------
IMPLICIT NONE
REAL,INTENT (IN)			::	Elcor(:,:)
INTEGER, INTENT (IN)	::	nodes, ldim
REAL,INTENT (OUT)			::	L
REAL 									::	B(ldim+1), distB3, distB2, p, a, c
INTEGER								::	Cdim
Cdim=ldim+1
SELECT CASE (Nodes)
	CASE (2)
		L=Dist(Elcor(:,1),Elcor(:,2),Cdim)
		RETURN
	CASE (3)
		B=(Elcor(:,1)+Elcor(:,2))/2.0
		distB3=Dist(B(:),Elcor(:,3),Cdim)
		distB2=Dist(B,Elcor(:,2),Cdim)
		IF (distB3/distB2 < 0.01) THEN
			L=Dist(Elcor(:,1),Elcor(:,2),Cdim)
			RETURN
		END IF
		IF (distB3/distB2 < 0.1) THEN
			c=distB3/distB2
			L=2.0*distB2*(1.0+2.0/3.0*c**2-2.0/5.0*c**4)	!Length Parabel linearisiert
			RETURN
		END IF
		p=distB2**2/(2*DistB3)							!Parabel Parameter p=y**2/(2*x)
		a=SQRT(p**2+distB2**2)
		L=2.0*(distB2/(2.0*p)*a + p/2.0*LOG((distB2 + a)/p))	!Length Parabel exakt
		RETURN
	CASE DEFAULT
END SELECT
END SUBROUTINE Elength


END MODULE Geometry_lib

