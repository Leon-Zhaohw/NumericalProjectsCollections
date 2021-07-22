!     Last change:  G     6 Dec 99    7:15 pm
MODULE Elast_lib
USE Laplace_lib
!    Library for elasticity problems
IMPLICIT NONE
CONTAINS

FUNCTION UK(dxr,r,E,ny,Cdim)
!--------------------------------------------
!   FUNDAMENTAL SOLUTION FOR DISPLACEMENTS
!   isotropic material (Kelvin solution)
!--------------------------------------------
IMPLICIT NONE
REAL,INTENT(IN)			:: dxr(:)        !   rx/r etc.
REAL,INTENT(IN)			:: r             !   r
REAL,INTENT(IN)			:: E             !   Young's modulus
REAL,INTENT(IN)			:: ny            !   Poisson's ratio
INTEGER,INTENT(IN)	:: Cdim          !   Cartesian dimension
REAL								:: UK(Cdim,Cdim)  !   Function returns array of same dim as dxr
REAL								:: G,c,c1,onr,clog,conr  !   Temps
G= E/(2.0*(1+ny))
c1= 3.0 - 4.0*ny
SELECT CASE (Cdim)
	CASE (2)											 !     Two-dimensional solution
		c= 1.0/(8.0*Pi*G*(1.0 - ny))
		clog= -c1*LOG(r)
		UK(1,1)= c*(clog + dxr(1)*dxr(1))
		UK(1,2)= c*dxr(1)*dxr(2)
		UK(2,2)= c*(clog + dxr(2)*dxr(2))
		UK(2,1)= UK(1,2)
	CASE(3)												!      Three-dimensional solution
		c= 1.0/(16.0*Pi*G*(1.0 - ny))
		conr=c/r
		UK(1,1)= conr*(c1 + dxr(1)*dxr(1))
		UK(1,2)= conr*dxr(1)*dxr(2)
		UK(1,3)= conr*dxr(1)*dxr(3)
		UK(2,1)= UK(1,2)
		UK(2,2)= conr*(c1 + dxr(2)*dxr(2))
		UK(2,3)= conr*dxr(2)*dxr(3)
		UK(3,1)= UK(1,3)
		UK(3,2)= UK(2,3)
		UK(3,3)= conr*(c1 + dxr(3)*dxr(3))
	CASE DEFAULT
END SELECT
RETURN
END FUNCTION UK

FUNCTION TK(dxr,r,Vnor,ny,Cdim)
!--------------------------------------------
!   FUNDAMENTAL SOLUTION FOR TRACTIONS
!   isotropic material (Kelvin solution)
!--------------------------------------------
IMPLICIT NONE
REAL,INTENT(IN)			:: dxr(:)							!   rx/r etc.
REAL,INTENT(IN)			:: r									!   r
REAL,INTENT(IN)			:: Vnor(:)						!   normal vector
REAL,INTENT(IN)			:: ny									!   Poisson's ratio
INTEGER,INTENT(IN)	:: Cdim								!   Cartesian dimension
REAL								:: TK(Cdim,Cdim)			!   Function returns array of same dim as dxr
REAL		            :: c2,c3,costh,Conr   !   Temps
c3= 1.0 - 2.0*ny
Costh= DOT_PRODUCT (Vnor,dxr)
SELECT CASE (Cdim)
	CASE (2)        
		c2= 1.0/(4.0*Pi*(1.0 - ny))
		Conr= c2/r
		TK(1,1)= -(Conr*(c3 + 2.0*dxr(1)*dxr(1))*Costh)
		TK(1,2)= -(Conr*(2.0*dxr(1)*dxr(2)*Costh + c3*(Vnor(1)*dxr(2) - Vnor(2)*dxr(1))))
		TK(2,2)= -(Conr*(c3 + 2.0*dxr(2)*dxr(2))*Costh)
		TK(2,1)= -(Conr*(2.0*dxr(1)*dxr(2)*Costh + c3*(Vnor(2)*dxr(1) - Vnor(1)*dxr(2))))
	CASE(3)           !    Three-dimensional
		c2= 1.0/(8.0*Pi*(1.0 - ny))
		Conr= c2/r**2
		TK(1,1)= -Conr*(c3 + 3.0*dxr(1)*dxr(1))*Costh
		TK(1,2)= -Conr*(3.0*dxr(1)*dxr(2)*Costh - c3*(Vnor(2)*dxr(1) - Vnor(1)*dxr(2)))	
		TK(1,3)= -Conr*(3.0*dxr(1)*dxr(3)*Costh - c3*(Vnor(3)*dxr(1) - Vnor(1)*dxr(3)))
		TK(2,1)= -Conr*(3.0*dxr(1)*dxr(2)*Costh - c3*(Vnor(1)*dxr(2) - Vnor(2)*dxr(1)))
		TK(2,2)= -Conr*(c3 + 3.0*dxr(2)*dxr(2))*Costh
		TK(2,3)= -Conr*(3.0*dxr(2)*dxr(3)*Costh - c3*(Vnor(3)*dxr(2) - Vnor(2)*dxr(3)))
		TK(3,1)= -Conr*(3.0*dxr(1)*dxr(3)*Costh - c3*(Vnor(1)*dxr(3) - Vnor(3)*dxr(1)))
		TK(3,2)= -Conr*(3.0*dxr(2)*dxr(3)*Costh - c3*(Vnor(2)*dxr(3) - Vnor(3)*dxr(2)))
		TK(3,3)= -Conr*(c3 + 3.0*dxr(3)*dxr(3))*Costh
	CASE DEFAULT
END SELECT
END FUNCTION TK

SUBROUTINE SK(TS,DXR,R,C2,C3)
!------------------------------------------------------------
!   KELVIN SOLUTION FOR STRESS  
!    TO BE MULTIPLIED WITH T
!------------------------------------------------------------
REAL, INTENT(OUT) :: TS(:,:)   ! Fundamental solution
REAL, INTENT(IN)  :: DXR(:)    ! rx , ry, rz
REAL, INTENT(IN)  :: R         ! r
REAL, INTENT(IN)  :: C2,C3     ! Elastic constants
INTEGER ::  Cdim      !   Cartesian dimension
INTEGER :: NSTRES  !  No. of stress components
INTEGER :: II(6), JJ(6) !  sequence of stresses in pseudo-vector
REAL    :: A
INTEGER :: I,N,J,K
Cdim= UBOUND(DXR,1)
IF(CDIM == 2) THEN
   NSTRES= 3
   II(1:3)= (/1,2,1/)
   JJ(1:3)= (/1,2,2/)   
ELSE
   NSTRES= 6
   II= (/1,2,3,1,2,3/)
   JJ= (/1,2,3,2,3,1/)
END IF
Coor_directions:&
DO K=1,Cdim
	Stress_components:&
	DO N=1,NSTRES
				I= II(N)
				J= JJ(N)
				A= 0.
				IF(K .EQ. J) A= A + DXR(I)
				IF(I .EQ. J) A= A - DXR(K)
				IF(K .EQ. I) A= A + DXR(J)
				A= A*C3
				TS(N,K)= C2/R*(A + Cdim*DXR(I)*DXR(J)*DXR(K))
				IF(Cdim .EQ. 3) TS(N,K)= TS(N,K)/R
	END DO &
	Stress_components
END DO &
Coor_directions
RETURN
END SUBROUTINE SK

SUBROUTINE RK(US,DXR,R,VNORM,C3,C5,C6,C7,ny)
!------------------------------------------------------------
!    KELVIN SOLUTION FOR STRESS COMPUTATION
!    TO BE MULTIPLIED WITH U
!------------------------------------------------------------
REAL, INTENT(OUT) :: US(:,:)        ! Fundamental solution
REAL, INTENT(IN)  :: DXR(:)         ! rx , ry, rz
REAL, INTENT(IN)  :: R              ! r
REAL, INTENT(IN)  :: VNORM(:)       ! nx , ny , nz
REAL, INTENT(IN)  :: C3,C5,C6,C7,ny ! Elastic constants
INTEGER ::  Cdim   !   Cartesian dimension
INTEGER :: NSTRES  !   No. of stress components
INTEGER :: II(6), JJ(6) !  sequence of stresses in pseudo-vector
REAL    :: costh, B,C,cny
INTEGER :: I,N,J,K
Cdim= UBOUND(DXR,1) 
IF(CDIM == 2) THEN
   NSTRES= 3
   II(1:3)= (/1,2,1/)
   JJ(1:3)= (/1,2,2/)   
ELSE
	NSTRES= 6                   
	II= (/1,2,3,1,2,3/)
	JJ= (/1,2,3,2,3,1/)
END IF  
COSTH= DOT_Product(dxr,vnorm) 
Cny= Cdim*ny
Coor_directions:&
DO K=1,Cdim
	Stress_components:&
	DO N=1,NSTRES 
				I= II(N)
				J= JJ(N)
				B= 0.
				C= 0.
				IF(I .EQ. J) B= Cdim*C3*DXR(K)
				IF(I .EQ. K) B= B + Cny*DXR(J)
				IF(J .EQ. K) B= B + Cny*DXR(I)
				B= COSTH *(B - C6*DXR(I)*DXR(J)*DXR(K))
				C= DXR(J)*DXR(K)*Cny
				IF(J .EQ.K) C= C + C3
				C= C*VNORM(I)
				B= B+C
				C= DXR(I)*DXR(K)*Cny
				IF(I .EQ. K) C=C + C3
				C= C*VNORM(J)
				B= B+C
				C= DXR(I)*DXR(J)*Cdim*C3
				IF(I .EQ. J) C= C - C7
				C= C*VNORM(K)
				US(N,K)= (B + C)*C5/R/R
				IF(Cdim .EQ. 3) US(N,K)= US(N,K)/R
	END DO &
	Stress_components
END DO &
Coor_directions 
RETURN
END SUBROUTINE RK

SUBROUTINE Trans_mat(v1,v2,v3, T)
!-----------------------------------------------
!  Computes Stress Transformation Matrix in 3-D
!----------------------------------------------
IMPLICIT NONE
REAL, INTENT(IN)     :: v1(3),v2(3),v3(3)       !  unit vectors in orthogonal directions
REAL, INTENT(OUT)    :: T(6,6)                  !  transformation matrix
REAL                 :: v1x,v1y,v1z,v2x,v2y,v2z,v3x,v3y,v3z   ! temps
v1x= v1(1) ; v1y= v1(2) ; v1z= v1(3)
v2x= v2(1) ; v2y= v2(2) ; v2z= v2(3)
v3x= v3(1) ; v3y= v3(2) ; v3z= V3(3)
!   T?11
T(1,1)= v1x**2  ; T(1,2)= v2x**2 ; T(1,3)= v3x**2
T(2,1)= v1y**2 ;  T(2,2)= v2y**2 ; T(2,3)= v3y**2
T(3,1)= v1z**2 ;  T(3,2)= v2z**2 ; T(3,3)= v3z**2
!   T?21
T(1,4)= 2.0*v1y*v1x ; T(1,5)= 2.0*v2y*v2x ; T(1,6)= 2.0*v3y*v3x
T(2,4)= 2.0*v1y*v1z ; T(2,5)= 2.0*v2y*v2z ;  T(2,6)= 2.0*v3y*v3z
T(3,4)= 2.0*v1x*v1z ; T(3,5)= 2.0*v2x*v2z ;  T(3,6)= 2.0*v3x*v3z
!   T?12
T(4,1)= v1x*v2x ; T(4,2)= v2x*v3x ; T(4,3)= v1x*v3x
T(5,1)= v1y*v2y ; T(5,2)= v2y*v3y ; T(5,3)= v1y*v3y
T(6,1)= v1z*v2z ; T(6,2)= v2z*v3z ; T(6,3)= v1z*v3z
!   T?22
T(4,4)= v1x*v2y+v1y*v2x ; T(4,5)= v2x*v3y+v2y*v3x ; T(4,6)= v1x*v3y+v1y*v3x
T(5,4)= v1y*v2z+v1z*v2y ; T(5,5)= v2y*v3z+v2z*v3y ; T(5,6)= v1y*v3z+v1z*v3y
T(6,4)= v1x*v2z+v1z*v2x ; T(6,5)= v2x*v3z+v2z*v3x ; T(6,6)= v1x*v3z+v1z*v3x
RETURN
END  SUBROUTINE Trans_mat
 
SUBROUTINE D_mat(E,ny,D,Cdim)
!-----------------------------------
!   Computes isotropic D-matrix
!   Plane-strain (Cdim= 2)
!   or 3-D       (Cdim= 3)
!-----------------------------------
IMPLICIT NONE
REAL, INTENT(IN)   :: E      !  Young's modulus
REAL, INTENT(IN)   :: ny     !  Poisson's ratio
INTEGER,INTENT(IN) :: Cdim   !  Cartesian Dimension
REAL, INTENT(OUT)  :: D(:,:) !  D-matrix
REAL               :: c1,c2,G
c1= E*(1.0-ny)/( (1.0+ny)*(1.0-2.0*ny) )
c2= ny/(1.0-ny)
G = E/(2.0*(1.0+ny))
D = 0.0
SELECT CASE (Cdim)
	CASE (2)
		D(1,1)= 1.0  ; D(2,2)= 1.0
		D(2,1)= c2   ; D(1,2)= c2
		D(3,3)= G/c1
	CASE (3)         !    3-D
		D(1,1)= 1.0  ;  D(2,2)= 1.0   ;  D(3,3)= 1.0
		D(2,1)= c2   ;  D(1,3)= c2   ;  D(2,3)= c2
		D(1,2)= c2   ;  D(3,1)= c2   ;  D(3,2)= c2
		D(4,4)= G/c1 ;  D(5,5)= G/c1 ;  D(6,6)= G/c1
	CASE DEFAULT
END SELECT
D= c1*D
RETURN
END SUBROUTINE D_mat

SUBROUTINE D_mat_anis(D,E1,G1,E2,G2,ny2,Cdim)
!-----------------------------------
!   Computes an-isotropic D-matrix
!   Plane-strain (Cdim= 2)
!   or 3-D       (Cdim= 3)
!-----------------------------------
IMPLICIT NONE
REAL, INTENT(OUT) :: D(:,:)!  D-matrix
REAL, INTENT(IN)  :: E1    ! Young's modulus, dir 1
REAL, INTENT(IN)  :: G1    ! Shear modulus , dir 1
REAL, INTENT(IN)  :: E2    ! Young's modulus, dir 2
REAL, INTENT(IN)  :: G2    ! Shear Modulus , dir 2
REAL, INTENT(IN)  :: ny2   ! Poisson's ratio, dir 2
INTEGER,INTENT(IN):: Cdim  !  Cartesian Dimension
REAL              :: n     !  ratio E1/E2
REAL :: cc,c1,c2,c3,c4,ny1   !   temps
ny1= 0.5*E1/G1 -1.0
n= E1/E2
cc= E2/(1.+ny1)/(1.-ny1-2.*n*ny2**2)
c1= n*(1.-n*ny2**2)*cc
c3= n*ny2*(1.0+ny1)*cc
c4= (1 - ny1**2)*cc
D= 0. !  only nonzero components of D are assigned
SELECT CASE (Cdim)
	CASE (2)          !    plane strain
		D(1,1)= c1 ; D(2,2)= c4
		D(1,2)= c3 ; D(2,1)= c3
		D(3,3)= G2
	CASE (3)          !    3-D
		c2= n*(ny1+n*ny2**2)*cc
		D(1,1)= C1 ; D(2,2)= c1 ; D(3,3)= c4
		D(1,2)= C2 ; D(1,3)= c3 ; D(2,3)= C3
		D(2,1)= C2 ; D(3,1)= c3 ; D(3,2)= C3
		D(4,4)= G1 ; D(5,5)= G2 ; D(6,6)= G2
	CASE DEFAULT
END SELECT
RETURN
END SUBROUTINE D_mat_anis

END MODULE Elast_lib
