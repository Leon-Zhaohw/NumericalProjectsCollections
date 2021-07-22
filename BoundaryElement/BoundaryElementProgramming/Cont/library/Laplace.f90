!     Last change:  D     1 Dec 1999    3:32 pm
MODULE Laplace_lib
REAL :: Pi= 3.14159265359
CONTAINS
  COMPLEX FUNCTION UW(r,k,Cdim)
  !-------------------------------
  !   Fundamental solution for scalar wave equation
  !   Pressure
  !------------------------------
  REAL,INTENT(IN)     ::  r  	!   Distance between source and field point
  REAL,INTENT(IN)     ::  k    	!   wave number
  INTEGER,INTENT(IN)  :: Cdim   !   Cartesian dimension (2-D,3-D)
  COMPLEX :: C0, Hankel0
  SELECT CASE (CDIM)
     CASE (2)          	!  Two-dimensional solution
        C0= CMPLX(0,0.25)
 !       UW= C0*Hankel0(k*r)
     CASE (3)          	!  Three-dimensional solution
        C0= CMPLX(0,k*r)
        UW= 1/(4.0*k*r)*EXP(C0)
     CASE DEFAULT
        UW=0.0
        WRITE (11,*)'Cdim not equal 2 or 3 in Function U(...)'
  END SELECT
  END FUNCTION UW
  REAL FUNCTION U(r,k,Cdim)
  !-------------------------------
  !   Fundamental solution for Potential problems
  !   Temperature/Potential
  !------------------------------
  REAL,INTENT(IN)     ::  r  	!   Distance between source and field point
  REAL,INTENT(IN)     ::  k    	!   Conducivity
  INTEGER,INTENT(IN)  :: Cdim   !   Cartesian dimension (2-D,3-D)
  SELECT CASE (CDIM)
     CASE (2)          	!  Two-dimensional solution
        U= 1.0/(2.0*Pi*k)*LOG(1/r)
     CASE (3)          	!  Three-dimensional solution
        U= 1.0/(4.0*Pi*r*k)
     CASE DEFAULT
        U=0.0
        WRITE (11,*)'Cdim not equal 2 or 3 in Function U(...)'
  END SELECT
  END FUNCTION U
 REAL FUNCTION T(r,dxr,Vnorm,Cdim)
 !-------------------------------
 !   Fundamental solution for Potential problems
 !   Normal gradient
 !------------------------------
 REAL,INTENT(IN)::            r	!   Distance between source and field point
 REAL,INTENT(IN)::       dxr(:)	!   rx/r , ry/r , rz/r
 REAL,INTENT(IN)::     Vnorm(:)	!   Normal vector
 INTEGER,INTENT(IN) ::    Cdim	!   Cartesian dimension
 SELECT CASE (Cdim)
    CASE (2)           	!  Two-dimensional solution
      T= -DOT_PRODUCT (Vnorm,dxr)/(2.0*Pi*r)
    CASE (3)           	!  Three-dimensional solution
      T= -DOT_PRODUCT (Vnorm,dxr)/(4.0*Pi*r*r)
    CASE DEFAULT
      T=0.0
      WRITE (11,*)'Cdim not equal 2 or 3 in Function U(...)'
 END SELECT
 END FUNCTION T
 FUNCTION dU(r,dxr,Cdim)
  !-------------------------------
  !   Derivatives of Fundamental solution for Potential problems
  !   Temperature/Potential
  !------------------------------
  REAL,INTENT(IN)::       r    !   Distance between source and field point
  REAL,INTENT(IN)::  dxr(:)    !   Distances in Cartesian directions divided by r
  REAL :: dU(UBOUND(dxr,1))    !   dU is array of same dim as dxr
  INTEGER ,INTENT(IN)            :: Cdim   !   Cartesian dimension (2-D,3-D)
  REAL :: C
  SELECT CASE (CDIM)
     CASE (2)           !  Two-dimensional solution
      C=1/(2.0*Pi*r)
      dU(1)= -C*dxr(1)
      dU(2)= -C*dxr(2)
     CASE (3)           !  Three-dimensional solution
      C=1/(4.0*Pi*r**2)
      dU(1)= C*dxr(1)
      dU(2)= C*dxr(2)
      dU(3)= C*dxr(3)
     CASE DEFAULT
 END SELECT
 END FUNCTION dU
 FUNCTION dT(r,dxr,Vnorm,Cdim)
 !-------------------------------
 !   derivatives of the Fundamental solution for Potential problems
 !   Normal gradient
 !------------------------------
 INTEGER,INTENT(IN) :: Cdim     !   Cartesian dimension
 REAL,INTENT(IN)::        r     !   Distance between source and field point
 REAL,INTENT(IN)::    dxr(:)    !   Distances in Cartesian directions divided by R
 REAL,INTENT(IN)::  Vnorm(:)    !   Normal vector
 REAL :: dT(UBOUND(dxr,1))      !   dT is array of same dim as dxr 
 REAL :: C,COSTH
 
COSTH= DOT_PRODUCT (Vnorm,dxr)
SELECT CASE (Cdim)
    CASE (2)           !  Two-dimensional solution
     C= 1/(2.0*Pi*r**2)
			dT(1)= C*(2.*dxr(1)*COSTH - Vnorm(1))	
			dT(2)= C*(2.*dxr(2)*COSTH - Vnorm(2))	
!     dT(1)= C*COSTH*dxr(1)
!     dT(2)= C*COSTH*dxr(2)
    CASE (3)           !  Three-dimensional solution
     C= 3/(4.0*Pi*r**3)
     dT(1)= C*COSTH*dxr(1)
     dT(2)= C*COSTH*dxr(2)
     dT(3)= C*COSTH*dxr(3)
    CASE DEFAULT
 END SELECT
 END FUNCTION dT
END MODULE Laplace_lib

