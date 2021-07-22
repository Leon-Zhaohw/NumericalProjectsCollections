!     Last change:  G     3 May 100    6:28 pm

PROGRAM Trefftz
!---------------------------------
!   Program to compute the heat flow past a cylindrical isolator
!   in a 2-D infinite domain using the Trefftz method
!---------------------------------
USE Laplace_lib ; USE Geometry_lib;   ! use Laplace and Geometry libs
IMPLICIT NONE                         ! declare all variables
INTEGER(4)			 ::  response
REAL             ::  q                ! inflow/outflow
REAL             ::  k                ! Thermal conductivity
INTEGER          ::  npnts            ! Number of points P,Q
REAL             ::  rq               ! radius of isolator
REAL             ::  rp               ! radius of source points
REAL(KIND=8),ALLOCATABLE ::  Lhs(:,:) ! left hand side of equations (coeff. matrix)
REAL(KIND=8),ALLOCATABLE ::  Rhs(:)   ! right hand side of system of equations
REAL(KIND=8),ALLOCATABLE ::  F(:)     ! fictitious source intensity
REAL             ::  dxr(2)           ! rx/r ,ry/r
REAL             ::  dx(2)            ! delta x(P,Q)
REAL             ::  vnorm(2)         ! normal vector
REAL             ::  Delth,Thetq,Thetp,xq,yq,xp,yp,xi,yi,r,uq,qx,qy,dUxy(2)
INTEGER          ::  npq,npp,ninpts,nin

 OPEN(UNIT=10,FILE='INPUT.DAT',STATUS='OLD',ACTION='READ')
 OPEN(UNIT=11,FILE='OUTPUT.DAT',STATUS='UNKNOWN',ACTION='WRITE')
 READ(10,*) q,k,npnts,rq,rp
 WRITE(11,*) ' Program 2 : heat flow past a cylinder with Trefftz method'
 WRITE(11,*) '  Heat inflow/outflow=  ',q
 WRITE(11,*) '  Thermal conductivity= ',k
 WRITE(11,*) '  Number of Points P,Q= ',npnts
 WRITE(11,*) '  Radius of Isolator=   ',rq
 WRITE(11,*) '  Radius of Sources =   ',rp

 ALLOCATE (Lhs(npnts,npnts),Rhs(npnts),F(npnts))     ! allocate arrays for equations
 Delth= 2*Pi/npnts                                   ! increment in angle theta between points (radians)
 Thetq= Pi/2.0                                       ! angle theta to first field point Q1

 Field_points: DO npq= 1,npnts
                  Rhs(npq)= q * SIN(Thetq)           ! right hand side
                  xq= rq*COS(Thetq)                  ! x-coordinate of field point
                  yq= rq*SIN(Thetq)                  ! y-coordinate of field point
                  vnorm(1)= -COS(Thetq)              ! normal vector to Q
                  vnorm(2)= -SIN(Thetq)
                  Thetq= Thetq + Delth               ! angle to next field point Q
                  Thetp= Pi/2.0                      ! angle to first source point P1
 Source_points:   DO npp= 1,npnts
                     xp= rp*COS(Thetp)               ! x-coordinate of source point
                     yp= rp*SIN(Thetp)               ! y-coordinate of source point
											dxr(1)= xp-xq
											dxr(2)= yp-yq
											r= SQRT(dxr(1)**2 + dxr(2)**2)  ! distance between field and source point
                     dxr= dxr/r                      ! normalise vector dxr
											Lhs(npq,npp)= T(r,dxr,vnorm,2)  ! use function T from Laplace_lib
                     Thetp= Thetp + Delth            ! angle to next source point P
	          END DO  Source_points
	       END DO  Field_points

 Lhs= - Lhs   !Multiplication with minus because of negative pivots in Lhs
 Rhs= - Rhs   !Multiplication with minus because of negative pivots in Lhs

 !  Solve system of equations: calculate F out of Lhs and Rhs
 CALL Solve(Lhs,Rhs,F)

 !   Postprocessing - Boundary values of temperature
 WRITE(11,*)  ''
 WRITE(11,*)  'Temperatures at Boundary points:'
 Thetq= Pi/2.0                                       ! angle to first field point Q1
 Field_points1: DO npq= 1,npnts
                   uq= 0.0
                   xq= rq*COS(Thetq)                 ! x-coordinate of field point
                   yq= rq*SIN(Thetq)                 ! y-coordinate of field point
                   Thetq= Thetq + Delth              ! angle to next field point Q
                   Thetp= Pi/2.0                     ! angle to first source point P1
 Source_points1:   DO npp= 1,npnts
                      xp= rp*COS(Thetp)              ! x-coordinate of source point
                      yp= rp*SIN(Thetp)              ! y-coordinate of source point
                      dxr(1)= xp-xq
                      dxr(2)= yp-yq
                      r= SQRT(dxr(1)**2 + dxr(2)**2) ! distance between field and source point
                      uq= uq + U(r,k,2)*F(npp)       ! use function U from Laplace_lib
                      Thetp= Thetp + Delth           ! angle to next source point P
                   END DO  Source_points1
                   uq=uq-q/k*yq
                   WRITE(11,*) 'Temperature at field point',npq,' =',uq
                END DO  Field_points1

 !   Postprocessing - Interior points
 WRITE(11,*)  ''
 WRITE(11,*)  'Temperatures at interior points:'
 READ(10,*) ninpts                                   ! read number of interior points
 Int_points:    DO nin= 1,ninpts
                   READ(10,*) xi,yi                  ! coordinates of interior points
                   uq= 0.0
									 qx= 0.0
									 qy= 0.0
                   Thetp= Pi/2.0                     ! angle to first source point P1
 Source_points2:   DO npp= 1,npnts
                      xp= rp*COS(Thetp)              ! x-coordinate of source point
                      yp= rp*SIN(Thetp)              ! y-coordinate of source point
											dx(1)= (-xp+xi)
                      dx(2)= (-yp+yi)                      
                      r= SQRT(dx(1)**2 + dx(2)**2)
											dxr=dx/r	
                      uq= uq + U(r,k,2)*F(npp)       ! use function U from Laplace_lib
											dUxy= dU(r,dxr,2)	
											qx= qx + dUxy(1)*F(npp)					
											qy= qy + dUxy(2)*F(npp)					
                      Thetp= Thetp + Delth           ! angle to next source point P
                   END DO  Source_points2
                   uq=uq-q/k*yi
										qy=qy+1.	
                   WRITE(11,*) 'Temperature at x=',xi,', y=',yi,' =',uq
									 WRITE(11,*) 'Flow qx at     x=',xi,', y=',yi,' =',qx
									 WRITE(11,*) 'Flow qy at     x=',xi,', y=',yi,' =',qy
									 WRITE(11,*) ' '
                END DO  Int_points
END PROGRAM Trefftz

