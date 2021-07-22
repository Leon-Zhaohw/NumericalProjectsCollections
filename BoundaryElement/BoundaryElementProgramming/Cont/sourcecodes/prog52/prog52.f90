!     Last change:  G     3 May 100    8:24 am
PROGRAM Direct_Method
!---------------------------------
!   Program to compute the heat flow  past a cylindrical isolator
!   in an 2-D infinite domain using the direct BE method
!   with constant line segments
!---------------------------------
USE Utility_lib ;USE Geometry_lib            !   subroutine to solve equations
INTEGER(4)			 ::  response
REAL          ::  q         !   inflow/outflow
REAL          ::  k         !   Thermal conductivity
INTEGER       ::  nseg      !   Number of segments
REAL          ::  rq        !   radius of isolator (inner)
REAL          ::  rqo       !   radius of isolator (outer)
REAL(KIND=8),ALLOCATABLE :: Lhs(:,:) !   left hand side coeff. matrix (DT)
REAL(KIND=8),ALLOCATABLE :: F(:)     !   right hand side vector {F}
REAL(KIND=8),ALLOCATABLE :: u(:)     !   Temperatures at segment centers
REAL,ALLOCATABLE :: Rhs(:,:)         !   right hand side coeff. matrix (DU)
REAL,ALLOCATABLE :: t0(:)            !   Applied flows
REAL,ALLOCATABLE :: xA(:,:),xB(:,:)  !   Start/end coordinates of segments
REAL,ALLOCATABLE :: xS(:,:)          !   Coords of points Pi
REAL,ALLOCATABLE :: Ve(:,:),Vn(:,:)  !  Vectors along and normal to segment
REAL :: PI=3.14159265359
REAL :: vrA(2),vrB(2)
REAL :: lens
C= 0.5/Pi
!    Read in job information
OPEN(UNIT=10,FILE='INPUT.DAT',STATUS='OLD',ACTION='READ')
OPEN(UNIT=11,FILE='OUTPUT.DAT',STATUS='UNKNOWN',ACTION='WRITE')
READ(10,*) q,k,nseg,rq
WRITE(11,*) 'Heat flow past a cylinder (direct BE method)'
WRITE(11,*) 'Input values:'
WRITE(11,*) ' Heat inflow/outflow= ',q
WRITE(11,*) ' Thermal conductivity=',k
WRITE(11,*) ' Radius of Isolator=  ',rq
WRITE(11,*) ' Number of segments=  ',nseg
!   allocate arrays
ALLOCATE (Lhs(nseg,nseg),Rhs(nseg,nseg),F(nseg))
ALLOCATE (xA(2,nseg),xB(2,nseg),t0(nseg),u(nseg))
ALLOCATE (xS(2,nseg),ve(2,nseg),vn(2,nseg))
C1=0.5/(Pi*k)
Delth= 2.0*Pi/nseg     !     increment in angle theta
rqo=rq/COS(Delth/2.0)  !     outer radius of isolator
Thet= (Pi-Delth)/2.0
!   Compute start/end  coordinates of segments
xA(1,1)= rqo*COS(Thet)
xA(2,1)= rqo*SIN(Thet)
Segments: &
DO ns= 1,nseg-1
	Thet= Thet + Delth
	xB(1,ns)= rqo*COS(Thet)
	xB(2,ns)= rqo*SIN(Thet)
	xA(1,ns+1)= xB(1,ns)
	xA(2,ns+1)= xB(2,ns)
END DO &
Segments
xB(1,nseg)= xA(1,1)
xB(2,nseg)= xA(2,1)
!   Compute center coordinates of segments (collocation point coords)
Segments1: &
DO ns= 1,nseg
	xS(1,ns)= (xA(1,ns) + xB(1,ns))/2.0
	xS(2,ns)= (xA(2,ns) + xB(2,ns))/2.0
END DO &
Segments1
!   Compute applied tractions at centers of elements
Thet= Pi/2.0
Segments2: &
DO ns= 1,nseg
	t0(ns)= q*SIN(Thet)
	Thet= Thet + Delth
END DO  &
Segments2
!    Assemble matrices DT and DU
Segments3: &
DO ns=1,nseg
	lens= dist(xA(:,ns),xB(:,ns),2)
!    Vector normal to segment A-B
        dx= xA(1,ns) - xB(1,ns)
        dy= xA(2,ns) - xB(2,ns)
        ve(1,ns)= dx/lens
        ve(2,ns)= dy/lens
        vn(1,ns)= ve(2,ns)
        vn(2,ns)=-ve(1,ns)
        Points_Pi: &
        DO np=1,nseg
!           WRITE(11,*) 'Point',np
           rA= Dist(xA(:,ns),xS(:,np),2)
           rB= Dist(xB(:,ns),xS(:,np),2)
           vrA(1)= xA(1,ns)- xS(1,np)
           vrA(2)= xA(2,ns)- xS(2,np)
           vrB(1)= xB(1,ns)- xS(1,np)
           vrB(2)= xB(2,ns)- xS(2,np)
           COSThA= -DOT_PRODUCT(vn(:,ns),vrA)/rA
           COSThB= -DOT_PRODUCT(vn(:,ns),vrB)/rB
           SINThA= -DOT_PRODUCT(ve(:,ns),vrA)/rA
           SINThB= -DOT_PRODUCT(ve(:,ns),vrB)/rB
           ThetA= ACOS(COSThA)*SIGN(1.0,SinThA)
           ThetB= ACOS(COSThB)*SIGN(1.0,SinThB)
	   IF(np == ns) THEN         ! Diagonal coefficients
		Lhs(np,np)=  0.5
		Rhs(np,np)= lens*C1*(LOG(lens/2.0)-1.0)
	   ELSE                      ! off-diagonal coeff.
		Lhs(np,ns)= C*(ThetB-ThetA)
		Rhs(np,ns)= C1*(rB*SINThB*(LOG(rB)-1)+ThetB*rB*COSThB &
		              - rA*SINThA*(LOG(rA)-1)-ThetA*rA*COSThA)
	   END IF
	END DO  &
	Points_Pi
END DO &
Segments3
F= MATMUL(Rhs,t0)    !    compute right hand side vector
CALL Solve(Lhs,F,u)  !    solve system of equations
!    output computed temperatures
WRITE(11,*) 'Temperatures at segment centers:'
Segments4: &
DO ns= 1,nseg
 WRITE(11,'(A,I5,A,F10.3)') ' Segment',ns,'  T=',u(ns)-q/k*xS(2,ns)
END DO &
Segments4
DEALLOCATE (xS)
!   Compute Temperatures and flows at interior points
READ(10,*,IOSTAT=IOS) NPoints
IF(NPoints == 0 .OR. IOS /= 0) THEN
! PAUSE 'program Finshed'
 STOP
END IF
ALLOCATE (xS(2,NPoints))
WRITE(11,*) 'Temperatures(T) and flow (q-x,q-y) at interior points:'
DO n=1,NPoints
READ(10,*) xS(1,n),xS(2,n)
END DO
Interior_points: &
DO np=1,Npoints
    up= 0.0
    qx= 0.0
    qy= 0.0
    Segments5 : &
    DO ns=1,nseg
       rA= Dist(xA(:,ns),xS(:,np),2)
       rB= Dist(xB(:,ns),xS(:,np),2)
       vrA(1)= xA(1,ns)- xS(1,np)
       vrA(2)= xA(2,ns)- xS(2,np)
       vrB(1)= xB(1,ns)- xS(1,np)
       vrB(2)= xB(2,ns)- xS(2,np)
       COSThA= -DOT_PRODUCT(vn(:,ns),vrA)/rA
       COSThB= -DOT_PRODUCT(vn(:,ns),vrB)/rB
       SINThA= -DOT_PRODUCT(ve(:,ns),vrA)/rA
       SINThB= -DOT_PRODUCT(ve(:,ns),vrB)/rB
       H= RA*CosThA
       ThetA= ACOS(COSThA)*SIGN(1.0,SinThA)
       ThetB= ACOS(COSThB)*SIGN(1.0,SinThB)
       IF(ThetB-ThetA > Pi) ThetA= 2.0*Pi + ThetA
       dT= C*(ThetB-ThetA)
       dU= C1*(rB*SINThB*(LOG(rB)-1)+ThetB*rB*COSThB &
	     - rA*SINThA*(LOG(rA)-1)-ThetA*rA*COSThA)
       dSx= C/k*(ThetB-ThetA)
       Fact= CosthB/CosthA
       IF(Fact > 0.0) THEN
        dSy= -C/k*LOG(Fact)
       ELSE
        dSy= 0.
       END IF
       dRx= -C/H*(costhB*SINThB - cosThA*sinThA)
       dRy= C/H*(costhB**2 - cosThA**2)
       up= up + dU*t0(ns) - dT*u(ns)
       qxp= -k*(dSx*t0(ns)-dRx*u(ns))         !   q-x'
       qyp= -k*(dSy*t0(ns)-dRy*u(ns))         !   q-y'
       qx= qx + qxp*vn(1,ns) - qyp*vn(2,ns)
       qy= qy + qxp*vn(2,ns) + qyp*vn(1,ns)
    END DO &
    Segments5
    Up= Up - q/k*xS(2,np)
    qy= qy + q
    WRITE(11,'(5(A,F10.3))') &
    'x=',xS(1,np),', y=',xS(2,np),', T=',up,',  q-x=',qx,',  q-y=',qy
END DO  &
Interior_points

END PROGRAM Direct_Method

