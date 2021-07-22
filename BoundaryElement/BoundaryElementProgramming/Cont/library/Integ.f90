!   Last change: CD  27 Mar 2000  2:10 pm
MODULE Integration_lib
USE Geometry_lib; USE Laplace_lib; USE Elast_lib
CONTAINS
SUBROUTINE Gauss_coor(Cor,Wi,Intord)
!------------------------------------
!   Returns Gauss coordinates and Weights for up to 8 Gauss points
!------------------------------------
IMPLICIT NONE
REAL, INTENT(OUT)  :: Cor(8)  !  Gauss point coordinates
REAL, INTENT(OUT)  :: Wi(8)   !  weigths
INTEGER,INTENT(IN) :: Intord  !  integration order
SELECT CASE (Intord)
 CASE (1)
  Cor(1)= 0.
  Wi(1) = 2.0
 CASE(2)
  Cor(1)= .577350269 ; Cor(2)= -Cor(1)
  Wi(1) = 1.0 ;  Wi(2) = Wi(1)
 CASE(3)
  Cor(1)= .774596669 ; Cor(2)= 0.0 ; Cor(3)= -Cor(1)
  Wi(1) = .555555555 ; Wi(2) = .888888888 ; Wi(3) = Wi(1)
 CASE(4)
  Cor(1)= .861136311 ; Cor(2)= .339981043 ; Cor(3)= -Cor(2) ; Cor(4)= -Cor(1)
  Wi(1) = .347854845 ; Wi(2) = .652145154 ; Wi(3) = Wi(2) ; Wi(4) = Wi(1)
 CASE(5)
  Cor(1)= .9061798459 ; Cor(2)= .5384693101 ; Cor(3)= .0 ; Cor(4)= -Cor(2)
  Cor(5)= -Cor(1)
  Wi(1) = .236926885 ; Wi(2) = .478628670 ; Wi(3) = .568888888 ; Wi(4) = Wi(2)
  Wi(5) = Wi(1)
 CASE(6)
  Cor(1)= .932469514 ; Cor(2)= .661209386 ; Cor(3)= .238619186
  Cor(4)= -Cor(3) ;  Cor(5)= -Cor(2) ; Cor(6)= -Cor(1)
  Wi(1) = .171324492 ; Wi(2) = .360761573 ; Wi(3) = .467913934
  Wi(4) = Wi(3) ; Wi(5) = Wi(2) ; Wi(6) = Wi(1)
 CASE(7)
  Cor(1)= .949107912 ; Cor(2)= .741531185 ; Cor(3)= .405845151
  Cor(4)= 0.
  Cor(5)= -Cor(3) ;Cor(6)= -Cor(2) ;Cor(7)= -Cor(1)
  Wi(1) = .129484966 ; Wi(2) = .279705391 ; Wi(3) = .381830050
  Wi(4) = .417959183
  Wi(5) = Wi(3) ; Wi(6) = Wi(2) ; Wi(7) = Wi(1)
 CASE(8)
  Cor(1)= .960289856 ; Cor(2)= .796666477 ; Cor(3)= .525532409 ; Cor(4)= .183434642
  Cor(5)= -Cor(4) ; Cor(6)= -Cor(3) ; Cor(7)= -Cor(2) ; Cor(8)= -Cor(1)
  Wi(1) = .101228536 ; Wi(2) = .222381034 ; Wi(3) = .313706645 ;Wi(4) = .362683783
  Wi(5) = Wi(4) ; Wi(6) = Wi(3) ; Wi(7) = Wi(2) ; Wi(8) = Wi(1)
 CASE DEFAULT
CALL Error_Message('Gauss points not in range 1-8')
END SELECT
END SUBROUTINE Gauss_coor

SUBROUTINE Gauss_Laguerre_coor(Cor,Wi,Intord)
!------------------------------------
!  Returns Gauss_Laguerre coordinates and Weights for up to 8 Gauss points
!------------------------------------
IMPLICIT NONE
REAL, INTENT(OUT)  :: Cor(8)  !  Gauss point coordinates
REAL, INTENT(OUT)  :: Wi(8)   !  weigths
INTEGER,INTENT(IN) :: Intord  !  integration order
SELECT CASE (Intord)
 CASE (1)
  Cor(1)= 0.5
  Wi(1) = 1.0
 CASE(2)
  Cor(1)= .112008806 ; Cor(2)=.602276908
  Wi(1) = .718539319 ; Wi(2) =.281460680
 CASE(3)
  Cor(1)= .063890793 ; Cor(2)= .368997063 ; Cor(3)= .766880303
  Wi(1) = .513404552 ; Wi(2) = .391980041 ; Wi(3) = .0946154065
 CASE(4)
  Cor(1)= .0414484801 ; Cor(2)=.245274914 ; Cor(3)=.556165453 ; Cor(4)= .848982394
  Wi(1) = .383464068 ; Wi(2) =.386875317 ; Wi(3) =.190435126 ; Wi(4) = .0392254871
 CASE(5)
  Cor(1)= .0291344721 ; Cor(2)= .173977213 ; Cor(3)= .411702520; Cor(4)=.677314174
  Cor(5)= .894771361
  Wi(1) = .297893471 ; Wi(2) = .349776226 ; Wi(3) =.234488290 ; Wi(4) = .0989304595
  Wi(5) = .0189115521
 CASE(6)
  Cor(1)= .0216340058 ; Cor(2)= .129583391 ; Cor(3)= .314020449
  Cor(4)= .538657217 ; Cor(5)= .756915337 ; Cor(6)=.922668851
  Wi(1) =  .238763662 ; Wi(2) =.308286573 ; Wi(3) =.245317426
  Wi(4) = .142008756 ; Wi(5) =.0554546223 ; Wi(6) =.0101689586
 CASE(7)
  Cor(1)= .0167193554 ; Cor(2)= .100185677 ; Cor(3)= .246294246
  Cor(4)= .433463493
  Cor(5)= .632350988 ; Cor(6)= .811118626 ; Cor(7)= .940848166
  Wi(1) = .196169389 ; Wi(2) = .270302644 ; Wi(3) = .239681873
  Wi(4) = .165775774
  Wi(5) = .0889432271 ; Wi(6) =.0331943043 ; Wi(7) = .00593278701
 CASE(8)
  Cor(1)= .0133202441 ; Cor(2)=.0797504290 ; Cor(3)= .197871029 ; Cor(4)= .354153994
  Cor(5)= .529458575 ; Cor(6)= .701814529 ; Cor(7)= .849379320 ; Cor(8)= .953326450
  Wi(1) = .164416604 ; Wi(2) = .237525610 ; Wi(3) = .226841984 ;Wi(4) = .175754079
  Wi(5) = .112924030 ; Wi(6) =.0578722107 ; Wi(7) =.0209790737 ;Wi(8) =.00368640710
 CASE DEFAULT
 CALL Error_Message('Gauss points not in range 1-8')
END SELECT
END SUBROUTINE Gauss_Laguerre_coor

INTEGER FUNCTION Ngaus(RonL,ne,Rlim)
!-----------------------------------------------------
!   Function returns number of Gauss points needed 
!   to integrate a function 1/rn
!------------------------------------------------------
IMPLICIT NONE
REAL , INTENT(IN)    :: RonL   !  R/L
INTEGER , INTENT(IN) :: ne    !  exponent (1,2,3)
REAL , INTENT(OUT)   :: Rlim(2)  !  array to store values of table
INTEGER         :: n
SELECT CASE(ne)
 CASE(1)
  Rlim= (/1.4025, 0.6736/)
 CASE(2)
  Rlim= (/2.3187, 0.9709/)
 CASE(3)
  Rlim= (/3.4170, 1.2908/)
 CASE DEFAULT
END SELECT
Ngaus=0
DO  N=1,2
 IF(RonL >= Rlim(N)) THEN
  Ngaus= N+2
  EXIT
 END IF
END DO 
IF (Ngaus == 0)THEN      ! Point is to near the surface
 Ngaus=5
END IF
RETURN
END FUNCTION Ngaus

SUBROUTINE Integ2P (Elcor,Inci,Nodel,Ncol,xP,k,dUe,dTe,Ndest,Isym)
!--------------------------------------------------
!   Computes  [dT]e and [dU]e for 2-D potential problems
!   by numerical integration
!-------------------------------------------------
IMPLICIT NONE
REAL, INTENT(IN)  :: Elcor(:,:)     ! Element coordinates
INTEGER, INTENT(IN) :: Ndest(:,:)     !   Node destination vector
INTEGER, INTENT(IN) :: Inci(:)      ! Element Incidences
INTEGER, INTENT(IN) :: Nodel       ! No. of Element Nodes
INTEGER , INTENT(IN):: Ncol        ! Number of points Pi (coll. points)
INTEGER , INTENT(IN):: Isym  
REAL , INTENT(IN)  :: xP(:,:)      ! Array with coll. points coords.
REAL , INTENT(IN)  :: k         ! Permeability
REAL(KIND=8) , INTENT(OUT) :: dUe(:,:),dTe(:,:) ! arrays
REAL :: epsi= 1.0E-10           !   Small value for comparing coords
REAL :: Eleng,Rmin,RonL,Glcor(8),Wi(8),Ni(Nodel),Vnorm(2),GCcor(2)
REAL :: UP,Jac,dxr(2),TP,r,pi,c1,c2,xsi,eta,dxdxb,Rlim(2),Xsi1,Xsi2,RJacB
INTEGER :: i,m,n,Mi,nr,ldim,cdim,nreg,id,NDIV,NDIVS,MAXDIVS
pi=3.14159265359
ldim= 1
cdim=ldim+1
CALL Elength(Eleng,Elcor,nodel,ldim)   ! Element Length
dUe= 0.0 ; dTe= 0.0            ! Clear arrays for summation
!-----------------------------------------------------------------
!    Integration off-diagonal coeff.  -> normal Gauss Quadrature
!-----------------------------------------------------------------
dUe= 0.0 ; dTe= 0.0                        ! Clear arrays for summation
MAXDIVS=1
Colloc_points: DO i=1,Ncol
         Rmin= Min_dist1(Elcor,xP(:,i),Nodel,inci,ELeng,Eleng,ldim)  ! Distance coll. point and element
         RonL= Rmin/Eleng                ! R/L
         Mi= Ngaus(RonL,1,Rlim)                ! Number of Gauss points for (1/r) singularity
         NDIVS= 1
         RJacB=1.0
         IF(Mi == 5) THEN
          IF(RonL > epsi) NDIVS= INT(RLim(2)/RonL) + 1   ! Number of subdvisions
          IF(NDIVS > MAXDIVS) MAXDIVS= NDIVS
          RJacB= 1.0/NDIVS
          Mi=4
         END IF
!         write(2,*) i,Rmin,RonL,NDIVS,RJacB
         Call Gauss_coor(Glcor,Wi,Mi)          ! Assign coords/Weights
         Xsi1=-1
  Subdivisions: DO NDIV=1,NDIVS
          Xsi2= Xsi1 + 2.0/NDIVS
   Gauss_points: DO m=1,Mi
           xsi= Glcor(m)
!           write(2,*) NDIV,m,xsi1,xsi2
!           write(2,*) xsi
           IF(NDIVS > 1) xsi= 0.5*(Xsi1+Xsi2)+xsi/NDIVS
!           write(2,*) xsi
           CALL Serendip_func(Ni,xsi,eta,ldim,nodel,Inci)      ! Shape function value
           Call Normal_Jac(Vnorm,Jac,xsi,eta,ldim,nodel,Inci,elcor) ! Jacobian and normal
           CALL Cartesian(GCcor,Ni,ldim,elcor)            ! Cart. coords of Gauss pt
           r= Dist(GCcor,xP(:,i),cdim)                ! Dist. P,Q
           dxr= (GCcor-xP(:,i))/r                  ! rx/r , ry/r
           UP= U(r,k,cdim) ; TP= T(r,dxr,Vnorm,cdim)         ! Kernels
     Node_points: DO n=1,Nodel
             IF(Isym == 0)THEN
              iD= i
             ELSE
              iD= Ndest(i,1)              !  line number in array
             END IF
             IF (id == 0) CYCLE
             IF(Dist(Elcor(:,n),xP(:,i),cdim) > epsi) THEN    ! Only case where coords of n and Pi not same
              dUe(id,n)= dUe(id,n) + Ni(n)*UP*Jac*Wi(m)*RJacB
              dTe(id,n)= dTe(id,n) + Ni(n)*TP*Jac*Wi(m)*RJacB
             END IF
            END DO Node_points
          END DO Gauss_points
          Xsi1= Xsi2
         END DO Subdivisions
        END DO Colloc_points
!------------------------------
!    Diagonal terms of dUe
!------------------------------
!write(2,*) 'Max. subdivisions=',MAXDIVS
c1= 1/(2.0*pi*k)
Colloc_points1: DO i=1,Ncol
  Node_points1: DO n=1,Nodel
          IF(Isym == 0)THEN
           iD= i
          ELSE
           iD= Ndest(i,1)              !  line number in array
          END IF
          IF (id == 0) CYCLE
          IF(Dist(Elcor(:,n),xP(:,i),cdim) > Epsi) CYCLE ! only do this when Pi is node n
          Nreg=1
          IF(n == 3) nreg= 2
    Subregions: DO nr=1,Nreg
           Mi= 4
           Call Gauss_Laguerre_coor(Glcor,Wi,Mi)
   Gauss_points1: DO m=1,Mi
            SELECT CASE (n)
             CASE (1)
              xsi= 2.0*Glcor(m)-1.0
              dxdxb= 2.0
             CASE (2)
              xsi= 1.0 -2.0*Glcor(m)
              dxdxb= 2.0
             CASE (3)
              dxdxb= 1.0
               IF(nr == 1) THEN
                xsi= -Glcor(m)
                ELSE
                xsi= Glcor(m)
               END IF
             CASE DEFAULT
            END SELECT
            CALL Serendip_func(Ni,xsi,eta,ldim,nodel,Inci)      ! Shape function value
            Call Normal_Jac(Vnorm,Jac,xsi,eta,ldim,nodel,Inci,elcor) ! Jacobian
            dUe(id,n)= dUe(id,n) + Ni(n)*c1*Jac*dxdxb*Wi(m)
           END DO Gauss_points1
          END DO Subregions
          Mi= 2
          Call Gauss_coor(Glcor,Wi,Mi)          ! Assign coords/Weights
  Gauss_points2: DO m=1,Mi
           SELECT CASE (n)
            CASE (1)
             c2=-LOG(Eleng)*c1
            CASE (2)
             c2=-LOG(Eleng)*c1
            CASE (3)
             c2=LOG(2/Eleng)*c1
            CASE DEFAULT
           END SELECT
           xsi= Glcor(m)
           CALL Serendip_func(Ni,xsi,eta,ldim,nodel,Inci)      ! Shape function value
           Call Normal_Jac(Vnorm,Jac,xsi,eta,ldim,nodel,Inci,elcor) ! Jacobian and normal
           dUe(id,n)= dUe(id,n) + Ni(n)*c2*Jac*Wi(m)
          END DO Gauss_points2
         END DO Node_points1
        END DO Colloc_points1
        RETURN
END SUBROUTINE Integ2P

SUBROUTINE Integ2E(Elcor,Inci,Nodel,Ncol,xP,E,ny,dUe,dTe,Ndest,Isym)
!--------------------------------------------------
!    Computes  [dT]e and [dU]e for 2-D elasticity problems
!    by numerical integration
!-------------------------------------------------
IMPLICIT NONE
REAL, INTENT(IN)    :: Elcor(:,:)     !   Element coordinates
INTEGER, INTENT(IN) :: Ndest(:,:)     !   Node destination vector
INTEGER, INTENT(IN) :: Inci(:)        !   Element Incidences
INTEGER, INTENT(IN) :: Nodel          !   No. of Element Nodes
INTEGER , INTENT(IN):: Ncol           !   Number of points Pi (coll. points)
INTEGER , INTENT(IN):: Isym  
REAL, INTENT(IN)  :: E,ny           !   Elastic constants
REAL, INTENT(IN)    :: xP(:,:)        !   Array with coll. points coords.
REAL(KIND=8), INTENT(OUT)   :: dUe(:,:),dTe(:,:) !  arrays for storing element coefficients
REAL        :: epsi= 1.0E-10  !    Small value for comparing coords
REAL        :: Eleng,Rmin,RonL,Glcor(8),Wi(8),Ni(Nodel),Vnorm(2),GCcor(2)
REAL        :: Jac,dxr(2),UP(2,2),TP(2,2), xsi, eta, r, dxdxb,Pi,C,C1,Rlim(2),Xsi1,Xsi2,RJacB
INTEGER       :: i,j,k,m,n,Mi,nr,ldim,cdim,iD,nD,Nreg,NDIV,NDIVS,MAXDIVS
Pi=3.14159265359
C=(1.0+ny)/(4*Pi*E*(1.0-ny))   
ldim= 1                             ! Element dimension
cdim=ldim+1
MAXDIVS=1
CALL Elength(Eleng,Elcor,nodel,ldim)  ! Element Length
dUe= 0.0 ; dTe= 0.0                 ! Clear arrays for summation
 Colloc_points: DO i=1,Ncol
         Rmin= Min_dist1(Elcor,xP(:,i),Nodel,inci,ELeng,Eleng,ldim) !  Distance coll. point and element
         RonL= Rmin/Eleng                   !  R/L
 !     Integration off-diagonal coeff.  -> normal Gauss Quadrature
         Mi= Ngaus(RonL,1,Rlim)            !  Number of Gauss points for (1/r) singularity
         NDIVS= 1
         RJacB=1.0
         IF(Mi == 5) THEN
          IF(RonL > epsi) NDIVS= INT(RLim(2)/RonL) + 1   ! Number of subdvisions
          IF(NDIVS > MAXDIVS) MAXDIVS= NDIVS
          RJacB= 1.0/NDIVS
          Mi=4
         END IF
!         write(2,*) i,Rmin,RonL,NDIVS
         Call Gauss_coor(Glcor,Wi,Mi)  ! Assign coords/Weights
          Xsi1=-1
Subdivisions: DO NDIV=1,NDIVS
          Xsi2= Xsi1 + 2.0/NDIVS
  Gauss_points: DO m=1,Mi
          xsi= Glcor(m)
          IF(NDIVS > 1) xsi= 0.5*(Xsi1+Xsi2)+xsi/NDIVS
          CALL Serendip_func(Ni,xsi,eta,ldim,nodel,Inci)  !   Shape function value
          Call Normal_Jac(Vnorm,Jac,xsi,eta,ldim,nodel,Inci,elcor) ! Jacobian and normal
          CALL Cartesian(GCcor,Ni,ldim,elcor)          ! Cart. coords of Gauss pt
          r= Dist(GCcor,xP(:,i),cdim)            !  Dist. P,Q
          dxr= (GCcor-xP(:,i))/r         !  rx/r , ry/r
          UP= UK(dxr,r,E,ny,Cdim) ; TP= TK(dxr,r,Vnorm,ny,Cdim) !  Kernels
   Node_points: DO n=1,Nodel
    Direction_P: DO j=1,2
            IF(Isym == 0)THEN
             iD= 2*(i-1) + j
            ELSE
             iD= Ndest(i,j)              !  line number in array
            END IF
            IF (id == 0) CYCLE
     Direction_Q: DO k= 1,2
             nD= 2*(n-1) + k             !  column number in array
             IF(Dist(Elcor(:,n),xP(:,i),cdim) > epsi) THEN  ! n and Pi not same
              dUe(iD,nD)= dUe(iD,nD) + Ni(n)*UP(j,k)*Jac*Wi(m)*RJacB
              dTe(iD,nD)= dTe(iD,nD) + Ni(n)*TP(j,k)*Jac*Wi(m)*RJacB
             ELSE
              dUe(iD,nD)= dUe(iD,nD) + Ni(n)*C*dxr(j)*dxr(k)*Jac*Wi(m)*RJacB  !  non-log part of U
             END IF
            END DO Direction_Q
           END DO Direction_P
          END DO Node_points
         END DO Gauss_points
         Xsi1= Xsi2
        END DO Subdivisions
       END DO Colloc_points
!      Diagonal terms of dUe
 C= C*(3.0-4.0*ny)        
 Colloc_points1: DO i=1,Ncol
  Node_points1:  DO n=1,Nodel
           IF(Dist(Elcor(:,n),xP(:,i),cdim) > Epsi) CYCLE ! only do when Pi is node n
            Nreg=1
           IF (n == 3) nreg= 2
     Subregions: DO nr=1,Nreg
                Mi= 4
            Call Gauss_Laguerre_coor(Glcor,Wi,Mi)
   Gauss_points1: DO m=1,Mi
            SELECT CASE (n)
             CASE (1)
              xsi= 2.0*Glcor(m)-1.0
              dxdxb= 2.0
             CASE (2)
              xsi= 1.0 -2.0*Glcor(m)
              dxdxb= 2.0
             CASE (3)
              dxdxb= 1.0
              IF(nr == 1) THEN
               xsi= -Glcor(m)
              ELSE
               xsi= Glcor(m)
              END IF
             CASE DEFAULT
            END SELECT
            CALL Serendip_func(Ni,xsi,eta,ldim,nodel,Inci)  !   Shape function value
            Call Normal_Jac(Vnorm,Jac,xsi,eta,ldim,nodel,Inci,elcor) ! Jacobian
      Direction1: DO j=1,2
             IF(Isym == 0)THEN
              iD= 2*(i-1) + j
             ELSE
              iD= Ndest(i,j)              !  line number in array
             END IF
             IF (id == 0) CYCLE           
             nD= 2*(n-1) + j              !  column number in array
             dUe(iD,nD)= dUe(iD,nD) + Ni(n)*C*Jac*dxdxb*Wi(m)
            END DO Direction1
           END DO Gauss_points1
          END DO Subregions
          Mi= 2
          Call Gauss_coor(Glcor,Wi,Mi)          ! Assign coords/Weights
  Gauss_points2: DO m=1,Mi
           SELECT CASE (n)
            CASE (1)
             C1=-LOG(Eleng)*C
            CASE (2)
             C1=-LOG(Eleng)*C
            CASE (3)
             C1=LOG(2/Eleng)*C
            CASE DEFAULT
           END SELECT
           xsi= Glcor(m)
           CALL Serendip_func(Ni,xsi,eta,ldim,nodel,Inci)      ! Shape function value
           Call Normal_Jac(Vnorm,Jac,xsi,eta,ldim,nodel,Inci,elcor) ! Jacobian and normal
     Direction2: DO j=1,2
            IF(Isym == 0)THEN
             iD= 2*(i-1) + j
            ELSE
             iD= Ndest(i,j)              !  line number in array
            END IF
            IF (id == 0) CYCLE           
            nD= 2*(n-1) + j              !  column number in array
            dUe(iD,nD)= dUe(iD,nD) + Ni(n)*C1*Jac*Wi(m)
           END DO Direction2
          END DO Gauss_points2
         END DO Node_points1
        END DO Colloc_points1
        RETURN
        END SUBROUTINE Integ2E

SUBROUTINE Integ3(Elcor,Inci,Nodel,Ncol,xPi,Ndof,E,ny,ko,dUe,dTe,Ndest,Isym)
!--------------------------------------------------
!    Computes  [dT]e and [dU]e for 3-D problems
!    by numerical integration
!-------------------------------------------------
IMPLICIT NONE
REAL, INTENT(IN)    :: Elcor(:,:)     !   Element coordinates
INTEGER, INTENT(IN) :: Ndest(:,:)     !   Node destination vector
INTEGER, INTENT(IN) :: Inci(:)        !   Element Incidences
INTEGER, INTENT(IN) :: Nodel          !   No. of Element Nodes
INTEGER , INTENT(IN):: Ncol           !   Number of points Pi (coll. points)
REAL , INTENT(IN)   :: xPi(:,:)       !   Array with coll. points coords.
INTEGER , INTENT(IN):: Ndof           !   Number dgrees of freedom /node (1 or 3)
INTEGER , INTENT(IN):: Isym          
REAL , INTENT(IN)   :: E,ny           !   Elastic constants (for elasticity problems)
REAL , INTENT(IN)   :: ko            
REAL(KIND=8) , INTENT(OUT)  :: dUe(:,:),dTe(:,:) !  arrays for storing coefficients
REAL :: Elengx,Elenge,Rmin,RLx,RLe,Glcorx(8),Wix(8),Glcore(8),Wie(8),Weit,r
REAL :: Ni(Nodel),Vnorm(3),GCcor(3),dxr(3),Jac,Jacb,xsi,eta,xsib,etab,Rlim(2)
REAL :: Xsi1,Xsi2,Eta1,Eta2,RJacB,RonL
REAL :: UP(Ndof,Ndof),TP(Ndof,Ndof)   !   Arrays for storing kernels
INTEGER :: i,m,n,k,ii,jj,ntr,Mi,Ki,id,nd,lnod,Ntri,NDIVX,NDIVSX,NDIVE,NDIVSE,MAXDIVS
INTEGER :: ldim= 2       !   Element dimension
INTEGER :: Cdim= 3       !   Cartesian dimension
ELengx= Dist((Elcor(:,3)+Elcor(:,2))/2.,(Elcor(:,4)+Elcor(:,1))/2.,Cdim)  ! Lxsi
ELenge= Dist((Elcor(:,2)+Elcor(:,1))/2.,(Elcor(:,3)+Elcor(:,4))/2.,Cdim)  ! Leta
dUe= 0.0 ; dTe= 0.0                 ! Clear arrays for summation
!---------------------------------------------------------------
!     Part 1 : Pi is not one of the element nodes
!---------------------------------------------------------------
Colloc_points: DO i=1,Ncol
         IF(.NOT. ALL(Inci /= i)) CYCLE     !  Check if incidence array contains i
         Rmin= Min_dist1(Elcor,xPi(:,i),Nodel,inci,ELengx,Elenge,ldim) !  Distance coll. point and element
         Mi= Ngaus(Rmin/Elengx,2,Rlim)           !  Number of G.P. in xsi dir. for (1/r2) sing.
         RonL= Rmin/Elengx 
         NDIVSX= 1 ; NDIVSE= 1
         RJacB=1.0
         IF(Mi == 5) THEN
          IF(RonL > 0.0) NDIVSX= INT(RLim(2)/RonL) + 1   ! Number of subdvisions
          IF(NDIVSX > MAXDIVS) MAXDIVS= NDIVSX
          Mi=4
         END IF
!         write(2,*) i,Rmin,RonL,NDIVSX
         Call Gauss_coor(Glcorx,Wix,Mi)     !  Assign coords/Weights xsi-direction
         Ki= Ngaus(Rmin/Elenge,2,Rlim)           !  Number of G.P. in eta dir. for (1/r2) sing.
         RonL= Rmin/Elenge 
         IF(Ki == 5) THEN
          IF(RonL > 0.0) NDIVSE= INT(RLim(2)/RonL) + 1   ! Number of subdvisions
          IF(NDIVSE > MAXDIVS) MAXDIVS= NDIVSE
          Ki=4
         END IF
         IF(NDIVSX > 1 .OR. NDIVSE>1) RJacB= 1.0/(NDIVSX*NDIVSE)
!         write(2,*) i,RonL,NDIVSE
         Call Gauss_coor(Glcore,Wie,Ki)     !  Assign coords/Weights eta-direction
          Xsi1=-1.0
Subdivisions_xsi: DO NDIVX=1,NDIVSX
          Xsi2= Xsi1 + 2.0/NDIVSX
          Eta1=-1.0
Subdivisions_eta: DO NDIVE=1,NDIVSE
          Eta2= Eta1 + 2.0/NDIVSE
Gauss_points_xsi: DO m=1,Mi
          xsi= Glcorx(m)
          IF(NDIVSX > 1) xsi= 0.5*(Xsi1+Xsi2)+xsi/NDIVSX
  Gauss_points_eta: DO k=1,Ki
           eta= Glcore(k)
           IF(NDIVSE > 1) eta= 0.5*(Eta1+Eta2)+eta/NDIVSE
           Weit= Wix(m)*Wie(k)*RJacB
           CALL Serendip_func(Ni,xsi,eta,ldim,nodel,Inci)  !   Shape function value
           Call Normal_Jac(Vnorm,Jac,xsi,eta,ldim,nodel,Inci,elcor) ! Jacobian and normal
           CALL Cartesian(GCcor,Ni,ldim,elcor)          ! Cart. coords of Gauss pt
           r= Dist(GCcor,xPi(:,i),Cdim)             !  Dist. P,Q
           dxr= (GCcor-xPi(:,i))/r             !  rx/r , ry/r
            IF(Ndof .EQ. 1) THEN
             UP= U(r,ko,Cdim) ; TP= T(r,dxr,Vnorm,Cdim) !  Kernels Potential problem
            ELSE
             UP= UK(dxr,r,E,ny,Cdim) ; TP= TK(dxr,r,Vnorm,ny,Cdim) !  Kernels elasticity
            END IF
!      Node_points: DO n=1,Nodes
       Direction_P: DO ii=1,Ndof
               IF(Isym == 0)THEN
                iD= Ndof*(i-1) + ii     !  line number in array
               ELSE
                iD= Ndest(i,ii)              !  line number in array
               END IF
               IF (id == 0) CYCLE
        Direction_Q: DO jj=1,Ndof
         Node_points: DO n=1,Nodel
                nD= Ndof*(n-1) + jj             !  column number in array
                dUe(iD,nD)= dUe(iD,nD) + Ni(n)*UP(ii,jj)*Jac*Weit
                dTe(iD,nD)= dTe(iD,nD) + Ni(n)*TP(ii,jj)*Jac*Weit
                END DO Node_points
               END DO Direction_Q
              END DO Direction_P
!             END DO Node_points
          END DO Gauss_points_eta
         END DO Gauss_points_xsi
         Eta1= Eta2
        END DO Subdivisions_eta
        Xsi1= Xsi2
       END DO Subdivisions_xsi
      END DO Colloc_points
!---------------------------------------------------------
!     Part 1 : Pi is one of the element nodes
!---------------------------------------------------------
Colloc_points1: DO i=1,Ncol
         lnod= 0
         DO n= 1,Nodel                      !   Determine which local node is Pi
          IF(Inci(n) .EQ. i) THEN
           lnod=n
          END IF
         END DO
          IF(lnod .EQ. 0) CYCLE          !  None -> next Pi
           Ntri= 2
           IF(lnod > 4) Ntri=3          !  Number of triangles
      Triangles: DO ntr=1,Ntri
            CALL Tri_RL(RLx,RLe,Elengx,Elenge,lnod,ntr)
            Mi= Ngaus(RLx,2,Rlim)          !  Number of G.P. in xsi dir. for (1/r2) sing.
! write(2,*) i, Mi
            IF(Mi == 5) Mi=4
!            Mi= 8
            Call Gauss_coor(Glcorx,Wix,Mi)   !  Assign coords/Weights xsi-direction
            Ki= Ngaus(RLe,2,Rlim)          !  Number of G.P. in eta dir. for (1/r2) sing.
! write(2,*) i, Mi
            IF(Ki == 5) Ki=4
!            Ki= 8
            Call Gauss_coor(Glcore,Wie,Ki)   !  Assign coords/Weights eta-direction
   Gauss_points_xsi1: DO m=1,Mi
              xsib= Glcorx(m)
    Gauss_points_eta1: DO k=1,Ki
               etab= Glcore(k)
               Weit= Wix(m)*Wie(k)
               CALL Trans_Tri(ntr,lnod,xsib,etab,xsi,eta,Jacb) !  Coord transf from triang coords
               CALL Serendip_func(Ni,xsi,eta,ldim,nodel,Inci)  !   Shape function value
               Call Normal_Jac(Vnorm,Jac,xsi,eta,ldim,nodel,Inci,elcor) ! Jacobian and normal
               Jac= Jac*Jacb
               CALL Cartesian(GCcor,Ni,ldim,elcor)             ! Cart. coords of Gauss pt
               r= Dist(GCcor,xPi(:,i),Cdim)                          !  Dist. P,Q
               dxr= (GCcor-xPi(:,i))/r                          !  rx/r , ry/r
               IF(Ndof .EQ. 1) THEN
                UP= U(r,ko,Cdim) ; TP= T(r,dxr,Vnorm,Cdim) !  Kernels Potential problem
               ELSE
                UP= UK(dxr,r,E,ny,Cdim) ; TP= TK(dxr,r,Vnorm,ny,Cdim) !  Kernels elasticity
               END IF
         Direction_P1: DO ii=1,Ndof
                 IF(Isym == 0)THEN
                  iD= Ndof*(i-1) + ii     !  line number in array
                 ELSE
                  iD= Ndest(i,ii)             !  line number in array
                 END IF                 
                 IF (id == 0) CYCLE
         Direction_Q1:  DO jj=1,Ndof
           Node_points1: DO n=1,Nodel
                    nD= Ndof*(n-1) + jj    !  column number in array
                   dUe(iD,nD)= dUe(iD,nD) + Ni(n)*UP(ii,jj)*Jac*Weit
                   IF(Inci(n) /= i) THEN  !   diagonal elements of dTe not computed
                    dTe(iD,nD)= dTe(iD,nD) + Ni(n)*TP(ii,jj)*Jac*Weit
                   END IF
                  END DO Node_points1
                 END DO Direction_Q1
                END DO Direction_P1
              END DO Gauss_points_eta1
             END DO Gauss_points_xsi1
            END DO Triangles
        END DO Colloc_points1
        RETURN
END SUBROUTINE Integ3

SUBROUTINE Triangel_Coord(cor_tri,lnod,ntr)
!---------------------------------------------
!   Assigns local coordinates of triangular
!   subelements
!---------------------------------------------
IMPLICIT NONE
INTEGER, INTENT(IN) :: lnod,ntr      !  node, subelement no. 
REAL , INTENT(OUT)  :: cor_tri(2,3)  !  xsi,eta of triangle nodes
REAL        :: xsii(8),etai(8)
INTEGER       :: nod_tri(3),n
SELECT CASE (ntr)
 CASE (1)
 SELECT CASE(lnod)
  CASE (1)
   nod_tri=(/2,3,1/)
  CASE (2)
   nod_tri=(/3,4,2/)
  CASE (3)
   nod_tri=(/1,2,3/)
  CASE (4)
   nod_tri=(/1,2,4/)
  CASE (5)
   nod_tri=(/4,1,5/)
  CASE (6)
   nod_tri=(/1,2,6/)
  CASE (7)
   nod_tri=(/4,1,7/)
  CASE (8)
   nod_tri=(/1,2,8/)
  CASE DEFAULT 
 END SELECT
 CASE (2)
 SELECT CASE(lnod)
  CASE (1)
   nod_tri=(/3,4,1/)
  CASE (2)
   nod_tri=(/4,1,2/)
  CASE (3)
   nod_tri=(/4,1,3/)
  CASE (4)
   nod_tri=(/2,3,4/)
  CASE (5)
   nod_tri=(/2,3,5/)
  CASE (6)
   nod_tri=(/3,4,6/)
  CASE (7)
   nod_tri=(/2,3,7/)
  CASE (8)
   nod_tri=(/3,4,8/)
  CASE DEFAULT
 END SELECT  
 CASE (3)
 SELECT CASE(lnod)
  CASE (5)
   nod_tri=(/3,4,5/)
  CASE (6)
   nod_tri=(/4,1,6/)
  CASE (7)
   nod_tri=(/1,2,7/)
  CASE (8)
   nod_tri=(/2,3,8/)
  CASE DEFAULT
 END SELECT  
 CASE DEFAULT
END SELECT
xsii=(/-1.0,1.0,1.0,-1.0,0.0,1.0,0.0,-1.0/)
etai=(/-1.0,-1.0,1.0,1.0,-1.0,0.0,1.0,0.0/)
cor_tri=0
DO n=1, 3
 cor_tri(1,n)= xsii(nod_tri(n))
 cor_tri(2,n)= etai(nod_tri(n))
END DO
END SUBROUTINE Triangel_Coord

SUBROUTINE Trans_Tri(ntr,lnod,xsib,etab,xsi,eta,Jacb)
!--------------------------------------------
!  Transforms from local triangle coordinates
!  to xsi,eta coordinates and computes 
!  the Jacobean of the transformation
!---------------------------------------------
IMPLICIT NONE
REAL, INTENT (IN)    :: xsib,etab     !  local coordinates
INTEGER, INTENT (IN) :: ntr,lnod      !   subelement-no, local node no.
REAL, INTENT (OUT)   :: Jacb,xsi,eta  !   Jacobean and xsi,eta coords
REAL         :: Nb(3),dNbdxb(3),dNbdeb(3),dxdxb,dxdeb,dedxb,dedeb,cor_tri(2,3)
INTEGER      :: n,i
CALL Triangel_Coord(cor_tri,lnod,ntr)
!
!    Transform xsi-bar and eta-bar to xsi and eta
!
Nb(1)=0.25*(1.0+xsib)*(1.0-etab)
Nb(2)=0.25*(1.0+xsib)*(1.0+etab)
Nb(3)=0.5*(1.0-xsib)
xsi = 0.0
eta = 0.0
DO n=1,3
 xsi = xsi+Nb(n)*cor_tri(1,n) 
 eta = eta+Nb(n)*cor_tri(2,n) 
END DO
!
!    Jacobian of Transformation xsi-bar and eta-bar to xsi and eta
!
dNbdxb(1)=0.25*(1.0-etab)
dNbdxb(2)=0.25*(1.0+etab)
dNbdxb(3)=-0.5
dNbdeb(1)=-0.25*(1.0+xsib)
dNbdeb(2)=0.25*(1.0+xsib)
dNbdeb(3)=0.0
dxdxb=0.0
dedxb=0.0
dxdeb=0.0
dedeb=0.0
DO i=1,3
 dxdxb=dxdxb+dNbdxb(i)*cor_tri(1,i)
 dedxb=dedxb+dNbdxb(i)*cor_tri(2,i)
 dedeb=dedeb+dNbdeb(i)*cor_tri(2,i)
 dxdeb=dxdeb+dNbdeb(i)*cor_tri(1,i)
END DO
Jacb=dxdxb*dedeb-dedxb*dxdeb
END SUBROUTINE Trans_Tri

SUBROUTINE Tri_RL(RLx,RLe,Elengx,Elenge,lnod,ntr)
!---------------------------------------  
!  Computes ize of triangular sub-element
!  and Rmin
!---------------------------------------
IMPLICIT NONE
REAL, INTENT (IN)     :: Elengx,Elenge  !  length of element in xsi,eta dir
INTEGER, INTENT (IN)  :: lnod,ntr       !  local node, subelement no.
REAL, INTENT (OUT)    :: RLx,RLe        !  Lenghts of subelement
IF(lnod <= 4) THEN
 RLx=Elengx/Elenge
 RLe=Elenge/Elengx
END IF
IF (lnod == 5 .or. lnod == 7) THEN
 SELECT CASE (ntr)
  CASE (1)
   RLx=(Elengx/2.0)/Elenge
   RLe=Elenge/(Elengx/2.0)
  CASE (2)
   RLx=(Elengx/2.0)/Elenge
   RLe=Elenge/(Elengx/2.0)
  CASE (3)
   RLx=Elengx/Elenge
   RLe=Elenge/Elengx
 END SELECT
END IF
IF (lnod == 6 .or. lnod == 8) THEN
 SELECT CASE (ntr)
  CASE (1)
   RLx=Elengx/(Elenge/2)
   RLe=(Elenge/2)/Elengx
  CASE (2)
   RLx=Elengx/(Elenge/2)
   RLe=(Elenge/2)/Elengx
  CASE (3)
   RLx=Elengx/Elenge
   RLe=Elenge/Elengx
 END SELECT
END IF
END SUBROUTINE Tri_RL
END MODULE Integration_lib
