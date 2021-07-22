PROGRAM General_purpose_BEM
!------------------------------------------------------
!     General purpose BEM program
!     for solving elasticity and potential problems
!------------------------------------------------------
!USE DFLIB;
USE Utility_lib ; USE Elast_lib ; USE Laplace_lib; USE Integration_lib; USE BodyForce_lib


IMPLICIT NONE
INTEGER(4)                                       :: response
INTEGER, ALLOCATABLE :: Inci(:,:)  !  Element Incidences
INTEGER, ALLOCATABLE :: BCode(:,:), NCode(:) !  Element BC�s
INTEGER, ALLOCATABLE :: Ldest(:,:) !  Element destination vector
INTEGER, ALLOCATABLE :: Ndest(:,:) !  Node destination vector
REAL, ALLOCATABLE :: Elres_u(:,:)  !  Element results , u
REAL, ALLOCATABLE :: Elres_t(:,:)  !  Element results , t
REAL, ALLOCATABLE :: Elcor(:,:)    !  Element coordinates
REAL, ALLOCATABLE :: xP(:,:)       !  Node co-ordinates
REAL(KIND=8), ALLOCATABLE :: dUe(:,:),dTe(:,:),Diag(:,:)
REAL(KIND=8), ALLOCATABLE :: Lhs(:,:),F(:)
REAL(KIND=8), ALLOCATABLE :: u1(:) !  global vector of unknown
CHARACTER (LEN=80) :: Title
INTEGER :: Cdim,Node,m,n,Istat,Nodel,Nel,Ndof,Toa
INTEGER :: Nreg,Ltyp,Nodes,Maxe,Ndofe,Ndofs,ndg,NE_u,NE_t               
INTEGER :: nod,nd,i,j,k,l,DoF,Pos,Isym,nsym,nsy
REAL,ALLOCATABLE    :: Fac(:)     !  Factors for symmetry
REAL,ALLOCATABLE    :: Elres_te(:),Elres_ue(:)   
INTEGER,ALLOCATABLE :: Incie(:)   !  Incidences for one element
INTEGER,ALLOCATABLE :: Ldeste(:)  !  Destination vector 1 elem
REAL :: Con,E,ny,Scat,Scad
!-----------------------------------------------------
!   Read job information
!-----------------------------------------------------
OPEN (UNIT=1,FILE='INPUT',FORM='FORMATTED') !  Input
OPEN (UNIT=2,FILE='OUTPUT',FORM='FORMATTED')!  Output
Call Jobin(Title,Cdim,Ndof,Toa,Nreg,Ltyp,Con,E,ny,&
           Isym,nodel,nodes,maxe)
Nsym= 2**Isym   !   number of symmetry loops
ALLOCATE(xP(Cdim,Nodes))   !  Array for node coordinates
ALLOCATE(Inci(Maxe,Nodel)) !  Array for incidences
CALL Geomin(Nodes,Maxe,xp,Inci,Nodel,Cdim)
Ndofe= Nodel*Ndof   !    Total degrees of freedom of element
ALLOCATE(BCode(Maxe,Ndofe))      !    Element Boundary codes
ALLOCATE(Elres_u(Maxe,Ndofe),Elres_t(Maxe,Ndofe))       
CALL BCinput(Elres_u,Elres_t,Bcode,nodel,ndofe,ndof) 
ALLOCATE(Ldest(maxe,Ndofe))  ! Elem. destination vector
ALLOCATE(Ndest(Nodes,Ndof))

!---------------------------------------------------------------------
!     Determine Node destination vector and Element destination vector 
!---------------------------------------------------------------------
CALL Destination(Isym,Ndest,Ldest,xP,Inci,Ndofs,nodes,Ndof,Nodel,Maxe)
!------------------------------------------
!     Determine global Boundary code vector
!---------------------------------------------
ALLOCATE(NCode(Ndofs))            

NCode=0
DoF_o_System: &
DO      nd=1,Ndofs
        DO Nel=1,Maxe
                DO m=1,Ndofe
                        IF (nd == Ldest(Nel,m) .and. NCode(nd) == 0) THEN 
                                NCode(nd)= NCode(nd)+BCode(Nel,m)
                        END IF
                END DO
        END DO
END DO &
DoF_o_System
IF(Ndof ==1)E= Con
CALL Scal(E,xP(:,:),Elres_u(:,:),Elres_t(:,:),Cdim,Scad,Scat)
ALLOCATE(dTe(Ndofs,Ndofe),dUe(Ndofs,Ndofe))             ! Elem. coef. matrices
ALLOCATE(Diag(Ndofs,Ndof))                                                                              ! Diagonal coefficients
ALLOCATE(Lhs(Ndofs,Ndofs),F(Ndofs),u1(Ndofs)) ! global arrays
ALLOCATE(Elcor(Cdim,Nodel))                                                                             !  Elem. Coordinates
ALLOCATE(Fac(Ndofe))                                                                                                    !  Factor for symmetric elements
ALLOCATE(Incie(Nodel))                                                                                          !  Element incidences
ALLOCATE(Ldeste(Ndofe))                                                                                         !  Element destination
ALLOCATE(Elres_te(Ndofe),Elres_ue(Ndofe))                                                                                       ! Tractions of Element  

!----------------------------------------------------------------
!  Compute element coefficient matrices
!----------------------------------------------------------------
Lhs(:,:) = 0.0; F(:) = 0.0; u1(:) = 0.0; Diag(:,:) = 0.0
Elements_1:&
DO Nel=1,Maxe
! write(2,*) 'Element=',Nel
        Symmetry_loop:&
        DO nsy= 1,Nsym
                Elcor(:,:)= xP(:,Inci(Nel,:))  !        gather element coordinates
                Incie= Inci(nel,:)             !        incidences
                Ldeste= Ldest(nel,:)           !        and destinations
                Fac(1:nodel*ndof)= 1.0
                Elres_te(:)=Elres_t(Nel,:)
                IF(Isym > 0) THEN
                        CALL Mirror(Isym,nsy,Nodes,Elcor,Fac,Incie,Ldeste,Elres_te,Elres_ue & 
                                                                        ,nodel,ndof,Cdim) 
                END IF
                IF(Cdim == 2) THEN
                        IF(Ndof == 1) THEN
                                CALL Integ2P(Elcor,Incie,Nodel,Nodes,xP,Con,dUe,dTe,Ndest,Isym)
                        ELSE
                                CALL Integ2E(Elcor,Incie,Nodel,Nodes,xP,E,ny,dUe,dTe,Ndest,Isym)  
                        END IF
                ELSE
                        CALL Integ3(Elcor,Incie,Nodel,Nodes,xP,Ndof &
                                                                        ,E,ny,Con,dUe,dTe,Ndest,Isym)    
                END IF
                CALL Assembly(Lhs,F,DTe,DUe,Ldeste,BCode(Nel,:),Ncode &
                                                                        ,Elres_u(Nel,:),Elres_te,Diag,Ndofe,Ndof,Nodel,Fac)     
        END DO &
        Symmetry_loop
END DO &
Elements_1
!WRITE(2,*)'LHS'
!DO n=1, 16
! 	WRITE(2,'(16F10.6)')Lhs(n,:)
! END DO
! WRITE(2,*)'DIAG'
! WRITE(2,'(16F10.6)')Diag(:,1)

CALL Body_force(F,CDim,xP,Nodes,Isym,E,ny)   !  add any Bodyforce effect
!------------------------------------------------------------
!  Add azimuthal integral for infinite regions
!------------------------------------------------------------
IF(Nreg == 2) THEN
        DO m=1, Nodes
                DO n=1, Ndof
                        IF(Ndest(m,n) == 0)CYCLE
                        k=Ndest(m,n)
!                       k=Ndof*(m-1)+n
                        Diag(k,n) = Diag(k,n) + 1.0
                END DO
        END DO           
END IF
!-------------------------------------------------------------
!  Add Diagonal coefficients
!-------------------------------------------------------------
DO m=1,Ndofs            ! Loop over collocationpoints
        Nod=0
        DO n=1, Nodes
                DO l=1,Ndof
                        IF (m == Ndest(n,l))THEN
                                Nod=n
                                EXIT
                        END IF  
                END DO
                IF (Nod /= 0)EXIT
        END DO
        DO k=1,Ndof
                DoF=Ndest(Nod,k)
                IF(DoF /= 0) THEN
                        IF(NCode(DoF) == 1) THEN
                                Nel=0
                                Pos=0
                                DO i=1,Maxe
                                        DO j=1,Ndofe
                                                IF(DoF == Ldest(i,j))THEN
                                                        Nel=i
                                                        Pos=j
                                                        EXIT
                                                END IF
                                        END DO
                                        IF(Nel /= 0)EXIT
                                END DO
                                F(m) = F(m) - Diag(m,k) * Elres_u(Nel,Pos)
                        ELSE
                                Lhs(m,DoF)= Lhs(m,DoF) + Diag(m,k)
                        END IF
                END IF
        END DO
END DO
! WRITE(2,*)'LHS'
! DO n=1, 16
! 	WRITE(2,'(16F10.6)')Lhs(n,:)
! END DO

!---------------------------------------------------------
!   Solve system of equations
!---------------------------------------------------------
CALL Solve(Lhs,F,u1)
! CLOSE(UNIT=2)
OPEN (UNIT=3,FILE='BERESULTS',FORM='FORMATTED')
!   Gather Element results from global result vector u1

WRITE (2,*) ' '
Elements_2:     &
DO nel=1,maxe
        D_o_F1:         &
        DO nd=1,Ndofe
                IF(Ldest(nel,nd) /= 0)THEN
                        IF(NCode(Ldest(nel,nd)) == 0) THEN
                                Elres_u(nel,nd) = Elres_u(nel,nd) + u1(Ldest(nel,nd))
                        ELSE
                                Elres_t(nel,nd) = Elres_t(nel,nd) + u1(Ldest(nel,nd))
                        END IF
                END IF
        END DO &
        D_o_F1  
        Elres_u(nel,:)= Elres_u(nel,:) * Scad
        Elres_t(nel,:)= Elres_t(nel,:) / Scat
 !       WRITE (2,*) 'Results, Element',nel
 !       WRITE(2,'(A,6E10.3/3x,6E10.3/3x,6E10.3/3x,6E10.3)') ' u=',(Elres_u(nel,m), m=1,Ndofe)
 !       WRITE(2,'(A,6F10.3/3x,6F10.3/3x,6F10.3/3x,6F10.3)') ' t=',(Elres_t(nel,m), m=1,Ndofe)
        WRITE(3,'(24F12.5)') (Elres_u(nel,m), m=1,Ndofe)
        WRITE(3,'(24F12.5)') (Elres_t(nel,m), m=1,Ndofe)
END DO &
Elements_2
CLOSE(UNIT=2)
CLOSE(UNIT=3)
!response = MESSAGEBOXQQ('Program finished'C,'General Purpose Program'C, MB$ICONASTERISK)
END PROGRAM

