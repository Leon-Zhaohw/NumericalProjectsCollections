!     Last change:  CD   21 Mar 2000   11:15 pm
MODULE Utility_lib
!    Utility programs
CONTAINS
	SUBROUTINE Error_Message(TEXT)
  !---------------------------------------
  ! Writes an error message onto an error file
  ! and the console and terminates the program
  !--------------------------------------
   implicit none
   CHARACTER (LEN=*) TEXT
   LOGICAL :: EXST
   INQUIRE(FILE='ERR.DAT', EXIST= EXST)
   IF(EXST) THEN
    OPEN (UNIT=99,FILE='ERR.DAT',STATUS='OLD',FORM='FORMATTED',POSITION='APPEND')
   ELSE
    OPEN (UNIT=99,FILE='ERR.DAT',STATUS='NEW',FORM='FORMATTED')
   END IF
    WRITE(99,'(A)') TEXT
   CALL PERROR('Fatal Error, see file ERR.DAT')
   CALL EXIT(1)
  END SUBROUTINE Error_Message

  SUBROUTINE Solve(Lhs,Rhs,F)
  !---------------------------------------------
  !    Solution of system of equations
  !    by Gauss Elimination
  !---------------------------------------------
  IMPLICIT NONE
  REAL(KIND=8)	::    Lhs(:,:)    !    Equation Left hand side
  REAL(KIND=8)	::    Rhs(:)      !    Equation right hand side
  REAL(KIND=8)	::    F(:)        !    Unknowns
  REAL(KIND=8)	::    FAC
  INTEGER				::		M           !    Size of system
  INTEGER				::		i,n
  M= UBOUND(Lhs,1)
  !  Reduction
  Equation_n:DO n=1,M-1
                IF(Lhs(n,n) < 1.0E-14 .and. Lhs(n,n) > -1.0E-14) THEN
     CALL Error_Message('Singular Matrix')
                END IF
    Equation_i: DO i=n+1,M
                 FAC= Lhs(i,n)/Lhs(n,n)
                 Lhs(i,n+1:M)= Lhs(i,n+1:M) - Lhs(n,n+1:M)*FAC
                 Rhs(i)= Rhs(i) - Rhs(n)*FAC
                END DO   Equation_i
             END DO Equation_n
  !     Backsubstitution
  Unknown_n: DO n= M,1,-1
                F(n)= -1.0/Lhs(n,n)*(SUM(Lhs(n,n+1:M)*F(n+1:M)) - Rhs(n))
        END DO Unknown_n
  RETURN
  END SUBROUTINE Solve

 SUBROUTINE Assembly(Lhs,Rhs,DTe,DUe,Ldest,BCode,Ncode &
           ,Elres_u,Elres_te,Diag,Ndofe,Ndof,Nodel,Fac)
!---------------------------------------------
!  Assembles Element contributions DTe , DUe
!  into global matrix Lhs and vector Rhs
!  Also sums off-diagonal coefficients 
!  for the computation of diagonal coefficients
!---------------------------------------------
REAL(KIND=8)            :: Lhs(:,:),Rhs(:)    ! Global arrays
REAL(KIND=8), INTENT(IN):: DTe(:,:),DUe(:,:)  ! Element arrays
INTEGER , INTENT(IN)    :: LDest(:)						! Element destination vector
INTEGER , INTENT(IN)		:: BCode(:)						! Boundary code(local)
INTEGER , INTENT(IN)		:: NCode(:)						! Boundary code (global) 
INTEGER , INTENT(IN)		:: Ndofe							! D.o.F큦 / Elem
INTEGER , INTENT(IN)		:: Ndof								! D.o.F큦 / Node
INTEGER , INTENT(IN)		:: Nodel							! Nodes/Element
REAL , INTENT(IN)				:: Elres_u(:)					! vector u for element
REAL , INTENT(IN)				:: Elres_te(:)				! vector t for element
REAL , INTENT(IN)				:: Fac(:)							! Mult. factors for symmetry  
REAL(KIND=8)						:: Diag(:,:)					! Array containing diagonal coeff of DT
INTEGER									:: n,Ncol
DoF_per_Element:&
DO m=1,Ndofe  
	Ncol=Ldest(m)      !   Column number 
	IF(BCode(m) == 0) THEN    !   Neumann BC
		Rhs(:) = Rhs(:) + DUe(:,m)*Elres_te(m)*Fac(m)
	!     The assembly of dTe depends on the global BC
		IF (NCode(Ldest(m)) == 0 .and. Ncol /= 0) THEN	
			Lhs(:,Ncol)=  Lhs(:,Ncol) + DTe(:,m)*Fac(m)
		ELSE
			Rhs(:) = Rhs(:) - DTe(:,m) * Elres_u(m)*Fac(m)
		END IF
	END IF
	IF(BCode(m) == 1) THEN   !   Dirichlet BC
		Lhs(:,Ncol) = Lhs(:,Ncol) - DUe(:,m)*Fac(m)
		Rhs(:)= Rhs(:) - DTe(:,m) * Elres_u(m)*Fac(m)
	END IF
END DO &
DoF_per_Element
!				Sum of off-diagonal coefficients
DO n=1,Nodel
	DO k=1,Ndof
		l=(n-1)*Ndof+k
		Diag(:,k)= Diag(:,k) - DTe(:,l)
	END DO
END DO
RETURN
END SUBROUTINE Assembly

SUBROUTINE Mirror(Isym,nsy,Nodes,Elcor,Fac,Incie,Ldeste,Elres_te,Elres_ue &
                 ,Nodel,Ndof,Cdim)
!--------------------------------------------
!     Creates mirror image of element
!--------------------------------------------
IMPLICIT NONE
INTEGER, INTENT(IN)			::	Isym       ! symmetry indicator 
INTEGER, INTENT(IN)			::  nsy        ! symmetry count
INTEGER, INTENT(IN)			::  nodes      ! highest node no
REAL, INTENT(IN OUT)		::  Elcor(:,:) ! Coords (will be modified)
REAL, INTENT(IN OUT)		::	Elres_te(:)! Tractions of Element	
REAL, INTENT(IN OUT)		::	Elres_ue(:)! Displacements of Element	
REAL, INTENT(OUT)				::  Fac(:)     ! Multiplication factors 
INTEGER, INTENT(IN OUT)	::  Incie(:)	 ! Incidences     (will be
INTEGER, INTENT(IN OUT)	::	Ldeste(:)	 ! Destinations    modified)
INTEGER, INTENT(IN)			::  Nodel      ! Nodes per element
INTEGER, INTENT(IN)			::  Ndof       ! d.o.F. per Node 
INTEGER, INTENT(IN)			::  Cdim       ! Cartesian dimension
REAL  :: TD(3) ! Transformation vector (diagonal elements of T)
INTEGER :: n,m,Ison1,Ison2,Ison3,i
IF(nsy == 1)RETURN
!     Assign coefficients of TD 
SELECT CASE (nsy-1)
	CASE(1)
		TD=(/-1.0,1.0,1.0/)
	CASE(2)
		TD=(/-1.0,-1.0,1.0/) 
	CASE(3)
		TD=(/1.0,-1.0,1.0/) 
	CASE(4)
		TD=(/1.0,1.0,-1.0/) 
	CASE(5)
		TD=(/-1.0,1.0,-1.0/) 
	CASE(6)
		TD=(/-1.0,-1.0,-1.0/) 
	CASE(7)
		TD=(/1.0,-1.0,-1.0/) 
END SELECT
!     generate coordinates and incidences
Nodes1: &
DO n=1,nodel
	Direction: &
	DO m=1,Cdim
		Elcor(m,n)= Elcor(m,n)*TD(m)
	END DO & 
	Direction 
	!	 Check if point is on any symmetry plane
	Ison1= 0 
	Ison2= 0 
	Ison3= 0
	IF(Elcor(1,n)==0.0) Ison1=1  
	IF(Elcor(2,n)==0.0) Ison2=1  
	IF(Cdim > 2 .AND. Elcor(3,n)==0.0) Ison3=1  
	!   only change incidences for unprimed nodes
	IF(ison1==1 .AND. nsy-1==1)CYCLE
	IF(ison2==1 .AND. nsy-1==3) CYCLE
	IF(ison1+ison2+ison3 > 1 .AND. nsy-1<4) CYCLE
	Incie(n)= Incie(n) + Nodes
END DO &
Nodes1
!     generate multiplication factors elast. Problems only
IF(Ndof > 1) THEN
	I=0
	Nodes2: &
	DO n=1,nodel
		Degrees_of_freedom1: &
		DO m=1,Ndof
			I=I+1
			Fac(I)= TD(m)		! Multiplication factor for symmetry
		END DO & 
		Degrees_of_freedom1
	END DO &
	Nodes2
END IF
!   Reverse destination vector for selected elem
SELECT CASE (nsy-1)
CASE (1,3,4,6)
	CALL Reverse(Incie,elcor,ldeste,Elres_te,Elres_ue,Ndof,Cdim,nodel)
CASE DEFAULT
END SELECT
RETURN
END SUBROUTINE Mirror

SUBROUTINE Reverse(Inci,elcor,ldest,Elres_te,Elres_ue,Ndof,Cdim,nodel)
!--------------------------------------
!  reverses incidences, destination vector 
!  and co-ordinates
!  so that outward normal is reversed
!--------------------------------------
IMPLICIT NONE
INTEGER, INTENT (INOUT) :: Inci(:)			!   Incidences
REAL, INTENT (INOUT)    :: Elcor(:,:)	  !   Coordinates
REAL, INTENT (INOUT)    :: Elres_te(:)  !   Tractions of Element
REAL, INTENT (INOUT)    :: Elres_ue(:)  !   Displacements of Element	
INTEGER, INTENT (INOUT) :: Ldest(:)	    !   Destination vector
INTEGER, INTENT (IN)		:: Ndof 				!   No of degrees of freedom per node
INTEGER, INTENT (IN)		:: Cdim					!   Cartesian dimension
INTEGER, INTENT (IN)		:: Nodel				!   No of nodes per element
REAL, ALLOCATABLE				:: Elcort(:,:)					!   Temps
REAL, ALLOCATABLE				:: Elres_tet(:)					!   Temps
REAL, ALLOCATABLE				:: Elres_uet(:)					!   Temps
INTEGER, ALLOCATABLE		:: Incit(:),Ldestt(:)   !   Temps
INTEGER									:: Node(8)          !    reversing sequence
INTEGER									:: n,nc,Nchanges
ALLOCATE (Incit(nodel),Elcort(Cdim,nodel),Ldestt(Nodel*ndof),Elres_tet(Nodel*ndof),Elres_uet(Nodel*ndof))
Incit= Inci
Elcort= Elcor
Ldestt= Ldest
Elres_tet= Elres_te
Elres_uet= Elres_ue
SELECT CASE (Cdim)
	CASE (2)    !   2-D problem
		Node(1:2)= (/2,1/) ;  Nchanges= 2
	CASE (3)    !    3-D  problem
		Node= (/1,4,3,2,8,7,6,5/) ; Nchanges= nodel !-1
END SELECT
Number_changes: &
DO n=1,Nchanges
	nc= Node(n)
	inci(n)= Incit(nc) ; Elcor(:,n)= Elcort(:,nc)
	Ldest(Ndof*(n-1)+1:Ndof*n)= Ldestt(Ndof*(nc-1)+1:Ndof*nc)
	Elres_te(Ndof*(n-1)+1:Ndof*n)=Elres_tet(Ndof*(nc-1)+1:Ndof*nc) 
	Elres_ue(Ndof*(n-1)+1:Ndof*n)=Elres_uet(Ndof*(nc-1)+1:Ndof*nc) 
END DO &
Number_changes
DEALLOCATE(Incit,Elcort,Ldestt,Elres_tet,Elres_uet)
RETURN
END SUBROUTINE Reverse

SUBROUTINE Jobin(Title,Cdim,Ndof,Toa,Nreg,Ltyp,Con,E,ny &
                ,Isym,nodel,nodes,maxe)
!------------------------------------------------
!    Subroutine to read in basic job information
!------------------------------------------------
CHARACTER(LEN=80), INTENT(OUT):: Title
INTEGER, INTENT(OUT) :: Cdim,Ndof,Toa,Nreg,Ltyp,Isym,nodel
INTEGER, INTENT(OUT) :: Nodes,Maxe
REAL, INTENT(OUT)    :: Con,E,ny
READ(1,'(A80)') Title
WRITE(2,*)'Project:',Title
READ(1,*) Cdim
WRITE(2,*)'Cartesian_dimension:',Cdim
READ(1,*) Ndof        !    Degrees of freedom per node
IF(NDof == 1) THEN
	WRITE(2,*)'Potential Problem'
ELSE
	WRITE(2,*)'Elasticity Problem'
END IF
IF(Ndof == 2)THEN 
	READ(1,*) Toa				! Toa ....Type of analysis (solid plane strain = 1,solid plane stress = 2)
	IF(Toa == 1)THEN
		WRITE(2,*)'Type of Analysis: Solid Plane Strain'
	ELSE
		WRITE(2,*)'Type of Analysis: Solid Plane Stress'
	END IF
END IF	
READ(1,*) Nreg       !   Type of region
IF(NReg == 1) THEN
	WRITE(2,*)'Finite Region'
ELSE
	WRITE(2,*)'Infinite Region'
END IF
READ(1,*) Isym       !   Symmetry code
SELECT CASE (isym)
CASE(0)
WRITE(2,*)'No symmetry'
CASE(1) 
WRITE(2,*)'Symmetry about y-z plane' 
CASE(2)
WRITE(2,*)'Symmetry about y-z and x-z planes'
CASE(3) 
WRITE(2,*)'Symmetry about all planes'
END SELECT
READ(1,*) Ltyp        !   Element type
IF(Ltyp == 1) THEN
WRITE(2,*)'Linear Elements'
ELSE
WRITE(2,*)'Quadratic Elements'
END IF
!     Determine number of nodes per element
IF(Cdim == 2) THEN    !    Line elements
 IF(Ltyp == 1) THEN
  Nodel= 2
 ELSE
  Nodel= 3
 END IF
ELSE                  !    Surface elements
 IF(Ltyp == 1) THEN
  Nodel= 4
 ELSE
  Nodel= 8
 END IF
END IF
!   Read properties
IF(Ndof == 1) THEN
	READ(1,*) Con
	WRITE(2,*)'Conductivity=',Con
ELSE	
	READ(1,*) E,ny
	IF(ToA == 2)THEN		!	Solid Plane Stress
		E=E*((1+2*ny)/(1+ny)**2)
		ny =ny/(1+ny)
	END IF			
	WRITE(2,*)'Modulus:',E
	WRITE(2,*)'Poissons ratio:',ny
END IF

READ(1,*) Nodes
WRITE(2,*)'Number of Nodes of System:',Nodes
READ(1,*) Maxe     
WRITE(2,*)'Number of Elements of System:', Maxe
RETURN
END SUBROUTINE Jobin

SUBROUTINE JobinMR(Title,Cdim,Ndof,Toa,Ltyp,Isym,nodel,nodes,maxe)
!------------------------------------------------
!    Subroutine to read in basic job information
!------------------------------------------------
IMPLICIT NONE
INTEGER, INTENT(OUT)			:: Cdim							!   Cartesian dimension
INTEGER										:: ldim							!   Dimension of Element
INTEGER, INTENT(OUT)			:: Ndof							!   No. of degeres of freedom per node
INTEGER, INTENT(OUT)			:: Toa							!   Type of analysis (plane strain = 1, plane stress = 2)
INTEGER										:: Ltyp							!   Element type(linear = 1, quadratic = 2)
INTEGER, INTENT(OUT)			:: Nodel						!   No. of nodes per element
INTEGER, INTENT(OUT)			:: Nodes						!   No. of nodes of System
INTEGER, INTENT(OUT)			:: Maxe							!   Number of Elements of System
INTEGER										:: Isym							!   Symmetry code
CHARACTER(LEN=80)					:: Title						!		Title of calculation
INTEGER										:: nr,nb

READ(1,'(A80)') Title
WRITE(2,*)'Project:',Title
READ(1,*) Cdim
WRITE(2,*)'Cartesian_dimension:',Cdim
ldim= Cdim - 1
READ(1,*) Ndof        !    Degrees of freedom per node
IF(NDof == 1) THEN
	WRITE(2,*)'Potential Problem'
ELSE
	WRITE(2,*)'Elasticity Problem'
END IF
IF(Ndof == 2)THEN 
	READ(1,*) Toa				! Toa ....Type of analysis (solid plane strain = 1,solid plane stress = 2)
	IF(Toa == 1)THEN
		WRITE(2,*)'Type of Analysis: Solid Plane Strain'
	ELSE
		WRITE(2,*)'Type of Analysis: Solid Plane Stress'
	END IF
END IF	
READ(1,*) Ltyp        !   Element type
IF(Ltyp == 1) THEN
WRITE(2,*)'Linear Elements'
ELSE
WRITE(2,*)'Quadratic Elements'
END IF
!     Determine number of nodes per element
IF(Cdim == 2) THEN    !    Line elements
 IF(Ltyp == 1) THEN
	Nodel= 2
 ELSE
	Nodel= 3
 END IF
ELSE                  !    Surface elements
 IF(Ltyp == 1) THEN
	Nodel= 4
 ELSE
	Nodel= 8
 END IF
END IF
READ(1,*) Nodes
WRITE(2,*)'Number of Nodes of System:',Nodes
READ(1,*) Maxe     
WRITE(2,*)'Number of Elements of System:', Maxe
END SUBROUTINE JobinMR

SUBROUTINE Reg_Info(Nregs,ToA,Ndof,TypeR,ConR,ER,nyR,Nbel,ListR)
!----------------------------------------------------------------
!    Subroutine to read in basic job information for each region
!----------------------------------------------------------------
IMPLICIT NONE
INTEGER,INTENT(IN)				::	Nregs					! Number of Regions	
INTEGER,INTENT(IN)				::	ToA						! Type of analysis (solid plane strain = 1,solid plane stress = 2)
INTEGER,INTENT(INOUT)			::	TypeR(:)			! Type of BE-regions (1 == finite, 2 == Infinite)
INTEGER,INTENT(INOUT)			::	Nbel(:)				! Number of Boundary Elements each region
INTEGER,INTENT(INOUT)			::	ListR(:,:)		! List of Elementnumbers each region	
INTEGER,INTENT(IN)				::	Ndof					! No. of degeres of freedom per node
REAL,INTENT(INOUT)				::	ConR(:)				! Conductivity of each region
REAL,INTENT(INOUT)				::	ER(:)					! Youngsmodulus of regions
REAL,INTENT(INOUT)				::	nyR(:)				! Poissons ratio of regions
INTEGER										::	Isym					! Symmetrycode
INTEGER										::	nr,nb
ListR= 0
Region_loop: &
DO nr=1,Nregs
	WRITE(2,*)' Region ',nr
	READ(1,*) TypeR(nr)       !   Type of region
	IF(TypeR(nr) == 1) THEN
		WRITE(2,*)'Finite Region'
	ELSE
		WRITE(2,*)'Infinite Region'
	END IF
	READ(1,*) Isym       !   Symmetry code
	SELECT CASE (Isym)
	CASE(0)
	WRITE(2,*)'No symmetry'
	CASE(1) 
	WRITE(2,*)'Symmetry about y-z plane' 
	CASE(2)
	WRITE(2,*)'Symmetry about y-z and x-z planes'
	CASE(3) 
	WRITE(2,*)'Symmetry about all planes'
	END SELECT
	!   Read properties
	IF(Ndof == 1) THEN
		READ(1,*) ConR(nr)
		WRITE(2,*)'Conductivity=',ConR(nr)
	ELSE	
		READ(1,*) ER(nr),nyR(nr)
		IF(ToA == 2)THEN			!	Solid Plane Stress
			ER(nr)=ER(nr)*((1+2*nyR(nr))/(1+nyR(nr))**2)
			nyR(nr) =nyR(nr)/(1+nyR(nr))
		END IF			
		WRITE(2,*)'Youngsmodulus:',ER(nr)
		WRITE(2,*)'Poissons ratio:',nyR(nr)
	END IF
	READ(1,*)Nbel(nr)
!	IF(Nbel(nr) > MaxeR)MaxeR= Nbel(nr)
	READ(1,*)(ListR(nr,nb),nb=1,Nbel(nr))
	WRITE(2,*) ' List of Boundary Elements: '
	WRITE(2,*)(ListR(nr,nb),nb=1,Nbel(nr))
END DO &
Region_loop
RETURN
END SUBROUTINE Reg_info

SUBROUTINE BCInput(Elres_u,Elres_t,Bcode,nodel,ndofe,ndof) 
!------------------------------------------
!			Reads boundary conditions
!-------------------------------------------
REAL,INTENT(OUT)    :: Elres_u(:,:)  !  Element results , u
REAL,INTENT(OUT)    :: Elres_t(:,:)  !  Element results , t 
INTEGER,INTENT(OUT) :: BCode(:,:)    !  Element BC큦
INTEGER,INTENT(IN)  :: nodel         !  Nodes per element
INTEGER,INTENT(IN)  :: ndofe         !  D.o.F. per Element
INTEGER,INTENT(IN)  :: ndof          !  D.o.F per Node
INTEGER :: NE_u,NE_t  
WRITE(2,*)''
WRITE(2,*)'Elements with Dirichlet BC큦: '
WRITE(2,*)''
Elres_u(:,:)=0  ! Default prescribed values for u = 0.0
BCode = 0				! Default BC= Neumann Condition			
READ(1,*)NE_u		
IF(NE_u > 0) THEN
Elem_presc_displ: &
DO n=1,NE_u
	READ(1,*) Nel,(Elres_u(Nel,m),m=1,Ndofe)			
!	READ(1,*) Nel,(BCode(Nel,m),m=1,Ndofe)
	BCode(Nel,:)=1
	WRITE(2,*)'Element ',Nel,'  Prescribed values: '
	Na= 1
	Nodes: &
	DO M= 1,Nodel
		WRITE(2,*) Elres_u(Nel,na:na+ndof-1)
		Na= na+Ndof
	END DO &
	Nodes
END DO &
Elem_presc_displ
END IF
WRITE(2,*)''
WRITE(2,*)'Elements with Neuman BC큦: '
WRITE(2,*)''
Elres_t(:,:)=0   !   Default prescribed values = 0.0
READ(1,*)NE_t			
Elem_presc_trac:  &
DO n=1,NE_t
	READ(1,*) Nel,(Elres_t(Nel,m),m=1,Ndofe)						
	WRITE(2,*)'Element ',Nel,'  Prescribed values: '
	Na= 1
	Nodes1: &
	DO M= 1,Nodel
		WRITE(2,*) Elres_t(Nel,na:na+ndof-1)
		Na= na+Ndof
	END DO &
	Nodes1
END DO &
Elem_presc_trac
RETURN
END SUBROUTINE BCInput

SUBROUTINE Geomin(Nodes,Maxe,xp,Inci,Nodel,Cdim)
!------------------------------------
!   Inputs mesh geometry 
!-------------------------------------
IMPLICIT NONE
INTEGER, INTENT(IN) ::	Nodes			!   Number of nodes
INTEGER, INTENT(IN) ::	Maxe			!   Number of elements
INTEGER, INTENT(IN) ::	Nodel			!   Number of Nodes of elements
INTEGER, INTENT(IN) ::	Cdim			!   Cartesian Dimension
REAL, INTENT(OUT)   ::	xP(:,:)		!   Node co-ordinates
REAL								::	xmax(Cdim),xmin(Cdim),delta_x(Cdim)
INTEGER, INTENT(OUT)::	Inci(:,:) !   Element incidences
INTEGER							::	Node,Nel,M,n
!-------------------------------------------------------
!		Read Node Co-ordinates from Inputfile
!-------------------------------------------------------
DO Node=1,Nodes
 READ(1,*) (xP(M,Node),M=1,Cdim)
 WRITE(2,'(A5,I5,A8,3F8.2)') 'Node ',Node,&
         '  Coor  ',(xP(M,Node),M=1,Cdim)
END DO

!-------------------------------------------------------
!		Read Incidences from Inputfile
!-------------------------------------------------------
WRITE(2,*)''
WRITE(2,*)'Incidences: '
WRITE(2,*)''
Elements_1:&
	DO Nel=1,Maxe
READ(1,*) (Inci(Nel,n),n=1,Nodel)
WRITE(2,'(A3,I5,A8,24I5)')'EL ',Nel,'  Inci  ',Inci(Nel,:)
END DO &
Elements_1
RETURN
END SUBROUTINE Geomin

LOGICAL FUNCTION Match(Inci1,Inci2)
!-----------------------------------
!    Returns a value of TRUE if the incidences 
!    Inci1 and Inci2 match
!------------------------------------
IMPLICIT NONE
INTEGER, INTENT (IN) :: Inci1(:) !  1. incidence array
INTEGER, INTENT (IN) :: Inci2(:) !  2. incidence array
INTEGER :: Nodes,Node,N1,N2,Ncount
Nodes= UBOUND(Inci1,1)
Ncount= 0
Node_loop1: &
DO N1=1,Nodes
	Node= Inci1(n1)
	Node_loop2: &
	DO N2=1,Nodes
		IF(Node == Inci2(n2)) Ncount= Ncount+1
	END DO &
	Node_loop2
END DO &
Node_loop1
IF(Ncount == Nodes) THEN
 Match= .TRUE.
ELSE
 Match= .FALSE.
END IF
END FUNCTION Match

SUBROUTINE Destination(Isym,Ndest,Ldest,xP,Inci,Ndofs,nodes,Ndof,Nodel,Maxe)
!-------------------------------------------------------------------
! Determine Node destination vector and Element destination vector  
!-------------------------------------------------------------------
IMPLICIT NONE
REAL, INTENT (IN)		  					::	xP(:,:)					!	Node co-ordinates
INTEGER, INTENT (IN OUT)				::	Ndest(:,:)			! Node destination vector
INTEGER, INTENT (IN OUT)				::	Ldest(:,:)			! Element destination vector
INTEGER, INTENT (IN OUT)				::	Ndofs						! DoF's of System
INTEGER, INTENT (IN)						::	Inci(:,:)				!	Element Incidences
INTEGER, INTENT (IN)						::	Isym,nodes,Ndof,Nodel,Maxe 
INTEGER													::	k,m,n,Nel,l
!---------------------------------------------------
!     Determine Node destination vector
!			Set Ndest == 0 if Point is on a symmetry plane
!---------------------------------------------------
!		no symmetry
IF(Isym == 0) THEN
	k=1
	Nodes0:	DO m=1, nodes
						DO n=1, Ndof 				
							Ndest(m,n)= k
							k=k+1
						END DO
					END DO Nodes0
!		y-z symmetry
ELSE IF(Isym == 1) THEN
	k=1
	Nodes1:	DO m=1, nodes
						IF(xP(1,m) == 0.0)THEN
							Ndest(m,1)= 0
						ELSE
							Ndest(m,1)= k
							k=k+1
						END IF
						DO n=2, Ndof
							Ndest(m,n)= k
							k=k+1
						END DO
					END DO Nodes1
!		x-z and y-z symmetry
ELSE IF(Isym == 2) THEN
	k=1
	Nodes2:	DO m=1, nodes
						IF(xP(1,m) == 0.0)THEN
							Ndest(m,1)= 0
						ELSE
							Ndest(m,1)= k
							k=k+1
						END IF
						IF(xP(2,m) == 0.0)THEN
							Ndest(m,2)= 0
						ELSE
							Ndest(m,2)= k
							k=k+1
						END IF
						IF (Ndof == 3) THEN
							Ndest(m,3)= k
							k=k+1
						END IF
					END DO Nodes2
!		x-y, x-z and y-z symmetry
ELSE
	k=1
	Nodes3: DO m=1, nodes
						IF(xP(1,m) == 0.0)THEN
							Ndest(m,1)= 0
						ELSE
							Ndest(m,1)= k
							k=k+1
						END IF
						IF(xP(2,m) == 0.0)THEN
							Ndest(m,2)= 0
						ELSE
							Ndest(m,2)= k
							k=k+1
						END IF
						IF(xP(3,m) == 0.0)THEN
							Ndest(m,3)= 0
						ELSE
							Ndest(m,3)= k
							k=k+1
						END IF
					END DO Nodes3
END IF						
Ndofs= k-1							! DoF's of System
!------------------------------------------
!     Determine Element destination vector
!---------------------------------------------
Elements:&
DO Nel=1,Maxe
	DO n=1,Nodel
		k= (n-1)*Ndof+1
		l= n*Ndof
		Ldest(Nel,k:l)=Ndest(Inci(Nel,n),:)	
	END DO		   
END DO &
Elements
END SUBROUTINE Destination

SUBROUTINE Scal(E,xP,Elres_u,Elres_t,Cdim,Scad,Scat)
IMPLICIT NONE
REAL, INTENT (INOUT)  					::	E								!	Youngsmodulus
REAL, INTENT (INOUT)  					::	xP(:,:)				!	Node co-ordinates
REAL, INTENT (INOUT)  					::	Elres_u(:,:)		! Element results , u
REAL, INTENT (INOUT)  					::	Elres_t(:,:)		! Element results , t
REAL, INTENT (OUT)							::	Scad
REAL, INTENT (OUT)							::	Scat
INTEGER, INTENT(IN)							::	Cdim						!   Cartesian Dimension
REAL														::	xmax(Cdim),xmin(Cdim),delta_x(Cdim)

!-------------------------------------------------------
!		Determine Scalefactor for Tractions
!		Scat ... 1/E
!-------------------------------------------------------
Scat= 1./E													! Scalefactor for Tractions
E=1.0																! Scaled Youngsmodulus
Elres_t=Elres_t*Scat								! Scaled prescribed Tractions by Scat
!----------------------------------------------------------
!		Determine Scalefactor for Displacements
!		Scad ... max. Distance in any co-ordinate direction
!----------------------------------------------------------
xmax(1)= MAXVAL(xp(1,:))
xmax(2)= MAXVAL(xp(2,:))
IF(Cdim == 3)xmax(3)= MAXVAL(xp(3,:))
xmin(1)= MINVAL(xp(1,:))
xmin(2)= MINVAL(xp(2,:))
IF(Cdim == 3)xmin(3)= MINVAL(xp(3,:))
delta_x= xmax - xmin
Scad= MAXVAL(delta_x)										!	Scad ... Scalefactor for Displacements	
xP=xP/Scad															! Scaled Node co-ordinates
Elres_u=Elres_u/Scad										! Scaled prescribed Displacements by Scad
END SUBROUTINE Scal

SUBROUTINE Scal_Discont(E,xP,xP_Disc,Elres_u,Elres_t,Cdim,Scad,Scat)
IMPLICIT NONE
REAL, INTENT (INOUT)  					::	E								!	Youngsmodulus
REAL, INTENT (INOUT)  					::	xP(:,:), xP_Disc(:,:)				!	Node co-ordinates
REAL, INTENT (INOUT)  					::	Elres_u(:,:)		! Element results , u
REAL, INTENT (INOUT)  					::	Elres_t(:,:)		! Element results , t
REAL, INTENT (OUT)							::	Scad
REAL, INTENT (OUT)							::	Scat
INTEGER, INTENT(IN)							::	Cdim						!   Cartesian Dimension
REAL														::	xmax(Cdim),xmin(Cdim),delta_x(Cdim)

!-------------------------------------------------------
!		Determine Scalefactor for Tractions
!		Scat ... 1/E
!-------------------------------------------------------
Scat= 1./E													! Scalefactor for Tractions
E=1.0																! Scaled Youngsmodulus
Elres_t=Elres_t*Scat								! Scaled prescribed Tractions by Scat
!----------------------------------------------------------
!		Determine Scalefactor for Displacements
!		Scad ... max. Distance in any co-ordinate direction
!----------------------------------------------------------
xmax(1)= MAXVAL(xp(1,:))
xmax(2)= MAXVAL(xp(2,:))
IF(Cdim == 3)xmax(3)= MAXVAL(xp(3,:))
xmin(1)= MINVAL(xp(1,:))
xmin(2)= MINVAL(xp(2,:))
IF(Cdim == 3)xmin(3)= MINVAL(xp(3,:))
delta_x= xmax - xmin
Scad= MAXVAL(delta_x)										!	Scad ... Scalefactor for Displacements	
xP=xP/Scad															! Scaled Node co-ordinates
xP_Disc=xP_Disc/Scad
Elres_u=Elres_u/Scad										! Scaled prescribed Displacements by Scad
END SUBROUTINE Scal_Discont

LOGICAL FUNCTION HasEntry(M, D)
!-------------------------------------------------------------------
! Cecks if an integer entry is in an integer matrix
!-------------------------------------------------------------------
IMPLICIT NONE
INTEGER, INTENT (IN)  					::	M(:)									!	matrix
INTEGER, INTENT (IN)  					::	D										!	possible matrix entry
INTEGER												:: k,  a, b


a= UBOUND(M,1)

hasEntry=.false.
DO k=1, a
		IF (M(k) == D)THEN
			HasEntry=.true.
		END IF
END DO
RETURN
END FUNCTION


END MODULE Utility_lib
