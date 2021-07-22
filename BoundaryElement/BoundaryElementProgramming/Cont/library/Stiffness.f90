MODULE Stiffness_lib
USE Utility_lib ; USE Integration_lib ; USE Geometry_lib
CONTAINS
SUBROUTINE Stiffness_BEM(nr,xP,Nodel,Ndof,Ndofe,NodeR,Ncode,NdofR,Ndofc,KBE,A,tc,Cdim,Elres_u,Elres_t,&
												IncieR,LdesteR,Nbel,ListR,TypeR,Bcode,Con,E,ny,Ndest,Isym)
!---------------------------------------------
!    Computes the stiffness matrix
!    of a boundary element region
!    no symmetry
!--------------------------------------------
IMPLICIT NONE
REAL, INTENT(INOUT):: xP(:,:)    !  Array of node coordinates
INTEGER, INTENT(IN):: nr  
INTEGER, INTENT(IN):: Ncode(:)   !  Global restraint code
INTEGER, INTENT(IN):: NdofR      
INTEGER, INTENT(IN):: Ndofc      !  No of interface degrees of freedom
INTEGER, INTENT(IN):: NodeR(:)    
INTEGER, INTENT(IN):: TypeR(:)    
INTEGER, INTENT(IN):: Cdim    
INTEGER, INTENT(IN):: IncieR(:,:) 
INTEGER, INTENT(IN):: LdesteR(:,:) 
INTEGER, INTENT(INOUT):: Bcode(:,:) 
INTEGER, INTENT(IN):: Nbel(:) 
INTEGER, INTENT(IN):: ListR(:,:) 
INTEGER, INTENT(IN):: Nodel 
INTEGER, INTENT(IN):: Ndof
INTEGER, INTENT(IN):: Ndofe 
INTEGER, INTENT(IN):: Isym 
INTEGER, INTENT(IN):: Ndest(:,:) 
REAL, INTENT(INOUT):: Elres_u(:,:),Elres_t(:,:)  
REAL, INTENT(INOUT)	 :: E,ny,Con    
REAL(KIND=8), INTENT(OUT)  :: KBE(:,:) !  Stiffness matrix
REAL(KIND=8), INTENT(OUT)  :: A(:,:)   ! u due to únit values ui
REAL(KIND=8), INTENT(OUT)  :: tc(:)    ! interface tractions
!   temporal arrays :
REAL(KIND=8), ALLOCATABLE :: dUe(:,:),dTe(:,:),Diag(:,:)
REAL(KIND=8), ALLOCATABLE :: Lhs(:,:)
REAL(KIND=8), ALLOCATABLE :: Rhs(:),RhsM(:,:) ! right hand sides
REAL(KIND=8), ALLOCATABLE :: u1(:),u2(:,:)    ! results
REAL, ALLOCATABLE :: Elcor(:,:) 
REAL		:: Scat,Scad
INTEGER					:: NdofF
INTEGER :: Dof,k,l,nel
INTEGER :: n,m,Pos,i,j,nd,ne
ALLOCATE(dTe(NdofR,Ndofe),dUe(NdofR,Ndofe))   
ALLOCATE(Diag(NdofR,Ndof))                   
ALLOCATE(Lhs(NdofR,NdofR),Rhs(NdofR),RhsM(NdofR,NdofR))
ALLOCATE(u1(NdofR),u2(NdofR,NdofR)) 
ALLOCATE(Elcor(Cdim,Nodel))        
!------------------------------------------
!     Scaling
!------------------------------------------
CALL Scal(E,xP,Elres_u,Elres_t,Cdim,Scad,Scat)
!----------------------------------------------------------------
!  Compute and assemble element coefficient matrices
!----------------------------------------------------------------
Lhs= 0.0
Diag= 0.0
Rhs= 0.0
RhsM= 0.0
Elements_1:&
DO Nel=1,Nbel(nr)
		ne= ListR(nr,Nel)
		Elcor(:,:)= xP(:,IncieR(ne,:))      !    gather element coords
		IF(Cdim == 2) THEN
			IF(Ndof == 1) THEN
				CALL Integ2P(Elcor,IncieR(ne,:),Nodel,NodeR(nr),xP,Con,dUe,dTe,Ndest,Isym)
			ELSE
				CALL Integ2E(Elcor,IncieR(ne,:),Nodel,NodeR(nr),xP,E,ny,dUe,dTe,Ndest,Isym)  
			END IF
		ELSE
			CALL Integ3(Elcor,IncieR(ne,:),Nodel,NodeR(nr),xP,Ndof,E,ny,Con,dUe,dTe,Ndest,Isym)    
		END IF
		CALL AssemblyMR(ne,Ndof,Ndofe,Nodel,Lhs,Rhs,RhsM,DTe,DUe,LdesteR(ne,:),Ncode,Bcode,Diag,Elres_u,Elres_t,Scad)
END DO &
Elements_1
!------------------------------------------------------------
!  Add azimuthal integral for infinite regions
!------------------------------------------------------------
IF(TypeR(nr) == 2) THEN
	DO m=1, NodeR(nr)
		DO n=1, Ndof
			k=Ndof*(m-1)+n
			Diag(k,n) = Diag(k,n) + 1.0
		END DO
	END DO		 
END IF
!-------------------------------------------------------------
!  Add Diagonal coefficients
!-------------------------------------------------------------
Nodes_global: &
DO m=1,NodeR(nr)
	Degrees_of_Freedoms_node: &
	DO n=1, Ndof
		DoF = (m-1)*Ndof + n   !  global degree of freedom no.
		k = (m-1)*Ndof + 1     !  address in coeff. matrix (row)
		l = k + Ndof - 1       !  address in coeff. matrix (column)
		IF (NCode(DoF) == 1 .or. NCode(DoF) == 2) THEN		! Dirichlet - Add Diagonal to Rhs 
			Pos = 0
			Nel = 0
			!   get local degree of freedom no corresponding to global one
			Elements_all: &
			DO i=1,Nbel(nr)	
				ne= ListR(nr,i)
				Degrees_of_freedom_elem: &
				DO j=1,Ndofe
					IF (DoF == LdesteR(ne,j)) THEN
						Nel = ne
						Pos = j
						EXIT				 
					END IF
				END DO &
				Degrees_of_freedom_elem
			IF (Nel /= 0) EXIT 
			END DO &
			Elements_all
			Rhs(k:l) = Rhs(k:l) - Diag(k:l,n)*Elres_u (Nel,Pos)
			IF(NCode(DoF) == 2)THEN
				RhsM(k:l,DoF) = RhsM(k:l,DoF) - Diag(k:l,n) / Scad 		
			END IF
		ELSE
			Lhs(k:l,Dof)= Lhs(k:l,Dof) + Diag(k:l,n)	! Neuman - Add to Lhs 
		END IF	
	END DO &
	Degrees_of_Freedoms_node
END DO	&
Nodes_global	
!    Solve problem 
CALL Solve_Multi(Lhs,Rhs,RhsM,u1,u2)
!------------------------------------------
!		Back - Scaling
!------------------------------------------
DO N=1,NdofC
	u1(N)= u1(N) / Scat
	u2(N,:)= u2(N,:) / Scat
END DO
M=NdofC
NdofF= NdofR-NdofC
DO N=1,NdofF
	M=M+1
	IF(NCode(M) == 0) THEN
		u1(M)= u1(M) * Scad
		u2(M,:)= u2(M,:) * Scad
	ELSE
		u1(M)= u1(M) / Scat
		u2(M,:)= u2(M,:) / Scat
	END IF
END DO
	Elres_u(:,:)= Elres_u(:,:) * Scad
	Elres_t(:,:)= Elres_t(:,:) / Scat

!--------------------------------------
!  Gather element results due to 
!  zero Dirichlet conditions at the interface
!--------------------------------------
Elements2:	&
DO nel=1,Nbel(nr)
	ne= ListR(nr,nel)
	D_o_F1:		&
	DO nd=1,Ndofe
		IF(Ncode(LdesteR(ne,nd)) == 0) THEN
			Elres_u(ne,nd) =  u1(LdesteR(ne,nd))
		ELSE IF(Bcode(ne,nd) == 1 .or. Bcode(ne,nd) == 2) THEN
			Elres_t(ne,nd) =  u1(LdesteR(ne,nd))
		END IF
	END DO &
	D_o_F1
END DO &
Elements2

!------------------------------------
!   Gather stiffness matrix KBE and matrix A
!------------------------------------
Interface_DoFs: &
DO N=1,Ndofc
	KBE(N,:)= u2(N,:)
	tc(N)= u1(N)
END DO &
Interface_DoFs
A= 0.0
M=NdofC
Free_DoFs: &
DO N=1,NdofF
	M= M+1
	A(N,1:NdofC)= u2(M,:)
END DO &
Free_DoFs
DEALLOCATE (dUe,dTe,Diag,Lhs,Rhs,RhsM,u1,u2,Elcor)
RETURN
END SUBROUTINE Stiffness_BEM

SUBROUTINE Solve_Multi(Lhs,Rhs,RhsM,u,uM)
!---------------------------------------------
!    Solution of system of equations
!    by Gauss Elimination
!    for multple right hand sides
!---------------------------------------------
REAL(KIND=8) ::    Lhs(:,:)    !    Equation Left hand side
REAL(KIND=8) ::    Rhs(:)      !    Equation right hand side 1
REAL(KIND=8) ::    RhsM(:,:)   !    Equation right hand sides 2
REAL(KIND=8) ::    u(:)        !    Unknowns 1
REAL(KIND=8) ::    uM(:,:)     !    Unknowns 2
REAL(KIND=8) ::    FAC
INTEGER  M,Nrhs            !    Size of system
INTEGER  i,n,nr
M= UBOUND(RhsM,1) ; Nrhs= UBOUND(RhsM,2)
!  Reduction
Equation_n: &
DO n=1,M-1
   IF(ABS(Lhs(n,n)) < 1.0E-10) THEN
     CALL Error_Message('Singular Matrix')
   END IF
   Equation_i: &
	 DO i=n+1,M
     FAC= Lhs(i,n)/Lhs(n,n)
     Lhs(i,n+1:M)= Lhs(i,n+1:M) - Lhs(n,n+1:M)*FAC
     Rhs(i)= Rhs(i) - Rhs(n)*FAC
		 RhsM(i,:)= RhsM(i,:) - RhsM(n,:)*FAC
   END DO  & 
	 Equation_i
END DO &
Equation_n
!     Backsubstitution 
Unknown_1: &
DO n= M,1,-1	  
	 u(n)= -1.0/Lhs(n,n)*(SUM(Lhs(n , n+1:M)*u(n+1:M)) - Rhs(n))
END DO &
Unknown_1
Load_case: &
DO Nr=1,Nrhs
  Unknown_2: &
	DO n= M,1,-1	  
	 uM(n,nr)= -1.0/Lhs(n,n)*(SUM(Lhs(n , n+1:M)*uM(n+1:M , nr)) - RhsM(n,nr))
  END DO &
	Unknown_2
END DO &
Load_case
RETURN
END SUBROUTINE Solve_Multi

SUBROUTINE AssemblyMR(Nel,Ndof,Ndofe,Nodel,Lhs,Rhs,RhsM,DTe,DUe,Ldest,Ncode,Bcode,Diag,Elres_u,Elres_t,Scad)
!---------------------------------------------
!  Assembles Element contributions DTe , DUe
!  into global matrix Lhs, vector Rhs
!  and matrix RhsM 
!---------------------------------------------
INTEGER,INTENT(IN)      :: NEL
REAL(KIND=8)            :: Lhs(:,:) !  Eq.left hand side
REAL(KIND=8)            :: Rhs(:)   !  Right hand side
REAL(KIND=8)            :: RhsM(:,:) ! Matrix of right hand sides
REAL(KIND=8), INTENT(IN):: DTe(:,:),DUe(:,:)   !  Element arrays
REAL, INTENT(INOUT)			:: Elres_u(:,:),Elres_t(:,:)  
INTEGER , INTENT(IN)    :: LDest(:) ! Element destination vector
INTEGER , INTENT(IN) :: NCode(:) ! Boundary code (global) 
INTEGER , INTENT(IN) :: BCode(:,:) ! Boundary code (global) 
INTEGER , INTENT(IN) :: Ndof
INTEGER , INTENT(IN) :: Ndofe
INTEGER , INTENT(IN) :: Nodel
REAL(KIND=8) :: Diag(:,:) ! Array containing diagonal coeff of DT
INTEGER :: n,Ncol,m,k,l
DoF_per_Element:&
DO m=1,Ndofe  
	Ncol=Ldest(m)      !   Column number 
	IF(BCode(nel,m) == 0) THEN    !   Neumann BC
		Rhs(:) = Rhs(:) + DUe(:,m)*Elres_t(nel,m)
!     The assembly of dTe depends on the global BC
		IF (NCode(Ldest(m)) == 0) THEN	
			Lhs(:,Ncol)=  Lhs(:,Ncol) + DTe(:,m)
		ELSE
			Rhs(:) = Rhs(:) - DTe(:,m) * Elres_u(nel,m)
		END IF
	ELSE IF(BCode(nel,m) == 1) THEN   !   Dirichlet BC
		Lhs(:,Ncol) = Lhs(:,Ncol) - DUe(:,m)
		Rhs(:)= Rhs(:) - DTe(:,m) * Elres_u(nel,m)
	END IF
	IF(BCode(nel,m) == 2) THEN   !   Interface
		Lhs(:,Ncol) = Lhs(:,Ncol) - DUe(:,m)
	END IF
	IF(NCode(Ldest(m)) == 2) THEN   !   Interface
		RhsM(:,Ncol)= RhsM(:,Ncol) - DTe(:,m) / Scad
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
END SUBROUTINE AssemblyMR

END MODULE Stiffness_lib
