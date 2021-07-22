MODULE Postproc_lib
USE Geometry_lib ; USE Utility_lib 
IMPLICIT NONE
CONTAINS
SUBROUTINE BFLOW(Flow,xsi,eta,u,Inci,Elcor,k)
!----------------------------------------------
!    Computes flow vectors in direction tangential to the
!    Boundary 
!-----------------------------------------------
REAL , INTENT(OUT)   :: Flow(:)  !  Flow vector
REAL , INTENT(IN)    :: xsi,eta  !  intrinsic coordinates of point
REAL , INTENT(IN)    :: u(:,:)     !  Nodal temperatures/potentials
INTEGER, INTENT (IN) :: Inci(:)  !  Element Incidences
REAL, INTENT (IN)    :: Elcor(:,:)  !  Element coordinates
REAL, INTENT (IN)    :: k        !  Conductivity   
REAL, ALLOCATABLE    :: Vxsi(:),Veta(:),DNi(:,:),V3(:)	
INTEGER :: Nodes,Cdim,Ldim
REAL :: Jxsi,Jeta,Flows(2),v1(3),v2(3),CosA,CosB,CosG,CosT
Nodes= UBOUND(ELCOR,2)  !  Number of nodes
Cdim= UBOUND(ELCOR,1)   !  Cartesian Dimension
Ldim= Cdim-1            !  Local (element) dimension
ALLOCATE (Vxsi(cdim),Dni(Nodes,Ldim),v3(cdim))
IF(ldim > 1)	ALLOCATE (Veta(cdim))
!   Compute Vector(s) tangential to boundary surface
CALL Serendip_deriv(DNi,xsi,eta,ldim,nodes,inci)
Vxsi(1)= Dot_Product(Dni(:,1),Elcor(1,:))
Vxsi(2)= Dot_Product(Dni(:,1),Elcor(2,:))
IF(Cdim == 2) THEN
	CALL Vector_norm(Vxsi,Jxsi)
	Flow(1)= -k*Dot_product(Dni(:,1),u(:,1))/Jxsi
ELSE
	Vxsi(3)= Dot_Product(Dni(:,1),Elcor(3,:))
	CALL Vector_norm(Vxsi,Jxsi)
	Veta(1)= Dot_Product(Dni(:,2),Elcor(1,:))
	Veta(2)= Dot_Product(Dni(:,2),Elcor(2,:))
	Veta(3)= Dot_Product(Dni(:,2),Elcor(3,:))
	CALL Vector_norm(Veta,Jeta)
	!   Flows in skew coordinate system
	Flows(1)= -k*Dot_product(Dni(:,1),u(:,1))/Jxsi
	Flows(2)= -k*Dot_product(Dni(:,2),u(:,1))/Jeta
	!   Orthoginal system
	v3= Vector_ex(Vxsi,Veta)
	Call Ortho(v3,v1,v2)
	CosA= DOT_Product(Vxsi,v1)
	CosB= DOT_Product(Veta,v2)
	CosG= DOT_Product(Veta,v1)
	CosT= DOT_Product(Vxsi,v2)
	Flow(1)= Flows(1)*CosA**2 + Flows(2)* CosG**2
	Flow(2)= Flows(1)*CosT**2 + Flows(2)* CosB**2
END IF
RETURN
END SUBROUTINE BFLOW

SUBROUTINE BStress(Stress,xsi,eta,u,t,Inci,Elcor,E,ny,IPS)
!----------------------------------------------
!    Computes stresses in a plane tangential to the
!    Boundary Element
!-----------------------------------------------
REAL , INTENT(OUT)   :: Stress(:)!  Stress vector
REAL , INTENT(IN)    :: xsi,eta  !  intrinsic coordinates of point
REAL , INTENT(IN)    :: u(:,:)   !  Nodal displacements
REAL , INTENT(IN)    :: t(:,:)   !  Nodal Tractions
INTEGER, INTENT (IN) :: Inci(:)  !  Element Incidences
REAL, INTENT (IN)    :: Elcor(:,:)  !  Element coordinates
REAL, INTENT (IN)    :: E,ny     !  Elastic constants   
INTEGER, INTENT (IN) ::  IPS ! IPS= 0 plane strain; =1 plane stress
REAL, ALLOCATABLE    :: Vxsi(:),Veta(:),DNi(:,:),Ni(:),trac(:)
INTEGER :: Nodes, Cdim, Ldim
REAL :: Jxsi,Jeta,v1(3),v2(3),CosA, CosB, CosG, CosT,v3(3)
REAL :: C1,C2,G,tn,ts,ts1,ts2
REAL , ALLOCATABLE :: Dudxsi(:),Dudeta(:),Strain(:),Strains(:)
Nodes= UBOUND(ELCOR,2)  !  Number of nodes
Cdim = UBOUND(ELCOR,1)  !  Cartesian Dimension
Ldim= Cdim-1            !  Local (element) dimension
ALLOCATE (Vxsi(cdim),Veta(cdim),Dni(Nodes,Ldim),Ni(Nodes))
ALLOCATE (Dudxsi(Cdim),Dudeta(Cdim),trac(Cdim))
!   Compute Vector(s) tangential to boundary surface
CALL Serendip_deriv(DNi,xsi,eta,ldim,nodes,inci)
CALL Serendip_func(Ni,xsi,eta,ldim,nodes,inci)
trac(1)= Dot_Product(Ni,t(:,1))
trac(2)= Dot_Product(Ni,t(:,2))
Vxsi(1)= Dot_Product(Dni(:,1),Elcor(1,:))
Vxsi(2)= Dot_Product(Dni(:,1),Elcor(2,:))
IF(Cdim == 2) THEN
		ALLOCATE (Strain(1))
		CALL Vector_norm(Vxsi,Jxsi)
		V3(1)=   Vxsi(2)
		V3(2)= - Vxsi(1)
		tn= Dot_Product(v3(1:2),trac)
		ts= Dot_Product(vxsi(1:2),trac)
		DuDxsi(1)= Dot_Product(Dni(:,1),u(:,1))
		DuDxsi(2)= Dot_Product(Dni(:,1),u(:,2))
		Strain(1)= Dot_Product(DuDxsi,Vxsi)/Jxsi
		IF(IPS == 2) THEN
		 Stress(1)= E*Strain(1) + ny*tn     !     plane stress
		ELSE
		 Stress(1)= 1/(1.0-ny)*(E/(1.0+ny)*Strain(1) - ny*tn)		! Plane strain
		END IF
		Stress(2)= tn
		Stress(3)= ts
ELSE
		ALLOCATE (Strain(3),Strains(3))
		trac(3)= Dot_Product(Ni,t(:,3))
		Vxsi(3)= Dot_Product(Dni(:,1),Elcor(3,:))
		CALL Vector_norm(Vxsi,Jxsi)
		Veta(1)= Dot_Product(Dni(:,2),Elcor(1,:))
		Veta(2)= Dot_Product(Dni(:,2),Elcor(2,:))
		Veta(3)= Dot_Product(Dni(:,2),Elcor(3,:))
		CALL Vector_norm(Veta,Jeta)
		V3= Vector_ex(Vxsi,veta)
		tn= Dot_Product(v3,trac)
		DuDxsi(1)= Dot_Product(Dni(:,1),u(:,1))
		DuDxsi(2)= Dot_Product(Dni(:,1),u(:,2))
		DuDxsi(3)= Dot_Product(Dni(:,1),u(:,3))
		DuDeta(1)= Dot_Product(Dni(:,2),u(:,1))
		DuDeta(2)= Dot_Product(Dni(:,2),u(:,2))
		DuDeta(3)= Dot_Product(Dni(:,2),u(:,3))
!   Strains in skew coordinate system
		Strains(1)= Dot_product(DuDxsi,Vxsi)/Jxsi
		Strains(2)= Dot_product(DuDeta,Veta)/Jeta
		Strains(3)= Dot_product(DuDeta,Vxsi)/Jeta + &
								Dot_product(DuDxsi,Veta)/Jxsi
!   Orthogonal system
		v3= Vector_ex(Vxsi,Veta)
		CALL Ortho(v3,v1,v2)
		ts1= DOT_Product(v1,trac)
		ts2= DOT_Product(v2,trac)
		CosA= DOT_Product(Vxsi,v1)
		CosB= DOT_Product(Veta,v2)
		CosG= DOT_Product(Veta,v1)
		CosT= DOT_Product(Vxsi,v2)
!   Compute Strains
		Strain(1)= Strains(1)*CosA**2 + Strains(2)*CosG**2 + &
							 Strains(3)*CosA*CosG
		Strain(2)= Strains(1)*CosT**2 + Strains(2)*CosB**2 + &
							 Strains(3)*CosT*CosB
		Strain(3)= Strains(1)*CosG*CosT + Strains(2)*CosG*CosB + &
							 Strain(3)*(CosA*CosB+CosG*CosT)
!   Compute stresses
		C1= E/(1.0-ny**2)  ;  C2= ny/(1.0-ny) ; G=E/(1.0-2*ny)
		Stress(1)= C1*(Strain(1)+ny*strain(2))+ C2*Tn
		Stress(2)= C1*(Strain(2)+ny*strain(1))+ C2*Tn
		Stress(3)= tn
		Stress(4)= G*Strain(3)
		Stress(5)= ts1
		Stress(6)= ts2
END IF
RETURN
End SUBROUTINE BStress
END MODULE Postproc_lib
