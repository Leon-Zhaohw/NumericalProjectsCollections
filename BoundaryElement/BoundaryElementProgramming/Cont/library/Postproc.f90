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
REAL, ALLOCATABLE    :: Vxsi(:),Veta(:),DNi(:,:),V3(:),V1(:),V2(:)
INTEGER :: Nodes,Cdim,Ldim
REAL :: Jxsi,Jeta,CosT,SinT,DuDxsi,DuDeta,V3_L
REAL :: DxsiDx,DxsiDy,DetaDx,DetaDy
Nodes= UBOUND(ELCOR,2)  !  Number of nodes
Cdim= UBOUND(ELCOR,1)   !  Cartesian Dimension
Ldim= Cdim-1            !  Local (element) dimension
ALLOCATE (Vxsi(cdim),Dni(Nodes,Ldim),v3(cdim))
ALLOCATE(v1(cdim), v2(cdim))
IF(ldim > 1)	ALLOCATE (Veta(cdim))
!   Compute Vector(s) tangential to boundary surface
CALL Serendip_deriv(DNi,xsi,eta,ldim,nodes,inci)
Vxsi(1)= Dot_Product(Dni(:,1),Elcor(1,:))
Vxsi(2)= Dot_Product(Dni(:,1),Elcor(2,:))
IF(Cdim == 2) THEN
	CALL Vector_norm(Vxsi,Jxsi)
	V1= Vxsi
	V2(1)=v1(2)
	V2(2)=-v1(1)
	Flow(1)= -k*Dot_product(Dni(:,1),u(:,1))/Jxsi
ELSE
	Vxsi(3)= Dot_Product(Dni(:,1),Elcor(3,:))
	CALL Vector_norm(Vxsi,Jxsi)
	Veta(1)= Dot_Product(Dni(:,2),Elcor(1,:))
	Veta(2)= Dot_Product(Dni(:,2),Elcor(2,:))
	Veta(3)= Dot_Product(Dni(:,2),Elcor(3,:))
	CALL Vector_norm(Veta,Jeta)
 v3= Vector_ex(Vxsi,Veta)
 Call Vector_norm(v2,v3_L)
 	v1=Vxsi
		v2= Vector_ex(v3,v1)		
		DuDxsi= Dot_Product(Dni(:,1),u(:,1))
		DuDeta= Dot_Product(Dni(:,2),u(:,1))
		CosT= DOT_Product(Vxsi,Veta)
		SinT= ABS(DOT_Product(V2,Veta))
		DxsiDx= 1/Jxsi
		DxsiDy= -CosT /(Jxsi*SinT)
		DetaDx= 0.0
		DetaDy= 1/(Jeta*SinT)
!   Flow in locall coordinate directions 
		Flow(1)= -k*DuDxsi*DxsiDx
  Flow(2)= -k*(DuDxsi*DxsiDy+DuDeta*DetaDy)	
END IF
!		Transformation of flux
CALL Flow_Transformation(v1,v2,v3,Flow,Cdim)
RETURN
END SUBROUTINE BFLOW

SUBROUTINE BStress(Stress,xsi,eta,u,t,Inci,Elcor,E,Ny,IPS)
!----------------------------------------------
!    Computes stresses in a plane tangential to the
!    Boundary Element
!-----------------------------------------------
REAL , INTENT(OUT)   :: Stress(:)!  Stress vector
REAL , INTENT(IN)    :: xsi,eta  !  intrinsic coordinates of point
REAL , INTENT(IN)    :: u(:,:)   !  Nodal displacements
REAL , INTENT(IN)    :: t(:,:)   !  Nodal Tractions
INTEGER, INTENT (IN) :: Inci(:)  !  Element Incidences
! INTEGER, INTENT (IN) :: nr			 !  Number of region
REAL, INTENT (IN)    :: Elcor(:,:)  !  Element coordinates
REAL, INTENT (IN)    :: E,Ny
INTEGER , INTENT (IN) :: IPS
REAL, ALLOCATABLE    :: Vxsi(:),Veta(:),DNi(:,:),Ni(:),trac_GP(:)
REAL :: Jxsi,Jeta,v1(3),v2(3),v3(3),v3_L,CosT,SinT
REAL :: DxsiDx, DxsiDy, DetaDx, DetaDy
REAL :: C1,C2,G,tn,ts,ts1,ts2
REAL , ALLOCATABLE :: Dudxsi(:),Dudeta(:),Strain(:)
INTEGER :: Nodes, Cdim, Ldim
Nodes= UBOUND(Elcor,2)
Cdim= UBOUND(Elcor,1)
ldim= Cdim-1
ALLOCATE (Vxsi(cdim),Veta(cdim),Dni(Nodes,Ldim),Ni(Nodes))
ALLOCATE (Dudxsi(Cdim),Dudeta(Cdim),trac_GP(Cdim))
!   Compute Vector(s) tangential to boundary surface
CALL Serendip_deriv(DNi,xsi,eta,ldim,nodes,inci)
CALL Serendip_func(Ni,xsi,eta,ldim,nodes,inci)
trac_GP(1)= Dot_Product(Ni,t(:,1))
trac_GP(2)= Dot_Product(Ni,t(:,2))
Vxsi(1)= Dot_Product(Dni(:,1),Elcor(1,:))
Vxsi(2)= Dot_Product(Dni(:,1),Elcor(2,:))
IF(Cdim == 2) THEN
!		2-D Calculation
		ALLOCATE (Strain(1))
		CALL Vector_norm(Vxsi,Jxsi)
		V1(1:2)= Vxsi
		V3(1)=   Vxsi(2)
		V3(2)= - Vxsi(1)
		tn= Dot_Product(v3(1:2),trac_GP)
		ts=	Dot_Product(v1(1:2),trac_GP)
		DuDxsi(1)= Dot_Product(Dni(:,1),u(:,1))
		DuDxsi(2)= Dot_Product(Dni(:,1),u(:,2))
		Strain(1)= Dot_Product(DuDxsi,V1(1:2))/Jxsi
!   Compute stresses in local directions, tangential and normal to the boundary
		IF(IPS == 2) THEN															!     plane stress
			Stress(1)= E*Strain(1) + ny*tn     
			Stress(2)= tn
			Stress(4)= ts
		ELSE
			Stress(1)= 1/(1.0-ny)*(E/(1.0+ny)*Strain(1) + ny*tn)		! Plane strain
			Stress(2)= tn
			Stress(4)= ts
		END IF
		DEALLOCATE (Strain)
ELSE
!		3-D Calculation
		ALLOCATE (Strain(3))
		trac_GP(3)= Dot_Product(Ni,t(:,3))
		Vxsi(3)= Dot_Product(Dni(:,1),Elcor(3,:))
		CALL Vector_norm(Vxsi,Jxsi)
		Veta(1)= Dot_Product(Dni(:,2),Elcor(1,:))
		Veta(2)= Dot_Product(Dni(:,2),Elcor(2,:))
		Veta(3)= Dot_Product(Dni(:,2),Elcor(3,:))
		CALL Vector_norm(Veta,Jeta)
		v3= Vector_ex(Vxsi,veta)
		CALL Vector_norm(v3,v3_L)
		v1=Vxsi
		v2= Vector_ex(v3,v1)		
		DuDxsi(1)= Dot_Product(Dni(:,1),u(:,1))
		DuDxsi(2)= Dot_Product(Dni(:,1),u(:,2))
		DuDxsi(3)= Dot_Product(Dni(:,1),u(:,3))
		DuDeta(1)= Dot_Product(Dni(:,2),u(:,1))
		DuDeta(2)= Dot_Product(Dni(:,2),u(:,2))
		DuDeta(3)= Dot_Product(Dni(:,2),u(:,3))
		CosT= DOT_Product(Vxsi,Veta)
		SinT= ABS(DOT_Product(V2,Veta))
		DxsiDx= 1/Jxsi
		DxsiDy= -CosT /(Jxsi*SinT)
		DetaDx= 0.0
		DetaDy= 1/(Jeta*SinT)
!   Strains 
		Strain(1)= Dot_product(DuDxsi,v1)*DxsiDx   
		Strain(2)= Dot_product(DuDxsi,v2)*DxsiDy + Dot_product(DuDeta,v2)*DetaDy
		Strain(3)= Dot_product(DuDxsi,v1)*DxsiDy + Dot_product(DuDeta,v1)*DetaDy &
							 +	Dot_product(DuDxsi,v2)*DxsiDx
		tn= Dot_Product(v3,trac_GP)
		ts1= Dot_Product(v1,trac_GP)
		ts2= Dot_Product(v2,trac_GP)
!   Compute stresses in local directions, tangential and normal to the boundary
		C1= E/(1.0-ny**2)  ;  C2= ny/(1.0-ny) ; G=E/(2*(1.0+ny))
		Stress(1)= C1*(Strain(1)+ny*strain(2))+ C2*tn
		Stress(2)= C1*(Strain(2)+ny*strain(1))+ C2*tn
		Stress(3)= tn
		Stress(4)= G*Strain(3)
		Stress(5)= ts2
		Stress(6)= ts1
 	DEALLOCATE (Strain)
END IF
!		Transformation of local stresses in global stresses
CALL Stress_Transformation(v1,v2,v3,Stress,Cdim)
RETURN
End SUBROUTINE BStress

SUBROUTINE Stress_Transformation(v1,v2,v3,Stress,Cdim)
!---------------------------------------------------------
!    Transforms stresses from local into global directions
!---------------------------------------------------------
REAL , INTENT(INOUT)		:: Stress(:)					!  Stress vector
REAL , INTENT(IN)				 :: v1(:),v2(:),v3(:)	!  Orthogonal vectors of local coordinate system
INTEGER , INTENT (IN) :: Cdim
REAL, ALLOCATABLE				 :: S(:)
IF(Cdim == 2)THEN
	ALLOCATE(S(3))
ELSE
	ALLOCATE(S(6))
END IF
S=0.0
IF(Cdim == 2) THEN
	S(1)=Stress(1)
	S(2)=Stress(2)
	S(3)=Stress(4)
	Stress=0.0
	Stress(1)= S(1)*v1(1)**2 + S(2)*v3(1)**2 + 2*S(3)*v1(1)*v3(1)
	Stress(2)= S(1)*v1(2)**2 + S(2)*v3(2)**2 + 2*S(3)*v1(2)*v3(2)
	Stress(4)= S(1)*v1(1)*v1(2) + S(2)*v3(1)*v3(2) + S(3)*(v1(1)*v3(2)+v1(2)*v3(1))
ELSE IF(Cdim == 3)THEN
	S=Stress
	Stress=0.0
	Stress(1)= S(1)*v1(1)**2 + S(2)*v2(1)**2 + S(3)*v3(1)**2 + 2*S(4)*v1(1)*v2(1) + 2*S(5)*v2(1)*v3(1)+ 2*S(6)*v1(1)*v3(1)
	Stress(2)= S(1)*v1(2)**2 + S(2)*v2(2)**2 + S(3)*v3(2)**2 + 2*S(4)*v1(2)*v2(2) + 2*S(5)*v2(2)*v3(2)+ 2*S(6)*v1(2)*v3(2)
	Stress(3)= S(1)*v1(3)**2 + S(2)*v2(3)**2 + S(3)*v3(3)**2 + 2*S(4)*v1(3)*v2(3) + 2*S(5)*v2(3)*v3(3)+ 2*S(6)*v1(3)*v3(3)
	Stress(4)= v1(1)*(S(1)*v1(2)+S(4)*v2(2)+S(6)*v3(2)) + v2(1)*(S(2)*v2(2)+S(4)*v1(2)+S(5)*v3(2)) + v3(1)*(S(3)*v3(2)&
                  +S(5)*v2(2)+S(6)*v1(2))
	Stress(5)= v1(2)*(S(1)*v1(3)+S(4)*v2(3)+S(6)*v3(3)) + v2(2)*(S(2)*v2(3)+S(4)*v1(3)+S(5)*v3(3)) + v3(2)*(S(3)*v3(3)&
                  +S(5)*v2(3)+S(6)*v1(3))
	Stress(6)= v1(1)*(S(1)*v1(3)+S(4)*v2(3)+S(6)*v3(3)) + v2(1)*(S(2)*v2(3)+S(4)*v1(3)+S(5)*v3(3)) + v3(1)*(S(3)*v3(3)&
                  +S(5)*v2(3)+S(6)*v1(3))
END IF	
DEALLOCATE(S)
END SUBROUTINE Stress_Transformation

SUBROUTINE Flow_Transformation(v1,v2,v3,Flux,Cdim)
!---------------------------------------------------------
!    Transforms stresses from local into global directions
!---------------------------------------------------------
REAL , INTENT(INOUT)		:: Flux(:)					!  Stress vector
REAL , INTENT(IN)				 :: v1(:),v2(:),v3(:)	!  Orthogonal vectors of local coordinate system
INTEGER , INTENT (IN)		::Cdim
REAL, ALLOCATABLE		 :: F(:)
IF(Cdim == 2)THEN
	ALLOCATE(F(2))
ELSE
	ALLOCATE(F(3))
END IF
F=0.0
IF(Cdim == 2) THEN
	F(1)=Flux(1)
	F(2)=Flux(2)
	Flux=0.0
	Flux(1)= F(1)*v1(1) + F(2)*v2(1)
	Flux(2)= F(1)*v1(2) + F(2)*v2(2)
ELSE IF(Cdim == 3)THEN
!	S=Stress
!	Stress=0.0
!	Stress(1)= S(1)*v1(1)**2 + S(2)*v2(1)**2 + S(3)*v3(1)**2 + 2*S(4)*v1(1)*v2(1) + 2*S(5)*v2(1)*v3(1)+ 2*S(6)*v1(1)*v3(1)
!	Stress(2)= S(1)*v1(2)**2 + S(2)*v2(2)**2 + S(3)*v3(2)**2 + 2*S(4)*v1(2)*v2(2) + 2*S(5)*v2(2)*v3(2)+ 2*S(6)*v1(2)*v3(2)
!	Stress(3)= S(1)*v1(3)**2 + S(2)*v2(3)**2 + S(3)*v3(3)**2 + 2*S(4)*v1(3)*v2(3) + 2*S(5)*v2(3)*v3(3)+ 2*S(6)*v1(3)*v3(3)
!	Stress(4)= v1(1)*(S(1)*v1(2)+S(4)*v2(2)+S(6)*v3(2)) + v2(1)*(S(2)*v2(2)+S(4)*v1(2)+S(5)*v3(2)) + v3(1)*(S(3)*v3(2)&
!                 +S(5)*v2(2)+S(6)*v1(2))
!	Stress(5)= v1(2)*(S(1)*v1(3)+S(4)*v2(3)+S(6)*v3(3)) + v2(2)*(S(2)*v2(3)+S(4)*v1(3)+S(5)*v3(3)) + v3(2)*(S(3)*v3(3)&
!                  +S(5)*v2(3)+S(6)*v1(3))
!	Stress(6)= v1(1)*(S(1)*v1(3)+S(4)*v2(3)+S(6)*v3(3)) + v2(1)*(S(2)*v2(3)+S(4)*v1(3)+S(5)*v3(3)) + v3(1)*(S(3)*v3(3)&
 !                 +S(5)*v2(3)+S(6)*v1(3))
END IF	
DEALLOCATE(F)
END SUBROUTINE Flow_Transformation

END MODULE Postproc_lib
