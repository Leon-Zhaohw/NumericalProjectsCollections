PROGRAM General_purpose_MRBEM
!------------------------------------------------
!     General purpose BEM program
!     for solving elasticity and potential problems
!     with multiple regions
!------------------------------------------------------
!USE DFLIB;
USE Utility_lib; USE Elast_lib; USE Laplace_lib 
USE Integration_lib; USE Stiffness_lib
IMPLICIT NONE
INTEGER(4)								:: response
INTEGER, ALLOCATABLE			:: NCode(:,:)				! Element BC´s
INTEGER, ALLOCATABLE			:: Ldest_KBE(:)			! Interface destination vector for region assembly
INTEGER, ALLOCATABLE			:: TypeR(:)					! Type of BE-regions (1 == finite, 2 == Infinite)
REAL, ALLOCATABLE					:: Elcor(:,:)				! Element coordinates
REAL, ALLOCATABLE					:: xP(:,:)					! Node co-ordinates
REAL, ALLOCATABLE					:: xP_discont(:,:)					! Node co-ordinates
REAL, ALLOCATABLE					:: Elres_u(:,:)			!	Element results
REAL, ALLOCATABLE					:: Elres_t(:,:)			!	Element results
REAL, ALLOCATABLE					:: Elres_ueNew(:)			!	Element results
REAL, ALLOCATABLE					:: Elres_teNew(:)			!	Element results
REAL(KIND=8), ALLOCATABLE :: KBE(:,:,:)				! Region stiffness 
REAL(KIND=8), ALLOCATABLE :: A(:,:,:)					! Results due to ui=1 
REAL(KIND=8), ALLOCATABLE :: Lhs(:,:),Rhs(:)	! global matrices
REAL(KIND=8), ALLOCATABLE :: uc(:)						! interface unknown
REAL(KIND=8), ALLOCATABLE :: ucr(:)						! interface unknown (region)
REAL(KIND=8), ALLOCATABLE :: tc(:)						! interface tractions
REAL(KIND=8), ALLOCATABLE :: xf(:)						! free unknown
REAL(KIND=8), ALLOCATABLE :: tcxf(:)					!	unknowns of region
REAL, ALLOCATABLE					:: XpR(:,:)					! Region node coordinates
REAL, ALLOCATABLE 				:: ConR(:)					! Conductivity of regions
REAL, ALLOCATABLE					:: ER(:)						! Youngsmodulus of regions
REAL, ALLOCATABLE					:: nyR(:)						! Poissons ratio of regions
REAL											:: E,ny,Con					! 
INTEGER,ALLOCATABLE				:: InciR(:,:)				! Global node numbers / region / local sequence
INTEGER,ALLOCATABLE				:: Incie(:,:)				! Element incidences (global) of system
INTEGER,ALLOCATABLE				:: Incie_Discont(:,:)				! Element incidences (global) of system of discontinuous elements
INTEGER, ALLOCATABLE			:: IncieR(:,:)			!	Element incidences (local) of region
INTEGER,ALLOCATABLE				:: ListC(:)					! List of interface nodes
INTEGER,ALLOCATABLE				:: ListEC(:,:)			! List of interface Elements / region
INTEGER,ALLOCATABLE				:: ListEF(:,:)			! List of free Elements / region
INTEGER,ALLOCATABLE				:: LdestR(:,:)			! Global D.o.F. numbers / region / local sequence
INTEGER, ALLOCATABLE			:: Nbel(:)					! Number of Boundary Elements each region
INTEGER,ALLOCATABLE				:: NbelC(:)					! Number of Interfaceelements / region
INTEGER,ALLOCATABLE				:: NbelF(:)					!	Number of free elements / region
INTEGER, ALLOCATABLE			:: Bcode(:,:)				!	Boundary code for all elements
INTEGER, ALLOCATABLE			:: Ldeste(:,:)			!	Destination Vector global of System		
INTEGER, ALLOCATABLE			:: LdesteR(:,:)			!	Destination Vector local of all regions	
INTEGER, ALLOCATABLE			:: NodeR(:)					! No. of nodes of Region
INTEGER, ALLOCATABLE			:: NodeC(:)					! No. of nodes of Interface / region
INTEGER, ALLOCATABLE			:: ListR(:,:)				! List of Elementnumbers each region
INTEGER, ALLOCATABLE			:: Ndest(:,:)
INTEGER										:: Cdim							! Cartesian dimension
INTEGER										:: Nodes						! No. of nodes of System
INTEGER										:: Nodel						! No. of nodes per element
INTEGER										:: Ndofe						! D.o.F´s of Element
INTEGER										:: Ndof							! No. of degeres of freedom per node
INTEGER										:: Ndofs						! D.o.F´s of System
INTEGER										:: NdofR						! Number of D.o.F. of region
INTEGER										:: NdofC						! Number of interface D.o.F. of region
INTEGER										:: NdofF						! Number D.o.F. of free nodes of region
INTEGER										:: NodeF						! Number of free Nodes of region
INTEGER										:: NodesC						! Total number of interface nodes of System
INTEGER										:: NdofsC						! Total number of interface D.o.F. of System
INTEGER										:: Toa							! Type of analysis (plane strain = 1, plane stress = 2)
INTEGER										:: Nregs						! Number of regions
INTEGER										:: Ltyp							! Element type(linear = 1, quadratic = 2)
INTEGER										:: Isym							! Symmetry code
INTEGER										:: Maxe							! Number of Elements of System
INTEGER										:: nr,nb,ne,ne1,nel
INTEGER										:: n,node,is,nc,no,ro,co
INTEGER										:: k,m,nd,nrow,ncln,DoF_KBE,DoF
INTEGER										::	NodeNr, found, nod, nod1, ldim
CHARACTER(LEN=80)					:: Title

!-----------------------------------------------------
!   Read job information
!-----------------------------------------------------
OPEN (UNIT=1,FILE='INPUT',FORM='FORMATTED') !  Input
OPEN (UNIT=2,FILE='OUTPUT',FORM='FORMATTED')!  Output
Call JobinMR(Title,Cdim,Ndof,Toa,Ltyp,Isym,nodel,nodes,maxe)
Ndofs= Nodes * Ndof							!	D.O.F's of System
Ndofe= Nodel * Ndof							!	D.O.F's of Element
Isym= 0  !   no symmetry considered here
ALLOCATE(Ndest(Nodes,Ndof))
Ndest= 0
READ(1,*)Nregs
ALLOCATE(TypeR(Nregs),Nbel(Nregs),ListR(Nregs,Maxe))
IF(Ndof == 1)THEN
	ALLOCATE(ConR(Nregs))
ELSE
	ALLOCATE(ER(Nregs),nyR(Nregs))
END IF
CALL Reg_Info(Nregs,ToA,Ndof,TypeR,ConR,ER,nyR,Nbel,ListR)
ALLOCATE(xP(Cdim,Nodes))  !  Array for node coordinates
ALLOCATE(Incie(Maxe,Nodel)) !  Array for incidences
CALL Geomin(Nodes,Maxe,xp,Incie,Nodel,Cdim)
ALLOCATE(BCode(Maxe,Ndofe))      
ALLOCATE(Elres_u(Maxe,Ndofe),Elres_t(Maxe,Ndofe))	
CALL BCinput(Elres_u,Elres_t,Bcode,nodel,ndofe,ndof) 

!--------------------------------------------------------------------------------------
! Transform continuous element incidencies into discontinuous
!--------------------------------------------------------------------------------------
ALLOCATE(Incie_Discont(Maxe,Nodel)) !  Array for incidences
Incie_Discont=0
NodeNr=1
DO ne=1,Maxe
	found=0
	IF(Incie_Discont(ne,1) /= 0)CYCLE
	DO ne1=ne+1,Maxe
		IF(Match(Incie(ne1,:),Incie(ne,:))) THEN
			found=1
			DO nod=1, Nodel
				DO nod1=1, Nodel
					IF(Incie(ne,nod) == Incie(ne1,nod1))THEN
						Incie_Discont(ne,nod)= NodeNr
						Incie_Discont(ne1,nod1)= NodeNr
						NodeNr= NodeNr + 1
					END IF
				END DO	
			END DO	
		END IF	
	END DO
	IF(found == 0)THEN
		DO nod=1, nodel
			Incie_Discont(ne,nod)= NodeNr
			NodeNr= NodeNr + 1
		END DO		
	END IF
END DO	

nodes= NodeNr -1
ndofs= nodes * Ndof	
WRITE(2,*)''
WRITE(2,*)'Incidences for discontinuous elements : '
WRITE(2,*)''
DO Nel=1,Maxe
	WRITE(2,'(A3,I5,A16,24I5)')'EL   ',Nel,'  Incidencies    ',Incie_Discont(Nel,:)
END DO 

!---------------------------------------------------------------------------------------------
! Evaluate coordinates of discontinuous nodes -> collocation nodes
!---------------------------------------------------------------------------------------------
ALLOCATE(xP_Discont(Cdim,nodes))  !  Array for discontinuous node coordinates
xP_Discont=0.0
CALL ContToDiscont(Incie_Discont, Incie, xP_Discont, xP, Nodel, Maxe, Cdim)



!------------------------------------------
!     Determine Element destination vector for assembly
!------------------------------------------
ALLOCATE(Ldeste(Maxe,Ndofe))
Elements_of_region2:&
DO Nel=1,Maxe
	k=0
	DO n=1,Nodel
		DO m=1,Ndof		
			k=k+1								
			IF(Ndof > 1) THEN
				Ldeste(Nel,k)= ((Incie_Discont(Nel,n)-1)*Ndof + m)
			ELSE
				Ldeste(Nel,k)= Incie(Nel,n)
			END IF
		END DO 	   
	END DO		   
END DO &
Elements_of_region2
!-------------------------------------------
!    Detect interface elements,
!    assign interface boundary conditions 
!    Determine number of interface nodes
!-------------------------------------------
ALLOCATE(ListC(Nodes))
NodesC=0
ListC=0
Elements_loop: &
DO ne=1,Maxe
	Elements_loop1: &
	DO ne1=ne+1,Maxe
		IF(Match(Incie_Discont(ne1,:),Incie_Discont(ne,:))) THEN
			BCode(ne,:)= 2 ; BCode(ne1,:)= 2   !  assign interface BC	
			Element_nodes: &
			DO n=1,nodel
				Node= Incie_Discont(ne,n)
				is= 0
				Interface_nodes: &
				DO nc=1,NodesC
					IF(Node == ListC(nc)) is= 1
				END DO &
				Interface_nodes
				IF(is == 0) THEN
					NodesC= NodesC + 1
					ListC(NodesC)= Node
				END IF
			END DO &
			Element_nodes
			EXIT
		END IF
	END DO &
	Elements_loop1
END DO &
Elements_loop



NdofsC= NodesC*Ndof
ALLOCATE(InciR(Nregs,Nodes),IncieR(Maxe,Nodel))
ALLOCATE(KBE(Nregs,NdofsC,NdofsC),A(Nregs,Ndofs,Ndofs))
ALLOCATE(Lhs(NdofsC,NdofsC),Rhs(NdofsC),uc(NdofsC),tc(NdofsC))
ALLOCATE(NodeR(Nregs),NodeC(Nregs))
ALLOCATE(ListEC(Nregs,maxe))
ALLOCATE(ListEF(Nregs,maxe))
ALLOCATE(LdesteR(Maxe,Ndofe))	! Elem. destination vector
ALLOCATE(Ldest_KBE(Ndofs))
ALLOCATE(NCode(Nregs,Ndofs))            
ALLOCATE(LdestR(Nregs,Ndofs))
ALLOCATE(NbelC(Nregs))
ALLOCATE(NbelF(Nregs))

LdesteR= 0
Ncode= 0
NbelF= 0
NbelC= 0
!-------------------------------------------
!    Assign local (region) numbering
!    and incidences of BE in local numbering
!--------------------------------------------  
ListEC= 0
ListEF= 0
DoF_KBE= 0
Regions_loop_1: &
DO nr=1,Nregs
	node= 0
	Elements_of_region: &
	DO nb=1,Nbel(nr)
		ne= ListR(nr,nb)
		Interface_elements: &
		IF(Bcode(ne,1) == 2) THEN    
			NbelC(nr)= NbelC(nr) + 1
			ListEC(nr,NbelC(nr))= ne
			Nodes_of_Elem: &
			DO n=1,Nodel
!   check if node has allready been entered
        is=0
				DO no=1,node
					IF(InciR(nr,no) == Incie_Discont(ne,n)) THEN
						is= 1
						EXIT
					END IF
				END DO
				IF(is == 0) THEN
				  node=node+1
					InciR(nr,node)= Incie_Discont(ne,n)
					IncieR(ne,n)= node					  
				ELSE
					IncieR(ne,n)= no
				END IF
			END DO &
			Nodes_of_Elem
		END IF &
		Interface_elements
	END DO &
	Elements_of_region
	NodeC(nr)= Node				! No of interface nodes of Region nr
	NdofC= NodeC(nr)*Ndof  ! D.o.F. at interface of Region nr  
	Elements_of_region1: &
	DO nb=1,Nbel(nr)
		ne= ListR(nr,nb)
		Free_elements: &
		IF(Bcode(ne,1) /= 2) THEN    
			NbelF(nr)= NbelF(nr) + 1
			ListEF(nr,NbelF(nr))= ne
			Nodes_of_Elem1: &
			DO n=1,Nodel
        is=0
				DO no=1,node
					IF(InciR(nr,no) == Incie_Discont(ne,n)) THEN
						is= 1
						EXIT
					END IF
				END DO
				IF(is == 0) THEN
					node=node+1
					InciR(nr,node)= Incie_Discont(ne,n)
					IncieR(ne,n)= node
				ELSE
					IncieR(ne,n)= no
				END IF
			END DO &
			Nodes_of_Elem1
		END IF &
		Free_elements
	END DO &
	Elements_of_region1
	NodeR(nr)= node               !   number of nodes per region
	!------------------------------------------
	!     Determine Local Element destination vector
	!------------------------------------------
	Elements:&
	DO Nel=1,Nbel(nr)
		k=0
		ne= ListR(nr,Nel)
		DO n=1,Nodel
			DO m=1,Ndof		
				k=k+1								
				IF(Ndof > 1) THEN
					LdesteR(ne,k)= ((IncieR(ne,n)-1)*Ndof + m)
				ELSE
					LdesteR(ne,k)= IncieR(ne,n)
				END IF
			END DO 	   
		END DO		   
	END DO &
	Elements
	!------------------------------------------
	!     Determine Local Node destination vector
	!------------------------------------------
	n= 0
	DO no=1, NodeR(nr)
		DO m=1, Ndof
			n= n + 1
			LdestR(nr,n)= (InciR(nr,no)-1) * Ndof + m
		END DO
	END DO


!------------------------------------------
!     Determine global Boundary code vector for assembly
!------------------------------------------
	NdofR= NodeR(nr)*Ndof							! Total degrees of freedom of region
	DoF_o_System: &
	DO	nd=1,NdofR
		DO Nel=1,Nbel(nr)
			ne=ListR(nr,Nel)
			DO m=1,Ndofe
				IF (nd == LdesteR(ne,m) .and. NCode(nr,nd) == 0) THEN 
					NCode(nr,nd)= NCode(nr,nd)+BCode(ne,m)
				END IF
			END DO
		END DO
	END DO &
	DoF_o_System
END DO &
Regions_loop_1

Regions_loop_2: &
DO nr=1,Nregs
!-----------------------------------
!   allocate coordinates in local(region) numbering
!-----------------------------------
	ALLOCATE(XpR(Cdim,NodeR(nr)))
	Region_nodes: &
	DO Node=1,NodeR(nr)
		XpR(:,Node)= Xp_Discont(:,InciR(nr,node))
	END DO &
	Region_nodes
!-----------------------------------------------------------------
!    Determine interface destination vector for region assembly
!-----------------------------------------------------------------
	No_o_Interfaceelements:&
	DO n=1, NbelC(nr)
		ne= ListEC(nr,n)
		DoF_o_Element:&
		DO m=1, Ndofe
			DoF= Ldeste(ne,m)
			IF(Ldest_KBE(DoF) == 0)THEN
				DoF_KBE= DoF_KBE + 1
				Ldest_KBE(DoF)= DoF_KBE
			END IF
		END DO DoF_o_Element
	END DO No_o_Interfaceelements
	NdofR= NodeR(nr)*Ndof							! Total degrees of freedom of region
	NdofC= NodeC(nr)*Ndof							! D.o.F. of interface of Region nr  
	E=ER(nr)
	ny=nyR(nr)
	CALL Stiffness_BEM(nr,XpR,Xp,Nodel,Ndof,Ndofe,NodeR,Ncode(nr,:),NdofR,NdofC,KBE(nr,:,:),A(nr,:,:),tc,Cdim,Elres_u,Elres_t,&
										IncieR,Incie,LdesteR,Nbel,ListR,TypeR,Bcode,Con,E,ny,Ndest,Isym)
	DO ro=1,NdofC
		DoF= LdestR(nr,ro)
		Nrow= Ldest_KBE(DoF)
		Rhs(Nrow)= Rhs(Nrow) + tc(ro)
		DO co=1, NdofC  
			DoF= LdestR(nr,co)
			Ncln= Ldest_KBE(DoF)
			Lhs(Nrow,Ncln)= Lhs(Nrow,Ncln) - KBE(nr,ro,co)
		END DO
	END DO
	DEALLOCATE (XPR)
END DO &
Regions_loop_2
DEALLOCATE(tc)
!------------------------------
!   solve for interface unknown
!------------------------------
CALL Solve(Lhs,Rhs,uc)
!-----------------------------
!  compute and add effect of interface displ.
!-----------------------------
Regions_loop_3: &
DO nr=1,Nregs
!  gather region interface displacements
  NdofC= NodeC(nr)*Ndof
	ALLOCATE(ucr(NdofC))
  Interface_dof: &
	DO n=1,NdofC
		DoF= LdestR(nr,n)
		ucr(n)= uc(Ldest_KBE(DoF))
	END DO &
	Interface_dof
!-------------------------------------------------------------
! Store Interfacedisplacements into Elres_u
!-------------------------------------------------------------
	Interface_DoF1:&
	DO nd=1, NdofC
		DO n=1, Nbel(nr)
			ne=ListR(nr,n)
			DO m=1,Ndofe
				IF(nd == LdesteR(ne,m))THEN
					Elres_u(ne,m)= Elres_u(ne,m) + ucr(nd)
				END IF
			END DO
		END DO
	END DO Interface_DoF1
!   effects of interface displacement in local (region) numbering
	NdofR= NodeR(nr)*Ndof
	NdofF= (NodeR(nr) - NodeC(nr))*Ndof          !   d.o.F , free nodes
	ALLOCATE(tc(NdofC),xf(NdofF),tcxf(NdofR))
	tc= 0.0;	xf= 0.0;	tcxf= 0.0
	tc= Matmul(KBE(nr,1:NdofC,1:NdofC),ucr)
	xf= Matmul(A(nr,1:NdofF,1:NdofC),ucr)
	tcxf(1:NdofC)= tc
	tcxf(NdofC+1:NdofR)= xf
	!-------------------------------------------------------------
	! Store Interfacetractions into Elres_t
	!-------------------------------------------------------------
	DO nd=1, NdofC
		DO n=1, NbelC(nr)
			ne=ListEC(nr,n)
			DO m=1, Ndofe
				IF(nd == LdesteR(ne,m))THEN
					Elres_t(ne,m)= Elres_t(ne,m) + tcxf(nd)
				END IF
			END DO
		END DO
	END DO
	!-------------------------------------------------------------
	! Store Results of free nodes into Elres_u or Elres_t
	!-------------------------------------------------------------
	DO nd=NdofC+1, NdofR
		DO n=1, NbelF(nr)
			ne=ListEF(nr,n)
			DO m=1, Ndofe
				IF(nd == LdesteR(ne,m))THEN
					IF(Ncode(nr,nd) == 0)THEN
						Elres_u(ne,m)= Elres_u(ne,m) + tcxf(nd)
					ELSE IF(Bcode(ne,m) == 1)THEN
						Elres_t(ne,m)= Elres_t(ne,m) + tcxf(nd)
					END IF
				END IF
			END DO
		END DO
	END DO					
!Elements_of_region3: &
!DO nb=1,Nbel(nr)
!	ne= Listr(nr,nb)
!	Element_Dofr: &
!	DO nd=1,Ndofe
!		IF(Ncode(nr,LdesteR(ne,nd)) == 2) THEN
!		 Elres_t(ne,nd)= Elres_t(ne,nd) + tcxf(LdesteR(ne,nd))
!		ELSE IF(Ncode(nr,LdesteR(ne,nd)) == 1) THEN
!		 Elres_t(ne,nd)= Elres_t(ne,nd) + tcxf(LdesteR(ne,nd))
!		ELSE IF(Ncode(nr,LdesteR(ne,nd)) == 0) THEN
!		 Elres_u(ne,nd)= Elres_u(ne,nd) + tcxf(LdesteR(ne,nd))
!		END IF
!	END DO &
!	Element_Dofr
!END DO &
!Elements_of_region3
	DEALLOCATE(tc,xf,tcxf,ucr)
END DO &
Regions_loop_3
!--------------------------
!    Print out results
!--------------------------
CLOSE(UNIT=2)
OPEN(UNIT=2,FILE= 'BERESULTS_AT_DISCONTINUOUS_NODES',FORM='FORMATTED' )
Elements_all:	&
DO nel=1,Maxe	
WRITE(2,*) ' Results, Element ',nel
 WRITE(2,*) 'u=' , (Elres_u(nel,m), m=1,Ndofe)
 WRITE(2,*) 't=' , (Elres_t(nel,m), m=1,Ndofe)
END DO &
Elements_all

CLOSE(UNIT=2)
OPEN (UNIT=2,FILE='BERESULTS_AT_CONTINUOUS_NODES',FORM='FORMATTED')
ldim=Cdim-1
ALLOCATE(Elres_teNew(Ndofe),Elres_ueNew(Ndofe))										

DO nel=1,maxe
!	Inci=Incie(nel,:)													
	
	CALL Transform_DiscToCont(Elres_u(nel,:),Elres_ueNew,Incie(nel,:),Nodel,Cdim)
	CALL Transform_DiscToCont(Elres_t(nel,:),Elres_teNew,Incie(nel,:), Nodel,Cdim)



	WRITE(2,'(24F12.5)') (Elres_ueNew(m), m=1,Ndofe)
	WRITE(2,'(24F12.5)') (Elres_teNew(m), m=1,Ndofe)
END DO



!response = MESSAGEBOXQQ('Program finished'C,'General Purpose MRBEM Program'C, MB$ICONASTERISK)

END PROGRAM General_purpose_MRBEM
