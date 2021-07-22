PROGRAM Post_processor
!------------------------------------------------
!     General purpose Postprocessor
!     for computing results at boundary and interior points
!     for elasticity and potential problems
!------------------------------------------------------
USE Utility_lib;USE Elast_lib;USE Laplace_lib;USE Integration_lib
USE Postproc_lib
IMPLICIT NONE
INTEGER(4) :: response
INTEGER, ALLOCATABLE :: Inci(:)                 !  Incidences (one elem.)
INTEGER, ALLOCATABLE :: Incie(:,:)              !  Incidences (all elem.)
INTEGER, ALLOCATABLE :: Ldest(:)                !  Destinations (one elem.)
REAL, ALLOCATABLE    :: Elcor(:,:)              !  Element coordinates 
REAL, ALLOCATABLE    :: El_u(:,:,:)             !
REAL, ALLOCATABLE    :: El_t(:,:,:)             !  Results of System
REAL, ALLOCATABLE    :: El_ue(:,:)              !  Diplacement results of Element
REAL, ALLOCATABLE    :: El_te(:,:)              !  Traction results of Element
REAL, ALLOCATABLE    :: Disp(:)                 !  Diplacement results of Node
REAL, ALLOCATABLE    :: Trac(:)                 !  Traction results of Node
REAL, ALLOCATABLE    :: El_trac(:)              !  Traction results of Element
REAL, ALLOCATABLE    :: El_disp(:)              !  Displacement results of Element
REAL, ALLOCATABLE    :: xP(:,:)                 !  Node co-ordinates of BE
REAL, ALLOCATABLE    :: xPnt(:)                 !  Co-ordinates of int. point
REAL, ALLOCATABLE    :: Ni(:),GCcor(:),dxr(:),Vnorm(:)  
CHARACTER (LEN=80)      :: Title 
REAL :: Elengx,Elenge,Rmin,Glcorx(8),Wix(8),Glcore(8),Wie(8)
REAL :: Jac,Step
REAL :: Xsi1,Xsi2,Eta1,Eta2,RJacB,RonL
REAL, ALLOCATABLE :: Flow(:),Stress(:)!  Results for bound.Point
REAL, ALLOCATABLE :: uPnt(:),SPnt(:) !  Results for int Point
REAL, ALLOCATABLE :: TU(:,:),UU(:,:) !  Kernels for u
REAL, ALLOCATABLE  :: TS(:,:),US(:,:) !  Kernels for q,s
REAL, ALLOCATABLE  :: Fac(:),Fac_nod(:,:)  !   Factors for symmetry 
INTEGER :: Cdim,Node,M,N,Istat,Nodel,Nel,Ndof,Cod,Nreg
INTEGER :: Ltyp,Nodes,Maxe,Ndofe,Ndofs,Ncol,ndg,ldim
INTEGER :: nod,nd,Nstres,Nsym,Isym,nsy,IPS,Nan,Nen,Ios,dofa,dofe
INTEGER :: Mi,Ki,K,I,NDIVX,NDIVSX,NDIVE,NDIVSE,MAXDIVS
INTEGER :: Nxs,Net,NXSI,NETA,j
REAL    :: Con,E,ny,Fact,G,C2,C3,C5,C6,C7
REAL    :: xsi,eta,Weit,R,Rlim(2)
REAL		:: epsi= 1.0E-10           !   Small value for comparing coords

!   Read job information
OPEN (UNIT=1,FILE='INPUT',FORM='FORMATTED')
OPEN (UNIT=2,FILE='OUTPUT',FORM='FORMATTED')
Call Jobin(Title,Cdim,Ndof,IPS,Nreg,Ltyp,Con,E,ny,&
           Isym,nodel,nodes,maxe)
Ndofe= nodel*ndof
ldim= Cdim-1
Nsym= 2**Isym   !   number of symmetry loops
ALLOCATE(xP(Cdim,Nodes))   !  Array for node coordinates
ALLOCATE(Incie(Maxe,Nodel),Inci(Nodel),Ldest(Ndofe)) !  Array for incidences
ALLOCATE(Ni(Nodel),GCcor(Cdim),dxr(Cdim),Vnorm(Cdim))
CALL Geomin(Nodes,Maxe,xp,Incie,Nodel,Cdim)
!   Compute constants
IF(Ndof == 1) THEN
 Nstres= Cdim
ELSE
 G= E/(2.0*(1+ny))
 C2= 1/(8*Pi*(1-ny))
 C3= 1.0-2.0*ny
 C5= G/(4.0*Pi*(1-ny))
 C6= 15
 C7= 1.0-4.0*ny
 Nstres= 6
 IF(Cdim == 2) THEN
  IF(IPS == 1) THEN                                                  ! Plane Strain
     C2= 1/(4*Pi*(1-ny))
     C5= G/(2.0*Pi*(1-ny))
     C6= 8
     Nstres= 4
  ELSE
     C2= (1+ny)/(4*Pi                )                               ! Plane Stress
     C3= (1.0-ny)/(1.0+ny)
     C5= (1.0+ny)*G/(2.0*Pi)
     C6= 8
     C7= (1.0-3.0*ny)/(1.0+ny)
     Nstres= 4 
  END IF 
 END IF 
END IF      
ALLOCATE(El_u(Maxe,Nodel,ndof),El_t(Maxe,Nodel,ndof),El_te(Nodel,ndof),El_ue(Nodel,ndof),Fac_nod(Nodel,ndof))
ALLOCATE(El_trac(Ndofe),El_disp(Ndofe))
CLOSE(UNIT=1)
OPEN (UNIT=1,FILE='BERESULTS',FORM='FORMATTED')
WRITE(2,*) ' '
WRITE(2,*) 'Post-processed Results'
WRITE(2,*) ' '
Elements1:&
DO Nel=1,Maxe
                READ(1,*) ((El_u(nel,n,m),m=1,ndof),n=1,Nodel)
                READ(1,*) ((El_t(nel,n,m),m=1,ndof),n=1,Nodel)
END DO &
Elements1
ALLOCATE(Elcor(Cdim,Nodel))
CLOSE(UNIT=1)
OPEN (UNIT=1,FILE='INPUT2',FORM='FORMATTED')
ALLOCATE(Flow(Cdim),Stress(Nstres))
!------------------------------------------------------------
!     Computation of boundary fluxes/stresses
!-------------------------------------------------------------
WRITE(2,*) 'Results at Boundary Elements:'
READ(1,*) Nan,Nen
IF(Nan > 0) THEN
        IF(LTYP == 1) THEN    ! linear elements
          NXSI=2
          NETA=2
          Step=2.0
        ELSE                  !  quadratic elements
          NXSI=3
          NETA=3
          Step= 1.0
        ENDIF
        IF(LDIM == 1) NETA=1    ! 1-D element
        Element_loop: &
        DO NEL= Nan,Nen
         Inci= Incie(nel,:)
         Elcor= xp(:,Inci(:))
         Eta= -1.0
         Eta_loop: &
         DO Net=1,NETA          
           Xsi= -1.0
           Xsi_loop: &
           DO Nxs=1,NXSI
                IF(Ndof == 1) THEN
                  Flow= 0.0
				  IF( Cdim == 2 .AND. Xsi == -1.0  )THEN
					Flow(2)= El_t(NEL, 1 , 1)
				  ELSE IF( Cdim == 2 .AND. Xsi == 0.0  )THEN
					Flow(2)= El_t(NEL, 3 , 1)
				  ELSE IF( Cdim == 2 .AND. Xsi == 1.0  )THEN
					Flow(2)= El_t(NEL, 2 , 1)
				  END IF		
                  Call BFLOW(Flow,xsi,eta,El_u(Nel,:,:),Inci,Elcor,Con)
                  WRITE(2,'(A,I5,A,F6.2,A,F6.2)') 'Element',Nel,' xsi=',xsi,' eta=',eta
                  WRITE(2,'(A,2F10.3)') 'Flux: ',Flow
                ELSE
                  Stress= 0.0
                  Call BStress(Stress,xsi,eta,El_u(Nel,:,:),El_t(Nel,:,:),Inci,Elcor,E,ny,IPS)    
				  IF( Cdim ==2)THEN
					WRITE(2,'(A,I5,A,F6.2)') 'Element',Nel,' xsi=',xsi
 				  ELSE IF( Cdim ==3)THEN
					WRITE(2,'(A,I5,A,F6.2,A,F6.2)') 'Element',Nel,' xsi=',xsi,' eta=',eta
				  END IF	
					WRITE(2,'(A,6F9.2)') ' Stress: ',Stress               
				END IF
                Xsi= Xsi + Step
            END DO Xsi_loop
            Eta=Eta + Step
           END DO Eta_loop
        END DO &
        Element_loop
END IF
ALLOCATE(uPnt(NDOF),SPnt(NSTRES),UU(NDOF,NDOF),TU(NDOF,NDOF))
ALLOCATE(TS(Nstres,Ndof),US(Nstres,Ndof))
ALLOCATE(xPnt(Cdim),Fac(Ndofe))
ALLOCATE(Disp(NDOF),Trac(NDOF))
WRITE(2,*)''
WRITE(2,*) 'Internal Results:'
WRITE(2,*)''
Internal_points: &
DO
        READ(1,*,IOSTAT=IOS) xPnt
        IF(IOS /= 0) EXIT
        Write(2,'(A,3F10.2)') 'Coordinates: ',xPnt
        !------------------------------------------------------------------------
        !    Computation of Temperatures/Displacements at Points inside a region
        !------------------------------------------------------------------------
        uPnt= 0.0
        Element_loop1: &
        DO NEL= 1,MAXE
                Symmetry_loop1:&
                DO nsy=1,Nsym
                        Inci= Incie(nel,:)
                        Elcor= xp(:,Inci(:))
                        IF(ldim == 2) THEN
                          ELengx= Dist((Elcor(:,3)+Elcor(:,2))/2.,(Elcor(:,4)+Elcor(:,1))/2.,Cdim)  ! Lxsi
                          ELenge= Dist((Elcor(:,2)+Elcor(:,1))/2.,(Elcor(:,3)+Elcor(:,4))/2.,Cdim)  ! Leta
                        ELSE
                          Call Elength(Elengx,Elcor,nodel,ldim)
                        END IF
                        Ldest= 1
                        Fac= 1.0
                        Fac_nod=1.0
                        El_ue(:,:)=El_u(Nel,:,:)
                        El_te(:,:)=El_t(Nel,:,:)
                        IF(Isym > 0) THEN
                                DO Nod=1,Nodel
                                        dofa= (nod-1)*Ndof+1
                                        dofe= dofa+Ndof-1
                                        El_trac(dofa:dofe)= El_te(Nod,:)
                                        EL_disp(dofa:dofe)= El_ue(Nod,:)
                                END DO
                                CALL Mirror(Isym,nsy,Nodes,Elcor,Fac,Inci,Ldest,El_trac,EL_disp & 
           ,nodel,ndof,Cdim) 
                                DO Nod=1,Nodel
                                        dofa= (nod-1)*Ndof+1
                                        dofe= dofa+Ndof-1                                       
                                        El_te(Nod,:)= El_trac(dofa:dofe)
                                        El_ue(Nod,:)= El_disp(dofa:dofe)
                                        Fac_nod(Nod,:)= Fac(dofa:dofe)
                                END DO
                        END IF
                        Rmin= Min_dist(Elcor,xPnt,Nodel,ldim,Inci) !  Distance to Pa 
                        Mi= Ngaus(Rmin/Elengx,Cdim-1,Rlim)   ! Gauss Pts. in x dir. 
                        NDIVSX= 1 ; NDIVSE= 1
                        RJacB=1.0
                        RonL= Rmin/Elengx
                        IF(Mi == 5) THEN
                         IF(RonL > epsi) NDIVSX= INT(RLim(2)/RonL) + 1   ! Number of subdvisions
                         IF(NDIVSX > MAXDIVS) MAXDIVS= NDIVSX
                         Mi=4
                        END IF
!         write(2,*) i,Rmin,RonL,NDIVSX
                        CALL Gauss_coor(Glcorx,Wix,Mi)  ! Coords/Wghts x dir
                        Ki= 1 ; Wie(:)= 1.0 ; Glcore(:)= 0.0
                        IF(Cdim == 3) THEN
                         Ki= Ngaus(Rmin/Elenge,Cdim-1,Rlim)  ! Gauss Pts. in h dir. 
                         RonL= Rmin/Elenge 
                         IF(Ki == 5) THEN
                          IF(RonL > epsi) NDIVSE= INT(RLim(2)/RonL) + 1   ! Number of subdvisions
                          IF(NDIVSE > MAXDIVS) MAXDIVS= NDIVSE
                          Ki=4
                         END IF
                        ! write(2,*) i,RonL,NDIVSX
                         CALL Gauss_coor(Glcore,Wie,Ki) ! Coords/Wghts h dir
                        END IF
                        IF(NDIVSX > 1 .OR. NDIVSE>1) RJacB= 1.0/(NDIVSX*NDIVSE)
!                        write(2,*) RJacB
                        Xsi1=-1.0                        
Subdivisions_xsi: DO NDIVX=1,NDIVSX
          Xsi2= Xsi1 + 2.0/NDIVSX
          Eta1=-1.0
Subdivisions_eta: DO NDIVE=1,NDIVSE
          Eta2= Eta1 + 2.0/NDIVSE
                        Gauss_points_xsi: &
                        DO m=1,Mi
                                xsi= Glcorx(m)
                                IF(NDIVSX > 1) xsi= 0.5*(Xsi1+Xsi2)+xsi/NDIVSX
                             !   write(2,*) Ndive,m,xsi
                                Gauss_points_eta: &
                                DO k=1,Ki                                        
                                        eta= Glcore(k)
                                        IF(NDIVSE > 1) eta= 0.5*(Eta1+Eta2)+eta/NDIVSE
                                        Weit= Wix(m)*Wie(k)
                                        CALL Serendip_func(Ni,xsi,eta,ldim,nodel,Inci) 
                                        CALL Normal_Jac(Vnorm,Jac,xsi,eta,ldim,nodel,Inci,elcor) 
                                        Fact= Weit*Jac*RJacB
                                        CALL Cartesian(GCcor,Ni,ldim,elcor)  ! x,y,z of Gauss pnt
                                        r= Dist(GCcor,xPnt,Cdim)      !  Dist. P,Q
                                        dxr= (GCcor-xPnt)/r         !  rx/r , ry/r  etc
                                        IF(Ndof .EQ. 1) THEN
                                                TU= U(r,Con,Cdim)  ; UU= T(r,dxr,Vnorm,Cdim)
                                        ELSE
                                                TU= UK(dxr,r,E,ny,Cdim) ; UU= TK(dxr,r,Vnorm,ny,Cdim)   
                                        END IF
                                        Node_loop1:&
                                        DO Node=1,Nodel
                                                Disp= El_ue(Node,:)* Fac_nod(Node,:)
                                                Trac= El_te(Node,:)* Fac_nod(Node,:)
                                                uPnt= uPnt + (MATMUL(TU,Trac)-&
                                                MATMUL(UU,Disp))* Ni(Node)* Fact
                                        END DO &
                                        Node_loop1
                                END DO &
                                Gauss_points_eta
                        END DO &
                        Gauss_points_xsi
                        Eta1= Eta2
                       END DO Subdivisions_eta
                       Xsi1= Xsi2
                   END DO Subdivisions_xsi
                END DO &
                Symmetry_loop1
        END DO &
        Element_loop1
        WRITE(2,'(A,3E10.3)') '          u: ',uPnt
        !------------------------------------------------------------
        !    Computation of Fluxes/Stresses at Points inside a region
        !------------------------------------------------------------
        SPnt= 0.0
        Element_loop2: &
        DO NEL= 1,MAXE
                Symmetry_loop2: &
                DO nsy=1,Nsym
                        Inci= Incie(nel,:)
                        Elcor= xp(:,Inci(:))
                        IF(ldim == 2) THEN
                                ELengx= Dist((Elcor(:,3)+Elcor(:,2))/2.,(Elcor(:,4)+Elcor(:,1))/2.,Cdim)  ! Lxsi
                                ELenge= Dist((Elcor(:,2)+Elcor(:,1))/2.,(Elcor(:,3)+Elcor(:,4))/2.,Cdim)  ! Leta
                        ELSE
                          Call Elength(Elengx,Elcor,nodel,ldim)
                        END IF
                        Ldest= 1
                        Fac= 1.0
                        El_ue(:,:)=El_u(Nel,:,:)
                        El_te(:,:)=El_t(Nel,:,:)
                        IF(Isym > 0) THEN
                                DO Nod=1,Nodel
                                        dofa= (nod-1)*Ndof+1
                                        dofe= dofa+Ndof-1
                                        El_trac(dofa:dofe)= El_te(Nod,:)
                                        EL_disp(dofa:dofe)= El_ue(Nod,:)
                                END DO
                                CALL Mirror(Isym,nsy,Nodes,Elcor,Fac,Inci,Ldest,El_trac,El_disp& 
           ,nodel,ndof,Cdim) 
                                DO Nod=1,Nodel
                                        dofa= (nod-1)*Ndof+1
                                        dofe= dofa+Ndof-1                                       
                                        El_te(Nod,:)= El_trac(dofa:dofe)
                                        El_ue(Nod,:)= El_disp(dofa:dofe)
                                        Fac_nod(Nod,:)= Fac(dofa:dofe)
                                END DO
                        END IF
                        Rmin= Min_dist(Elcor,xPnt,Nodel,ldim,Inci) !  Distance to Pa 
                        Mi= Ngaus(Rmin/Elengx,Cdim,Rlim)     ! Gauss Pts. in x dir.
                        NDIVSX= 1 ; NDIVSE= 1
                        RJacB=1.0
                        RonL= Rmin/Elengx
                        IF(Mi == 5) THEN
                         IF(RonL > epsi) NDIVSX= INT(RLim(2)/RonL) + 1   ! Number of subdvisions
                         IF(NDIVSX > MAXDIVS) MAXDIVS= NDIVSX
                         Mi=4
                        END IF
!         write(2,*) Nel,Nsy,Rmin,RonL,NDIVSX 
!                       Mi = 8
                        CALL Gauss_coor(Glcorx,Wix,Mi)  ! Coords/Wghts x dir
                        Ki= 1 ; Wie(:)= 1.0 ; Glcore(:)= 0.0
                        IF(Cdim == 3) THEN
                         Ki= Ngaus(Rmin/Elenge,Cdim,Rlim)    ! Gauss Pts. in h dir. 
                        ! Ki= Ngaus(Rmin/Elenge,Cdim-1,Rlim)  ! Gauss Pts. in h dir. 
                         RonL= Rmin/Elenge 
                         IF(Ki == 5) THEN
                          IF(RonL > epsi) NDIVSE= INT(RLim(2)/RonL) + 1   ! Number of subdvisions
                          IF(NDIVSE > MAXDIVS) MAXDIVS= NDIVSE
                          Ki=4
                         END IF                         
 !                        write(2,*) Rmin,RonL,NDIVSE
                         CALL Gauss_coor(Glcore,Wie,Ki) ! Coords/Wghts h dir
                        END IF
                        IF(NDIVSX > 1 .OR. NDIVSE>1) RJacB= 1.0/(NDIVSX*NDIVSE)
                        Xsi1=-1.0                        
Subdivisions_xsi1: DO NDIVX=1,NDIVSX
						Xsi2= Xsi1 + 2.0/NDIVSX
          Eta1=-1.0
Subdivisions_eta1: DO NDIVE=1,NDIVSE
						Eta2= Eta1 + 2.0/NDIVSE
                        Gauss_points_xsi2: &
                        DO m=1,Mi
                                xsi= Glcorx(m)
                                IF(NDIVSX > 1) xsi= 0.5*(Xsi1+Xsi2)+xsi/NDIVSX
                                Gauss_points_eta2: &
                                DO k=1,Ki
                                        eta= Glcore(k)
                                        IF(NDIVSE > 1) eta= 0.5*(Eta1+Eta2)+eta/NDIVSE
                                        Weit= Wix(m)*Wie(k)
                                        CALL Serendip_func(Ni,xsi,eta,ldim,nodel,Inci) 
                                        CALL Normal_Jac(Vnorm,Jac,xsi,eta,ldim,nodel,Inci,elcor) 
                                        Fact= Weit*Jac*RJacB
                                        CALL Cartesian(GCcor,Ni,ldim,elcor)  ! x,y,z of Gauss pnt
                                        r= Dist(GCcor,xPnt,Cdim)      !  Dist. P,Q
                                        dxr= (GCcor-xPnt)/r         !  rx/r , ry/r  etc
                                        IF(Ndof .EQ. 1) THEN
                                         TS(:,1)= dU(r,dxr,Cdim); US(:,1)= dT(r,dxr,Vnorm,Cdim)
                                        ELSE
                                         CALL SK(TS,DXR,R,C2,C3) 
                                         CALL RK(US,DXR,R,VNORM,C3,C5,c6,C7,ny)
                                        END IF
                                        Node_loop2:&
                                        DO Node=1,Nodel
                                                Disp= El_ue(Node,:)* Fac_nod(Node,:)
                                                Trac= El_te(Node,:)* Fac_nod(Node,:)
                                                SPnt= SPnt + (MATMUL(TS,Trac)- &
                                                MATMUL(US,Disp))* Ni(Node)* Fact 
                                        END DO &
                                        Node_loop2
                                END DO &
                                Gauss_points_eta2
                        END DO &
                        Gauss_points_xsi2
                        Eta1= Eta2
                       END DO Subdivisions_eta1
                       Xsi1= Xsi2
                   END DO Subdivisions_xsi1
                END DO &
                Symmetry_loop2
        END DO &
        Element_loop2 
        IF(Ndof == 1) THEN
        WRITE(2,'(A,6F10.3)') '       Flux: ',SPnt
        ELSE
           IF(CDIM == 2) THEN
            IF(IPS==1) THEN
             SPnt(4)=SPnt(3)
             SPnt(3)= ny*(SPnt(1)+SPnt(2))
            ELSE
             SPnt(4)=SPnt(3)
             SPnt(3)= 0
            END IF
           END IF
                WRITE(2,'(A,6F10.3)') '     Stress: ',SPnt
        END IF
END DO &
Internal_points

END PROGRAM Post_processor



