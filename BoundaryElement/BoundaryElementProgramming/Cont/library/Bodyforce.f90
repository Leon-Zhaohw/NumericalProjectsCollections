MODULE BodyForce_lib
USE Utility_Lib; USE Elast_lib; USE Integration_lib
!    Library for Body_Force problems
IMPLICIT NONE
CONTAINS

SUBROUTINE Body_force(F,CDim,xP,NCol,Isym,E,ny)
!--------------------------------------------------------
!  Adds contribution of body force terms (initial strain)
!  to the right hand side vector F
!  This implementation is only for linear cells and plane 
!  problems
!--------------------------------------------------------

IMPLICIT NONE

INTEGER , INTENT(IN)        :: CDim
REAL , INTENT(IN)           :: E
REAL , INTENT(IN)           :: ny
INTEGER , INTENT(IN)        :: NCol   !   Number of collocation points
INTEGER , INTENT(IN)        :: Isym
REAL , INTENT (IN)          :: xP(Cdim,Ncol)  !  Collocation point coords
REAL(KIND=8), INTENT(INOUT) :: F(CDim*Ncol)   !  right hand side vector
!
INTEGER, ALLOCATABLE :: InciC(:,:)  !  Cell Incidences
INTEGER, ALLOCATABLE :: Inci(:)
REAL, ALLOCATABLE    :: xPC(:,:)    !  Cell Node co-ordinates
REAL, ALLOCATABLE    :: Ni(:),Elcor(:,:)
!
INTEGER   :: NodesC,Ncells,NodelC,ldimC,IOS
INTEGER   :: m,n,k,Node,Nc,i,ii,jj,Mi,Ki,iD
REAL      :: Eps0(2),SigK(Cdim,Cdim),GCcor(3)
REAL      :: Glcorx(8),Glcore(8),Wix(8),Wie(8),Vnorm(3)
REAL      :: Jac,Weit,xsi,eta,r,dxr(Cdim)
IF(Cdim > 2) RETURN  ! This coding is for plane problems only
IF(ISym > 0)  RETURN  ! Symmetry not considered
READ(1,*,IOSTAT=IOS) NodesC
IF(IOS /= 0) RETURN    !   No body force effects
Write(2,*) 'Number of cell nodes=',NodesC
READ(1,*) Ncells
Write(2,*) 'Number of cells=',Ncells
READ(1,*) Eps0
Write(2,*) 'Eps0=',Eps0
NodelC=4       !  only linear elements considered
ldimC= 2       !  plane cells
ALLOCATE(xPC(Cdim,NodesC))   !  Array for node coordinates
ALLOCATE(InciC(Ncells,NodelC),Inci(NodelC)) !  Array for incidences
ALLOCATE(Ni(nodelC),ELCOR(Cdim+1,nodelC))
!-------------------------------------------------------
!		Read Cell Node Co-ordinates 
!-------------------------------------------------------
DO Node=1,NodesC
 READ(1,*) (xPC(M,Node),M=1,Cdim)
 WRITE(2,'(A5,I5,A8,3F8.2)') 'Node ',Node,&
         '  Coor  ',(xPC(M,Node),M=1,Cdim)
END DO
!-------------------------------------------------------
!		Read Cell Incidences 
!-------------------------------------------------------
WRITE(2,*)''
WRITE(2,*)'Incidences: '
WRITE(2,*)''
DO Nc=1,Ncells
  READ(1,*) (InciC(Nc,n),n=1,NodelC)
  WRITE(2,'(A3,I5,A8,4I5)')'EL ',Nc,'  Inci  ',InciC(Nc,:)
END DO
!--------------------------------------------------------
!   compute contribution to right hand side
!------------------------------------------------------
Colloc_points: DO i=1,Ncol
  Cells: DO nc=1,Ncells
         Mi=4 ;  Ki=4  !  Assumes that cell is sufficiently far away from Pi
         Call Gauss_coor(Glcorx,Wix,Mi)     !  Assign coords/Weights xsi-direction        
         Call Gauss_coor(Glcore,Wie,Ki)     !  Assign coords/Weights eta-direction
         Elcor(1:2,:)= xPC(1:2,InciC(nc,:))
         Elcor(3,:)= 0.0 ! Need to do this because we are using 2-D BE element for cell
         Inci=InciC(nc,:)
     Gauss_points_xsi: DO m=1,Mi
          xsi= Glcorx(m)
       Gauss_points_eta: DO k=1,Ki
           eta= Glcore(k)
           CALL Serendip_func(Ni,xsi,eta,ldimC,nodelC,Inci)  !   Shape function value
           Call Normal_Jac(Vnorm,Jac,xsi,eta,ldimC,nodelC,Inci,elcor) ! Jacobian and normal
           Weit= Wix(m)*Wie(k)*Jac
           CALL Cartesian(GCcor,Ni,ldimC,elcor)      ! Cart. coords of Gauss pt
           r= Dist(GCcor(1:2),xP(:,i),Cdim)             !  Dist. P,Q
           dxr= (GCcor(1:2)-xP(1:2,i))/r             !  rx/r , ry/r
           SigK= SigmaK(dxr,r,E,ny,Cdim)
         Direction_P: DO ii=1,Cdim
            iD= Cdim*(i-1) + ii     !  Position in F array
            Direction_Q: DO jj=1,Cdim
               Node_points: DO n=1,NodelC
                 F(iD)= F(iD) + Ni(n)*SigK(ii,jj)*Weit*Eps0(jj)
               END DO Node_points
           END DO Direction_Q
         END DO Direction_P
       END DO Gauss_points_eta
      END DO Gauss_points_xsi
     END DO Cells
  END DO Colloc_points
DEALLOCATE(xPC)   !  Array for node coordinates
DEALLOCATE(InciC) !  Array for incidences
DEALLOCATE(Ni,Elcor)
RETURN
END SUbroutine Body_force


END MODULE BodyForce_lib
                