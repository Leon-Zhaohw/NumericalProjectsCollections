!     Last change:  G     3 May 100    9:14 am
PROGRAM Compute_Area
!---------------------------------
!   Program to compute the length/surface area
!   of a line/surface modelled by boundary elements
!---------------------------------
USE Geometry_lib ; USE Utility_lib ; USE Integration_lib; 
IMPLICIT NONE
INTEGER(4) :: response
INTEGER :: ldim,noelem,nelem,lnodes,maxnod,node,Cdim, intord, n, i, j
INTEGER,ALLOCATABLE :: inciG(:,:)  ! Incidences (all elements)
INTEGER,ALLOCATABLE :: inci(:)     ! Incidences one element
REAL,ALLOCATABLE :: corG(:,:)      ! Coordinates (all nodes)
REAL,ALLOCATABLE :: cor(:,:)       ! Coordinates one element
REAL,ALLOCATABLE :: v3(:)          ! Normal vector
REAL             :: Gcor(8),Wi(8)  ! Gauss point coords and weights
REAL             :: Jac,xsi,eta, Area
 OPEN(UNIT=10,FILE='INPUT.DAT',STATUS='OLD',ACTION='READ')
 OPEN(UNIT=11,FILE='OUTPUT.DAT',STATUS='UNKNOWN',ACTION='WRITE')
 READ(10,*) ldim,lnodes,noelem,intord
 WRITE(11,*) ' Element dimension=',ldim
 WRITE(11,*) ' No. of elementnodes=',lnodes
 WRITE(11,*) ' Number of elements=',noelem
 WRITE(11,*) ' Integration order =',intord
 Cdim= ldim+1    !  Cartesian dimension
 ALLOCATE(v3(Cdim))
 ALLOCATE(inciG(lnodes,noelem))  !  Allocate global incidence array
 DO nelem=1,noelem
  READ(10,*) (inciG(n,nelem),n=1,lnodes)
 END DO
 maxnod= MAXVAL(inciG)
 ALLOCATE(corG(Cdim,0:maxnod))   !   Allocate global coordinate array
 corG(:,0)= 0.0                  !   Node # 0 means node is missing (coord. value of 0.0)
 DO node=1,maxnod
  READ(10,*) (corG(i,node),i=1,Cdim)
 END DO
 ALLOCATE(inci(lnodes),cor(Cdim,lnodes))  !   Element incidence and coordinate array
 CALL Gauss_coor(Gcor,Wi,Intord)          !   Get Gauss point coordinates and weights
 Area= 0.0                                !   Start sum for area/length
 Element_loop: DO nelem=1,noelem
                inci=  inciG(:,nelem)     !   Store incidences locally
                cor= corG(:,inci)         !   gather element coordinates
                SELECT CASE (ldim)
                    CASE (1)              !   One-dim. problem determine length
           Gauss_loop:DO I=1,INTORD
                       xsi= Gcor(i)
                       CALL Normal_Jac(v3,Jac,xsi,eta,ldim,lnodes,inci,cor)
                       Area= Area + Jac*Wi(i)
                      END DO Gauss_loop
                    CASE (2)              !   Two-dim. problem determine area
          Gauss_loop1:DO I=1,INTORD
                       DO j=1,INTORD
                        xsi= Gcor(i)
                        eta= Gcor(j)
                        CALL Normal_Jac(v3,Jac,xsi,eta,ldim,lnodes,inci,cor)
!                        WRITE(11, '(3F10.6)') Jac,Wi(i),Wi(j)
			Area= Area + Jac*Wi(i)*Wi(j)
                       END DO
                      END DO Gauss_loop1
                    CASE DEFAULT
                END SELECT
               END DO Element_loop
 IF (ldim .EQ. 1) THEN
 WRITE(11,*) ' Length =',Area
 Else
 WRITE(11,*) ' Area =',Area
 END IF
END PROGRAM Compute_Area
