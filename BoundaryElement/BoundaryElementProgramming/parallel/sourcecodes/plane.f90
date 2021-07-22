PROGRAM plane
! Incidences for planes 1 to 6
IMPLICIT NONE
INTEGER::i,j,k,m,n,p,   count  
INTEGER,ALLOCATABLE::plane1(:,:),plane2(:,:),plane3(:,:),inci(:,:),   &
                     plane4(:,:),plane5(:,:),plane6(:,:)
n = 100; m = n + 1  ; p = m*m 
ALLOCATE(plane1(m,m),plane2(m,m),plane3(m,m),inci(1:6*n*n,1:4),       &
         plane4(m,m),plane5(m,m),plane6(m,m))
OPEN(11,FILE='plane_3.dat')

! Plane 1 z = .0
plane1(1,1) = 1; plane1(2,1) = 2
DO i = 4 , 2*m , 2
   plane1(1,i/2) = i
END DO
DO i = 3 , 2*m - 1 , 2
   plane1(2,(i+1)/2) = i
END DO
DO j = 2*m+1 , n*m+1 , m
   DO i = 1 , m
      plane1((j-1)/m+1, i ) = j + i - 1
   END DO
END DO
count = 0
DO i = 1 , n
   DO j = 1 , n
      count = count + 1
      inci(count,1) = plane1(i,j) ; inci(count,2) = plane1(i+1,j)
      inci(count,3) = plane1(i+1,j+1); inci(count,4) = plane1(i,j+1)
   END DO
END DO

! Plane 2 z = 10.0
k = p + 1
plane2(1,1) = k; plane2(1,2) = k + 1
j = 1
DO i = k+2,k+2*m-2,2
   plane2(j+1,2) = i   ; j = j + 1
END DO
j = 1
DO i = k+3,k+2*m-1,2
   plane2(j+1,1) = i   ; j = j + 1
END DO
DO j = k + 2*m, k + n*m , m
   DO i = 1 , m
      plane2(i, (j-k)/m+1 )  =  i + j - 1
   END DO
END DO
DO j = 1 , n
   DO i = 1 , n
      count = count + 1
      inci(count,1) = plane2(i,j) ; inci(count,2) = plane2(i,j+1)
      inci(count,3) = plane2(i+1,j+1); inci(count,4) = plane2(i+1,j)
   END DO
END DO

! Plane 3 y = 10.0
k = 2*p + 1
DO j = k , k + (n-2)*m  , m
   DO i = 1 , m 
      plane3( (j-k)/m+2 , i  )  =  i + j - 1
   END DO
END DO
plane3(1,:) = plane1(m,:) ; plane3(m,:) = plane2(m,:)

DO i = 1 , n
   DO j = 1 , n
      count = count + 1
      inci(count,1) = plane3(i,j) ; inci(count,2) = plane3(i+1,j)
      inci(count,3) = plane3(i+1,j+1); inci(count,4) = plane3(i,j+1)
   END DO
END DO

! Plane 4 y = 0.0
k = 3*p -2*m + 1
j = 1
DO i = k , k+2*(n-2) , 2
   plane4(j+1,2) = i   ; j = j + 1
END DO
j = 1
DO i = k+1, k+2*(n-2)+1 , 2
   plane4(j+1,1) = i   ; j = j + 1
END DO
DO j = k + 2*(n-1), k + n*(n-1) , n-1
   DO i = 2 , n
      plane4(i, (j-k)/(n-1)+1 )  =  i + j - 2
   END DO
END DO
plane4(m,:) = plane2(1,:) ; plane4(1,:) = plane1(1,:)
DO j = 1 , n
   DO i = 1 , n
      count = count + 1
      inci(count,1) = plane4(i,j) ; inci(count,2) = plane4(i,j+1)
      inci(count,3) = plane4(i+1,j+1); inci(count,4) = plane4(i+1,j)
   END DO
END DO

! Plane 5 x = 0.0
k = 4*p -4*m + 1
DO j = k , k + (n-1)*(n-2)  , n-1
   DO i = 2 , n   
      plane5( i, (j-k)/(n-1) +2 )  =  i + j - 2
   END DO
END DO
plane5(2:n,1) = plane4(2:n,m) ; plane5(2:n,m) = plane3(2:n,m)
plane5(1,:) = plane1(:,m) ; plane5(m,:) = plane2(:,m)
DO j = 1 , n
   DO i = 1 , n
      count = count + 1
      inci(count,1) = plane5(i,j) ; inci(count,2) = plane5(i,j+1)
      inci(count,3) = plane5(i+1,j+1); inci(count,4) = plane5(i+1,j)
   END DO
END DO

! Plane 6 x = 10.0
k = 5*p -6*m -2*(n-1) + 1
DO j = k , k + (n-1)*(n-2)   ,  n-1 
   DO i = 2 , n 
      plane6( (j-k)/(n-1) + 2 , i )  =  i + j - 2
   END DO
END DO
plane6(1,:) = plane1(:,1) ; plane6(m,:) = plane2(:,1)
plane6(2:n,1) = plane4(2:n,1) ; plane6(2:n,m) = plane3(2:n,1)
DO i = 1 , n
   DO j = 1 , n
      count = count + 1
      inci(count,1) = plane6(i,j) ; inci(count,2) = plane6(i+1,j)
      inci(count,3) = plane6(i+1,j+1); inci(count,4) = plane6(i,j+1)
   END DO
END DO
DO i = 1 , count
   WRITE(11,'(4I7)') inci(i,:)
END DO
END PROGRAM plane

