PROGRAM coordi
IMPLICIT NONE
INTEGER:: i,j,k,m,n,p,s  ; REAL:: aa
REAL,ALLOCATABLE::coord(:,:)
OPEN(11,FILE='coord_3.dat')
n = 100; m = n+1; p = m*m     ; aa = 10.0/n
ALLOCATE(coord(1:6*n*n+2,1:3))

! Plane 1 z = .0
DO i = 1 , p      ;   coord(i,3) = .0  ; END DO
coord(1,1) = 10.0; coord(1,2) = 0.0
coord(2,1) = 10.0; coord(2,2) = aa
j = 1
DO i = 4 , 2*m , 2
 coord(i,1) = 10.0 - j*aa; coord(i,2) = .0
 j = j + 1
END DO
j = 1
DO i = 3 , 2*m - 1 , 2
 coord(i,1) = 10.0 - j*aa; coord(i,2) = aa
 j = j + 1
END DO
DO j = 2*m+1,n*m+1,m  
 DO i = 1 , m
    coord(i+j-1,1) = (m-i)*aa ; coord(i+j-1,2) = (j-1)/m*aa
 END DO
END DO

! Plane 2 z = 10.0
k = p + 1
DO i = k , 2*p   ;   coord(i,3) = 10.0;  END DO
coord(k,1) = 10.0; coord(k,2) = 0.0
coord(k+1,1) = 10.0 - aa; coord(k+1,2) = 0.0
j = 1
DO i = k+2 ,k+2*m-2 , 2
 coord(i,1) = 10.0 - aa; coord(i,2) = j*aa
 j = j + 1
END DO
j = 1
DO i = k+3 ,k + 2*m - 1 , 2
 coord(i,1) = 10.0 ; coord(i,2) =j * aa
 j = j + 1
END DO
s = 2
DO j = k + 2*m ,k + n*m , m  
 DO i = 1 , m
    coord(i+j-1,1) = 10.0 - s * aa ; coord(i+j-1,2) = (i - 1) * aa   
 END DO
 s = s + 1
END DO

! Plane 3 y = 10.0
k = 2*p + 1
DO i = k , 3*p - 2 * m ; coord(i,2) = 10.0;  END DO
s = 1
DO j = k , k + (n-2)*m  , m
   DO i = 1 , m
      coord(i + j -1,1) = 10.0 - (i-1)*aa ; coord(i + j - 1, 3) = s * aa
   END DO
   s = s + 1
END DO

! Plane 4 y = .0
k = 3*p -2*m + 1
DO i = k , 4*p - 4*m ; coord(i,2)  =  0.0;  END DO
j = 1
DO i = k , k + 2*(n-2) , 2
 coord(i,1) = 10.0 - aa;  coord(i,3) = j*aa
 j = j + 1
END DO
j = 1
DO i = k + 1 , k + 2*(n-2) + 1 , 2
 coord(i,1) =  10.0 ; coord(i,3) =  j * aa
 j = j + 1
END DO
s = 1
DO j = k + 2*(n-1) , k + n * (n - 1)  , n - 1   
 DO i = 1 , n - 1
    coord(i+j-1,1) = 9.0 - s * aa ; coord(i+j-1,3) = i * aa  
 END DO
 s = s + 1
END DO

! Plane 5 x = .0
k =  4*p -4*m + 1 
DO i = k , 5*p -6*m -2*(n-1)  ; coord(i,1) = .0  ; END DO
s = 1
DO j = k , k + (n-1)*(n-2) , n-1
   DO i = 1 , n-1
      coord( i + j - 1 ,2) =  s * aa
      coord( i + j - 1 ,3) =  i * aa
   END DO
   s = s + 1
END DO

! Plane 6 x = 10.0
k = 5*p -6*m -2*(n-1) +  1
DO i = k , 6*n*n+2 ; coord(i,1) = 10.0 ;  END DO
s = 1
DO j  = k , k + (n-1)*(n-2), n-1
   DO i = 1 , n-1
      coord( i + j - 1 ,2) =  i * aa
      coord( i + j - 1 ,3) =  s * aa
   END DO
   s = s + 1
END DO
DO i = 1 , 6*n*n+2   ;   WRITE(11,*)coord(i,:)   ; END DO
END PROGRAM coordi



