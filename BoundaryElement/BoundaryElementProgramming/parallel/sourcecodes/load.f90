PROGRAM load
! Loading conditions for simple boxes
IMPLICIT NONE
INTEGER::i,n,count,pl,nels,bc 
REAL::p(1:6,1:3)
OPEN(11,FILE='load_3.dat')    
n = 100; nels = n*n ; bc = 0  
WRITE(11,*) bc;  WRITE(11,*) 6*nels
p(1,:) = (/.0,.0,10./) ;  p(2,:) = (/.0,.0,-10./); p(3,:) = (/.0,-10.,.0/)
p(4,:) = (/.0,10.,.0/) ; p(5,:) = (/10.,.0,.0/) ; p(6,:) = (/-10.,.0,.0/)
count = 0
DO pl = 1 , 6
 DO i = 1 , nels
   count = count + 1
   WRITE(11,'(I7,12F7.2)') count,p(pl,:),p(pl,:),p(pl,:),p(pl,:)  
 END DO
END DO
END PROGRAM load
