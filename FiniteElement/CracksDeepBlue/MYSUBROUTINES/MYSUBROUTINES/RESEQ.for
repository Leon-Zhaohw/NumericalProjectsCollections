C Copyright (c) 2014 David Patrick Hodapp
C
C Permission is hereby granted, free of charge, to any person obtaining a copy
C of this software and associated documentation files (the "Software"), to deal
C in the Software without restriction, including without limitation the rights   
C to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
C copies of the Software, and to permit persons to whom the Software is
C furnished to do so, subject to the following conditions:
C
C The above copyright notice and this permission notice shall be included in
C all copies or substantial portions of the Software.
C
C THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
C IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
C FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
C AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
C LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
C OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
C THE SOFTWARE.
C
      SUBROUTINE RESEQ(TPARRAY,ARRAYDIM1,N,FRST)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC     
C     PURPOSE:
C     RESEQUENCES DATA SET 
C
C     PARAMETERS:
C     - TPARRAY:   ARRAY OF STRESS SEQUENCES
C     - ARRAYDIM1: TPARRAY LENGTH (DIMENSION 1)
C                  **CANNOT CREATE EXPLICIT INTERFACE WITH ABAQUS
C     - N:         LENGTH OF ARRAY
C     - FRST:      INDEX OF REVERSAL TO BEGIN RESEQUENCED DATA SET
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      IMPLICIT NONE
C      
C     DECLARE CALLING PARAMETERS
C     **TPARRAY IS PASSED-BY REFERENCE, ADDT'L MEMORY NOT ALLOCATED 
      INTEGER, INTENT(IN) :: ARRAYDIM1,N,FRST
      REAL(KIND=8), DIMENSION(ARRAYDIM1,2), INTENT(INOUT) :: TPARRAY
C
C     DECLARE LOCAL VARIABLES
C     **DUMARRAY MIGHT BE LARGE, ARRAY IS DYNAMICALLY ALLOCATED ON 
C       THE HEAP INSTEAD OF ON THE STACK
      INTEGER :: I,J
      REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: DUMARRAY
      ALLOCATE(DUMARRAY(ARRAYDIM1,2))
C
C     RESEQUENCE ARRAY TO START WITH REVERSAL AT INDEX FRST
      DO I=FRST,N
         J=I-FRST+1
         DUMARRAY(J,1)=TPARRAY(I,1)
         DUMARRAY(J,2)=TPARRAY(I,2)
      END DO
      DO I=1,(FRST-1)
         J=I+N-FRST+1
         DUMARRAY(J,1)=TPARRAY(I,1)
         DUMARRAY(J,2)=TPARRAY(I,2)
      END DO
      DO I=1,N
          TPARRAY(I,1)=DUMARRAY(I,1)
          TPARRAY(I,2)=DUMARRAY(I,2)
      END DO
C      
      DEALLOCATE(DUMARRAY)
      RETURN
      END SUBROUTINE