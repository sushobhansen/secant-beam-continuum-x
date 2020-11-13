      MODULE MATRIXFUNCTIONS
          CONTAINS
          
          FUNCTION INVERSEMAT(AIN,N) RESULT(C)
              
              !DOLITLLTE LU METHOD
              !https://ww2.odu.edu/~agodunov/computing/programs/book2/Ch06/Inverse.f90
              INCLUDE 'ABA_PARAM.INC'
              
              INTEGER, INTENT(IN) :: N
              DOUBLE PRECISION, INTENT(IN) :: AIN(N,N)
              DOUBLE PRECISION :: A(N,N), C(N,N), L(N,N), U(N,N), B(N), 
     & D(N), X(N)
              DOUBLE PRECISION :: COEFF
              
              !INITIALIZATION
              A=AIN
              L=0.0
              U=0.0
              B=0.0
              
              !STEP 1
              DO K=1,N-1
                  DO I=K+1,N
                      COEFF=A(I,K)/A(K,K)
                      L(I,K)=COEFF
                      DO J=K+1,N
                          A(I,J)=A(I,J)-COEFF*A(K,J)
                      END DO
                  END DO
              END DO
              
              !STEP 2
              DO I=1,N
                  L(I,I)=1.0
              END DO
              
              DO J=1,N
                  DO I=1,J
                      U(I,J)=A(I,J)
                  END DO
              END DO
              
              !STEP 3
              DO K=1,N
                  B(K)=1.0
                  D(1)=B(1)
                  
                  !STEP 3A: SOLVE LD=B
                  DO I=2,N
                      D(I)=B(I)
                      DO J=1,I-1
                          D(I)=D(I)-L(I,J)*D(J)
                      END DO
                  END DO
                  
                  !STEP 3B: SOLVE UX=D
                  X(N)=D(N)/U(N,N)
                  DO I=N-1,1,-1
                      X(I)=D(I)
                      DO J=N,I+1,-1
                          X(I)=X(I)-U(I,J)*X(J)
                      END DO
                      X(I)=X(I)/U(I,I)
                  END DO
                  
                  !STEP 3C: FILL IN SOLUTION INTO C
                  DO I=1,N
                      C(I,K)=X(I)
                  END DO
                  B(K)=0.0
              END DO
              
          END FUNCTION INVERSEMAT
          
          FUNCTION DIRACDELTA(X) RESULT(Y)
              DOUBLE PRECISION, INTENT(IN) :: X
              DOUBLE PRECISION :: Y
              
              IF(X .EQ. 0.0D0) THEN
                  Y = 1.0D0
              ELSE
                  Y = 0.0D0
              END IF
          END FUNCTION DIRACDELTA
          
      END MODULE MATRIXFUNCTIONS
      
      SUBROUTINE READINPUT(NSTEPS,NELEM,NIX,NIY,DBCOH,ELEMENTTYPE,
     & AX,BX,AY,BY,THICK,PROPS,YOUNG,POISSON,NCONSTRAINED,
     & CONSTRAINTS,NKNOWN,KNOWNDOFS,UKNOWN,FULLINT,FINPUT)
          INTEGER :: NSTEPS,NIX,NIY,DBCOH,NCONSTRAINED,
     & NKNOWN
          INTEGER:: ELEMENTTYPE(NELEM),CONSTRAINTS(100),
     & KNOWNDOFS(100)
          DOUBLE PRECISION :: AX(3),BX(3),AY(3),BY(3),THICK,
     & YOUNG(3),POISSON(3),UKNOWN(100)
          DOUBLE PRECISION :: PROPS(9)
          LOGICAL :: FULLINT
          CHARACTER(LEN=200) :: FINPUT
          
          OPEN(11,FILE=FINPUT)
          
          READ(11,*) NSTEPS

          DO I=1,NELEM
                  READ(11,*) ELEMENTTYPE(I)
          END DO
          
          DO I=1,NELEM
            READ(11,*) AX(I),BX(I),AY(I),BY(I)
          END DO
          
          DO I=1,NELEM
            READ(11,*) YOUNG(I),POISSON(I)
          END DO
          
          READ(11,*) THICK
          
          READ(11,*) FULLINT
          
          READ(11,*) NCONSTRAINED
          
          IF(NCONSTRAINED .GT. 100) THEN
              PRINT*, 'Too many constraints (max 100), 
     & please recompile!'
              CALL EXIT(0)
          END IF
          
          READ(11,*) (CONSTRAINTS(I),I=1,NCONSTRAINED)
          
          READ(11,*) NKNOWN
          
          IF(NKNOWN .GT. 100) THEN
              PRINT*, 'Too many known DOFs (max 100), 
     & please recompile!'
              CALL EXIT(0)
          END IF
          
          READ(11,*) (KNOWNDOFS(I),I=1,NKNOWN)
          
          READ(11,*) (UKNOWN(I),I=1,NKNOWN)
          
          READ(11,*) NIX,NIY
          
          READ(11,*) PROPS(1),PROPS(2)
          
          READ(11,*) PROPS(3),PROPS(4)
          
          READ(11,*) PROPS(5),PROPS(6)
          
          READ(11,*) PROPS(7),PROPS(8)
          
          READ(11,*) PROPS(9)
          
          READ(11,*) DBCOH
          
          CLOSE(11)
          
      END SUBROUTINE READINPUT