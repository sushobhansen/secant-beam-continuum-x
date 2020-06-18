      PROGRAM MAIN
          USE MATRIXFUNCTIONS
          INCLUDE 'ABA_PARAM.INC' 
          
          INTEGER :: NELEM=5,NIX,NIY,DBCOH,NDOFS,NCOMPDOFS,
     & NBCOMPDOFS,NICOMPDOFS,FSEP=11,FDISPB=12,FDISPI=13,
     & FFORCE=14,TSTEP,NSTEPS
          INTEGER :: LOCALDOFS(6)
          INTEGER, ALLOCATABLE :: ELEMENTDOFS(:,:),
     & ELEMENTTYPE(:)
          DOUBLE PRECISION :: AX(5),BX(5),AY(5),BY(5),THICK,
     & YOUNG(5),POISSON(5)
          DOUBLE PRECISION :: PROPS(9),KCOMPBEAM(6,6),KBEAM(6,6)
          DOUBLE PRECISION, ALLOCATABLE :: KGLOBAL(:,:),INVKII(:,:,:),
     & INVKIIE(:,:),KIB(:,:,:),KIBE(:,:)
          CHARACTER(LEN=200) :: FINPUT,FOUTPUT
          LOGICAL :: FULLINT
          
          INTEGER :: NCONSTRAINED,NKNOWN,NFREE,II
          INTEGER :: CONSTRAINTS(100),KNOWNDOFS(100)
          INTEGER, ALLOCATABLE :: FREEDOFS(:)
          DOUBLE PRECISION :: UB(6),UKNOWN(100)
          DOUBLE PRECISION, ALLOCATABLE :: UFREE(:),U(:),
     & DUFREE(:),KFF(:,:),KFK(:,:),UI(:,:),UIE(:),
     & T_d(:,:,:,:),FGLOBAL(:)      

          ALLOCATE(ELEMENTTYPE(NELEM))
          CALL GETARG(1,FINPUT)
          CALL READINPUT(NSTEPS,NELEM,NIX,NIY,DBCOH,ELEMENTTYPE,
     & AX,BX,AY,BY,THICK,PROPS,YOUNG,POISSON,NCONSTRAINED,
     & CONSTRAINTS,NKNOWN,KNOWNDOFS,UKNOWN,FULLINT,FINPUT)
          
          NDOFS = 3*(NELEM+1)
          NIY = 2*NIY+1
          NCOMPDOFS = (NIX+1)*(NIY+1)*2
          NBCOMPDOFS = 4*(NIX+1)
          NICOMPDOFS = NCOMPDOFS-NBCOMPDOFS
          NFREE = NDOFS-NCONSTRAINED-NKNOWN
          
          ALLOCATE(ELEMENTDOFS(NELEM,6),
     & KGLOBAL(NDOFS,NDOFS),U(NDOFS),UI(NELEM,NICOMPDOFS),
     & UIE(NICOMPDOFS),INVKII(NELEM,NICOMPDOFS,NICOMPDOFS),
     & INVKIIE(NICOMPDOFS,NICOMPDOFS),KIB(NELEM,NICOMPDOFS,NBCOMPDOFS),
     & KIBE(NICOMPDOFS,NBCOMPDOFS),T_d(NIX,2,2,2))
      
          DO I=1,NELEM
            ELEMENTDOFS(I,1:6) = (/(J,J=3*I-2,3*I+3)/)
          END DO


          ALLOCATE(FREEDOFS(NFREE),UFREE(NFREE),DUFREE(NFREE),
     & KFF(NFREE,NFREE),KFK(NFREE,NKNOWN),FGLOBAL(NDOFS))
          II = 1
          DO I=1,NDOFS
            IF(ALL(CONSTRAINTS(1:NCONSTRAINED) .NE. I) .AND. 
     & ALL(KNOWNDOFS(1:NKNOWN) .NE. I)) THEN
                FREEDOFS(II) = I
                II = II+1
            END IF
          END DO  
     
          UI = 0.0D0
          UIE = 0.0D0
          U = 0.0D0
          UB =0.0D0
          UFREE = 0.0D0
          FGLOBAL = 0.0D0
          KBEAM = 0.0D0
          KCOMPBEAM = 0.0D0
          INVKII = 0.0D0
          INVKIIE = 0.0D0
          KIB = 0.0D0
          KIBE = 0.0D0
          T_d = 0.0D0
          
          FOUTPUT = TRIM(FINPUT)//'_boundary-disp.txt'
          OPEN(UNIT=FDISPB,FILE=FOUTPUT) 
          FOUTPUT = TRIM(FINPUT)//'_internal-disp.txt'
          OPEN(UNIT=FDISPI,FILE=FOUTPUT)
          FOUTPUT = TRIM(FINPUT)//'_separation.txt'
          OPEN(UNIT=FSEP,FILE=FOUTPUT) 
          FOUTPUT = TRIM(FINPUT)//'_boundary-force.txt'
          OPEN(UNIT=FFORCE,FILE=FOUTPUT)
          
          DO TSTEP=1,NSTEPS
              KGLOBAL = 0.0D0
              PRINT*, 'TIME STEP ',TSTEP,' OF ',NSTEPS
          DO I=1,NELEM
                IF(ELEMENTTYPE(I) .EQ. 2) THEN
C Composite
                  LOCALDOFS(:) = ELEMENTDOFS(I,:)
                  DO II=1,6
                      UB(II) = U(LOCALDOFS(II))
                  END DO
                  INVKIIE(:,:) = INVKII(I,:,:)
                  KIBE = KIB(I,:,:)
                  CALL CRACKELEMENTX(UB,AX(I),BX(I),AY(I),BY(I),
     & THICK,YOUNG(I),POISSON(I),PROPS,TSTEP,INVKIIE,KIBE,KCOMPBEAM,
     & T_d,NIX,NIY,NCOMPDOFS,NBCOMPDOFS,NICOMPDOFS,FULLINT,FSEP,
     & FDISPI,DBCOH)
                    KIB(I,:,:) = KIBE(:,:)
                    INVKII(I,:,:) = INVKIIE(:,:)
                    CALL ASSEMBLEBEAM(KGLOBAL,KCOMPBEAM,NDOFS,
     & LOCALDOFS)  
                ELSE
C Elastic
                  LOCALDOFS(:) = ELEMENTDOFS(I,:)
                  CALL BEAMSTIFFNESSMAT(AX(I),BX(I),AY(I),BY(I),THICK,
     & YOUNG(I),KBEAM)
                  CALL ASSEMBLEBEAM(KGLOBAL,KBEAM,NDOFS,
     & LOCALDOFS)  
                END IF
          END DO
          
          DO I=1,NFREE
                DO J=1,NKNOWN
                    KFK(I,J)=KGLOBAL(FREEDOFS(I),KNOWNDOFS(J))
                END DO
            END DO

            DO I=1,NFREE
                DO J=1,NFREE
                      KFF(I,J)=KGLOBAL(FREEDOFS(I),FREEDOFS(J))
                END DO
            END DO
            
            UFREE(1:NFREE) = -MATMUL(INVERSEMAT(KFF(1:NFREE,1:NFREE),
     & NFREE),MATMUL(KFK(1:NFREE,1:NKNOWN),
     & UKNOWN(1:NKNOWN)*TSTEP/NSTEPS))
            
            DO I=1,NKNOWN
                U(KNOWNDOFS(I))=UKNOWN(I)*TSTEP/NSTEPS
            END DO

            DO I=1,NFREE
                U(FREEDOFS(I))=UFREE(I)
            END DO
            
            FGLOBAL = MATMUL(KGLOBAL,U)
            
            WRITE(FDISPB,'(I5,*(E15.6))') TSTEP,U
            WRITE(FFORCE,'(I5,*(E15.6))') TSTEP,FGLOBAL
          END DO
          
          CLOSE(FDISPB)
          CLOSE(FDISPI)
          CLOSE(FSEP) 
          
C          CALL WRITEEXTERNALMESH(NELEM,AX,BX,AY,BY,FINPUT)
C          CALL WRITEINTERNALMESH(NIX,NIY,AX(2),BX(2),BY(1),
C     & SUM(BY(1:2)),FINPUT)
          
      
      END PROGRAM MAIN