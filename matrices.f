      SUBROUTINE CONNECTIONMAT(NX,NBDOFS,AX,BX,AY,BY,SMAT,TMAT)
          
          INCLUDE 'ABA_PARAM.INC'
          
          INTEGER :: NX,NBDOFS
          DOUBLE PRECISION :: AX,BX,AY,BY
          DOUBLE PRECISION :: SMAT(NBDOFS,6),TMAT(6,NBDOFS)
          
          DOUBLE PRECISION :: XMID,X
          
          XMID = (BX-AX)/2.0D0
          SMAT = 0.0D0
          
          DO I=1,NX+1
C Lower edge
              X = AX+(BX-AX)*(I-1)/NX
              SMAT(2*I-1,1) = 1.0D0
              SMAT(2*I,2:3) = (/1.0D0,XMID-X/)
C Upper edge
              SMAT(2*I+2*NX+1,4) = 1.0D0
              SMAT(2*I+2*NX+2,5:6) = (/1.0D0,XMID-X/)
          END DO
          
          TMAT = TRANSPOSE(SMAT)
          
      END SUBROUTINE CONNECTIONMAT
      
      SUBROUTINE BEAMSTIFFNESSMAT(AX,BX,AY,BY,THICK,YOUNG,KBEAM)
          
          DOUBLE PRECISION :: AX,BX,AY,BY,THICK,YOUNG
          DOUBLE PRECISION :: KBEAM(6,6)
          
          DOUBLE PRECISION :: LX,LY,A,IZ
          
C Note: KBEAM takes inputs in the global coordinate system but 
C is defined in the local coordinate system. This is why
C LX and LY are flipped.
          LY = BX-AX
          LX = BY-AY
          A = LY*THICK
          IZ = THICK*(LY**3)/12.0D0
          
          KBEAM = 0.0D0
          
C Axial components
          KBEAM(2,2) = A*YOUNG/LX
          KBEAM(2,5) = -KBEAM(2,2)
          KBEAM(5,2) = -KBEAM(2,2)
          KBEAM(5,5) = KBEAM(2,2)
          
C Transverse and rotation components
          KBEAM(1,1) = 12.0D0*YOUNG*IZ/LX**3
          KBEAM(1,3) = 6.0D0*YOUNG*IZ/LX**2
          KBEAM(3,1) = KBEAM(1,3)
          KBEAM(3,3) = 4.0D0*YOUNG*IZ/LX
          
          KBEAM(1,4) = -KBEAM(1,1)
          KBEAM(1,6) = KBEAM(1,3)
          KBEAM(3,4) = -KBEAM(1,3)
          KBEAM(3,6) = 2.0D0*YOUNG*IZ/LX
          
          KBEAM(4,1) = -KBEAM(1,1)
          KBEAM(4,3) = -KBEAM(1,3)
          KBEAM(6,1) = KBEAM(1,3)
          KBEAM(6,3) = 2.0D0*YOUNG*IZ/LX
          
          KBEAM(4,4) = KBEAM(1,1)
          KBEAM(4,6) = -KBEAM(1,3)
          KBEAM(6,4) = -KBEAM(1,3)
          KBEAM(6,6) = 4.0D0*YOUNG*IZ/LX
          
      END SUBROUTINE BEAMSTIFFNESSMAT
      
      SUBROUTINE ASSEMBLEBEAM(KGLOBAL,KLOCAL,NDOFS,ELEMENTDOFS)
          INTEGER :: NDOFS
          INTEGER :: ELEMENTDOFS(6)
          DOUBLE PRECISION :: KGLOBAL(NDOFS,NDOFS),KLOCAL(6,6)
          
          DO I=1,6
            DO J=1,6
                KGLOBAL(ELEMENTDOFS(I),ELEMENTDOFS(J))=
     & KGLOBAL(ELEMENTDOFS(I),ELEMENTDOFS(J))+KLOCAL(I,J)
              END DO
          END DO
      END SUBROUTINE ASSEMBLEBEAM
      
            SUBROUTINE LINEARELASTIC(YOUNG,POISSON,DMAT,ISOTROPIC)
c This subroutine evaluates the plane stress elastic matrix (D)
c Young = Young's Modulus, Poisson = Poisson ratio
          
          INCLUDE 'ABA_PARAM.INC'
          
          DOUBLE PRECISION :: YOUNG, POISSON
          DOUBLE PRECISION :: DMAT(3,3)
          LOGICAL :: ISOTROPIC
          
          DOUBLE PRECISION :: STIFFNUM=1.0D8,YOUNGX,YOUNGY
          
          YOUNGX=YOUNG*STIFFNUM
          YOUNGY=YOUNG
          
          DMAT = 0.0D0
          
          IF(ISOTROPIC) THEN !Isotropic
            DMAT(1,1) = YOUNG/(1.0-POISSON**2)
            DMAT(1,2) = YOUNG*POISSON/(1.0-POISSON**2)
            DMAT(2,1) = YOUNG*POISSON/(1.0-POISSON**2)
            DMAT(2,2) = YOUNG/(1.0-POISSON**2)
            DMAT(3,3) = 0.5*YOUNG*(1-POISSON)/(1.0-POISSON**2)
          ELSE !Orthotropic
            DMAT(1,1) = YOUNGX/(1.0-POISSON**2)
            DMAT(1,2) = YOUNGY*POISSON/(1.0-POISSON**2)
            DMAT(2,1) = YOUNGY*POISSON/(1.0-POISSON**2)
            DMAT(2,2) = YOUNGY/(1.0-POISSON**2)
            
            DMAT(3,3) = 0.5*YOUNG*(1-POISSON)/(1.0-POISSON**2)
            DMAT(3,3) = DMAT(3,3)*STIFFNUM
          END IF
          
          RETURN
      END SUBROUTINE LINEARELASTIC
      
      SUBROUTINE MATRICES(X,Y,AX,AY,BX,BY,NMAT,BMAT)
c This subroutine evaluates the matrix (N) of linear shape functions
c and the matrix of their derivatives B at point (X,Y) within an element
c The x-coords of the rectangular elements range from ax to bx
c and the y-coords range from ay to by
          
          INCLUDE 'ABA_PARAM.INC'
          
          DOUBLE PRECISION X,Y,AX,AY,BX,BY
          DOUBLE PRECISION :: NMAT(2,8), NN(4), XX(4), YY(4), XX2(4), 
     & YY2(4), SCALE(4), BMAT(3,8), DNNX(4), DNNY(4)
          
          XX = (/AX,BX,BX,AX/)
          YY = (/AY,AY,BY,BY/)
          XX2 = (/BX,AX,AX,BX/)
          YY2 = (/BY,BY,AY,AY/)
          NMAT = 0.0D0
          BMAT = 0.0D0
          
          SCALE =(XX-XX2)*(YY-YY2)
          NN = (X-XX2)*(Y-YY2)/SCALE
          DNNX = (Y-YY2)/SCALE
          DNNY = (X-XX2)/SCALE
          
          DO I=1,4
              NMAT(1,2*I-1) = NN(I)
              NMAT(2,2*I) = NN(I)
              
              BMAT(1,2*I-1) = DNNX(I)
              BMAT(2,2*I) = DNNY(I)
              BMAT(3,2*I-1) = DNNY(I)
              BMAT(3,2*I) = DNNX(I)
          END DO
          
          RETURN
      END SUBROUTINE MATRICES
      
      SUBROUTINE STIFFNESSMATELASTIC(AX,BX,AY,BY,K,YOUNG,POISSON,
     & THICK,FULLINT,ISOTROPIC)
c Evaluates the local stiffness matrix K for a rectangular element
c If FULLINT is true, K = integral(B'DB) at Gauss points
c If FULLINT is false, K = B'DB*area of element at mid-point
c The x-coords go from ax to bx and y-coords from ay to by
c The thickness is THICK
      
        INCLUDE 'ABA_PARAM.INC'
      
        DOUBLE PRECISION :: AX, BX, AY, BY, THICK, YOUNG, POISSON
        DOUBLE PRECISION :: K(8,8)
        LOGICAL :: FULLINT
        INTEGER :: NGP=2
        INTEGER :: CONSTRAINTS(2)
        DOUBLE PRECISION :: DMAT(3,3), NMAT(2,8), BMAT(3,8), GP(2)
        DOUBLE PRECISION :: AREA, JACOBIAN, X, Y
        LOGICAL :: ISOTROPIC

        AREA = (BX-AX)*(BY-AY)
        JACOBIAN = 0.25*AREA
        GP = (/-0.57735026919, 0.57735026919/)
        K = 0.0D0
        NMAT = 0.0D0
        BMAT = 0.0D0
        
        CALL LINEARELASTIC(YOUNG,POISSON,DMAT,ISOTROPIC)

        IF(FULLINT) THEN !Full integration
          DO I=1,NGP
              DO J=1,NGP
                X = 0.5*(BX-AX)*GP(I)+0.5*(BX+AX)
                Y = 0.5*(BY-AY)*GP(J)+0.5*(BY+AY)
                CALL MATRICES(X,Y,AX,AY,BX,BY,NMAT,BMAT)
C K = sum(B'DB) at Gauss points
                K = K + MATMUL(MATMUL(TRANSPOSE(BMAT),DMAT),BMAT)
              END DO
          END DO
          K = K*JACOBIAN
        ELSE ! Reduced integration
            X = 0.5*(AX+BX)
            Y = 0.5*(AY+BY)
            CALL MATRICES(X,Y,AX,AY,BX,BY,NMAT,BMAT)
            K = K + MATMUL(MATMUL(TRANSPOSE(BMAT),DMAT),BMAT)
            K = K*AREA
        END IF
        
        K = K*THICK

        RETURN
      END SUBROUTINE STIFFNESSMATELASTIC
      
      SUBROUTINE ASSIGNNODEDOFS(GLOBALDOFS,NIX,NIY)
c Assigns global DOFs to each node (2 per node, ABAQUS convention)
c Boundary nodes are assigned 1-8 in ABAQUS' order first (CC-wise)
c Then, DOFs are assigned to other nodes by column
          INCLUDE 'ABA_PARAM.INC'
          
          INTEGER :: NIX,NIY
          INTEGER :: GLOBALDOFS((NIX+1),(NIY+1),2)
          INTEGER :: II
          
          II=1
          DO J=1,NIY+1
              DO I=1,NIX+1
                GLOBALDOFS(I,J,1:2) = (/II,II+1/)
                II=II+2
              END DO
          END DO
          

      END SUBROUTINE ASSIGNNODEDOFS
      
      SUBROUTINE ASSIGNELEMENTDOFS(ELEMENTDOFS,GLOBALDOFS,NIX,NIY)
C Accumulate DOFs corresponding to the nodes of each element
c in ABAQUS order (CC-wise)
c Essentially creates a connectivity matrix for each element
c but for DOFs instead of node numbers
          INCLUDE 'ABA_PARAM.INC'
          INTEGER :: NIX,NIY
          INTEGER :: ELEMENTDOFS(NIX,NIY,8),GLOBALDOFS(NIX+1,NIY+1,2)
          
          DO I=1,NIX
              DO J=1,NIY
                ELEMENTDOFS(I,J,1:2) = GLOBALDOFS(I,J,:)
                ELEMENTDOFS(I,J,3:4) = GLOBALDOFS(I+1,J,:)
                ELEMENTDOFS(I,J,5:6) = GLOBALDOFS(I+1,J+1,:)
                ELEMENTDOFS(I,J,7:8) = GLOBALDOFS(I,J+1,:)
              END DO
          END DO
          RETURN
      END SUBROUTINE ASSIGNELEMENTDOFS
      
      SUBROUTINE ASSEMBLEGLOBAL(KGLOBAL,KLOCAL,NDOFS,ELEMENTDOFS)
c Assemble local stiffness matrix KLOCAL consisting of DOFs
c specified in ELEMENTDOFs into the global matrix KGLOBAL
c KLOCAL is always (8,8) while KGLOBAL is (NDOFS,NDOFS)
c NDOFs is total # of DOFs of the system
          INTEGER :: NDOFS, ELEMENTDOFS(8)
          DOUBLE PRECISION :: KGLOBAL(NDOFS,NDOFS), KLOCAL(8,8)
          
          DO I=1,8
              DO J=1,8
                  KGLOBAL(ELEMENTDOFS(I),ELEMENTDOFS(J))=
     & KGLOBAL(ELEMENTDOFS(I),ELEMENTDOFS(J))+KLOCAL(I,J)
              END DO
          END DO
          RETURN
      END SUBROUTINE ASSEMBLEGLOBAL
      
            SUBROUTINE STIFFNESSMATCOHESIVE(PROPS,AMATRX,COORDS,
     & KSTEP,U,T_d,NDOFEL)
C Interface to the ABAQUE UEL for the PPR code from Park & Paulino, 2012
C AMATRX is the cohesive stiffness matrix generated by the UEL
C Global coordinates of the nodes of the element, in ABAQUS order (CC-wise)
C are stores in COORDS, corresponding displacements in U
C NDOFEL is # of DOFs of the element, always 8
C FGPSEP, FGPTRAC - unit handles for writing separation and traction to files
C KSTEP - time step, not used by the model but only for controlling writing
          INCLUDE 'ABA_PARAM.INC'
          
          INTEGER :: NDOFEL,KSTEP
          INTEGER :: NRHS=1, MCRD=2, NNODE=4, MLVARX=8, NPREDF=8,
     & MDLOAD=8
          DOUBLE PRECISION :: AMATRX(NDOFEL,NDOFEL), COORDS(2,4),
     & U(NDOFEL), PROPS(9), T_d(2,2,2)
          
C The next two allocate statements are to create an interface with UEL
C Most of the inputs are never used
          ALLOCATABLE :: RHS(:,:), SVARS(:), ENERGY(:), 
     & DU(:,:), V(:), A(:), TIME(:), PARAMS(:),
     & JDLTYP(:,:), ADLMAG(:,:), DDLMAG(:,:),
     & PREDEF(:, :, :), LFLAGS(:), JPROPS(:)
      
          ALLOCATE(RHS(MLVARX,NRHS), SVARS(4), ENERGY(8), 
     & DU(MLVARX,100), V(NDOFEL), A(NDOFEL), TIME(2), PARAMS(100),
     & JDLTYP(MDLOAD,100), ADLMAG(MDLOAD,100), DDLMAG(MDLOAD,100),
     & PREDEF(2, NPREDF, NNODE), LFLAGS(100), JPROPS(100))

C PROPS - 9 parameters of the PPR model
C PROPS(1),(2) are tangential and normal fracture energy respectively
c PROPS(3),(4) are tan and norm strength respectively
C PROPS(5),(6) are alpha and beta respectively (shape parameters)
C PROPS(7),(8) are lambda_t and lambda_n respectively (softening slopes)
C PROPS(9) is thickness, should agree with THICK variable for elastic K
      
          AMATRX = 0.0
C SVARS is used to collect displacement history, but ignored for now
          SVARS = 0.0

C Call UEL for PPR model (in file PPR_model_educational.f)
C Returns RHS vector (of internal forces: not used in this code)
C Returns AMATRX - local cohesive stiffness matrix
C Uses COORDS and U to determine separation, from where traction
C and stiffness are calculated
          CALL UEL(RHS, AMATRX, SVARS, ENERGY, NDOFEL, NRHS, NSVARS,
     & PROPS, NPROPS, COORDS, MCRD, NNODE, U, DU, V, A, JTYPE, TIME,
     & DTIME, KSTEP, KINC, JELEM, PARAMS, NDLOAD, JDLTYP, ADLMAG,
     & PREDEF, NPREDF, LFLAGS, MLVARX, DDLMAG, MDLOAD, PNEWDT, JPROPS,
     & NJPRO, PERIOD, T_d)
      
          DEALLOCATE(RHS, SVARS, ENERGY, DU, V, A, TIME, PARAMS,
     & JDLTYP, ADLMAG, DDLMAG, PREDEF, LFLAGS, JPROPS)
          
          RETURN
          
      END SUBROUTINE STIFFNESSMATCOHESIVE
      
      SUBROUTINE CONDENSEDMAT(K,KBBTILDE,INVKII,KIB,BDOFS,
     & IDOFS,NDOFS,NBDOFS,NIDOFS)
C Code for static condensation of the global stiffness mat K
C The problem Ku=f is partitioned as  [Kbb Kbi;Kib Kii](ub;ui)=(fb;fi)
C This partitioning can be done easily since ub and ui are ordered
C consecutively in the ASSIGNNODEDOFS subroutine
C Note that fb = 0 = fi
C The condensed stiffness matrix is Kbbtilde = Kbb-Kbi*Kii^(-1)*Kib
C The subroutine also stores Kii^(-1) and Kib because they are needed
C to calculate ui = -Kii^(-1)*Kib*ub
C IDOFS - list of internal DOFs to be condensed out
C BDOFS - list of boundary DOFs (always 1-8)
          
          USE MATRIXFUNCTIONS    
          INCLUDE 'ABA_PARAM.INC'
      
          INTEGER :: NDOFS, NBDOFS, NIDOFS
          INTEGER :: BDOFS(NBDOFS), IDOFS(NIDOFS)
          DOUBLE PRECISION :: K(NDOFS,NDOFS), KBBTILDE(NBDOFS,NBDOFS), 
     & KIB(NIDOFS,NBDOFS), INVKII(NIDOFS,NIDOFS)
          DOUBLE PRECISION :: KBB(NBDOFS,NBDOFS), KBI(NBDOFS,NIDOFS),
     & KII(NIDOFS,NIDOFS)
      
          KBB = 0.0
          KBI = 0.0
          KIB = 0.0
          KII = 0.0
          KBBTILDE = 0.0

C Partition K into Kbb,Kbi,Kib,Kii
          DO I=1,NBDOFS
              DO J=1,NBDOFS
                  KBB(I,J)=K(BDOFS(I),BDOFS(J))
              END DO
          END DO
          
          DO I=1,NBDOFS
              DO J=1,NIDOFS
                  KBI(I,J)=K(BDOFS(I),IDOFS(J))
              END DO
          END DO
          
          DO I=1,NIDOFS
              DO J=1,NBDOFS
                  KIB(I,J)=K(IDOFS(I),BDOFS(J))
              END DO
          END DO
          
          DO I=1,NIDOFS
              DO J=1,NIDOFS
                  KII(I,J)=K(IDOFS(I),IDOFS(J))
              END DO
          END DO
          
          INVKII = INVERSEMAT(KII,NIDOFS)
          KBBTILDE = KBB-MATMUL(KBI,MATMUL(INVKII,KIB))
      
          RETURN
      END SUBROUTINE CONDENSEDMAT
      
      SUBROUTINE WRITEEXTERNALMESH(NELEM,AX,BX,AY,BY,FINPUT)
          INTEGER :: NELEM
          DOUBLE PRECISION :: AX(NELEM),BX(NELEM),AY(NELEM),
     & BY(NELEM)
          CHARACTER(LEN=200) :: FINPUT,FOUTPUT
          
          FOUTPUT = TRIM(FINPUT)//'_extmesh.txt'
          OPEN(UNIT=20,FILE=FOUTPUT)
          
          WRITE(20,*) AX(2)+(BX(2)-AX(2))/2,AY(1)
          WRITE(20,*) AX(2)+(BX(2)-AX(2))/2,BY(1)
          WRITE(20,*) AX(2)+(BX(2)-AX(2))/2,SUM(BY(1:2))
          WRITE(20,*) AX(2)+(BX(2)-AX(2))/2,SUM(BY(1:3))
          
          CLOSE(20)
          
      END SUBROUTINE WRITEEXTERNALMESH
      
      SUBROUTINE WRITEINTERNALMESH(NX,NY,AX,BX,AY,BY,FINPUT)
          INTEGER :: NX,NY
          DOUBLE PRECISION :: AX,BX,AY,BY
          CHARACTER(LEN=200) :: FINPUT,FOUTPUT
          
          FOUTPUT = TRIM(FINPUT)//'_intmesh.txt'
          OPEN(UNIT=20,FILE=FOUTPUT)
          
          DO J=1,(NY+1)/2
              DO I=1,(NX+1)
                    WRITE(20,*) AX+(BX-AX)*(I-1)/NX,
     & AY+(BY-AY)*(J-1)/(NY-1)
              END DO
          END DO
                  
            J=(NY+3)/2
            DO I=1,(NX+1)
                  WRITE(20,*) AX+(BX-AX)*(I-1)/NX,
     & AY+(BY-AY)*(J-2)/(NY-1)
            END DO
              
      
          DO J=(NY+5)/2,NY+1
              DO I=1,(NX+1)
                    WRITE(20,*) AX+(BX-AX)*(I-1)/NX,
     & AY+(BY-AY)*(J-2)/(NY-1)
              END DO
          END DO
          
          CLOSE(20)
          
      END SUBROUTINE WRITEINTERNALMESH