#include <fintrf.h>
C
C laplaceDLPfmm.F 
C
C This subroutine is the main gateway to MATLAB.  When a MEX function
C  is executed MATLAB calls the MEXFUNCTION subroutine in the corresponding
C  MEX file.  
C
C
      SUBROUTINE MEXFUNCTION(NLHS, PLHS, NRHS, PRHS)
C-----------------------------------------------------------------------
C     (pointer) Replace integer by integer*8 on 64-bit platforms
C
      MWPOINTER PLHS(*), PRHS(*)

C-----------------------------------------------------------------------
C

      INTEGER NLHS, NRHS
      INTEGER ERROR
      
C
C-----------------------------------------------------------------------
C     (pointer) Replace integer by integer*8 on 64-bit platforms
C
      MWPOINTER MXCREATEDOUBLEMATRIX, MXGETPR

C-----------------------------------------------------------------------
C

      MWSIZE MXGETM, MXGETN
C
C KEEP THE ABOVE SUBROUTINE, ARGUMENT, AND FUNCTION DECLARATIONS FOR USE
C IN ALL YOUR FORTRAN MEX FILES.
C---------------------------------------------------------------------
C
C-----------------------------------------------------------------------
C     (pointer) Replace integer by integer*8 on 64-bit platforms
C
      MWPOINTER Q,XS,YS,DIR1,DIR2
      MWPOINTER POT

C-----------------------------------------------------------------------
C

      MWSIZE MS, NS, NEL
      REAL *8, ALLOCATABLE :: RQ(:),RXS(:),RYS(:)
      REAL *8, ALLOCATABLE :: RDIR1(:),RDIR2(:)
      REAL *8, ALLOCATABLE :: RPOT(:)

#if defined MSWIND
C     For Windows only!
C     This resets the floating point exception to allow divide by zero,
C     overflow and invalid numbers. 
C     
      INTEGER(2) CONTROL
      CALL GETCONTROLFPQQ(CONTROL)
      CONTROL = CONTROL .OR. FPCW$ZERODIVIDE
      CONTROL = CONTROL .OR. FPCW$INVALID
      CONTROL = CONTROL .OR. FPCW$OVERFLOW
      CALL SETCONTROLFPQQ(CONTROL)
#endif

C CHECK FOR PROPER NUMBER OF ARGUMENTS
C
      IF (NRHS .NE. 5) THEN
        CALL MEXERRMSGTXT('FMM requires five input arguments')
      ELSEIF (NLHS .NE. 1) THEN
        CALL MEXERRMSGTXT('FMM Requires one output argument')
      ENDIF
      

c     GET THE NUMBER OF CHARGES/SOURCES.
      MS = MXGETM(PRHS(1))
      NS = MXGETN(PRHS(1));
      MS = MAX(MS,NS)
c     MS is the number of sources

      ALLOCATE(RQ(MS),STAT=ERROR)
      ALLOCATE(RXS(MS),STAT=ERROR)
      ALLOCATE(RYS(MS),STAT=ERROR)
      ALLOCATE(RDIR1(MS),STAT=ERROR)
      ALLOCATE(RDIR2(MS),STAT=ERROR)
      ALLOCATE(RPOT(MS),STAT=ERROR)
C
C CREATE A MATRIX FOR RETURN ARGUMENT
      PLHS(1) = MXCREATEDOUBLEMATRIX(MS,1,0)
      PLHS(2) = MXCREATEDOUBLEMATRIX(MS,1,0)
C
C ASSIGN POINTERS TO THE VARIOUS PARAMETERS
C
      Q = MXGETPR(PRHS(1))
      XS = MXGETPR(PRHS(2))
      YS = MXGETPR(PRHS(3))
      DIR1 = MXGETPR(PRHS(4))
      DIR2 = MXGETPR(PRHS(5))
C
      POT = MXGETPR(PLHS(1))
C
C COPY RIGHT HAND ARGUMENTS TO LOCAL ARRAYS OR VARIABLES
      NEL = MS
      CALL MXCOPYPTRTOREAL8(Q, RQ, NEL)
      CALL MXCOPYPTRTOREAL8(XS, RXS, NEL)
      CALL MXCOPYPTRTOREAL8(YS, RYS, NEL)
      CALL MXCOPYPTRTOREAL8(DIR1, RDIR1, NEL)
      CALL MXCOPYPTRTOREAL8(DIR2, RDIR2, NEL)

C
C DO THE ACTUAL COMPUTATIONS IN A SUBROUTINE
C       CREATED ARRAYS.
C

      CALL laplaceDLP(RQ,RXS,RYS,MS,RDIR1,RDIR2,RPOT)

C
C COPY OUTPUT WHICH IS STORED IN LOCAL ARRAY TO MATRIX OUTPUT
      NEL = MS
      CALL MXCOPYREAL8TOPTR(RPOT, POT, NEL)
C

      DEALLOCATE(RQ)
      DEALLOCATE(RXS)
      DEALLOCATE(RYS)
      DEALLOCATE(RDIR1)
      DEALLOCATE(RDIR2)
      DEALLOCATE(RPOT)

      RETURN
      END
