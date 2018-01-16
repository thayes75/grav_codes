      PROGRAM MAKE_HIST_SUMA

      parameter (NF=800, NTSTP=11000)
      character fnout*30
      dimension tslipp(NF,NTSTP)
      dimension islip(NF),slip1(NF),slip2(NF)
      dimension dslpS(NF,NTSTP),dslpD(NF,NTSTP), dslpT(NF,NTSTP)
      dimension initeq(NF,NTSTP)

c     ******************************************************

      PRINT *,' What name for FAULT_SLP_HIST.d output file?'
      READ(*,'(A30)') fnout
      PRINT *,' '

      OPEN(UNIT=30,NAME=fnout,STATUS='new')
      CLOSE(30)

C     ----------------------------------------------------------
C     CREATE FAULT_SLP_HIST.d FILE
C     ----------------------------------------------------------
      

C     SET NUMBER OF FAULTS
      numfs     = 4
C     SET SLIP INTERVAL IN UNITS OF YEARS (e.g. 1. = 1 year)
      tstepp    = 0.1
      itime     = 61

C     INITIAL YEAR TO START COUNTING 
      tminn     = 2000.

C     SET THE TIME OF THE SLIP FOR EACH EVENT IN DECIMAL YEARS
      tsuma     = 2004.9
      tnias     = 2005.3
C     INITIALIZE THE SLIP VECTORS AND THE TIMES OF SLIP
C     NOTE: EACH SEGMENT MUST SLIP ONCE. WE DEAL WITH THIS
C     BY LETTING THE SYSTEM SLIP BY A NEGLIABLE AMOUNT IN 
C     THE FIRST YEAR (tone) BY THE AMOUNT slip_neg IF NEEDED
c      slip_neg  = 5.D0

      do i=1,numfs
         islip(i) = 1
         slip1(i) = 0.
         slip2(i) = 0.
      end do

C     SET SLIP *MAGNITUDE* IN [cm]
C     IS DETERMINED IN THE CALCULATION OF THE GREEN's FUNCTIONS
C     RECALL WE MUST MULTIPLY DIPPING COMPONENTS OF THE SLIP VECTOR
C     BY -1.0 TO RECOVER THE PROPER
C     VECTOR FOR THRUST MOTION.

C     ONLY MAGNITUDES REQUIRED FOR STRIKE-SLIP MOTION
      slip1(1)  = 80.
      slip1(2)  = 170.
      slip1(3)  = 280.
      slip1(4)  = 0.

C     MUST INCLUDE THE MINUS SIGN FOR DIP-SLIP VALUES
      slip2(1)  = -590.
      slip2(2)  = -1190.
      slip2(3)  = -1980.
      slip2(4)  = -1100.

C     CREATE SLIP MATRIX FOR EACH FAULT
      PRINT *, 'INITIALIZING THE SLIP AND INITIATOR VECTORS'
      PRINT *, ' '
      k = 0
      DO i = 1,numfs
         k = k + 1                     
         DO j = 1,islip(k)
            IF (i .LT. 4) then 
               tslipp(i,j)  = tsuma
               initeq(i,j)  = 0
               IF (i .EQ. 3) initeq(i,j) = 1
            ELSE
               tslipp(i,j)  = tnias
               initeq(i,j)  = 1
            END IF
         END DO
      END DO

C     SET THE TOTAL SLIP VECTORS
C     INITIALIZE THE TOTAL SLIP VECTORS
      DO j = 1,itime
         DO i = 1,numfs
            dslpS(i,j) = 0.
            dslpD(i,j) = 0.
            dslpT(i,j) = 0.            
         END DO
      END DO

C     NOW SET UP THE TOTAL SLIP
      tmp = tminn
      crt = 0.05
      k = 1 
      DO j = 2,itime
         DO i = 1,numfs
            tst1 = ABS(tsuma - tmp)
            tst2 = ABS(tnias - tmp)
            IF (tst1 .LE. crt) THEN
               IF (i .LT. 4) THEN
                  dslpS(i,j) =  slip1(i)
                  dslpD(i,j) =  slip2(i)
               END IF
            END IF
            IF (tst2 .LE. crt) THEN
               IF (i .EQ. 4) THEN
                  dslpS(i,j) =  slip1(i)
                  dslpD(i,j) =  slip2(i)
               END IF
            END IF
         END DO
         tmp = tmp + tstepp
      END DO

C     WRITE OUT HISTORY FILE
      PRINT *, ' '
      PRINT *, 'NOW WRITING OUT FAULT SEGMENT HISTORY DATA TO: ', fnout
      PRINT *, ' '
      OPEN(UNIT=30,NAME=fnout,STATUS='old')            
      WRITE(30,'(A30)') fnout
      WRITE(30,*) tmp, tstepp, itime
      WRITE(30,*) (islip(n), n=1,numfs)
      DO i=1,numfs
         WRITE(30,*) (tslipp(i,k), k=1,islip(i))
      END DO      
      WRITE(30,*) ((dslpS(i,j), i=1,numfs), j=1,itime)      
      WRITE(30,*) ((dslpD(i,j), i=1,numfs), j=1,itime)
      WRITE(30,*) ((dslpT(i,j), i=1,numfs), j=1,itime)
      DO i=1,numfs
         WRITE(30,*) (initeq(i,k), k=1,islip(i))
      END DO      
      CLOSE(30)            

C     END OF MAIN PROGRAM
      STOP
      END
