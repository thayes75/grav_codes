      PROGRAM MAKE_HIST_ALASKA

      parameter (NF=800, NTSTP=11000)
      character fnout*30
      dimension tslipp(NF,NTSTP)
      dimension islip(NF),slip(NF)
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
      numfs     = 29
C     SET SLIP INTERVAL IN UNITS OF YEARS (e.g. 1. = 1 year)
      tstepp    = 1.
C     81 TIME STEPS -> 8 YEARS IN 0.1 YEAR INTERVALS
C     THIS WILL TAKES US FROM 1992.0 TO 2000.0 -> HENCE 81
      itime     = 51
C     INITIAL YEAR TO START COUNTING 
      tminn     = 1934.D0

C     SET THE TIME OF THE SLIP FOR EACH EVENT IN DECIMAL YEARS
      tone      = 1935.D0
      talask    = 1964.D0
      slip_neg  = 5.D0
C     INITIALIZE THE SLIP VECTORS AND THE TIMES OF SLIP
C     NOTE: EACG SEGMENT MUST SLIP ONCE. WE DEAL WITH THIS
C     BY LETTING THE SYSTEM SLIP BY A NEGLIABLE AMOUNT IN 
C     THE FIRST YEAR (tone) BY THE AMOUNT slip_neg
      do i=1,numfs
         islip(i) = 1
         slip(i)  = 0.
      end do

C     SET NUMBER TO AMOUNT OF TIMES SEGMENT SLIPPED
      islip(12)  = 2
      islip(13)  = 2
      islip(14)  = 2
      islip(15)  = 2
      islip(16)  = 2
      islip(17)  = 2
      islip(18)  = 2
      islip(19)  = 2
      islip(20)  = 2

C     SET ALASKA SLIP *MAGNITUDE* IN [cm]
C     IS DETERMINED IN THE CALCULATION OF THE GREEN's FUNCTIONS
C     RECALL WE MUST MULTIPLY BY -1.0 TO RECOVER THE PROPER
C     VECTOR FOR THRUST MOTION.
      slip(12)  = -1000.
      slip(13)  = -1200.
      slip(14)  = -1100.
      slip(15)  = -1000.
      slip(16)  = -1000.
      slip(17)  = -900.
      slip(18)  = -900.
      slip(19)  = -800.
      slip(20)  = -800.

C     CREATE SLIP MATRIX FOR EACH FAULT
      PRINT *, 'INITIALIZING THE SLIP AND INITIATOR VECTORS'
      PRINT *, ' '
      k = 0
      DO i = 1,numfs
         k = k + 1                     
         DO j = 1,islip(k)
            tslipp(i,j)  = tone
            initeq(i,j)  = 0
            IF (j .EQ. 2) then 
               tslipp(i,j)  = talask
               initeq(i,j)  = 0
            END IF
            IF ((i .EQ. 13) .and. (j .eq. 2)) initeq(i,j) = 1
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
            IF (tmp .EQ. tone) then
               if (i .lt. 12) THEN
                  dslpS(i,j) = slip_neg
               else
                  dslpD(i,j) = slip_neg
               end if
            end if
            if (tmp .eq. talask) then
               IF ((i .GT. 11) .AND. (i .LT. 21)) THEN
                  dslpD(i,j) =  slip(i) + slip_neg
               END IF
            end if
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
