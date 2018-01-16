      PROGRAM MAKE_HIST_BEAR

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
      numfs     = 10
C     SET SLIP INTERVAL IN UNITS OF YEARS (e.g. 1. = 1 year)
      tstepp    = 0.1000000000
C     81 TIME STEPS -> 8 YEARS IN 0.1 YEAR INTERVALS
C     THIS WILL TAKES US FROM 1992.0 TO 2000.0 -> HENCE 81
      itime     = 81
C     INITIAL YEAR TO START COUNTING 
      tminn     = 1992.D0

C     SET THE TIME OF THE SLIP FOR EACH EVENT IN DECIMAL YEARS
      tjosh     = 1992.3D0
      tbigbear  = 1992.5D0
      thector   = 1999.8D0

C     INITIALIZE THE SLIP VECTORS AND THE TIMES OF SLIP
C     JOSHUA TREE -> SET NUMBER TO AMOUNT OF TIMES SEGMENT SLIPPED
      islip(1)  = 1
C     BIG BEAR
      islip(2)  = 1
C     LANDERS
      islip(3)  = 1
      islip(4)  = 1
      islip(5)  = 1
      islip(6)  = 1
      islip(7)  = 1
C     HECTOR MINE
      islip(8)  = 1
      islip(9)  = 1
      islip(10) = 1

C     SET THE *MAGNITUDE* OF SLIP IN [cm] -> SENSE OF SLIP ALREADY 
C     IS DETERMINED IN THE CALCULATION OF THE GREEN's FUNCTIONS
C     JOSHUA TREE / EUREKA --- Bennett et al. JGR [1995]
      slip(1)  = 80.      
C     BIG BEAR --- Hudnut BSSA [1994]
      slip(2)  = 44.
C     LANDERS --- Hudnut BSSA [1994] averaged values
      slip(3)  = 157. 
      slip(4)  = 610. 
      slip(5)  = 284.
      slip(6)  = 233.
      slip(7)  = 298.
C     HECTOR MINE --- Simons et al. BSSA [2002] averaged values
      slip(8)  = 480.
      slip(9)  = 350.
      slip(10) = 150.



C     CREATE SLIP MATRIX FOR EACH FAULT
      PRINT *, 'INITIALIZING THE SLIP AND INITIATOR VECTORS'
      PRINT *, ' '
      k = 0
      DO i = 1,numfs
C         IF (islip(i) .NE. 0) THEN
         k = k + 1
            DO j = 1,islip(k)
C               tslipp(i,j) = 0.
               IF (i .EQ. 1) then 
                  tslipp(i,j)  = tjosh
                  initeq(i,j)  = 1
                  PRINT*, 'JOSHUA  :',  tjosh, slip(i)
               ELSE IF (i .EQ. 2) THEN
                  tslipp(i,j) = tbigbear
                  initeq(i,j) = 0
                  PRINT*, 'BIG BEAR  :',  tbigbear, slip(i)
               ELSE IF ((i .GT.2) .AND. (i .LT. 8)) THEN
                  tslipp(i,j) = tbigbear
                  initeq(i,j) = 0
                  IF (i .EQ. 7) initeq(i,j) = 1
                  PRINT*, 'LANDERS :',  tbigbear, slip(i)
               ELSE IF (i .GT. 7) THEN
                  tslipp(i,j) = thector
                  initeq(i,j) = 0
                  IF (i .EQ. 9) initeq(i,j) = 1
                  PRINT*, 'HECTOR  :',  thector, slip(i)
               END IF
            END DO
c         END IF
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
         tmp1 = abs(tmp - tjosh)
         tmp2 = abs(tmp - tbigbear)
         tmp3 = abs(tmp - thector)  
         DO i = 1,numfs

C     JOSHUA TREE CASE
            IF ((tmp1 .LT. crt) .AND. (i .EQ. 1)) THEN
c               dslpS(i,j) = dslpS(i,j-1) + slip(i)
               dslpS(i,j) =  slip(i)

C     LANDERS / BIG BEAR CASE
            ELSE IF ((tmp2 .LT. crt) .AND. (i .GT. 1)) THEN
c               IF (i .LT. 9) dslpS(i,j) = dslpS(i,j-1) + slip(i)
               IF (i .LT. 8) dslpS(i,j) =  slip(i)

C     HECTOR MOINE CASE
            ELSE IF ((tmp3 .LT. crt) .AND. (i .GT. 7)) THEN
c               dslpS(i,j) = dslpS(i,j-1) + slip(i)
               dslpS(i,j) = slip(i)

C     ALL OTHER CASES
            ELSE
               dslpS(i,j) = dslpS(i,j)
            END IF
c            PRINT *, j, tmp , dslpS(i,j)
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
