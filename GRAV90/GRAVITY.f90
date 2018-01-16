! STARTED BY TYLER JOSEPH HAYES JULY 20, 2006...AND COUNTING
! 
! DESCRIPTION
! -----------
! GRAVITY.f90 TAKES THE OUTPUT OF THE CODE GRAV_GREEN.f90, TOGETHER
! WITH THE FAULT HISTORY DATA FILEE, AND COMPUTES THE
! GRAVITY VECTORS AT THE POINTS ON A SQUARE GRID CENTERED
! ON THE FAULT SYSTEM.
! 
! INPUTS
! ------
! FAULT_SLP_HIST.d   -- CONTAINS THE HISTORY OF THE SLIP ON EACH 
!                       SEGMENT.  THIS CAN BE GENERATED VIA 
!                       A DYNAMIC STRESS-EVOLUTION MODEL SUCH AS 
!                       Virtual California OR FROM A PRESCRIBED SUITE
!                       OF KNOWN SEGMENT SLIPS.
! GRAV_GREEN_DATA.d  -- THE EXACT NAME OF THE FILE WAS CHOSEN BY THE
!                       USER DURING RUN TIME FOR THE GREEN's FUNCTION
!                       LOOK UP TABLE. THIS FILE CONTAINS THE NAMES 
!                       OF THE "FAULT_SEG_INFO.d" FILE AND THE 
!                       "MEDIUM_INFO.d" FILE USED TO CREATE THE 
!                       LOOK-UP TABLE. AS SUCH, THESE FILE NEED TO
!                       BE IN THE CURRENT WORKING DIRECTORY AS WELL.
! MEDIUM_INFO.d      -- SEE ABOVE AND NOTES IN GRAV_GREEN.f90
! FAULT_SEG_INFO.d   -- SEE ABOVE AND NOTES IN GRAV_GREEN.f90
! 
! 
! SUBROUTINES
! -----------
! SUBROUTINE DSUM()  -- CALCULATES THE MULTIPLICATIVE SUMMED
!	 		TIME FACTORS FOR THE ENTIRE SLIP
!			HISTORY OF A GIVEN FAULT
! SUBROUTINE DEFOR() -- COMPUTES THE DISPLACEMENT AT THE TIME t
! 
! 
! MODULES
! -------
! DATA_DEFS    -- SETS THE ARRAY SIZES FOR THE MAIN PROGRAM
! CHAR_DEFS    -- SETS CHARACTER SIZES USED IN THE MAIN PROGRAM
! GRV_RD       -- CONTAINS THE GRAVITY VALUES READ IN AND CREATES A 
!                 COMMON STORAGE AREA FOR THEM
! GRV_WR       -- CONTAINS THE DEFINTIONS FOR THE OUTPUT GRAVITY VALUES
! SUBVALS      -- CREATES THE COMMON AREA FOR SEVERAL VALUES USED BY THE
!                 SUBROUTINES 
! SUBARRS      -- CREATES THE COMMON AREA FOR SEVERAL ARRAYS USED BY THE
!                 SUBROUTINES 
! REST_OF_DEFS -- DEFINES THE REST OF THE NECESSARY VARIABLES REQUIRED.
! 
! 
! OUTPUT FROM GRAVITY.f90 CAN THEN BE VISUALIZED USING GMT, IDL, MATLAB
! OR ANY OTHER SOFTWARE THE USER IS FAMILIAR WITH. 
!---------------------------------------------------------------------


! THIS IS THE MAIN MODULE WHICH DETERMINES THE ARRAY SIZES USED
! DATA DEFINTION MODULE ==> SHOULD BE THE SAME AS MAIN_PARS FOR EACH MODEL
! THAT WAS USED IN THE GREEN's FUNCTION CALCULATIONS WITH GRAV_GREEN.f90 
MODULE DATA_DEFS
  IMPLICIT NONE 
! THIS MODULE ASSIGNS THE SIZES OF THE SUBSEQUENT VARIABLES.  THIS
! IS THE FILE WHOSE VALUES CHANGE FOR DIFFERENT SEG. AND OBS. VALUES
! NPT   = max number of obsevation points
! NF    = number of fault segments used
! NTSTP = max number of time steps
  INTEGER,PARAMETER :: NPT=10000,NF=800,NPF=NPT*NF,NTSTP=11000 
END MODULE DATA_DEFS

! CAHRACTER DEFINTION MODULE
MODULE CHAR_DEFS
  IMPLICIT NONE 
! THIS MODULE DEFINES THE CHARACTER VARIABLES WHICH ARE EITHER
!      (1) A FILE NAME and/or,
!      (2) A USER RESPONSE VALUE
! fname_fault  = fault geometry file name 
! fhname       = fault history file name
! fdname       = fault segment location file name
! fname_gr     = gravity green's function output file
! fnstr        = garbage read-in value
! resp[XX] = user response variable
  CHARACTER (LEN=30) :: fdname,fnstr,fname
  CHARACTER (LEN=30) :: fhist_name,fmed_name,fgeo_name,fgreen_name
  CHARACTER (LEN=1) :: resp,resp1,respfb,respt,resp_start
END MODULE CHAR_DEFS


! SET UP THE "COMMON" AREA VARIABLES WITHIN MODULES TO AVOID ERRORS
! -------------------------------------------------------------------
! COMMON GRAVITY READ-IN VARIABLES
MODULE GRV_RD
  USE DATA_DEFS
  IMPLICIT NONE 
  REAL,DIMENSION(NPF) :: GZTES,GZDVES,GZFAES
  REAL,DIMENSION(NPF) :: GZTED,GZDVED,GZFAED
  REAL,DIMENSION(NPF) :: GZTET,GZDVET,GZFAET
  REAL,DIMENSION(NPT) :: GZTS,GZDVS,GZFAS
  COMMON/GRVRE/GZTES,GZDVES,GZFAES,&
       & GZTED,GZDVED,GZFAED,&
       & GZTET,GZDVET,GZFAET,&
       & GZTS,GZDVS,GZFAS
END MODULE GRV_RD

! COMMON MODULE FOR SUBROUTINE NON-ARRAY VALUES
MODULE SUBVALS
  IMPLICIT NONE 
  INTEGER :: ntime,npoint,nfault
  REAL :: tlater,tbase,tlast,tstep,timmin
  COMMON/SUBALL/ntime,npoint,tlater,tbase,tlast,tstep,timmin,nfault
END MODULE SUBVALS

! MODULE FOR COMMON SUBROUTINE ARRAYS
MODULE SUBARRS
  USE DATA_DEFS
  IMPLICIT NONE 
  REAL,DIMENSION(NF) :: islip,slpvlS,slpvlD,slpvlT
  REAL,DIMENSION(NF,NTSTP) :: dslpS,dslpD,dslpT
  COMMON/ARRSUB/islip,dslpS,dslpD,dslpT,slpvlS,slpvlD,slpvlT
END MODULE SUBARRS


! VARIABLES READ IN WHICH ARE NOT DEFINED IN THE ABOVE MAIN MODULES
! -------------------------------------------------------------------
! GRAVITY OUTPUT VARIABLES MODULE
MODULE GRV_WR
  USE DATA_DEFS
  IMPLICIT NONE 
  REAL,DIMENSION(NPT) :: GZT,GZDV,GZFA
END MODULE GRV_WR

! DATA NOT ALREADY DEFINED
MODULE REST_OF_DEFS
  USE DATA_DEFS
  IMPLICIT NONE 
  INTEGER :: ncycle,KMTST,ios
  INTEGER :: i,j,k,l,m,n,ig,ip,ix,mend,isum,isumm

! REAL VARIABLES
  REAL :: amu,alam,bulkm,vplatS,vplatD,vplatT,vplateX,vplateY,t_max
! REAL ARRAYS
  REAL,DIMENSION(NPT) :: x,y
  REAL,DIMENSION(NF) :: xe,ye,xw,yw,db,dt,segdp
  REAL,DIMENSION(NF,NTSTP) :: dslip,tslip
END MODULE REST_OF_DEFS


! MAIN PROGRAM
! -------------------------------------------------------------------
PROGRAM GRAVITY
  USE DATA_DEFS
  USE CHAR_DEFS
  USE REST_OF_DEFS
  USE SUBVALS
  USE SUBARRS
  USE GRV_RD
  USE GRV_WR
  IMPLICIT NONE 

! PRINT OUT REMINDER FOR FILE REQUIREMENTS IN THE DIRECTORY
  PRINT *, ' '
  PRINT *, '      ---->  THIS IS GRAVITY.f90  <----'
  PRINT *,' '
  PRINT *,'  AN INTERACTIVE PROGRAM FOR CALCUALTING GRAVITY'
  PRINT *,'  CHANGES OVER DIFFERENT TIME PERIODS.'
  PRINT *,' '
  PRINT *,' THIS PROGRAM REQUIRES THE FOLLOWING IN YOUR CURRENT'
  PRINT *,' WORKING DIRECTORY.'
  PRINT *,' ' 
  PRINT *,' (1) A *PRE-EXISTING* FILE CALLED gravity90.dat,'
  PRINT *,'     TO WRITE THE OUTPUT TO.'
  PRINT *,' (2) THE GRAVITY GREENs FUNCTIONS OUTPUT FROM'
  PRINT *,'     GRAV_GREEN.f90.'
  PRINT *,' (3) A SLIP HISTORY FILE, e.g., -> FAULT_SLP_HIST.d.'
  PRINT *,' (4) A FAULT GEOMETRY FILE, e.g., -> FAULT_SEG_INFO.d '
  PRINT *,' (5) A MEDIUM PARAMETRIZATION FILE, e.g., -> MEDIUM_INFO.d'
  PRINT *,' '

! UNCOMMENT BELOW IF YOU DO NOT WISH TO HAVE THIS WARNING
!   PRINT *,' MAKE SURE THESE FILE ARE IN YOUR WORKING DIRECTORY NOW AND'
!   PRINT *,' PRESS <ENTER> TO BEGIN.'
!   READ(*,'(A1)') resp_start

! QUERY FOR SLIP HISTORY FILE FROM USER
  DO 
!     PRINT *,' HEY! YOU THERE. PAY ATTENTION.'
     PRINT *,' ENTER THE NAME OF THE SLIP HISTORY FILE:'
     PRINT *, ' '
     READ(*,'(A30)') fhist_name
     OPEN(UNIT=10,FILE=fhist_name,STATUS="old",ACTION="read",IOSTAT=ios)
     IF (ios == 0) EXIT         ! EXIT LOOP IF OPEN IS SUCCESSFUL
! PRINT ERROR MESSAGE AND PROMPT AGAIN FOR NEW NAME
     PRINT *,' '
     PRINT *, 'ERROR OPENING FILE: ', fhist_name
     PRINT *,' '
     PRINT *, 'MAKE SURE THE FILE EXISTS IN YOUR CURRENT DIRECTORY'
     PRINT *, 'AND RE-ENTER THE NAME'
     PRINT *,' '
  END DO
  CLOSE (10)

! THE CURRENT INCARNATION HAS THESE TWO FILE NAMES READ IN FROM THE 
! GRAV_GREEN.f90 PROGRAM OUTPUT FILE
! QUERY FOR THE FAULT GEOMETRY FILE FROM USER
!  PRINT *,' '
!  PRINT *,' ENTER THE NAME OF THE FAULT GEOMETRY FILE:'
!  PRINT *, ' '
!  READ(*,'(A30)') fgeo_name

! QUERY FOR THE MEDIUM INFORMATION FILE FROM USER
!   PRINT *,' '
!   PRINT *,' ENTER THE NAME OF THE MEDIUM INFORMATION FILE:'
!   PRINT *, ' '
!   READ(*,'(A30)') fmed_name

! GET THE NAME OF THE GRAVITY GREEN's FUNCTION COEFFICIENTS FILE
  DO 
     PRINT *,' ENTER THE NAME OF THE INPUT FILE GRAVITY'
     PRINT *,' GREENS FUNCTIONS FOR OBSERVATION POINTS.'
     PRINT *, ' '
     READ(*,'(A30)') fgreen_name
     OPEN(UNIT=70,FILE=fgreen_name,STATUS="old",ACTION="read",IOSTAT=ios)
     IF (ios == 0) EXIT         ! EXIT LOOP IF OPEN IS SUCCESSFUL
! PRINT ERROR MESSAGE AND PROMPT AGAIN FOR NEW NAME
     PRINT *,' '
     PRINT *, 'ERROR OPENING FILE: ', fgreen_name
     PRINT *, 'MAKE SURE THE FILE EXISTS IN YOUR CURRENT DIRECTORY'
     PRINT *, 'AND RE-ENTER THE NAME'
     PRINT *,' '
  END DO

! READ THE OBSERVATION POINT GRAVITY GREEN's FUNCTION COEFFICIENTS
  READ(70,'(A30)') fgeo_name
  READ(70,'(A30)') fmed_name
  READ(70,*) nfault,npoint
  READ(70,*) (x(n), n=1,npoint)
  READ(70,*) (y(n), n=1,npoint)
  isum = nfault*npoint
  isumm=isum*6

! ELASTIC SLIP GRAVITY VALUES -> STRIKE SLIP
  READ(70,*) (GZTES(n), n=1,isum)
  READ(70,*) (GZDVES(n), n=1,isum)
  READ(70,*) (GZFAES(n), n=1,isum)

! ELASTIC SLIP GRAVITY VALUES -> DIP SLIP
  READ(70,*) (GZTED(n), n=1,isum)
  READ(70,*) (GZDVED(n), n=1,isum)
  READ(70,*) (GZFAED(n), n=1,isum)

! ELASTIC SLIP GRAVITY VALUES -> TENSILE MOTION
  READ(70,*) (GZTET(n), n=1,isum)
  READ(70,*) (GZDVET(n), n=1,isum)
  READ(70,*) (GZFAET(n), n=1,isum)

! THE STEADY STATE COMPONENTS
  READ(70,*) (GZTS(n), n=1,npoint)
  READ(70,*) (GZDVS(n), n=1,npoint)
  READ(70,*) (GZFAS(n), n=1,npoint)
  CLOSE(70)

! OPEN THE SEGMENT INFORMATION FILE --> FAULT_SEG_INFO.d
  OPEN(UNIT=40,FILE=fgeo_name,STATUS='old')
  READ(40,'(A30)') fname
  DO i=1,nfault
     READ(40,*) db(i),dt(i),xw(i),yw(i),xe(i),ye(i),&
          &  slpvlS(i),slpvlD(i),slpvlT(i),segdp(i)
  END DO
  CLOSE(40)

! OPEN AND READ MEDIUM_INFO.d FILE <- DON'T NEED
!  OPEN(UNIT=50,FILE=fmed_name,STATUS='old')      
!  READ(50,'(A30)') fname
!  READ(50,*) nfault,KMTST
!  READ(50,*) amu,alam,bulkm
!  READ(50,*) vplatS, vplatD, vplatT, vplateX, vplateY
! CLOSE FILE SINCE WE DON'T NEED TO READ IN THE REST OF THE VARIABLES
!  CLOSE(50)

! OPEN THE HISTORY_INFO.d FILE AND READ THE STUFF IN IT
  OPEN(UNIT=10,FILE=fhist_name,STATUS='old')
  READ(10,'(A30)') fnstr
  READ(10,*) tlast,tstep,ntime
  READ(10,*) (islip(n), n=1,nfault)
!   DO i=1,nfault
!      READ(10,*) (dslip(i,k), k=1,islip(i))
!   END DO
  DO i=1,nfault
     READ(10,*) (tslip(i,k), k=1,islip(i))
  END DO
  READ(10,*) ((dslpS(i,j), i=1,nfault), j=1,ntime)
  READ(10,*) ((dslpD(i,j), i=1,nfault), j=1,ntime)
  READ(10,*) ((dslpT(i,j), i=1,nfault), j=1,ntime)
  CLOSE(10)	

! State variables for each segment:
! islip(i): # of times segment i has slipped
! dslip(i,j): amount of slip of segment i when it slipped for the jth time
! tslip(i,j): time at which segment i slipped for the jth time 


! ENTER MAIN LOOP
  resp='y'
  DO WHILE(resp .EQ. 'y')

! GET THE TIMES FROM COMMON AREA
     t_max  = tlast
     timmin = t_max - float(ntime - 1)*tstep

! QUERY FOR TIME INTERVAL OF INTEREST FROM THE USER
     PRINT *, ' '
     PRINT  '(" MINIMUM time (yr): ",F7.2)', timmin
     PRINT *, ' '
     PRINT  '(" MAXIMUM time (yr): ",F7.2)', t_max
     PRINT *,' '
     PRINT *,' We need two dates, a later one'
     PRINT *,' and an earlier one to calculate'
     PRINT *,' the differential offsets.'
     PRINT *,' '
     PRINT *,' Enter EARLIER date (years)'
     READ *, tbase
     PRINT *,' '
     PRINT *,' Enter LATER date (years)'
     READ *, tlater
     
! CALCULATE THE GRAVITY VALUES FOR EACH OBS. POINT USING DEFOR
! AND STORE THE OUTPUT INTO THE COMMON AREA VARIABLES GZTp,GZDVp,GZFAp
     CALL DEFOR(GZT,GZDV,GZFA)

! WRITE OUT THE GRAVITY FIELD DATA FOR THE REQUESTED TIME PERIOD TO 
! gravity_data.d AND QUERY FOR ANOTHER CALCULATION. THIS IS DONE TO 
! FACILITATE INTERACTIVE PLOTTING.
     OPEN(UNIT=2,FILE='gravity90.dat',STATUS='old')
     WRITE(2,*) fmed_name
     WRITE(2,*) fhist_name
     WRITE(2,*) fgeo_name
     WRITE(2,*) tbase, tlater, tlast, npoint
     WRITE(2,*) (x(j), j=1,npoint)
     WRITE(2,*) (y(j), j=1,npoint)
     WRITE(2,*) (GZT(j), j=1,npoint)
     WRITE(2,*) (GZDV(j), j=1,npoint)
     WRITE(2,*) (GZFA(j), j=1,npoint)
     CLOSE(2)
     PRINT *, 'Compute another gravity field? (y/n)'
     READ(*,'(A1)') resp  
  END DO
END PROGRAM GRAVITY


! SUBROUTINES
! ----------------------------------------------------------------

SUBROUTINE DEFOR(dGZT,dGZDV,dGZFA)
  USE DATA_DEFS
  USE GRV_RD
  USE SUBVALS
  IMPLICIT NONE 
  INTEGER :: l,k,ig,n,i,ip
  REAL,DIMENSION(NPT),INTENT(OUT):: dGZT,dGZDV,dGZFA
  REAL,DIMENSION(NF) :: deltotS,deltotD,deltotT

!
! INTERFACE:
! dispx, dispy: x and y components of the displacement of a point
!               of the grid
! nfault: number of fault segments
! other input parameters passed through COMMON statement
!
! RETURN VALUE:
! dispx, dispy for all points on the grid
!
! DESCRIPTION:
! This subroutine computes the displacement at the time t
!

! COMPUTE TOTAL SLIP FOR THE REQUESTED TIME PERIOD
  CALL DSUM(deltotS,deltotD,deltotT)

! FOR DEBUGGING ONLY
!   DO i=1,nfault
!      PRINT *, deltotS(i)
!   END DO

! INITIALIZE
  DO i=1,npoint
     dGZT(i)  = 0.
     dGZDV(i) = 0.
     dGZFA(i) = 0.
  END DO

! FOR EACH SEGMENT FIND ITS CONTRIBUTION FROM STRIKE SLIP MOTION
  l=0
  k=0
  DO ig=1,nfault
     DO n=1,npoint
        l=l + 1
        dGZT(n)  = dGZT(n)  +  GZTES(l)*deltotS(ig)
        dGZDV(n) = dGZDV(n) + GZDVES(l)*deltotS(ig)
        dGZFA(n) = dGZFA(n) + GZFAES(l)*deltotS(ig)
     END DO
  END DO

! FOR EACH SEGMENT ADD ITS CONTRIBUTION FROM DIP SLIP MOTION
  l=0
  k=0
  DO ig=1,nfault
     DO n=1,npoint
        l=l + 1
        dGZT(n)  = dGZT(n)  +  GZTED(l)*deltotD(ig)
        dGZDV(n) = dGZDV(n) + GZDVED(l)*deltotD(ig)
        dGZFA(n) = dGZFA(n) + GZFAED(l)*deltotD(ig)
     END DO
  END DO

! FOR EACH SEGMENT ADD ITS CONTRIBUTION FROM TENSILE MOTION
  l=0
  k=0
  DO ig=1,nfault
     DO n=1,npoint
        l=l + 1
        dGZT(n)  = dGZT(n)  +  GZTET(l)*deltotT(ig)
        dGZDV(n) = dGZDV(n) + GZDVET(l)*deltotT(ig)
        dGZFA(n) = dGZFA(n) + GZFAET(l)*deltotT(ig)
     END DO
  END DO

! ADD THE STEADY STATE COMPONENT
  DO ip=1,npoint
!       dGZT(ip)  = dGZT(ip)  + 0.
!       dGZDV(ip) = dGZDV(ip) + 0.
!       dGZFA(ip) = dGZFA(ip) + 0.
     dGZT(ip)  = dGZT(ip)  +  GZTS(ip)*(tlater-tbase)
     dGZDV(ip) = dGZDV(ip) + GZDVS(ip)*(tlater-tbase)
     dGZFA(ip) = dGZFA(ip) + GZFAS(ip)*(tlater-tbase)
  END DO
END SUBROUTINE DEFOR


SUBROUTINE DSUM(deltotSS,deltotDD,deltotTT)
  USE DATA_DEFS
  USE CHAR_DEFS
  USE SUBVALS
  USE SUBARRS
  IMPLICIT NONE 
  INTEGER :: ix,j,jend 
  REAL :: time,tbarb,tbarl
  REAL,DIMENSION(NF),INTENT(OUT) :: deltotSS,deltotDD,deltotTT
  REAL,DIMENSION(NF) :: dtotlSS,dtotlDD,dtotlTT,dtotbSS,dtotbDD,dtotbTT
! INTERFACE:
! PARAMETERS ARE PASSED THROUGH THE COMMON AREA
!
! RETURN VALUE:
! deltot[SS,DD,TT](i): DIFFERENTIAL SLIP DEFICIT PART OF SEGEMENT i

  DO ix=1,nfault
! INITIALIZE EACH VECTOR FOR EVERY NEW SEGEMNT
     deltotSS(ix) = 0.
     deltotDD(ix) = 0.
     deltotTT(ix) = 0.
     dtotbSS(ix)  = 0.
     dtotlSS(ix)  = 0.
     dtotbDD(ix)  = 0.
     dtotlDD(ix)  = 0.
     dtotbTT(ix)  = 0.
     dtotlTT(ix)  = 0.

! SEE IF SEGMENT SLIPPED VIA jend
     jend=islip(ix)
     respt='s'
     time = timmin - tstep
     j=0

! SKIP WHILE LOOP IF SEGMENT ix DID NOT RECORD ANY SLIPS
     IF (jend .GT. 0) respt='g'     

! IF SEGMENT RECORDED SOME SLIP, SEE HOW MUCH OCCURED DURING 
! THE USER DEFINED TIME PERIOD
     DO WHILE(respt .EQ. 'g')
        j=j+1
        time  = time + tstep
        tbarl = (tlater - time + .001*tstep)
        tbarb = (tbase  - time + .001*tstep)
!        PRINT *, tbarl, tbarb
        IF (tbarl .LT. 0.) respt='s'
        IF (tbarl .GE. 0.) THEN     
           dtotlSS(ix) = dtotlSS(ix) + dslpS(ix,j)     
           dtotlDD(ix) = dtotlDD(ix) + dslpD(ix,j)     
           dtotlTT(ix) = dtotlTT(ix) + dslpT(ix,j)           
           IF (tbarb .GT. 0.) THEN
              dtotbSS(ix) = dtotbSS(ix) + dslpS(ix,j)
              dtotbDD(ix) = dtotbDD(ix) + dslpD(ix,j)
              dtotbTT(ix) = dtotbTT(ix) + dslpT(ix,j)
           END IF
        END IF
     END DO

! TOTAL TIME-INDEPENDENT DIFFERENTIAL SLIP DEFICIT
      deltotSS(ix) = dtotlSS(ix) - dtotbSS(ix)
      deltotDD(ix) = dtotlDD(ix) - dtotbDD(ix)
      deltotTT(ix) = dtotlTT(ix) - dtotbTT(ix)
  END DO  
END SUBROUTINE DSUM

