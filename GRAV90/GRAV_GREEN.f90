! STARTED BY TYLER JOSEPH HAYES JULY 20, 2006
! 
! DESCRIPTION
! -----------     
! THIS FORTRAN PROGRAM TAKES THE FAULT GEOMETRY DATA FROM A DATA
! THE FILE, MAKES THE NECESSARY COORDINATE TRANSFORMATIONS, AND THEN
! CALCULATES THE GRAVITY GREENS FUNCTIONS FOR A GRID OF 
! SURFACE OBSERVATION POINTS ABOUT THE FAULT SYSTEM.  
!
! INPUTS
! ------
! THE MAJOR FILES YOU REQUIRE TO HAVE IN YOUR WORKING DIRECTORY ARE 
! AS FOLLOWS:
! 
! FAULT_SEG_INFO.d   -- PROVIDES THE FAULT TOPOGRAPHY INFORMATION
!                       SUCH AS THE FAULT LENGTHS, ORIENTATION, 
!                       VELOCITY, & NUMBER OF SEGMENTS IN THE MODEL 
! MEDIUM_INFO.d      -- CONTAINS THE NECESSARY SUBSURFACE MEDIUM 
!                       PARAMETERS AND THE OBSERVATION GRID PARAMETERS
! 
! NOTE -- THE ABOVE ARE GENERIC NAMES FOR THE FILES, AND IT IS 
!         ASSUMED THAT YOU HAVE USED MORE MEANINGFUL NAMES. FOR EXAMPLES
!         PLEASE SEE THE SAMPLE FILES PROVIDED.
! 
! OUTPUTS
! -------
! THIS CODE OUTPUTS A GRAVITY GREEN'S FUNCTION LOOK UP TABLE OF
! THE GRAVITY VALUES FOR EACH POINT ON AN OBSERVATION GRID DUE TO
! UNIT [cm] SLIP. THE NAME OF THE DATA FILE WITH THESE VALUES IS 
! CHOSEN BY THE USER AND ONLY NEEDS TO BE CALCULATED ONCE FOR EACH 
! FAULT TOPOGRAPHY INPUT.
!                        
!     
! 
! MODULES
! -------
! MAIN_PARS  -- SETS THE ARRAY SIZES FOR THE MAIN PROGRAM
! GRAV_PARS  -- SETS CONSTANTS AND ARRAY SIZES USED 
!               FOR THE GRV_MAIN SUBROUTINE
! CONST      -- CONTAINS CONVERSION PARAMETERS USED THROUGHOUT 
! ELAS       -- CONTAINS THE "COMMON" AREA MEDIUM PARAMETERS
!               FOR THE GRV_MAIN SUBROUTINE
! 
!
! SUBROUTINES
! -----------
! SUBROUTINE GRV_MAIN()	-- CALLS THE GRAVITY GREEN'S 
!                                      FUNCTION COEFFICIENT CODES
!                                      OF OKUBO (1992)              
! SUBROUTINE GRVFLT()	-- CALCULATES THE GRAVITY
!                                      GREEN'S FUNCTION USING 
!                                      OKUBO (1992)
! SUBROUTINE ELVFLT()	-- CALCULATES THE VERTICAL 
!                                      DEFORMATION DUE THE SLIP 
!                                      ON THE FAULT --> USED FOR 
!                                      FREE-AIR CALCULATIONS
! SUBROUTINE CONVRT()   -- CONVERTS LAT/LON INTO KMs FROM A REFERENCE PT.  
!
!
! FUNCTIONS
! ---------
! FUNCTION SGN()  -- SIGNUM FUNCTION (X < 0, OR X > 0)
! FUNCTION HFLT   -- CALCUALTE UPLIFT DUE TO FAULTING
! FUNCTION GFLT0  -- CALCULATE GRAVITY CHANGE: STRIKE-SLIP 
! FUNCTION GFLT1  -- CALCULATE GRAVITY CHANGE: DIP-SLIP & TENSILE 
! FUNCTION GFLT2  -- CALCULATE GRAVITY CHANGE: CAVITY FORMATION
! 
!
! THE OUTPUT FILE FROM THIS CODE IS USED, TOGETHER
! WITH A SLIP HISTORY FILE (FAULT_SLP_HIST.d) IN GRAVITY.f90.
! THAT CODE THEN COMPUTES THE ACTUAL GRAVITY DIFFERENCE
! BETWEEN ANY TWO TIMES SELECTED BY THE USER.
! 
! OUTPUT FROM GRAVITY.f90 CAN THEN BE VISUALIZED USING GMT, IDL, MATLAB
! OR ANY OTHER SOFTWARE THE USER IS FAMILIAR WITH. AN IDL ROUTINE IS
! PROVIDED FOR THE ALASKA AND JOSHUA TREE-BIG BEAR-HECTOR MINE SEQUENCES.
! THE IDL FILE CAN BE EASILY HACKED TO GENERALIZE FOR YOUR OWN DATA SETS
! OR USED AS A PSEUDO-CODE TO CREATE A MATLAB, GMT, OR "R" VISUALIZATION
! CODE.
! ------------------------------------------------------------------------


! THIS IS THE MAIN MODULE WHICH DETERMINES THE ARRAY SIZES USED
! DATA DEFINTION MODULE
MODULE MAIN_PARS
  IMPLICIT NONE 
  SAVE
! THIS MODULE ASSIGNS THE SIZES OF THE SUBSEQUENT VARIABLES.  THIS
! IS THE FILE WHOSE VALUES CHANGE FOR DIFFERENT SEG. AND OBS. VALUES
! NPT   = max number of obsevation points
! NF    = number of fault segments used
! NTSTP = max number of time steps
  INTEGER,PARAMETER :: NPT=10000,NF=800,NPF=NPT*NF,NTSTP=11000
END MODULE MAIN_PARS

! GRAVITY SUBROUTINE PARAMETER DEFINTIONS
MODULE GRAV_PARS
  IMPLICIT NONE 
  INTEGER,PARAMETER :: NGR=200,NGR2=NGR**2
  REAL,PARAMETER :: GNEWTN=6.67D-11,GAMMA=3.086D-6,GCONVR=1.D2*1.D6
END MODULE GRAV_PARS

! CONSTANT PARAMETERS MODULE USED THROUGHOUT
MODULE CONST
  IMPLICIT NONE 
  REAL,PARAMETER :: M2KM=1.D-3,CM2M=1.D-2,KM2M=1.D3
  REAL,PARAMETER :: RADIUS=6371.0
  REAL,PARAMETER :: PI=3.14159265358979323846,RAD=PI/180.D0
END MODULE CONST         

! COMMON ELASTIC MEDIUM PARAMETERS USED FOR GRAVITY CALCULATIONS
MODULE ELAS
  IMPLICIT NONE 
  REAL :: XNU,MUETA
  COMMON/ELASTC/XNU,MUETA
END MODULE ELAS


! MAIN PROGRAM MODULE
! --------------------------------------------------------------------
MODULE MAIN_MOD
  USE MAIN_PARS
  USE CONST
  IMPLICIT NONE 
! FILE NAMES AND VARIABLES
  CHARACTER(LEN=30) :: fgeo_name,fmed_name,fname
  CHARACTER(LEN=1)  :: b1,b2

! INTEGER NAMES
  INTEGER :: i,nfault,KMTST,ntime,inorth,isouth,inew,ios,ik
  INTEGER :: indn,indnm,npoint,IPTS,JPTS,l,j,ITST,n,np,isum,nflt
  INTEGER :: ibcS,ibcD,ibcT

! REAL VARIABLES
  REAL :: vplatS,vplatD,vplatT,REFLON,REFLAT,vplateX,vplateY
  REAL :: tlast,tstepp,bulkm,AINCX,AINCY,TARG
  REAL :: delx,dely,ZNU,SGNDP,dip1,dip2,dbplat
  REAL :: xfe,yfe,xfw,yfw,dt,db,testy,xmax,xmin,xmax1,xmin1
  REAL :: ymax,ymin,ymax1,ymin1,xmidf,ymidf,XSTRT,YSTRT,OFF
  REAL :: st,ct,SGNW,SGNE,STEST,ATEST,xmid,ymid,sl,FLEN,amu,alam
  REAL :: USTR,UDIP,UTNS,xnorth,ynorth,xsouth,ysouth,vplateAVG
  REAL :: UDS,UDD,UDT,udsV,uddV,udtV,dipp,dips,SYSLAT
  REAL :: vS1,vD1,vT1,vS2,vD2,vT2,VTEST,dipcor

! ARRAYS
  REAL,DIMENSION(NF) :: dbn,dtn,xfwnt,yfwnt,xfent,yfent
  REAL,DIMENSION(NF) :: slpvlS,slpvlD,slpvlT,segdp1,segdp

! USED WITHIN THE MAIN PROGRAM
  REAL,DIMENSION(NF) :: xfwn,yfwn,xfen,yfen,ymidn,xmidn,thetan
  REAL,DIMENSION(NF) :: xmidt,ymidt,SUBPI
  REAL,DIMENSION(NPT) :: xpt,ypt,x,y,xxpp,yypp

! GRAVITY MATRICES
  REAL,DIMENSION(NPF) :: GZTES,GZDVES,GZFAES
  REAL,DIMENSION(NPF) :: GZTES1,GZDVES1,GZFAES1
  REAL,DIMENSION(NPF) :: GZTES2,GZDVES2,GZFAES2

  REAL,DIMENSION(NPF) :: GZTED,GZDVED,GZFAED
  REAL,DIMENSION(NPF) :: GZTED1,GZDVED1,GZFAED1
  REAL,DIMENSION(NPF) :: GZTED2,GZDVED2,GZFAED2

  REAL,DIMENSION(NPF) :: GZTET,GZDVET,GZFAET
  REAL,DIMENSION(NPF) :: GZTET1,GZDVET1,GZFAET1
  REAL,DIMENSION(NPF) :: GZTET2,GZDVET2,GZFAET2

  REAL,DIMENSION(NPT) :: GZTS,GZDVS,GZFAS
  REAL,DIMENSION(NPT) :: GTS,GDVS,GFAS

! EXTERNAL FUNCTION DECLARATIONS
  REAL,EXTERNAL :: SGN
END MODULE MAIN_MOD



! START THE MAIN PROGRAM GRAV_GREEN.f90
PROGRAM GRAV_GREEN
  USE MAIN_MOD
  IMPLICIT NONE 

! INTRODUCTORY PROGRAM SCREEN INFORMATION
  PRINT *,' '
  PRINT *,'              --===<  GRAV_GREEN.f90  >===--'
  PRINT *,' '
  PRINT *,'    ====> THIS IS THE GOOD FRIDAY VERSION (2006) <===='
  PRINT *,' '
  PRINT *,' THIS IS THE GRAVITY GREENS FUNCTION CODE FOR '
  PRINT *,' CREATING THE GRAVITY GREENs FUNCTIONS TO BE USED '
  PRINT *,' BY GRAVITY.f90'      
  PRINT *,' '
  PRINT *,' PLEASE MAKE NOTE OF THE FOLLOWING:'
  PRINT *,' '
  PRINT *,' (1) THIS VERSION HAS POSTIVE Y-AXIS AS: **NORTH** '
  PRINT *,' (2) THIS VERSION DEFINES POSITIVE SLIP AS: **LEFT LATERAL** '
  PRINT *,' (3) CONTAINS A CUSTOMIZABLE OBSERVATION GRID '
  PRINT *,' (4) ELASTIC GREENs FUNCTION OUTPUT VALUES:'
  PRINT *, '    ==> PRODUCES OUTPUT IN UNITS OF [microGals/cm]'
  PRINT *,' (5) STEADY STATE GREENs FUNCTION OUTPUT VALUES:'
  PRINT *,'     ==> PRODUCES OUTPUT IN UNITS OF [microGals/yr]'
  PRINT *,' '
  PRINT *,'                -----====<<<*>>>====------     '
  PRINT *,' ' 

! (1) ENTER AND READ IN DATA FROM FAULT/MEDIUM INFO FILES
! ----------------------------------------------------------
! GET NAME OF MEDIUM INFO FILE AND OPEN
  DO 
     PRINT *,' '
     PRINT *,' ENTER THE NAME OF THE MEDIUM INFORMATION FILE:'
     READ(*,'(A30)') fmed_name
     OPEN(UNIT=50,FILE=fmed_name,STATUS="old",ACTION="read",IOSTAT=ios)
     IF (ios == 0) EXIT         ! EXIT LOOP IF OPEN IS SUCCESSFUL
! PRINT ERROR MESSAGE AND PROMPT AGAIN FOR NEW NAME
     PRINT *, ' '
     PRINT *, 'ERROR OPENING FILE: ', fmed_name
     PRINT *, 'MAKE SURE THE FILE EXISTS IN YOUR CURRENT DIRECTORY'
     PRINT *, 'AND RE-ENTER THE NAME'
  END DO
! SEE MEDIUM_INFO.d FOR VARIABLE MEANINGS AND AN EXAMPLE FILE
! READ DATA FROM MEDIUM INFO FILE
  READ(50,'(A30)') fname
  READ(50,*) amu,alam,bulkm
  READ(50,*) vplatS, vplatD, vplatT, vplateX, vplateY
  READ(50,*) ibcS, ibcD, ibcT, b1, b2
  READ(50,*) dip1, dip2, dbplat
  READ(50,*) nfault,KMTST
  READ(50,*) REFLON, REFLAT
  READ(50,*) npoint,IPTS,JPTS,XSTRT,YSTRT
  CLOSE(50)

! GET THE FAULT GEOMETRY INFORMATION FILE NAME AND OPEN FILE
  DO 
     PRINT *, ' '
     PRINT *,' ENTER THE NAME OF THE SEGMENT GEOMETRY INFORMATION FILE:'
     READ(*,'(A30)') fgeo_name
     OPEN(UNIT=20,FILE=fgeo_name,STATUS="old",ACTION="read",IOSTAT=ios)
     IF (ios == 0) EXIT         ! EXIT LOOP IF OPEN IS SUCCESSFUL
! PRINT ERROR MESSAGE AND PROMPT AGAIN FOR NEW NAME
     PRINT *, ' '
     PRINT *, 'ERROR OPENING FILE: ', fgeo_name
     PRINT *, 'MAKE SURE THE FILE EXISTS IN YOUR CURRENT DIRECTORY'
     PRINT *, 'AND RE-ENTER THE FILE NAME'
  END DO
! SEE FAULT_SEG_INFO.d FOR VARIABLE MEANINGS AND AN EXAMPLE FILE
  READ(20,'(A30)') fname
  DO i=1,nfault
     READ(20,*) dbn(i),dtn(i),xfwnt(i),yfwnt(i),xfent(i),yfent(i),&
          &        slpvlS(i),slpvlD(i),slpvlT(i),segdp1(i)
  END DO
  
! OPEN THE OUTPUT FILE TO WRITE THE UNIT GREEN'S FUNCTIONS 
  DO 
     PRINT *,' '
     PRINT *,' ENTER THE NAME FOR THE OUTPUT FILE:'
     READ(*,'(A30)') fname
     OPEN(UNIT=90,FILE=fname,STATUS="new",ACTION="write",IOSTAT=ios)
     IF (ios == 0) EXIT         ! EXIT LOOP IF OPEN IS SUCCESSFUL
! PRINT ERROR MESSAGE AND PROMPT AGAIN FOR NEW NAME
     PRINT *, ' '
     PRINT *, 'ERROR OPENING FILE: ', fname
     PRINT *, 'MAKE SURE THE FILE DOES *NOT* ALREADY EXIST IN YOUR'
     PRINT *, 'CURRENT DIRECTORY AND RE-ENTER A NEW OUTPUT FILE NAME'
  END DO

! (2) CREATE THE OBSERVATION GRID FOR THE GREEN'S FUNCTION 
!     AND START THE FIRST LOOP FOR THE ELASTIC DEFORMATION
!     CALCULATION
! --------------------------------------------------------
! TEST FOR LAT/LON OR KM INPUT USING "KMTST"
  
! CONVERT IF NEED BE, OTHERWISE SET TO INPUT VALUES
  IF (KMTST == 1) THEN
     PRINT *, 'CONVERTING LATITUDE AND LONGITUDE VALUES TO [km]'
     DO i=1,nfault
        xfwn(i) = 0.
        yfwn(i) = 0.
        xfen(i) = 0.
        yfen(i) = 0.
     END DO
     CALL CONVRT(xfwn,yfwn,xfwnt,yfwnt,REFLON,REFLAT,nfault)
     CALL CONVRT(xfen,yfen,xfent,yfent,REFLON,REFLAT,nfault)
      PRINT *,' '
!      PRINT *, 'THESE ARE THE CONVERTED COORDINATES'
!      DO i =1,nfault
!         PRINT *, i,xfwn(i),yfwn(i),xfen(i),yfen(i)
!      END DO 
  ELSE
     PRINT *, 'COORDINATES ALREADY IN UNITS OF [km]...CONTINUING'
     DO i=1,nfault
        xfwn(i) = xfwnt(i)
        yfwn(i) = yfwnt(i)
        xfen(i) = xfent(i)      
        yfen(i) = yfent(i) 
     END DO 
  END IF

! EXAMINE THE BOUNDARY CONDITIONS TO SEE IF THEY ARE MIXED, i.e., 
! ONE BOUNDARY IS STRIKE-SLIP AND THE OTHER IS DIP-SLIP, SUCH AS IN 
! THE ALASKA EXAMPLE
! RECALL INTEGER FLAG MEANINGS:
!   1 = ONLY NORTH-WESTERN MOST BOUNDARY HAS THIS COMPONENT
!   2 = ONLY SOUTH-EASTERN MOST BOUNDARY HAS THIS COMPONENT
!   3 = BOTH BOUNDARIES HAVE THIS COMPONENT (DEFAULT)

! THE STRIKE-SLIP VECTOR BOUNDARY CONDITION 
  SELECT CASE (ibcS)
  CASE (1)
     vS1 = vplatS
     vS2 = 0.
  CASE (2)
     vS1 = 0.
     vS2 = vplatS
  CASE (3)
     vS1 = vplatS
     vS2 = vplatS
  CASE DEFAULT
     vS1 = vplatS
     vS2 = vplatS
  END SELECT

! THE DIP-SLIP VECTOR BOUNDARY CONDITION 
  SELECT CASE (ibcD)
  CASE (1)
     vD1 = vplatD
     vD2 = 0.
  CASE (2)
     vD1 = 0.
     vD2 = vplatD
  CASE (3)
     vD1 = vplatD
     vD2 = vplatD
  CASE DEFAULT
     vD1 = vplatD
     vD2 = vplatD
  END SELECT

! THE TENSILE-MOTION VECTOR BOUNDARY CONDITION 
  SELECT CASE (ibcT)
  CASE (1)
     vT1 = vplatT
     vT2 = 0.
  CASE (2)
     vT1 = 0.
     vT2 = vplatT
  CASE (3)
     vT1 = vplatT
     vT2 = vplatT
  CASE DEFAULT
     vT1 = vplatT
     vT2 = vplatT
  END SELECT

! DEFINE VARIABLES FOR ROTATION LOOP
  indn  = 0
  indnm = 0
  DO i=1,nfault
     xfe = xfen(i)
     yfe = yfen(i)
     xfw = xfwn(i)
     yfw = yfwn(i)
     dt  = dtn(i)
     db  = dbn(i)

! COMPUTE THE FAULT MIDPOINT
     xmidn(i) = .5*(xfw+xfe)
     ymidn(i) = .5*(yfw+yfe)

! CHECK TO SEE IF THE DIP HAS A "NEGATIVE" VALUE AND NEED TO SUB-PI
     IF (segdp1(i) < 0.) THEN
        segdp(i) = segdp1(i) * -1.0
        SUBPI(i) = -PI
     ELSE
        segdp(i) = segdp1(i)
        SUBPI(i) = 0.
     END IF

! COMPUTE THE ROTATION ANGLE AND ADD SUBPI RADIANS OF ROTATION
! IF DIP CONTAINS A "NEGATIVE" VALUE
     testy = yfw-yfe
     IF(ABS(testy) <= .01) THEN
        xmax = amax1(xfw,xfe)
        xmin = amin1(xfw,xfe)
        xfw  = xmax
        xfe  = xmin
        thetan(i) = 0.
     ELSE
        delx = xfe-xmidn(i)
        dely = yfe-ymidn(i)
        TARG = dely/delx
        thetan(i) = ATAN(TARG) + (SUBPI(i)*SGN(TARG))
!        PRINT *, i, ATAN(TARG)*180./PI, thetan(i)*180./PI
     END IF
  END DO

! FIND CENTER POINT OF FAULT GRID TO ADD TO EACH POINT
  xmax = -10000.
  xmin =  10000.
  ymax = -10000.
  ymin =  10000.
  DO i = 1,nfault
     xmax1 = amax1(xfen(i),xfwn(i))
     IF (xmax1 > xmax) xmax = xmax1
     xmin1 = amin1(xfen(i),xfwn(i))
     IF (xmin1 < xmin) xmin = xmin1
     ymax1 = amax1(yfen(i),yfwn(i))
     IF (ymax1 > ymax) ymax = ymax1
     ymin1 = amin1(yfen(i),yfwn(i))
     IF (ymin1 < ymin) ymin = ymin1
  END DO  
  xmidf = .5 *(xmax + xmin)
  ymidf = .5 *(ymax + ymin)

! DEFINE GRID OF "npoint" OBSERVATION POINTS
! MODIFY BELOW FOR MORE OR LESS OBSERVATION POINTS 
  AINCX   = (ABS(XSTRT)+ABS(YSTRT))/FLOAT(IPTS)
  AINCY   = (ABS(XSTRT)+ABS(YSTRT))/FLOAT(JPTS)

! DEFINE POINTS WITH RESPECT TO CENTER OF FAULT GRID
! THEN ADD MID POINT TO THE COORDINATES
! **i * j MUST BE AT LEAST THE SAME SIZE AS NGR2 IN GRV_MAIN**
! NOTE: TO MAKE THE IDL PLOTS WE START THE GRID FROM THE NW CORNER
  l = 0
  DO i = 1,IPTS
     DO j = 1,JPTS
        l = l + 1
        xpt(l) = XSTRT + AINCX*FLOAT(i-1) + xmidf
        ypt(l) = YSTRT - AINCY*FLOAT(j-1) + ymidf ! NW-CORNER START PT
!         ypt(l) = -YSTRT. + AINC*FLOAT(j-1) + ymidf ! SW-CORNER START PT
     END DO
  END DO

! WRITE OUT OBSERVATION GRID TO OUTPUT FILE -> UNIT 90
  WRITE(90,'(A30)') fgeo_name  
  WRITE(90,'(A30)') fmed_name
  WRITE(90,*) nfault,npoint
  WRITE(90,*) (xpt(n), n=1,npoint)
  WRITE(90,*) (ypt(n), n=1,npoint)
  PRINT *,' '
  PRINT *, ' SURFACE OBSERAVTION GRID HAS BEEN CREATED. '
  PRINT *, ' ===> NOW CALCULATING GREENs FUNCTIONS. '
  PRINT *,' '

! START THE MAIN LOOP FOR EACH FAULT SEGMENT -> *ELASTIC* SLIP
  DO i=1,nfault
     st=SIN(thetan(i))
     ct=COS(thetan(i))
     
! ROTATION ORIGIN ASSIGNMENT FOR OKUBO's GRAVITY ROUTINE
! CHECK THE TWO CASES FOR POSITIVE OR "NEGATIVE" DIP
! FIRST CASE IS THE STANDARD VERSION, SECOND FOR (-)VE DIP
! PROCEDURE:
! (1) SET THE ROTATION ANGLE THEN TEST        
! (2) TEST FOR PATCHES IN NEGATIVE X HALF OF OBS. GRID
! (3) TEST FOR PATCHES IN POSITIVE X HALF OF OBS. GRID
! (4) TEST FOR COMPONETS WHICH MAY CROSS FROM (+) TO (-)
! THEN,
! SECOND CASE: THE SYSTEM NEEDS TO BE ROTATED PI EXTRA 
! DEGREES TO ACCOUNT FOR DIP ORIENTATION W.R.T. THE
! ORIGINAL, POSTIVE (NORTH), AXIS (REPEAT ABOVE TEST)
     SGNW  = SGN(xfwn(i))
     SGNE  = SGN(xfen(i))
     STEST = SGNW + SGNE
     ATEST = ABS(xfwn(i)) - ABS(xfen(i))
     SGNDP = SGN(segdp1(i))

! TEST LOOPS
     IF (SGNDP >= 0.) THEN     
        xmid = xfwn(i)
        ymid = yfwn(i)     
        IF((STEST == -2.) .AND. (ATEST < 0.)) THEN
           xmid = xfen(i)
           ymid = yfen(i)
        END IF
        IF((STEST == 2.) .AND. (ATEST > 0.)) THEN
           xmid = xfen(i)
           ymid = yfen(i)
        END IF     
        IF((STEST == 0.) .AND. (SGNW > SGNE)) THEN
           xmid = xfen(i)
           ymid = yfen(i)
        END IF
     ELSE
        xmid = xfen(i)
        ymid = yfen(i)     
        IF((STEST == -2.) .AND. (ATEST < 0.)) THEN
           xmid = xfwn(i)
           ymid = yfwn(i)
        END IF
        IF((STEST == 2.) .AND. (ATEST > 0.)) THEN
           xmid = xfwn(i)
           ymid = yfwn(i)
        END IF
        IF((STEST == 0.) .AND. (SGNW > SGNE)) THEN
           xmid = xfwn(i)
           ymid = yfwn(i)
        END IF
     END IF

! ALL POINTS ARE AFFECTED BY ALL FAULTS
     DO n=1,npoint
        xpt(n) =  xpt(n) - xmid
        ypt(n) =  ypt(n) - ymid 
        x(n)   =  xpt(n)*ct + ypt(n)*st
        y(n)   = -xpt(n)*st + ypt(n)*ct
! THE DIPCOR SHIFTS THE TOP FAULT EDGE TO X=0
        dipcor = dbn(i) - dtn(i)
        xxpp(n)=  x(n) 
        IF (ABS(segdp1(i)) < 90.) THEN
           yypp(n)=  y(n) + dipcor*COS((ABS(segdp(i))*RAD))
        ELSE
           yypp(n)=  y(n)
        END IF
        xpt(n) =  xpt(n) + xmid
        ypt(n) =  ypt(n) + ymid 
     END DO

! HAVE TO ADD *FLEN* VARIABLE SINCE OKUBO USES FULL FAULT LENGTH
     FLEN = SQRT((xfwn(i)-xfen(i))**2 + (yfwn(i)-yfen(i))**2)
     dt   = dtn(i)
     db   = dbn(i)
     dipp = segdp(i)

! SET MEDIUM PARAMETERS, IF AMU = ALAM -> POISSON SOLID. 
! THESE ARE READ IN FROM THE MEDIUM_INFO.d FILE
!     amu = 1.             
!     alam = bulkm - 2./3. 
! BELOW FROM [Physics of the Earth] --- by STACEY (App. D)
     ZNU = alam/(2.D0*(amu + alam))

! HERE WE SET EACH SEGMENT's FAULT MECHANISM.  ASSUME ALL ARE 
! ACTIVE 
     OFF  = 0.
     USTR = 1.D0
     UDIP = 1.D0
     UTNS = 1.D0
     
!     IF (slpvlS(i) == 0.) USTR = OFF
!     IF (slpvlD(i) == 0.) UDIP = OFF
     IF (slpvlT(i) == 0.) UTNS = OFF
     
! OKUBO USES (+)VE SLIP AS LEFT LATERAL
! UDS = UNIT STRIKE-SLIP DISPLACEMENT 
! UDD = UNIT DIP-SLIP DISPLACEMENT
! UDT = UNIT TENSILE DISPLACEMENT
! NOTE:  SLIP HISTORY FILE DEFINES SLIP IN TERMS OF [cm]
! BUT OKUBO DEFINES SLIP [m] -> CHANGE TO UNITS OF [m] HERE
! IN ORDER TO OBTAIN UNIT-CM SLIP GREEN's FUNCTIONS
     UDS = CM2M*USTR*SGN(slpvlS(i))
     UDD = CM2M*UDIP*SGN(slpvlD(i))
     UDT = CM2M*UTNS*SGN(slpvlT(i))
     PRINT *, ' ELASTIC-SLIP CALCULATIONS FOR FAULT SEGMENT: ', i

! (3) CALCULATE THE GRAVITY GREEN'S FUNCTION COEEFFCIENTS
!     FOR THE NON-STEADY STATE COMPONENT OF THE SYSTEM
!     NOTE:  THE CODE IS VALID ONLY FOR FAULTS
!     *BURIED* WITHIN THE MEDIUM.  (SEE OKUBO 1990 AND 1992)
! ------------------------------------------------------------
!
! FIND ELASTIC SLIP CHANGES FOR EACH SLIP TYPE POSSIBLE IN CASE OF 
! MIXED SLIP MECHANISMS. 
!
! STRIKE SLIP MOTION VECTOR
     CALL GRV_MAIN(GZTES,GZDVES,GZFAES,xxpp,yypp,UDS,OFF,OFF,&
          &        dipp,dt,db,FLEN,npoint,ZNU)
! DIP SLIP MOTION VECTOR
     CALL GRV_MAIN(GZTED,GZDVED,GZFAED,xxpp,yypp,OFF,UDD,OFF,&
          &        dipp,dt,db,FLEN,npoint,ZNU)     
! TENSILE MOTION VECTOR
     CALL GRV_MAIN(GZTET,GZDVET,GZFAET,xxpp,yypp,OFF,OFF,UDT,&
          &        dipp,dt,db,FLEN,npoint,ZNU)
 
! STORE EACH SEGMENT INTO A TEMPORARY VECTOR
     DO np=1,npoint
        indn          = indn + 1
        GZTES1(indn)  = GZTES(np)
        GZDVES1(indn) = GZDVES(np)
        GZFAES1(indn) = GZFAES(np)

        GZTED1(indn)  = GZTED(np)
        GZDVED1(indn) = GZDVED(np)
        GZFAED1(indn) = GZFAED(np)

        GZTET1(indn)  = GZTET(np)
        GZDVET1(indn) = GZDVET(np)
        GZFAET1(indn) = GZFAET(np)
     END DO
! THE BELOW *END DO* GOES WITH THE LOOP OVER i (EACH FAULT SEG)
  END DO

! WRITE THE OUTPUT INTO FILE 90
  isum  = nfault*npoint
  DO ik = 1,isum

! INITIALIZE THE FINAL GRAVITY VALUES TO BE WRITTEN TO A FILE
     GZTES2(ik)  = 0.
     GZDVES2(ik) = 0.
     GZFAES2(ik) = 0.

     GZTED2(ik)  = 0.
     GZDVED2(ik) = 0.
     GZFAED2(ik) = 0.

     GZTET2(ik)  = 0.
     GZDVET2(ik) = 0.
     GZFAET2(ik) = 0.

! COMMENT OUT BELOW IF ONLY STEADY COMPONENTS WANTED
     IF(ABS(GZTES1(ik))  >= 1.e-10)  GZTES2(ik) = GZTES1(ik)
     IF(ABS(GZDVES1(ik)) >= 1.e-10) GZDVES2(ik) = GZDVES1(ik)
     IF(ABS(GZFAES1(ik)) >= 1.e-10) GZFAES2(ik) = GZFAES1(ik)

     IF(ABS(GZTED1(ik))  >= 1.e-10)  GZTED2(ik) = GZTED1(ik)
     IF(ABS(GZDVED1(ik)) >= 1.e-10) GZDVED2(ik) = GZDVED1(ik)
     IF(ABS(GZFAED1(ik)) >= 1.e-10) GZFAED2(ik) = GZFAED1(ik)

     IF(ABS(GZTET1(ik))  >= 1.e-10)  GZTET2(ik) = GZTET1(ik)
     IF(ABS(GZDVET1(ik)) >= 1.e-10) GZDVET2(ik) = GZDVET1(ik)
     IF(ABS(GZFAET1(ik)) >= 1.e-10) GZFAET2(ik) = GZFAET1(ik)
  END DO

! WRITE OUT ELASTIC COMPONENETS OF GRAVITY GREEN's FUNCTIONS
  PRINT *, ' '
  PRINT *, 'WRITING OUT UNIT-SLIP ELASTIC GREENs FUNCTIONS...'
  PRINT *, ' '

! STRIKE SLIP 
  WRITE(90,*) (GZTES2(n), n=1,isum)
  WRITE(90,*) (GZDVES2(n), n=1,isum)
  WRITE(90,*) (GZFAES2(n), n=1,isum)

! DIP SLIP
  WRITE(90,*) (GZTED2(n), n=1,isum)
  WRITE(90,*) (GZDVED2(n), n=1,isum)
  WRITE(90,*) (GZFAED2(n), n=1,isum)

! TENSILE MOTION
  WRITE(90,*) (GZTET2(n), n=1,isum)
  WRITE(90,*) (GZDVET2(n), n=1,isum)
  WRITE(90,*) (GZFAET2(n), n=1,isum)
! END OF THE ELASTIC SLIP GREEN's FUNCTION CALCULATIONS

! (4) CALCULATE THE STEADY STATE PART
! 
! STEADY STATE CALCULATIONS THIS USES A 
! *VELOCITY* FOR UDS, WHICH IS THEN CONVERTED TO
! SLIP OVER A CHOSEN EPOCH IN GRAVITY.f90 
! ------------------------------------------------------------
! FIND NORTH-WESTERN MOST END
  xnorth = xfwn(1)
  ynorth = yfwn(1)
  SELECT CASE (b1)
  CASE ("N")
     DO i=1,nfault              
        IF(yfwn(i) >= ynorth) THEN ! NORTHERN HEMISPHERE 
!      if(yfwn(i) <= ynorth) then 
           inorth = i
           ynorth = yfwn(i)
           xnorth = xfwn(i)
        END IF
     END DO
  CASE ("W")
     DO i=1,nfault              
        IF(xfwn(i) <= xnorth) THEN ! NORTHERN HEMISPHERE 
!      if(xfwn(i) <= xnorth) then 
           inorth = i
           ynorth = yfwn(i)
           xnorth = xfwn(i)
        END IF
     END DO
  CASE DEFAULT           
     DO i=1,nfault              
        IF(yfwn(i) >= ynorth) THEN ! NORTHERN HEMISPHERE 
!      if(yfwn(i) <= ynorth) then
           inorth = i
           ynorth = yfwn(i)
           xnorth = xfwn(i)
        END IF
     END DO
  END SELECT

! FIND SOUTH-EASTERN MOST END
  xsouth = xnorth
  ysouth = ynorth
  SELECT CASE (b2)
  CASE ("S")
     DO i=1,nfault
        IF(yfen(i) <= ysouth) THEN ! NORTH (+)VE AXIS 
!      if(yfen(i) >= ysouth) then ! SOUTH (+)VE AXIS 
           isouth = i
           ysouth = yfen(i)
           xsouth = xfen(i)
        END IF
     END DO
  CASE ("E")
     DO i=1,nfault
        IF(xfen(i) >= xsouth) THEN ! NORTH (+)VE AXIS 
!      if(yfen(i) >= ysouth) then ! SOUTH (+)VE AXIS 
           isouth = i
           ysouth = yfen(i)
           xsouth = xfen(i)
        END IF
     END DO
  CASE DEFAULT
     DO i=1,nfault
        IF(yfen(i) <= ysouth) THEN ! NORTH (+)VE AXIS 
!      if(yfen(i) >= ysouth) then ! SOUTH (+)VE AXIS 
           isouth = i
           ysouth = yfen(i)
           xsouth = xfen(i)
        END IF
     END DO
  END SELECT

! ADD TWO NEW SEGMENTS FOR THE UNDERLYING TECTONIC PLATE, i.e., THIS IS
! THE BOUNDARY-CONDITION PART OF THE CODE
  vplateAVG  = SQRT(vplatS**2 + vplatD**2 + vplatT**2)
  VTEST      = ABS(SGN(vplateX) + SGN(vplateY))

! BOUNDARY FOR NORTH-WESTERN SEGMENT
! IN THIS CASE THE SEGMENT IS *PARALLEL* TO THE AVERAGE PLATE VELOCITY
! VECTOR DEFINED BY vplatX AND vplatY
! NOTA BENE: WE DO NOT NOT TO WORRY ABOUT MIXED COMPONENTS OF THE 
! AVERAGE SURFACE PLATE VELOCITY AS WE ONLY USE THE X-dir AND Y-dir 
! COMPONENTS WHEN DESCRIBING THE VELOCITY VECTOR.
! THIS CASE LOOKS AT PURE STRIKE-SLIP AND OBLIQUE END-SEGMENTS
  inew       =  nfault + 1
  yfen(inew) =  ynorth
  xfen(inew) =  xnorth
  IF (vS1 /= 0.) THEN
     IF (VTEST > 0.) THEN
        yfwn(inew) =  ABS(vplateY*3000./vplateAVG) + ynorth
        xfwn(inew) =  ABS(vplateX*3000./vplateAVG) + xnorth  
     ELSE
        yfwn(inew) =  ABS(vplateY*3000./vplateAVG) + ynorth
        xfwn(inew) = -ABS(vplateX*3000./vplateAVG) + xnorth  
     END IF

! PURE DIP-SLIP/TENSILE MOTION BOUNDARY FOR NORTH-WESTERN SEGMENT
! IN THIS CASE THE SEGMENT IS *PERPENDICULAR* TO THE AVERAGE PLATE 
! VELOCITY VECTOR DEFINED BY vplatX AND vplatY
  ELSE IF (vS1 == 0.) THEN        
     IF (VTEST > 0.) THEN
        yfwn(inew) =  ABS(vplateX*3000./vplateAVG) + ynorth
        xfwn(inew) = -ABS(vplateY*3000./vplateAVG) + xnorth  
     ELSE
        yfwn(inew) = -ABS(vplateX*3000./vplateAVG) + ynorth
        xfwn(inew) = -ABS(vplateY*3000./vplateAVG) + xnorth 
     END IF
  END IF
  
! SET VELOCITIES FOR NORTH-WESTERN-MOST BOUNDARY SEGMENT
  segdp1(inew) = dip1
  IF ((SGN(slpvlS) < 0.) .AND. (SGN(segdp1) < 0.)) THEN
     slpvlS(inew) = vS1 * -1.
  ELSE
     slpvlS(inew) = vS1
  END IF
  slpvlD(inew) = vD1
  slpvlT(inew) = vT1

! SOUTH-EASTERN-MOST EXTENSION/BOUNDARY
  inew       =  nfault + 2
  yfwn(inew) =  ysouth
  xfwn(inew) =  xsouth
! STRIKE-SLIP BOUNDARY FOR SOUTH-EASTERN SEGMENT
! IN THIS CASE THE SEGMENT IS *PARALLEL* TO THE AVERAGE PLATE VELOCITY
! VECTOR DEFINED BY vplatX AND vplatY
! NOTA BENE: WE DO NOT NOT TO WORRY ABOUT MIXED COMPONENTS OF THE 
! AVERAGE SURFACE PLATE VELOCITY AS WE ONLY USE THE X-dir AND Y-dir 
! COMPONENTS WHEN DESCRIBING THE VELOCITY VECTOR.
! THIS CASE LOOKS AT PURE STRIKE-SLIP AND OBLIQUE END-SEGMENTS
  IF (vS2 /= 0.) THEN
     IF (VTEST > 0.) THEN
        yfen(inew) = -ABS(vplateY*3000./vplateAVG) + ysouth
        xfen(inew) = -ABS(vplateX*3000./vplateAVG) + xsouth  
     ELSE
        yfen(inew) = -ABS(vplateY*3000./vplateAVG) + ysouth
        xfen(inew) =  ABS(vplateX*3000./vplateAVG) + xsouth  
     END IF

! DIP-SLIP/TENSILE MOTION BOUNDARY FOR SOUTH-EASTERN SEGMENT
! IN THIS CASE THE SEGMENT IS *PERPENDICULAR* TO THE AVERAGE PLATE 
! VELOCITY VECTOR DEFINED BY vplatX AND vplatY
  ELSE IF (vS2 == 0.) THEN        
     IF (VTEST > 0.) THEN
        yfen(inew) = -ABS(vplateX*3000./vplateAVG) + ysouth
        xfen(inew) =  ABS(vplateY*3000./vplateAVG) + xsouth  
     ELSE
        yfen(inew) =  ABS(vplateX*3000./vplateAVG) + ysouth
        xfen(inew) =  ABS(vplateY*3000./vplateAVG) + xsouth 
     END IF
  END IF

! SET VELOCITIES FOR SOUTH-EASTERN-MOST BOUNDARY SEGMENT
  segdp1(inew) = dip2
  IF ((SGN(slpvlS) < 0.) .AND. (SGN(segdp1) < 0.)) THEN
     slpvlS(inew) = vS2 * -1.
  ELSE
     slpvlS(inew) = vS2
  END IF
  slpvlD(inew) = vD2
  slpvlT(inew) = vT2

! INITIALIZE THE GRAVITY MATRICES
  DO i=1,NPT
     GZTS(i)  = 0.
     GZDVS(i) = 0.
     GZFAS(i) = 0.
  END DO
  PRINT *, ' '

! DEFINE VARIABLES
  nflt=nfault + 2
  DO i=1,nflt
     xfe = xfen(i)
     yfe = yfen(i)
     xfw = xfwn(i)
     yfw = yfwn(i)

! WE HAVE TO REDEFINE WHERE WE ROTATE THE COORDINATE SYSTEM
! OKUBO USES OKADA'S COORDINATE SYSTEM AND ITS 
! ORIGIN IS THE FAULT END WHICH CORRESPONDS TO: (xfwn, yfwn)  
! BUT SOMETIMES, xfwn < xfen, THEREFORE xfen BECOMES THE ORIGIN.

! COMPUTE THE FAULT MIDPOINTS
     xmidn(i) = .5*(xfw+xfe) ! CHINNERY ORIGIN - USED FOR ANGLE CALC
     ymidn(i) = .5*(yfw+yfe) ! CHINNERY ORIGIN - USED FOR ANGLE CALC

! ROTATION ORIGIN ASSIGNMENT FOR OKUBO's GRAVITY ROUTINE
! CHECK THE TWO CASES FOR POSITIVE OR "NEGATIVE" DIP
! FIRST CASE IS THE STANDARD VERSION, SECOND FOR (-)VE DIP
! PROCEDURE:
! (1) SET THE ROTATION ANGLE THEN TEST        
! (2) TEST FOR PATCHES IN NEGATIVE X HALF OF OBS. GRID
! (3) TEST FOR PATCHES IN POSITIVE X HALF OF OBS. GRID
! (4) TEST FOR COMPONETS WHICH MAY CROSS FROM (+) TO (-)
! THEN,
! SECOND CASE: THE SYSTEM NEEDS TO BE ROTATED 90 EXTRA 
! DEGREES TO ACCOUNT FOR DIP ORIENTATION W.R.T. THE
! ORIGINAL POSTIVE (NORTH) AXIS (REPEAT ABOVE TEST)
     SGNW  = SGN(xfw)
     SGNE  = SGN(xfe)
     STEST = SGNW + SGNE
     ATEST = ABS(xfw) - ABS(xfe)
     SGNDP = SGN(segdp1(i))

! INITIALIZE THE ROTATION POINT TO THE *WEST* ENDS OF THE SEGMENTS
     xmidt(i) = xfw
     ymidt(i) = yfw

! TEST LOOPS
     IF (SGNDP >= 0.) THEN     
        xmidt(i) = xfw
        ymidt(i) = yfw     
        IF((STEST == -2.) .AND. (ATEST < 0.)) THEN
           xmidt(i) = xfe
           ymidt(i) = yfe
        END IF
        IF((STEST == 2.) .AND. (ATEST > 0.)) THEN
           xmidt(i) = xfe
           ymidt(i) = yfe
        END IF     
        IF((STEST == 0.) .AND. (SGNW > SGNE)) THEN
           xmidt(i) = xfe
           ymidt(i) = yfe
        END IF
     ELSE
        xmidt(i) = xfe
        ymidt(i) = yfe     
        IF((STEST == -2.) .AND. (ATEST < 0.)) THEN
           xmidt(i) = xfw
           ymidt(i) = yfw
        END IF
        IF((STEST == 2.) .AND. (ATEST > 0.)) THEN
           xmidt(i) = xfw
           ymidt(i) = yfw
        END IF
        IF((STEST == 0.) .AND. (SGNW > SGNE)) THEN
           xmidt(i) = xfw
           ymidt(i) = yfw
        END IF
     END IF

! CHECK TO SEE IF THE DIP HAS A "NEGATIVE" VALUE AND ADD PI RADIANS
     IF (segdp1(i) < 0.) THEN
        segdp(i) = segdp1(i) * -1.0
        SUBPI(i) = -PI
     ELSE
        segdp(i) = segdp1(i)
        SUBPI(i) = 0.
     END IF
     
! COMPUTE THE ROTATION ANGLE AND ADD PI RADIANS OF ROTATION IN 
! THE SAME SENSE AS ORIGINAL ROTATION IF DIP CONTAINS THE 
! "NEGATIVE" VALUE -> WE SUBTRACT PI TO KEEP THE ROTATION UNDER 180
     testy = yfw-yfe
     IF(ABS(testy) <= 1.) THEN
        xmax = amax1(xfw,xfe)
        xmin = amin1(xfw,xfe)
        xfw  = xmax
        xfe  = xmin
        thetan(i) = 0.
     ELSE
        delx = xfe-xmidn(i)
        dely = yfe-ymidn(i)
        TARG = dely/delx
        thetan(i) = ATAN(TARG) + (SUBPI(i)*SGN(TARG))
     END IF
  END DO

! ROTATE OBSERVATION POINTS ABOUT THE PROPER ORIGIN -> OKUBO'S
  DO i=1,nflt
     st = SIN(thetan(i))
     ct = COS(thetan(i))
!     xmid=xmidn(i)  ! CHINNERY ORIGIN
!     ymid=ymidn(i)  ! CHINNERY ORIGIN
     xmid = xmidt(i)   ! THIS IS OKUBO'S ORIGIN
     ymid = ymidt(i)   ! THIS IS OKUBO'S ORIGIN

! INTIALIZE THE DEPTHS FOR THE MAIN DRIVING PLATE    
! SET THE DEPTH FOR THE PLATE WHERE WE USE THE HANGING WALL FOR 
! REFERENCE, OTHERWISE COMMENT OUT AND SET MANUALLY
! LEFT LATERAL COMPONENT -> "NORTHERN"-BOUNDARY SEGMENT    
     IF (i == (nflt-1))  THEN 
        dt = dbn(inorth)     
! RIGHT LATERAL COMPONENT -> "SOUTHERN"-BOUNDARY SEGMENT
     ELSE IF (i == (nflt)) THEN
        dt = dbn(isouth)
! IF NOT A BOUNDARY SEGMENT USE THE TOP OF THE FAULT SEGMENT SPECIFIED
     ELSE
        dt = dtn(i)
     END IF
     db = dbplat

     DO n=1,npoint
        xpt(n) =  xpt(n) - xmid
        ypt(n) =  ypt(n) - ymid
! THE DIPCOR SHIFTS THE TOP FAULT EDGE TO X=0
        dipcor = db - dt
        x(n)   =  xpt(n)*ct + ypt(n)*st
        IF (ABS(segdp1(i)) < 90.) THEN
           y(n)   = -xpt(n)*st + ypt(n)*ct + dipcor*COS((ABS(segdp(i))*RAD))
        ELSE
           y(n)   = -xpt(n)*st + ypt(n)*ct
        END IF
        xpt(n) =  xpt(n) + xmid
        ypt(n) =  ypt(n) + ymid
     END DO

! BELOW WILL YIELD THE GREEN'S FUNCTIONS FOR UNIT [cm] SLIP PER YEAR
! WE ALSO ASSUME THAT THE MOTION PER YEAR IS A CONSTANT, THUS THERE
! IS NO NEED FOR THE THREE *SEPARATE* MOTION-TYPE CALCULATIONS.
! CONVERT UNIT [cm] SLIP TO [m] FOR USE IN OKUBO
     udsV = slpvlS(i)*CM2M
     uddV = slpvlD(i)*CM2M
     udtV = slpvlT(i)*CM2M
     dips = segdp(i)
     FLEN = SQRT((xfwn(i)-xfen(i))**2 + (yfwn(i)-yfen(i))**2)

! CALL THE GRAVITY GREEN'S FUCNTION SUBROUTINE FOR STEADY STATE MOTION
     CALL GRV_MAIN(GTS,GDVS,GFAS,x,y,udsV,uddV,udtV,&
          &      dips,dt,db,FLEN,npoint,ZNU)

! SUM THE STEADY-STATE COMPNONENTS FOR ALL OBSERVATION POINTS
     DO n=1,npoint
        GZTS(n)  = GZTS(n)  + GTS(n)
        GZDVS(n) = GZDVS(n) + GDVS(n)
        GZFAS(n) = GZFAS(n) + GFAS(n)
! TO SET STEADY STATE TO ZERO & OBSERVE ELASTIC COMPONENT ONLY, THEN
! TURN OFF THIS COMPONENT BY UNCOMMENTING BELOW
!        GZTS(n)  = 0.
!        GZDVS(n) = 0.
!        GZFAS(n) = 0.
     END DO

! PRINT OUT STATEMENT FOR UNDERLYING PLATE
     IF (i <= nfault) THEN
        PRINT *, ' STEADY-STATE CALCULATIONS FOR FAULT SEGMENT: ', i
     ELSE IF (i == (nflt-1)) THEN
        PRINT *, ' STEADY-STATE CALCULATIONS FOR NW-BOUNDARY  : ', i
     ELSE
        PRINT *, ' STEADY-STATE CALCULATIONS FOR SE-BOUNDARY  : ', i
     END IF

! THIS *END DO* GOES WITH THE LOOP OVER i FOR EACH SEGMENT
! IN THE NEW FAULT GEOMETRY WITH TWO EXTRA END SEGMENTS 
  END DO
! END OF STEADY STATE SLIP COMPONENT PART
! ----------------------------------------------------------

! WRITE OUT STEADY STATE COMPONENETS OF GRAVITY GREEN's FUNCTIONS
  PRINT *, ' '
  PRINT *, 'WRITING OUT STEADY STATE GREENs FUNCTIONS...'
  PRINT *, ' '
  
! WRITE THE OUTPUT INTO FILE 90
  WRITE(90,*) (GZTS(n), n=1,npoint)
  WRITE(90,*) (GZDVS(n), n=1,npoint)
  WRITE(90,*) (GZFAS(n), n=1,npoint)
  
! CLOSE FILES
  CLOSE(20)  ! FAULT TOPOLOGY FILE
  CLOSE(90)  ! OUTPUT FILE

! END OF THE MAIN ALGORITHM
END PROGRAM GRAV_GREEN
! ----------------------------------------------------------



! SUBROUTINES AND FUNCTIONS 
! ----------------------------------------------------------


! GRAVITY SUBROUTINE FOR FAULTING ON A FINITE RECT PLANE
! ----------------------------------------------------------
SUBROUTINE GRV_MAIN(GSUM,GDIV,GFB,XKM,YKM,U1,U2,U3,&
     &     DIP,DEPTOP,DEPBOT,LENKM,NOBS,ZNU)
  USE GRAV_PARS
  USE ELAS
  USE CONST
  IMPLICIT NONE 
  INTEGER :: I
  INTEGER, INTENT(IN) :: NOBS
  REAL,DIMENSION(NGR2),INTENT(IN) :: XKM,YKM
  REAL,DIMENSION(NGR2),INTENT(OUT) :: GSUM,GDIV,GFB
  REAL,INTENT(IN) :: U1,U2,U3,DIP,DEPTOP,DEPBOT,LENKM,ZNU
  REAL,DIMENSION(NGR2) :: DH,X,Y,GVOID
  REAL :: WIDKM,DEPBTM,DEPTPM,WID,LEN,RDIP,CD,SD,RHO,DELRHO
  REAL :: RHOCGS,DROCGS,DEPORG,FACFB,FACDIV,FACVID

! INPUT: GRV_MAIN NEEDS EVERY INPUT DATA IN MKS-UNIT FOR MICROGAL
!     OUTPUT [mu-gal]. UNITS INDICATED BELOW ARE THOSE UNITS *PASSED*
!     TO THE SUBROUTINE
!-----------------------------------------
! XKM    = OBS GRID X-OBSERVATION POINT [km]
! YKM    = OBS GRID Y-OBSERVATION POINT [km]
! U1     = STRIKE-SLIP VECTOR MAGNITUDE [m] <- LEFT-LATERAL IS (+)VE!
! U2     = DIP-SLIP VECTOR MAGNITUDE [m] 
! U3     = TENSILE MOTION VECTOR MAGNITUDE [m] 
! DEPTOP = DEPTH TO TOP OF FAULT [km]
! DEPBOT = DEPTH TO BOTTOM OF FAULT [km]
! LENKM  = RUNDLE FAULT LENGTH [km]
! NOBS   = NUMBER OF OBS. POINTS <-- MUST BE <= MIN OF "NGR2" IN GRV_PARS!
!----------------------------------------


! OUTTPUT: RETURNS GRAVITY VALUES
!----------------------------------------
! GSUM  = TOTAL GRAVITY CHANGES [mu-gal]
! GDIV  = DILATATIONAL EFFECTS ONLY [mu-gal]
! GFB   = FREE-AIR CORR./BOUGUER EFFECT OF GRAVITY [mu-gal]
! GVOID = CAVITATION COMPONENT OF GRAVITY [mu-gal]
! -> THE VOID CALCULATION IS NOT CURRENTLY USED 
!----------------------------------------

! NGR2 NEEDS TO BE AT LEAST AS LARGE AS i * j IN THE LOOP FOR
! THE OBSERVATION GRID...MAKE SURE THIS IS SO!
! -> TO BE SAFE, MAKE NGR2 LARGER THAN IPTS * JPTS IN MODULE ABOVE
  
! INPUTS NEEDED TO BE SPECIFIED BY THE USER OR INPUT FORM CODE
! ------------------------------------------------------------
! 
! FAULT GEOMETRY --- NEED [mks]-UNITS SO WE MUST CONVERT THOSE
! NOT ALREADY GIVEN IN [mks] NOW.
! ------------------------------------------------------------
! U1 = STIKE-SLIP VECTOR [m]
! U2 = DIP-SLIP/TENSILE-MOTION VECTOR [m]
! U3 = EXTENSIONAL VECTOR [m]
  WIDKM  = DEPBOT - DEPTOP
  DEPBTM = DEPBOT * KM2M
  DEPTPM = DEPTOP * KM2M
  WID    = WIDKM  * KM2M
  LEN    = LENKM  * KM2M
  RDIP   = RAD*DIP
  CD     = COS(RDIP)
  SD     = SIN(RDIP)

! MEDIUM VALUES
! ------------------------------------------------------------
  RHO    = 2.67D3
  DELRHO = 0.D0  ! NO CAVITY FILLING MATTER INCLUDED FOR NOW
  XNU    = ZNU
! IF XNU = 1/4 -> POISSON SOLID i.e., LAMDA = MU 
!  XNU = 0.25D0 ! <- ONLY USE IF NOT PASSED THROUGH INTERFACE
  MUETA  = 1.D0 - 2.D0*XNU
  RHOCGS = RHO*1.D-3
  DROCGS = DELRHO*1.D-3

! *** NOTE: USED OKADA'S NOTATION 'd' [SEE OKADA PAPER, 1985]***
! 
! SET UP THE X AND Y COORDINATES OF THE FAULT PLANE
! ------------------------------------------------------
! (EG) - USING A STRIKE SLIP FAULT WITH A DEPTPM OF 1000[m]
! THEN THE FAULT SOURCE ORIGIN IS LOCATED AS FOLLOWS:
! 
! EARTH'S -->  ================================   x_3 AXIS
! SURFACE                                             |
!                                                     |
!                                                     V
! 
!              |------------------- <-- DEPTPM
!              |                  | ^
!              |                  | |
!              |      FAULT       | |
!              |     SURFACE      | |== WID 
!              |   (FOOT WALL)    | |
!              |                  | |
!              |                  | v
!              0------------------- <-- DEPBTM  
!              ^
!              |
!              | 
!              '----  WHERE "O" IS THE LOCATION OF DEPORG
!                     DEPORG = (d - x_3) IN OKUBO(1992)
! 
! IF THE FAULT HAS A DIP ANGLE THEN WID*SD IS THE "HEIGHT"
! ON THE x_3 AXIS, OTHERWISE, WID == db AS IN A PURE STRIKE-SLIP SYSTEM

  IF(SD > 0.0) THEN
     DEPORG = DEPTPM + WID*SD  
  ELSE
     DEPORG = DEPTPM 
  ENDIF

! CONVERT OBSERVATION GRID INTO UNITS OF [m] FOR OKUBO
  DO I=1,NOBS
     X(I) = XKM(I)*KM2M 
     Y(I) = YKM(I)*KM2M
  END DO

! CALL MAIN SUBROUTINES TO CALCULATE GRAVITY CHANGES
! USES ELVFT FOR THE FREE-AIR AND GRVFLT FOR DIV*U
! -------------------------------------------------------
  DO I=1,NOBS
     CALL ELVFLT(X(I),Y(I),U1,U2,U3,DEPORG,LEN,WID,SD,CD,&
          &               DH(I))
     CALL GRVFLT(X(I),Y(I),U1,U2,U3,DEPORG,LEN,WID,SD,CD,&
          &               GDIV(I),GVOID(I))
  END DO

! MULTIPLY OUTPUTS FROM ELVFLT (DH) AND GRVFLT (GDIV,GVOID)
! BY THEIR RESPECTIVE FACTORS TO OBTAIN FINAL RESULTS.
! --------------------------------------------------------

! UNCOMMENT/COMMENT OUT WHICH CORRECTION OR GRAVITY VALUE 
! YOU WISH TO USE
!
! FREE-AIR AND BOUGUER EFFECT <-- THIS IS THE STANDARD ONE
  FACFB  = -(GAMMA - 2.D0*PI*RHO*GNEWTN)*GCONVR

! BOUGUER EFFECT ONLY
!  FACFB  = (2.D0*PI*RHO*GNEWTN)*GCONVR

! FREE-AIR EFFECT ONLY
!  FACFB  = -(GAMMA)*GCONVR

! ON TO THE CHANGES CAUSED BY "-DIV(U)"
  FACDIV = RHO*GNEWTN*GCONVR
  FACVID = DELRHO*GNEWTN*GCONVR
  
  DO I=1,NOBS
     GFB(I)   = FACFB*DH(I)
     GDIV(I)  = FACDIV*GDIV(I) 
! DON'T WANT THE VOID CALCULATION FOR NOW -- MAY 18, 2006
! ALSO MODIFIED THE *GSUM* COMPONENT TO REFLECT THIS
!
! *** THIS IMPLIES TENSILE MOTION NOT USED CURRENTLY      ***
! *** IF YOU WANT TENSILE MOTION UNCOMMENT BELOW          ***
! *** AND YOU MUST MAKE SURE THAT THE FAULT DOES NOT      ***
! *** BREAK THE SURFACE DUE TO SINGULARITIES (OKUBO 1992) ***
! *** MAKE SURE TO ALSO UNCOMMENT LINES FOR VOID CALC IN  ***
! *** GRVFLT SUBROUTINE BELOW!!!!                         ***
!        GVOID(I) = FACVID*GVOID(I) 
!        GSUM(I)  = GFB(I) + GDIV(I) + GVOID(I)
     GSUM(I)  = GFB(I) + GDIV(I)
  END DO
  RETURN
END SUBROUTINE GRV_MAIN

! END OF GRV_MAIN SUB
! --------------------------------------------------------
  
! BELOW ARE THE SUBROUTINES CALLED   
! BY THE MAIN GRV_MAIN SUBROUTINE    
! -------------------------------------------------------

SUBROUTINE ELVFLT(X,Y,U1,U2,U3,DEP,LEN,WID,SD,CD,UZ)
  USE CONST
  IMPLICIT NONE 
  REAL,INTENT(IN) :: X,Y,U1,U2,U3,DEP,LEN,WID,SD,CD
  REAL,INTENT(OUT) :: UZ
  REAL :: P,Q,U1PI,U2PI,U3PI,US,UD,UT
  REAL :: USTRK1,USTRK2,USTRK3,USTRK4
  REAL :: UDIP1,UDIP2,UDIP3,UDIP4
  REAL :: UTENS1,UTENS2,UTENS3,UTENS4

! FUNCTION DECLARATION
!  REAL,EXTERNAL :: HFLT

! SURFACE UPLIFT DUE TO RECTANGULAR SHEAR FAULTING
! INPUTS
! X,Y    : COOORDINATES OF THE OBSERVATION SITE
! U1,U2  : SLIP VECTOR
! U3     : TENSILE MOVEMENT
! WID    : FAULT WIDTH
! DEP    : OKADA'S NOTATION FOR 'd' (SEE BELOW FOR REF.)
! LEN    : FAULT LENGTH (OKUBO DOESN'T HAVE THIS LINE!)
! SD     : SIN (DELTA) [DELTA = DIP ANGLE]
! CD     : COS (DELTA)

! OUTPUT
! UZ     : SURFACE UPLIFT

! GEOMETRY OF THE FAULT IS FOUND IN OKADA 1985 BSSA, v.75 (1138)
!  PI  = ATAN(1.D0)*4.D0 <- PASSED THROUGH CONST MODULE NOW
  P    = Y*CD + DEP*SD
  Q    = Y*SD - DEP*CD
  U1PI = U1/(2.D0*PI) 
  U2PI = U2/(2.D0*PI)
  U3PI = U3/(2.D0*PI)
  
! INTEGRATE THE SOLUTION A LA CHINNERY
  CALL HFLT(X    ,P    ,Q,SD,CD,USTRK1,UDIP1,UTENS1)
  CALL HFLT(X    ,P-WID,Q,SD,CD,USTRK2,UDIP2,UTENS2)
  CALL HFLT(X-LEN,P    ,Q,SD,CD,USTRK3,UDIP3,UTENS3)
  CALL HFLT(X-LEN,P-WID,Q,SD,CD,USTRK4,UDIP4,UTENS4)

! STRIKE-SLIP MOTION COMPONENT
  IF (U1 == 0.D0) THEN
     US = 0.D0
  ELSE
     US = USTRK1 - USTRK2 - USTRK3 + USTRK4
  END IF
! DIP-SLIP AND TENSILE MOTION COMPONENT
  IF (U2 == 0.D0) THEN
     UD = 0.D0
  ELSE
     UD = UDIP1  - UDIP2  - UDIP3  + UDIP4
  END IF
! THE CAVITATION AND VOID FILLING MATTER COMPONENT
  IF (U3 == 0.D0) THEN
     UT = 0.D0
  ELSE
     UT = UTENS1 - UTENS2 - UTENS3 + UTENS4
  END IF

! OUTPUT
  UZ = -U1PI*US - U2PI*UD + U3PI*UT
  
  RETURN
END SUBROUTINE ELVFLT

! END OF ELVFLT SUBROUTINE
! -------------------------------------------------------------


SUBROUTINE HFLT(GZI,ETA,Q,SD,CD,ZSTRK,ZDIP,ZTENS)
  USE ELAS
  IMPLICIT NONE 
  REAL,INTENT(IN) :: GZI,ETA,Q,SD,CD
  REAL,INTENT(OUT) :: ZSTRK,ZDIP,ZTENS
  REAL :: I4,I5,CRT,DCRT,Q2,GZI2,R,XX,YBAR,DBAR,RD
  REAL :: RDLOG,RE,RELOG,FACI,RXX,SDCD
  REAL :: QR,QRRE,QRRG,GZIETA,QRCRT,ARG,ARGI1,ARGI2,ATA
! UPLIFT DUE TO FAULTING (CALLED BY SUBROUTINE "ELVFLT")

! INPUT     :
! GZI,ETA,Q : SEE OKADA 1985 BSSA, v.75 (1144)
! SD,CD     : SIN (DELTA), COS (DELTA)

! OUTPUT    :
! ZSTRK     : UZ FOR STRIKE-SLIP
! ZDIP      : UZ FOR DIP-SLIP
! ZTENS     : UZ FOR TENSILE FAULT 

! MUETA = MU/(LAMDA+MU)

! SIN(DIP)...COS(DIP)...TAN(DIP)...CRITERION
  CRT  = 1.D-7
  DCRT = 1.D0
  Q2   = Q*Q
  GZI2 = GZI**2
  R    = SQRT(GZI2 + ETA**2 + Q2)
  XX   = SQRT(GZI2 + Q2)
  YBAR = ETA*CD + Q*SD
  DBAR = ETA*SD - Q*CD
  RD   = R + DBAR
  RDLOG= LOG(RD)
  RE   = R + ETA
  QR   = Q*R
  
! DEFINE LOGARITHM PARAMETERS IN I_4
  IF(ABS(RE) > CRT) THEN
     RELOG = LOG(RE)
     QRRE  = Q/(R*RE)
  ELSE
     RELOG = -LOG(R - ETA)
     QRRE  = 0.
  END IF
  
! DEFINE THE ARC-TANGENT TERM
  GZIETA = GZI*ETA      
  QRCRT  = ABS(GZIETA)*CRT 
  IF(ABS(QR) > QRCRT) THEN
     ARG = GZIETA/QR
     ATA = ATAN(ARG)
  ELSE
     ATA = 0.
  END IF
  
! SET I_4 AND I_5 EQNS ARGS USING OKADA (1985) EQN (28) AND *NOT*
! OKUBO (1992)
  IF(ABS(CD) > CRT) THEN
     FACI  = MUETA/CD
     RXX   = R + XX
     ARGI1 = ETA*(XX + Q*CD) + XX*RXX*SD
     ARGI2 = GZI*RXX*CD
! MAKE SURE DENOMINATOR IN ARC-TANGENT IS NOT ZERO B/C OF GZI
     IF((ABS(GZI) < 1.D0) .OR. (ABS(RXX) < 1.D0)) THEN
        I5 = 0.
     ELSE
        I5 = 2.D0*FACI*ATAN(ARGI1/ARGI2) 
     END IF
     I4 = FACI*(RDLOG - SD*RELOG)
  ELSE
! IF COSDIP == 0 THEN USE OKUBO EQNs (63)--(64) FOR I_4 AND I_5 
     I5 = -MUETA*GZI*SD/RD
     I4 = -MUETA*Q/RD
  END IF  

  QRRG = Q/(R*(R+GZI))
  SDCD = SD*CD
  
! OUTPUT -> OKUBO (1992) EQUATIONS (58), (59), AND (60)
  ZSTRK = DBAR*QRRE + QRRE*R*SD + I4*SD
  ZDIP  = DBAR*QRRG + SD*ATA - I5*SDCD
  ZTENS = YBAR*QRRG + CD*(QRRE*GZI - ATA) - I5*(SD**2)
  
  RETURN
END SUBROUTINE HFLT

! END OF HFLT SUBROUTINE
! -------------------------------------------------------------


SUBROUTINE GRVFLT(X,Y,U1,U2,U3,DEP,LEN,WID,SD,CD,GDIVRG,GVOID)
  USE ELAS
  IMPLICIT NONE 
  REAL,INTENT(IN) :: X,Y,U1,U2,U3,DEP,LEN,WID,SD,CD
  REAL,INTENT(OUT) :: GDIVRG,GVOID
  REAL :: P,Q,XL,PW,FACT0,GDIV0,GDIV1

! FUNCTION DECLARATIONS
  REAL,EXTERNAL :: GFLT0,GFLT1,GFLT2

! GRAVITY CHANGE DUE TO FAULTING MOTION ON A FINITE RECT. PLANE

! INPUT :
! U1    : HORIZONTAL SLIP
! U2    : DIP-SLIP
! U3    : TENSILE MOTION
! LEN   : HORIZONTAL LENGTH OF FAULT
! WID   : FAULT WIDTH
! DEP   : OKADA'S NOTATION 'd' [BSSA 1985]
! SD    : SIN (DELTA) - [DELTA = DIP ANGLE]
! CD    : COS (DELTA)
! ---------------------------------------------------------------
! NOTE  : GZ = Xi IN OKUBO 1992   
!       : ET = Eta IN OKUBO 1992 -> FOUND IN GFLT-# CALCULATIONS
! ---------------------------------------------------------------

! OUTPUT:
! GDIVRG: GRAVITY CHANGE FROM (-DIV(U))
!         (MUST BE MULTIPLIED BY "DELTA RHO*G" LATER)
! GVOID : GRAVITY CHANGE FROM VOID FORMATION
!         (MUST BE MULTIPLIED BY "DELTA RHO*G" LATER)

! SEE OKADA 1985 FOR FAULT GEOMETRY [v.75 - p.1138]
  P  = Y*CD + DEP*SD
  Q  = Y*SD - DEP*CD
  XL = X - LEN
  PW = P - WID
  FACT0 = (1.D0 - 2.D0*XNU)*SD

! INTEGRATE OVER THE FULL RANGE TO YIELD GRAVITY OUTPUT
! VIA CHINNERY "DOUBLE BAR" NOTATION -> f(x,y)||
! STRIKE SLIP
  GDIV0 = FACT0*U1*&
       &        (GFLT0(XL,PW,Q,SD,CD) - GFLT0(XL,P,Q,SD,CD)&
       &        - GFLT0(X,PW,Q,SD,CD) + GFLT0(X,P,Q,SD,CD))

! DIP SLIP (NOTA BENE: USES TENSILE VECTOR AS WELL)
  IF (U2 == 0. .AND. U3 == 0.) THEN ! TO AVOID SINGULARITIES
     GDIV1 = 0.
  ELSE
     GDIV1 = FACT0*(U2*CD - U3*SD)*&
          &        (GFLT1(XL,PW,Q,SD,CD) - GFLT1(XL,P,Q,SD,CD)&
          &        - GFLT1(X,PW,Q,SD,CD) + GFLT1(X,P,Q,SD,CD))
  END IF
  GDIVRG = GDIV0 + GDIV1

! SET GVOID TO ZERO TO AVOID THE UNNECESSARY CALCULATION
! OTHERWISE UNCOMMENT BELOW AND SEE COMMENTS ABOVE
  GVOID = 0. 
!  GVOID = U3*&
!       &        (GFLT2(XL,PW,Q,SD,CD) - GFLT2(XL,P,Q,SD,CD)&
!       &        - GFLT2(X,PW,Q,SD,CD) + GFLT2(X,P,Q,SD,CD))
  RETURN
END SUBROUTINE GRVFLT

! END OF GRVFLT SUBROUTINE
! -------------------------------------------------------------


! FUNCTION TO BE USED BY GRVFLT
! -------------------------------------------------------------

FUNCTION GFLT0(GZ,ET,Q,SINDIP,COSDIP)
  IMPLICIT NONE 
  REAL :: GFLT0
  REAL,INTENT(IN) :: GZ,ET,Q,SINDIP,COSDIP
  REAL :: RET,RETLG,TANDIP,RETINV
  REAL :: CRT,RR,R
! FUNCTION TO GIVE THE GRAVITY CHANGE DUE TO HORIZONTAL
! SLIP ON A FINITE RECTANGULAR FAULT.  IT IS CALLED BY 
! "GRVFLT".

! ARGUMENTS:
! GZ,ET,Q, : FAULT COORDINATES - SEE OKADA BSSA 1985
! SINDIP   : SIN (DIP ANGLE) [DELTA]
! COSDIP   : COS (DIP ANGLE) 

! GFLT0 FACILITATES THE CALCULATION OF S_g* IN OKUBO 1992
! BY CALCULATING THE EQUATION I_4 -> SEE EQN (67) AND (63)
! IF COSDIP /= O THEN EQN (61), ELSE IF COSDIP = 0, EQN (63)
  CRT = 1.D-7 
  RR  = GZ**2 + ET**2 + Q**2
  R   = SQRT(RR)
  RET = R + ET
  
! EQUATION (61) CALCUALTION
  IF(ABS(COSDIP) > CRT) THEN
     IF(ABS(RET) > CRT) THEN
        RETLG =  LOG(R + ET)
     ELSE
        RETLG = -LOG(R - ET)
     END IF
     TANDIP = SINDIP/COSDIP
     GFLT0  = - RETLG*TANDIP&
          &           + LOG(R + ET*SINDIP - Q*COSDIP)/COSDIP
  ELSE
! EQUATION (63) CALCULATION
     IF(ABS(RET) > CRT) THEN
        RETINV = 1.D0/(R + ET*SINDIP)
     ELSE
        RETINV = 0
     END IF
     GFLT0 = -Q*SINDIP*RETINV
  END IF
  RETURN
END FUNCTION GFLT0


! FUNCTION TO BE USED BY GRVFLT
! -------------------------------------------------------------

FUNCTION GFLT1(GZ,ET,Q,SINDIP,COSDIP)
  IMPLICIT NONE 
  REAL :: GFLT1
  REAL,INTENT(IN) :: GZ,ET,Q,SINDIP,COSDIP
  REAL :: CRT,RR,R
  REAL :: RET0,SINFAC,SNFAC,GCOS,QCOS,FI1,RET,RETINV
! FUNCTION TO GIVE THE GRAVITY CHANGE DUE TO DIP
! SLIP AND TENSILE MOTION ON A FINITE RECTANGULAR FAULT.
! CALLED BY "GRVFLT"

! ARGUMENTS
! GZ,ET,Q : FAULT COORDINATES AS PER OKADA 1985
! SINDIP  : SIN (DELTA) [DIP ANGLE]
! COSDIP  : COS (DELTA)
! XNU     : POISSON'S RATIO


! GFLT1 FACILITATES THE CALCULATION OF D_g* AND T_g* BY
! CALCULATING THE EQUATION I_5 IN OKUBO (EQN 62)
! IF COSDIP /= 0, THEN ANOTHER CHECK MUST BE MADE FOR GZ IN
! EQN 31 (THAT IS I_4) -> IF GZ = ZERO THEN I_1 = 0.
! IF COSDIP = 0 THEN USE EQN 64 FOR I_5
  CRT    = 1.D-7
  RR     = GZ**2 + ET**2 + Q**2
  R      = SQRT(RR)
  RET0   = R + ET
  SINFAC = 1 + SINDIP
  SNFAC  = RET0*SINFAC
  GCOS   = GZ*COSDIP
  QCOS   = Q*COSDIP

! NOTE THE FI1 IT REPRESENTS THE FUNCTION I_1 FROM OKUBO (1992)
! -> THE "IF" STATEMENT IS REQUIRED AS MENTIONED IN OKUBO
! 1992 (RULE #1) FOR SINGULARITIES (GCOS WOULD BE ZERO HERE
! AND CAUSE A PROBLEM WHEN IN THE DENOMINATOR).
  IF(ABS(COSDIP) > CRT) THEN
     IF(ABS(GZ) < 1.D0) THEN   
        GFLT1 = 0.D0
     ELSE
! EQUATION (62) CALCULATION
        FI1   = ATAN((SNFAC - QCOS)/GCOS)
        GFLT1 = -2.D0*FI1/COSDIP
     END IF
  ELSE
! EQUATION (64) CALCULATION
     RET    = R + ET*SINDIP - Q*COSDIP
     RETINV = 1.D0/RET
     GFLT1  = GZ*SINDIP*RETINV
  END IF
  
  RETURN
END FUNCTION GFLT1


! FUNCTION TO BE USED BY GRVFLT
! -------------------------------------------------------------

FUNCTION GFLT2(GZ,ET,Q,SINDIP,COSDIP)
  IMPLICIT NONE
  REAL :: GFLT2
  REAL,INTENT(IN) :: GZ,ET,Q,SINDIP,COSDIP
  REAL :: CRT,RR,R
  REAL :: RGZ,RGE,RGECRT,FI2,RGZLOG
! FUNCTION TO GIVE THE GRAVITY CHANGE DUE TO THE 
! FORMATION OF A VOID ACROSS A FINITE RECT FAULT

! CALLED BY "GRVFLT"

! ARGUMENTS
! GZ,ET,Q : FAULT COORDINATES AS PER OKADA 1985
! SINDIP  : SIN (DELTA) [DIP ANGLE]
! COSDIP  : COS (DELTA)

! GFLT2 FACILITATES THE CALCULATION OF C_g AND IS THE
! CALCULATION OF EQN (55) IN OKUBO 1992. 
  CRT   = 1.D-20
  RR    = GZ**2 + ET**2 + Q**2
  R     = SQRT(RR)
  RGZ   = R + GZ
  RGE   = R + GZ + ET
  RGECRT= ABS(RGE)*CRT
  
! EQUATION (32) CALCULATION
  IF(ABS(Q) > RGECRT) THEN
     FI2 = ATAN(RGE/Q)
  ELSE
     FI2 = 0.
  END IF

! CHECK BELOW PUT IN BY TYLER HAYES. CORRESPONDS TO A VALUE AT THE
! ORIGIN -> ET = 0, GZ = 0, AND Q = 0.  
! GZ : FAULT PLANE "X" COORDINATE
! ET : FAULT PLANE "Y" COORDINATE
! Q  : DEPTH PARAMETER -> Q = x2*SIN(delta) - (d-x3)COS(delta)

  IF (ABS(RGZ) > CRT) THEN
     RGZLOG = LOG(RGZ)
  ELSE
! IF AT ORIGIN WE GET A SINGULARITY USING LOG(), SO SET TO ZERO
     RGZLOG = 0.
  END IF
  
! EQUATION (55) CALCULATION
  GFLT2 = 2.D0*FI2*COSDIP - RGZLOG*SINDIP
  RETURN
END FUNCTION GFLT2

! END OF GRAVITY SUBROUTIINES 
! -------------------------------------------------------------


SUBROUTINE CONVRT(XX,YY,XLON,YLAT,RLON,RLAT,NFLT)
  USE MAIN_PARS
  USE CONST
  IMPLICIT NONE 
  INTEGER :: i
  INTEGER,INTENT(IN) :: NFLT
  REAL :: FAC1
  
  REAL,INTENT(IN) :: RLON,RLAT
  REAL,DIMENSION(NF),INTENT(IN) :: XLON,YLAT
  REAL,DIMENSION(NF),INTENT(INOUT) :: XX,YY    
! THIS IS A SUBROUTINE FOR THE CONVERSION OF 
! CONVERSION FROM LATITUDE & LONGITUDE TO [km] VALUES
! 
! INPUT :
! YLAT  : LATITDE VALUES
! XLON  : LONGITUDE VALUES 
! RLAT  : REFERENCE LATITDE VALUE FROM WHICH km's WILL BE 
!         CALCULATED ==> USE BOTTOM LEFT COORDINATE VALUE
! RLON  : REFERENCE LONGITUDE VALUE FROM WHICH km's WILL BE 
!         CALCULATED ==> USE BOTTOM LEFT COORDINATE VALUE
! NOTE  : WEST LONGITUDE IS CONSIDERED *NEGATIVE*
! COS() : THIS TERM COMPUTES A SIMPLIFIED MERCATOR CORRECTION
! 
! OTPUT :
! XX    : LONGITUDE VALUES CONVERTED TO KILOMETRES
! YY    : LATITUDE VALUES CONVERTED TO KILOMETRES
! NOTE  : NSEGS SHOULD MATCH THE VALUE OF NF IN THE MAIN ROUTINE
!         SO AS TO KEEP THE DIMENSIONS CORRECT, OR AT LEAST AS LARGE  

! SET THE CONSTANT
  FAC1   = (RADIUS*PI)/180.0

! THE MAIN LOOP 
  DO i = 1,NFLT
     XX(i) = ((XLON(i) - RLON) * COS(YLAT(i)*RAD)) * FAC1
     YY(i) = (YLAT(i) - RLAT) * FAC1
  END DO
END SUBROUTINE CONVRT


! SIGNUM FUNCTION COMPUTES THE SIGN OF THE SUPPLIED VARIABLE
FUNCTION SGN(x)
  IMPLICIT NONE 
  REAL :: SGN
  REAL,INTENT(IN) :: x
  SGN = 1.
  IF(x < 0.)    SGN = -1.
  IF(x  ==  0.) SGN =  0.
END FUNCTION SGN
! -------------------------------------------------------------
