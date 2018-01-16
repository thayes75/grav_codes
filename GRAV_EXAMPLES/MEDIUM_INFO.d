!--------------------------------------------------------------------------------
! THIS IS AN EXAMPLE OF THE MEDIUM_INFO.d FILE WHICH
! CONTAINS THE SUBSURFACE MEDIUM PARAMETERS.  ALSO CONTAINS THE OBSERVATION
! GRID INFORMATION
!
! FOR SIGN CONVENTIONS FOR SLIP VECTORS SEE:
!
! SHUHEI OKUBO 
! /GRAVITY AND POTENTIAL CHANGES DUE TO SHEAR AND TENSILE FAULTS IN A HALF-SPACE/
! JGR, v.97, pp.7137--7144 [1992]
!
! THIS FILE IS USED AS INPUT FOR THE GREEN's FUNCTION CALCULATION AND THE GRAVITY
! CALCULATION IN GRAVITY.f90.  THE FILE NAME IS MANUALLY INPUT IN GRAV_GREEN.f90
! AND AUTOMATICALLY CALLED BY GRAVITY.f90
!
! CREATED BY TYLER JOSEPH HAYES, JULY 19, 2006
!--------------------------------------------------------------------------------


! EXAMPLE WITH VBARIABLE NAMES -> DO NOT USE THIS!!!!!!!

MEDIUM_INFO.d
amu alam bulkm
vplatS vplatD vplatT vplateX vplateY
ibcS ibcD ibcT b1 b2
dip1 dip2 dbplat
nfault KMTST
REFLON REFLAT
npoint IPTS JPTS XSTRT YSTRT


!--------------------------------------------------------------------------------
! PARAMETER DEFINITIONS
!--------------------------------------------------------------------------------
MEDIUM_INFO_FILE = header with the file name (NOTA BENE: use something 
	  descriptive, such as ALASKA_64.d, for the file name -> up to 30 chars)

Below are the underlying driving plate parameters
!--------------------------------------------------------------------------------
amu	= E.G.) 1. -> Rigidity (\mu) or Shear Modulus value for a Poisson solid
alam	= E.G.) 1. -> Lame Parameter (\lambda) value for a Poisson solid
bulkm	= E.G.) 1.66665995 -> Bulk Modulus (K) value for a Poisson solid
	  (NOTA BENE: The elastic parameters can set be to whatever medium
	   values you desire. The current incarnation uses a Poisson solid.)
vplatS	= Avg tectonic plate strike slip velocity in [cm/yr] -> LEFT-LATERAL (+) 
vplatD	= Avg tectonic plate dip slip velocity in [cm/yr] -> THRUSTING (+)
vplatT	= Avg tectonic plate tensile motion velocity in [cm/yr] -> OPENNING (+)
vplateX = Average SURFACE expression of the X-component of the tectonic 
	  plate velocity 
vplateY = Average SURFACE expression of the Y-component of the tectonic
	  plate velocity 
ibcS	= An INTEGER flag indicating which BOUNDARY the underlying 
	  plate strike-slip velocity (vplatS) is applied to:
		1=NORTH-WESTERN segment ONLY
		2=SOUTH-EASTERN segment ONLY
		3=BOTH NORTH-WESTERN and SOUTH-EASTERN boundaries
ibcD	= An INTEGER flag indicating which BOUNDARY the underlying 
	  plate dip-slip velocity (vplatD) is applied to:
		1=NORTH-WESTERN segment ONLY
		2=SOUTH-EASTERN segment ONLY
		3=BOTH NORTH-WESTERN and SOUTH-EASTERN boundaries
ibcT	= An INTEGER flag indicating which BOUNDARY the underlying 
	  plate tensile motion velocity (vplatT) is applied to:
		1=NORTH-WESTERN segment ONLY
		2=SOUTH-EASTERN segment ONLY
		3=BOTH NORTH-WESTERN and SOUTH-EASTERN boundaries
b1	= CHARACTER flag which indicates where the NORTH-WESTERN boundary
	  should go. MUST BE UPPER CASE!!!
		N=the first boundary is located on the NORTHERN most segment
		  in the fault system (e.g., San Andreas)
		W=the first boundary is located on the WESTERN most segment
		  in the fault system (e.g., Aleutians)

b2	= CHARACTER flag which indicates where the SOUTH-EASTERN boundary
	  should go. MUST BE UPPER CASE!!!
		S=the first boundary is located on the SOUTHERN most segment
		  in the fault system (e.g., San Andreas)
		E=the first boundary is located on the EASTERN most segment
		  in the fault system (e.g., Aleutians)
dip1	= Dip angle for the NORTH-WESTERN boundary
dip2	= Dip angle for the SOUTH-EASTERN boundary
dbplat 	= Depth to the bottom of the underlying plate

DIP ANGLE CONVENTION:
The observation grid is calculated assuming that X_2 is postive towards A (see figure below).
As such, for a dip-slip system such as the Aleutians, the magnitude of dip used is the same
as would be expected, but a NEGATIVE sign is attached to it. This triggers the intialization
of the code to rotate the observation grid by an extra 180 degrees so as to correctly assign 
the sense of motion with respect to the unrotated observation grid.


Below are the observation grid parameters
!--------------------------------------------------------------------------------
nfault 	= number of fault segments in the current model
KMTST  	= INTEGER switch to indicate if segment endpoint locations are 
	  provided in terms of kilometres(x,y) or latitude and longitude
		KMTST = 0 -> SEGMENT INFO IN KMs FROM A REFERENCE LAT/LON
		KMTST = 1 -> SEGMENT INFO IN LAT/LONs (preferred)
REFLON	= Reference longitude in terms of EAST [decimal degrees]
REFLAT	= Reference latitude in terms of NORTH [decimal degrees]
npoint 	= Number of observation points on the surface
IPTS	= Number of points for "i" in grid loop (i.e., X-direction)
JPTS	= Number of points for "j" in grid loop (i.e., Y-direction)
XSTRT	= Defines where to start the X values for the observation grid -> [km]
	  assuming that the origin is 0
YSTRT	= Defines where to start the Y values for the observation grid -> [km]
	  assuming that the origin is 0

	NOTA BENE: The loop to create the observation grid starts in the 
		NORTH-WEST corner, BUT the REFERENCE LONGITUDE and REFERENCE
		LATITUDE are in the SOUTH-WEST corner. See diagram below.


              A------------------| 	A = (XSTRT,YSTRT)   <- required units [km]
              |                  |	O = GRID ORIGIN = (0,0)  <- calculated
              |                  |      R = (REFLON,REFLAT) <- required units 
              |                  |                             [decimal degrees]
              |         O        | 
              |                  | 
              |                  | <- Observation grid, i.e, Earth's Surface
              |                  | 
              R------------------|   

 
!--------------------------------------------------------------------------------


! EXAMPLE WITH NUMEBERS -> DO NOT USE THIS!!!!!!!

MEDIUM_INFO_FILE.d
1. 1. 1.6666666598
-3.49657 0. 0. -2.36 2.58
3 3 3 N S
90. 90. 300.
50 1
-121 31.
10000 100 100 -600. 600.

!--------------------------------------------------------------------------------
! LINE-BY-LINE INTERPRETATION OF ABOVE FILE
!--------------------------------------------------------------------------------

(1) THE FILE IS CALLED MEDIUM_INFO_FILE.d

(2) THIS IS A POISSON SOLID WITH APPROPRIATE amu, alam, AND bulkm

(3) THE OVERALL UNDERLYING PLATE VELOCITY IS STRIKE SLIP ONLY WITH A MAGNITUDE OF 
    3.49657 [cm/yr] AND IS A RIGHT-LATERAL SYSTEM AS INDICATED BY THE 
    NEGATIVE SIGN.

(4) THE "SURFACE" EXPRESSION OF THE UNDERLYING PLATE STRIKE SLIP VELOCITY VECTOR 
    IN TERMS OF ITS X- AND Y-COMP. IS -2.36 [cm/yr] AND 2.58 [cm/yr] RESPECTIVELY.

(5) ALL THREE COMPONENTS (vplayS,vplatD,vplatT) ARE APPLIED TO BOTH BOUNDARIES.
    SINCE THE ONLY NON-ZERO COMPONENT IS vplatS, THE SYSTEM IS A PURE STRIKE-SLIP.
    IF THESE WERE DIFFERENT, ONE BOUNDARY COULD HAVE A DIP-SLIP VALUE, AND THE 
    STRIKE-SLIP.

(6) THE BOUNDARIES ARE LOCATED ON THE NORTHERN AND SOUTHERN MOST EXTENSIONS
    IN THE SYSTEM

(7) THE DIPS OF THE NORTH-WESTERN AND SOUTH-EASTERN BOUNDARY ARE BOTH 90. (PERFECT 
    STRIKE-SLIP SYSTEM.)

(8) THE DEPTH TO THE BOTTOM OF THE UNDERLYING PLATE IS SET AT 300[km].

(9) THERE ARE A TOTAL 50 FAULT SEGMENTS IN TERMS OF LAT/LONs GIVEN IN THE 
    SEGMENT INFORMATION FILE (TO BE READ IN LATER).

(10)THE OBSERVATION GRID WILL BE REFERENCED TO -121[E] LONGITUDE AND 31[N] LATITUDE.

(11)10000 POINT GRID THAT IS SQUARE (100 POINTS IN X-dir AND 100 IN Y-dir)
    AND THE SYSTEM STARTS IN THE TOP-LEFT CORNER WITH THE COORDINATES:
	X = -600 [km]
	Y = +600 [km]
	
	NOTA BENE: THESE COORDINATES ARE *ARBITRARY* AND USED FOR CALCULATING 
		THE GREEN's FUNCTION VALUES. THE ABSOLUTE CORRDINATES ARE 
		PRESERVED IN THE ROUTINE.
