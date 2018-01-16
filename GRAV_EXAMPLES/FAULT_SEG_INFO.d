!--------------------------------------------------------------------------------
! THIS IS AN EXAMPLE FILE FOR THE INPUT TO GRAV_GREEN.f
!
! THIS FILE IS USED AS INPUT FOR THE GREEN's FUNCTION CALCULATION AND THE GRAVITY
! CALCULATION IN GRAVITY.f90.  THE FILE NAME IS MANUALLY INPUT IN GRAV_GREEN.f90 
! AND AUTOMATICALLY CALLED BY GRAVITY.f90
!
! CREATED BY TYLER JOSEPH HAYES, JULY 19, 2006
!--------------------------------------------------------------------------------

! EXAMPLE WITH VARIABLE NAMES -> DO NOT USE THIS!!!!!!!

FAULT_SEG_INFO.d
dbn dtn xfwnt yfwnt xfent yfent slpvlS slpvlD slpvlT segdp
dbn dtn xfwnt yfwnt xfent yfent slpvlS slpvlD slpvlT segdp
dbn dtn xfwnt yfwnt xfent yfent slpvlS slpvlD slpvlT segdp
dbn dtn xfwnt yfwnt xfent yfent slpvlS slpvlD slpvlT segdp


!--------------------------------------------------------------------------------
! PARAMETER DEFINITIONS
!--------------------------------------------------------------------------------
FAULT_SEG_INFO.d = Header with file name
dbn	= depth to bottom [km] -> Z (+)ve DOWN
dtn	= fault top [km] -> Z (+)ve DOWN
xfwnt	= x-coordinates of "west" end (or north) [km]
yfwnt	= y-coordinates of "west" end (or north) [km]
xfent	= x-coordinates of "east" end (or south) [km]
yfent	= y-coordinates of "east" end (or south) [km]
	  NOTA BENE: WEST longitude is NEGATIVE
slpvlS	= STRIKE SLIP component of slip velocity [cm/yr]
		(+)ve indicates LEFT-LATERAL slip
slpvlD	= DIP SLIP component of slip velocity [cm/yr]
		(-)ve is used for a THRUSTING SLIP system
slpvlT	= TENSILE component of slip velcity [cm/yr]
		(+)ve indicates OPENNING 
segdp	= DIP ANGLE of the fault segment [degrees]
		90. -> PURE STRIKE SLIP fault

DIP ANGLE CONVENTION:
For a dip-slip system, such as the Aleutians, the magnitude of dip used is the same
as would be expected, but a NEGATIVE sign is attached to it. This triggers the intialization
of the code to rotate the observation grid by an extra 180 degrees so as to correctly assign 
the sense of motion with respect to the unrotated system. See MEDIUM_INFO.d for more
details.

CLARIFICATION OF slipvlD SIGN CONVENTION:
For a thrust fault system, the sense of motion during the STEADY-STATE (back-slip accumulation
phase) is the same as for a NORMAL FAULT. As such, for a THRUSTING fault system, a negative
value of slip is used to generate the initial unit gravity Green's functions. When the system
does finally slip elastically during an earthquake, the magnitude of the slip is ALSO ASSIGNED 
A NEGATIVE VALUE, thus recovering the POSTIVE SLIP VECTOR VALUE for the gravity calculations.
See FAULT_SLP_HIST.d for more details. 

FAULT LOCATION ASSIGNMENT:
One should note that for the fault coordinates, the values are defined by the surface
expression (or fault top) of the fault. This should be noted because the typical definition
found in Okada/Okubo use the fault bottom as the fault location reference. In the GREEN's
FUNCTION calculation, this discrepancy is accounted for by the DIPCOR variable. I use this
convention, because often what are most readily available are the fault surface expressions.


!--------------------------------------------------------------------------------


! EXAMPLE WITH NUMEBERS -> DO NOT USE THIS!!!!!!!

FAULT_SEG_INFO.d
15. 2. -20. 0. -10. 0. 1.2 0. 0. 90.
15. 2. -10. 0. 0. 0. 1.2 0. 0. 90.
15. 2. 0. 0. 10. 0. 1.6 0. 0. 90.
15. 2. 10. 0. 20. 0. 1.2 0. 0. 90.

!--------------------------------------------------------------------------------
! INTERPRETATION OF ABOVE FILE
!--------------------------------------------------------------------------------

(1) THERE ARE A TOTAL OF FOUR FAULT SEGMENTS IN THIS SYSTEM

(2) EACH SEGEMNT HAS A DEPTH OF 15.0[km] AND TOP AT 2.0[km]

(3) THE SEGMENTS SPAN THE LONGITUDE OF -20[E] TO 20[E] AND THE LATITUDES
    OF 0[N] -> ON THE EQUATOR, SPANNING 40 DEGREES OF LONGITUDE

(4) THE FIRST TWO SEGMENTS AND THE FOURTH SEGMENT  HAVE A STRIKE-SLIP 
    VELOCITY VECTOR OF 1.2[cm/yr] INDICATING A LEFT-LATERAL SYSTEM. THE 
    THIRD SEGMENT MOVES SLIGHTLY FASTER AT A RATE OF 1.6[cm/yr]. NOT ONE 
    OF THE SEGMENTS HAVE DIP-SLIP OR TENSILE VECTORS ASSOCIATED WITH THEM. 

(4) THE SYSTEM IS PURELY STRIKE SLIP FOR ALL FOUR SEGMENTS
