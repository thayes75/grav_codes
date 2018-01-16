!--------------------------------------------------------------------------------
! THIS IS AN EXAMPLE FILE FOR THE INPUT TO GRAV_GREEN.f
!
! THIS FILE IS USED AS INPUT FOR THE GRAVITY CALCULATION IN GRAVITY.f90 AND
! CONTAINS THE INFORMATION FOR THE SLIP HISTORY OF EACH SEGMENT DEFINED IN 
! FAULT_SEG_INFO.d. 
!
! CREATED BY TYLER JOSEPH HAYES, JULY 19, 2006
!--------------------------------------------------------------------------------


! EXAMPLE WITH VBARIABLE NAMES -> DO NOT USE THIS!!!!!!!

FAULT_SLP_HIST.d
tlast tstep ntime
islip(i)
tslip(i,j)
dslpS(i,j)
dslpD(i,j)
dslpT(i,j)

!--------------------------------------------------------------------------------
! PARAMETER DEFINITIONS
!--------------------------------------------------------------------------------

FAULT_SLP_HIST.d = header with the file name -> up to 30 char. long
tlast	= the year of the last time recorded for the history file
tstep	= the interval of the time step in years (NOTA BENE: can be decimal)
ntime	= number of time steps (INTEGER)
islip	= number of times segemnt i slipped
tslip	= time of the slip of segment i for the jth slip <- from islip
dslpS	= total accumulated slip on segment i at time step j -> STRIKE SLIP	
dslpD	= total accumulated slip on segment i at time step j -> DIP SLIP	
dslpT	= total accumulated slip on segment i at time step j -> TENSILE MOTION

NOTA BENE I: islip, tslip, dslpS, dslpD, and dslpT are all formatted as 
	FORTRAN LISTED data via WRITE(6,*) statements. 	Also, i=nfault and 
	j=ntime with each j-value corresponding to an interval of tstep in time.

NOTA BENE II: The values of ntime and nfault are typically smaller than the
	sizes set aside in the modules DATA_DEFS and MAIN_PARS. This is done on
	purpose so as to facilitate the following two features:
		 (i) the auotmatic boundary segment calculations which add two
		     additonal segments for the driver plate steady state 
	 	     calculations.
		(ii) to prevent continual modification of the code itself for
		     each new system modelled.  By keeping the values somewhat
		     large, numerous fault systems can be modelled with the same 
		     compiled code by only modifying the input files.

DIP-SLIP SLIP VECTOR CONVENTION:
In the FAULT_SEG_INFO.d file, it was explained that the STEADY-STATE component of slip for a 
REVERSE FAULT is set to a NEGATIVE VALUE. This was to facilitate the proper sense of 
accumulated backslip in the system. That is, the reverse fault system in the steady-state 
behaves as though it was a NORMAL FAULT. As such, we must then MULTIPLY the magnitude of
thrust faulting in the slip history file by -1.0 to recover the sense of motion for a thrust
(reverse) fault. See the MAKE_HIST.f file in the ALASKA example directory. 
