   ;; PRO gravity90_view
   ;; This procedure is used to plot the gravity data from the FORTRAN
   ;; 90 routine GRAVITY.f90.  It requires the following files:
   ;; 
   ;; (1) gravity90.dat
   ;; (2) FAULT_SLP_HIST.d
   ;; (3) FAULT_SEG_INFO.d
   ;; (4) MEDIUM_INFO.d
   ;; (5) LATLON_CONV.PRO
   ;; (6) COLORBAR.PRO <- by David Fanning
   ;; (7) NUMBER_FORMATER.PRO <- by David Fanning
   ;; (8) The TeXtoIDL routines 
   ;; 
   ;; The David Fanning files can be obtained at the following URL:
   ;; http://www.dfanning.com/programs/
   ;; 
   ;; The TeXtoIDL routines can be obtained at the following URL:
   ;; http://physweb.mnstate.edu/mcraig/TeXtoIDL/
   ;; 
   ;; This file was created on Tuesday, August 8th, 2006 by Tyler
   ;; Hayes
   ;; ----------------------------------------------------------------

   ;; Get the fault history file name
   fnhist = ''
   
   ;; PRINT, ' '
   ;; PRINT, ' Enter the name of the slip history file'
   ;; PRINT, ' '
   ;; READ, fnhist

   ;; Get the fault segment geometry file name
   fngeom = ''

   ;; PRINT, ' '
   ;; PRINT, ' Enter the name of the fault segment geometry file'
   ;; PRINT, ' '
   ;; READ, fngeom

   ;; Get the medium information file name
   fninfo = ''

   ;; PRINT, ' '
   ;; PRINT, ' Enter the name of the medium information file'
   ;; PRINT, ' '
   ;; READ, fninfo

   ;; New for this code, get file names from gravity90.dat file
   PRINT, ''
   PRINT, 'Make sure you have all of your data file in the same '
   PRINT, 'working directory as gravity90.dat. As well as the'
   PRINT, 'required IDL *PRO* files.'
   PRINT, ''   
   PRINT, 'The data files you need are:'
   PRINT, ''
   PRINT, '(1) FAULT_SLP_HIST.d'
   PRINT, '(2) FAULT_SEG_INFO.d'
   PRINT, '(3) MEDIUM_INFO.D'
   PRINT, ''

   ;; Get the file names
   grvname = 'gravity90.dat'
   ;; print, 'Enter Name of File with GRAVITY DATA: '
   ;; read, grvname
   OPENR, 44, grvname
   READF, 44, fninfo
   READF, 44, fnhist
   READF, 44, fngeom
   READF, 44, tbase, tlater, tmin, npint
   
   ;; These arrays are used for calculation purposes	
   npoint = LONG(npint)
   NX = LONG(SQRT(npoint))
   NY = LONG(SQRT(npoint))
   
   xpt		=	FLTARR(npoint,/NO)
   ypt		=	FLTARR(npoint,/NO)
   grt     	=	FLTARR(npoint,/NO) ; Total gravity
   grdv    	=	FLTARR(npoint,/NO) ; Dilatational gravity
   grfa    	=	FLTARR(npoint,/NO) ; Free-Air component
   ;; grvo    	=	FLTARR(npoint,/NO) ; Cavity gravity
   grv          =       FLTARR(npoint,/NO) ; User gravity choice
   
   ;; These are the variables which are read in from gravity_data.d
   xpt_read	=	FLTARR(npoint,/NO)
   ypt_read	=	FLTARR(npoint,/NO)
   grt_read	=	FLTARR(npoint,/NO) ; Total gravity
   grdv_read	=	FLTARR(npoint,/NO) ; Dilatational gravity
   grfa_read	=	FLTARR(npoint,/NO) ; Free-Air component
   ;; grvo_read	=	FLTARR(npoint,/NO) ; Cavity gravity
	
   ;; Read the FORTRAN formatted data
   READF, 44, xpt_read
   READF, 44, ypt_read
   READF, 44, grt_read
   READF, 44, grdv_read
   READF, 44, grfa_read
   ;; READF, 44, grvo_read

   ;; Close gravity90.dat
   CLOSE, 44
	
   PRINT, ' '
   PRINT, 'Finished reading in gravity data created using the files:'
   PRINT, ''
   PRINT, 'MEDIUM INFO file mame  : ', fninfo
   PRINT, 'FAULT HISTORY file name: ', fnhist 
   PRINT, 'SEGMENT INFO file name : ', fngeom 
   PRINT, ''
   PRINT, 'Make sure these files are in your working diretory to continue.'
   PRINT, ''

   ;; Get the medium information file name
   fninfo = ''

   PRINT, ' '
   PRINT, ' Enter the name of the medium information file'
   PRINT, ' '
   READ, fninfo
   ;; Get the fault history file name
   fnhist = ''

   PRINT, ' '
   PRINT, ' Enter the name of the slip history file'
   PRINT, ' '
   READ, fnhist

   ;; Get the fault segment geometry file name
   fngeom = ''

   PRINT, ' '
   PRINT, ' Enter the name of the fault segment geometry file'
   PRINT, ' '
   READ, fngeom

   ;; Open and read the stuff in the files
   ;; MEDIUM INFORMATION FILE   
   ;; OPENR, 1, 'big_bear_med.d'
   ;; b1 = ''
   ;; b2 = ''

;;    Open file demo.doc for reading.  
;;    OPENR, 1, fninfo  
  
;;    Create a variable of type string.  
;;    LINE = ''  
  
;;    Read and print each line until the end of the file is encountered.   
;; WHILE(~ EOF(1)) DO BEGIN READF,1,LINE & PRINT,LINE & END  
  
;;    Done with the file.  
;;       CLOSE, 1   



   OPENR, 1, fninfo
   READF, 1, fninfo
   READF, 1, amu, alam, bulkm
   READF, 1, vplatS, vplatD, vplatT, vplateX, vplateY
   READF, 1, ibcS, ibcD, ibcT ; Stop reading since we do not need  b1, b2
   READF, 1, dip1, dip2, dbplat
   READF, 1, nfault, KMTST
   READF, 1, REFLON, REFLAT 
   READF, 1, npoint, IPTS, JPTS, XSTRT, YSTRT   
   CLOSE, 1

   ;; Create vectors for history data
   NFT   = 11000
   NTSTP = 11000
   NF    = LONG(nfault)
   NFF   = 500000
   NMF   = 7500000
   NPT   = 500000
   
   ;; History file vectors
   islip  = INTARR(NF,/NO)
   tslip  = FLTARR(NF,NFT,/NO)
   initeq = FLTARR(NF, NFT, /NO)
   dslpS  = FLTARR(NF,NTSTP,/NO)
   dslpD  = FLTARR(NF,NTSTP,/NO)
   dslpT  = FLTARR(NF,NTSTP,/NO)

   ;; Now open and read in data from the history file
   ;; FAULT_SLP_HIST.d
   ;; OPENR, 20, 'big_bear_hist.d'
   OPENR, 20, fnhist
   READF, 20, fnhist
   READF, 20, tlast, tstep, ntime
   READF, 20, islip

   ;; Read in the tslip array and put into 2 column array
   FOR i=0L,NF-1 DO BEGIN
      NFI = islip(i)            
      IF (NFI GT 0) THEN BEGIN
         tslipi  =  FLTARR(NFI,/NO)
         READF, 20, tslipi
         FOR j=0L,NFI-1 DO BEGIN
            tslip(i,j) = tslipi(j)
         ENDFOR
      ENDIF
   ENDFOR

   ;; Now read in dslpS 
   PRINT, ' ==>  Now reading in dslpS'
   nthi = ntime + 1
   IF (nthi LT NTSTP) THEN BEGIN
      ndim = (nthi-1)*NF
      dslpS_vector = FLTARR(ndim,/NO)
      READF, 20, dslpS_vector
      k = long(0)
      FOR j = 0L,nthi-2 DO BEGIN
         FOR i = 0, NF-1 DO BEGIN
            dslpS(i,j) = dslpS_vector(k)
            k = k+1
         ENDFOR
      ENDFOR
   ENDIF ELSE BEGIN
      ndim = long(ntime*NF)
      dslpS_vector = FLTARR(ndim,/no)
      READF, 20, dslpT_vector
      k = long(0)
      FOR j = 0L,ntime-1 DO BEGIN
         FOR i = 0, NF-1 DO BEGIN
            dslpS(i,j) = dslpS_vector(k)
            k = k+1
         ENDFOR
      ENDFOR
   ENDELSE

   ;; Now read in dslpD 
   PRINT, ' ==>  Now reading in dslpD'
   nthi = ntime + 1
   IF (nthi LT NTSTP) THEN BEGIN
      ndim = (nthi-1)*NF
      dslpD_vector = FLTARR(ndim,/NO)
      READF, 20, dslpD_vector
      k = long(0)
      FOR j = 0L,nthi-2 DO BEGIN
         FOR i = 0, NF-1 DO BEGIN
            dslpD(i,j) = dslpD_vector(k)
            k = k+1
         ENDFOR
      ENDFOR
   ENDIF ELSE BEGIN
      ndim = long(ntime*NF)
      dslpD_vector = FLTARR(ndim,/no)
      READF, 20, dslpD_vector
      k = long(0)
      FOR j = 0L,ntime-1 DO BEGIN
         FOR i = 0, NF-1 DO BEGIN
            dslpD(i,j) = dslpD_vector(k)
            k = k+1
         ENDFOR
      ENDFOR
   ENDELSE

   ;; Now read in dslpT 
   PRINT, ' ==>  Now reading in dslpT'
   nthi = ntime + 1
   IF (nthi LT NTSTP) THEN BEGIN
      ndim = (nthi-1)*NF
      dslpT_vector = FLTARR(ndim,/NO)
      READF, 20, dslpT_vector
      k = long(0)
      FOR j = 0L,nthi-2 DO BEGIN
         FOR i = 0, NF-1 DO BEGIN
            dslpT(i,j) = dslpT_vector(k)
            k = k+1
         ENDFOR
      ENDFOR
   ENDIF ELSE BEGIN
      ndim = long(ntime*NF)
      dslpT_vector = FLTARR(ndim,/no)
      READF, 20, dslpT_vector
      k = long(0)
      FOR j = 0L,ntime-1 DO BEGIN
         FOR i = 0, NF-1 DO BEGIN
            dslpT(i,j) = dslpT_vector(k)
            k = k+1
         ENDFOR
      ENDFOR
   ENDELSE

   ;; Now read in initiatot data
   PRINT, ' ==>  Now reading in initeq'
   FOR i=0L,NF-1 DO BEGIN
      NFI = islip(i)	
      IF (NFI GT 0) THEN BEGIN
         initeqi  =  intarr(NFI,/NO)
         READF, 20, initeqi
         FOR j=0L,NFI-1 DO BEGIN
            initeq(i,j) = initeqi(j)
         ENDFOR
      ENDIF
   ENDFOR
   PRINT,' '
   PRINT,' Finished reading in data from history file:  ', fnhist
   PRINT,' '
   CLOSE, 20   
   
   ;; Now open and read in data from the segment information file
   ;; FAULT_SEG_INFO.d
   ;; Create segment information file vectors/matrix
   db	  = FLTARR(NF,/NO)
   dt	  = FLTARR(NF,/NO)
   xfe	  = FLTARR(NF,/NO)
   yfe	  = FLTARR(NF,/NO)
   xfw	  = FLTARR(NF,/NO)
   yfw	  = FLTARR(NF,/NO)
   slpvlS = FLTARR(NF,/NO)
   slpvlD = FLTARR(NF,/NO)
   slpvlT = FLTARR(NF,/NO)
   dip    = FLTARR(NF,/NO)
   fname  = ''
   
   fault_data = FLTARR(10,NF)

   ;; Open and read in the data
   OPENR, 3, fngeom
   READF, 3, fngeom
   READF, 3, fault_data
   CLOSE, 3

   ;; Convert fault data matrix to individual vectors
   db     = fault_data(0,*)
   dt     = fault_data(1,*)	
   xfw    = fault_data(2,*)
   yfw    = fault_data(3,*)
   xfe    = fault_data(4,*)
   yfe    = fault_data(5,*)
   slpvlS = fault_data(6,*)
   slpvlD = fault_data(7,*)
   slpvlT = fault_data(8,*)
   dip    = fault_data(9,*)

   PRINT, 'Finished reading fault model input data files.'

   ;; Supress annoying underflow errors
   !EXCEPT = 0

   ;; Set the driver for creating an EPS figure for LaTeX
   fname6 = ''
   PRINT,''
   PRINT,'Now enter the name for the EPS output file:'
   PRINT,'(Do not forget to include the .eps extension)'
   PRINT,''
   SET_PLOT,'ps'
   DEVICE,/COLOR,BITS=8
   READ, fname6

   ;; Change ggv to gv depending on the ghostview you use
   view = 'ggv ' + fname6 + ' &'
   DEVICE, FILENAME=fname6, /ENCAPSULATED, PREVIEW=2

   ;; Begin the gravity plotting routine
   ;; Create arrays used in plotting
   XL	= FLTARR(NF,/NO)     
   XR	= FLTARR(NF,/NO)     
   YL	= FLTARR(NF,/NO)     
   YR	= FLTARR(NF,/NO)     

   ;; These are nolonger used in the routine
   ;; XINT	= FLTARR(NF,/NO)
   ;; YINT	= FLTARR(NF,/NO)

   ;; NOTA BENE: IDL uses C-style increments starting from zero
   FOR i = 0L,npoint-1 DO BEGIN
      xpt(i)  = xpt_read(i)
      ypt(i)  = ypt_read(i)
      grt(i)  = grt_read(i)
      grdv(i) = grdv_read(i)
      grfa(i) = grfa_read(i)
      ;; grvo(i) = grvo_read(i)
   ENDFOR

   ;; Create new varaiables. NOTE: that's an 'L' not a one!
   xpt_center_1	=	FLTARR(npoint,/NO)
   xpt_center_2	=	FLTARR(npoint,/NO)
   ypt_center_1	=	FLTARR(npoint,/NO)
   ypt_center_2	=	FLTARR(npoint,/NO)
   grv          =       FLTARR(npoint,/NO)
   gx           =       FLTARR(npoint,/NO)
   gy           =       FLTARR(npoint,/NO)

   ;; Reset the grid values to the centre of the fault grid
   FOR i = 0L, npoint-1 DO BEGIN
      xpt_center_1(i) = xpt(i)  ; This is the FAULT centred grid
      ypt_center_1(i) = ypt(i)  ; This is the FAULT centred grid
      gx(i) = 0.0
      gy(i) = 0.0
   ENDFOR

   ;; Prompt user to choose gravity plot type
   resp_S = ''
   PRINT, ' '
   PRINT, ' Enter the gravity component you wish to plot: '
   PRINT, ' '
   PRINT, ' [T] = TOTAL gravity change '
   PRINT, ' [D] = DILATATIONAL gravity component '
   PRINT, ' [F] = FREE-AIR gravity component '
   PRINT, ' '
   ;; Void formation not considered at the moment
   ;; print, ' [V] = GRAVITY FROM VOID FORMATION '
   READ, FORMAT='(a1)', resp_S

   ;; Set the gravity values to plotted. Only allow one type at a
   ;; time for now.
   grvchoices = ['T','D','F']
   choice = WHERE(STRUPCASE(resp_S) eq grvchoices)
   CASE choice[0] OF
      0: grv = grt              ; Total gravity
      1: grv = grdv             ; Dilatational gravity
      2: grv = grfa             ; Free-Air gravity
      ;; 3: grv = grvo          ; Void formation
      ELSE: BEGIN     
         grv = grt              ; Default is Total gravity 
         PRINT, ' That is an Invalid choice guy... '
         PRINT, ' TOTAL GRAVITY being used as default '
         PRINT, ' '     
      END
   ENDCASE 

   ;; Create max and min values for annotating plots
   grv_max1 = MAX(grv)
   grv_min1 = MIN(grv)
   ;; Use number_formatter.pro to round to 3 decimal places 
   ;; for printing purposes
   grv_max  = NUMBER_FORMATTER(grv_max1,DECIMALS=3)
   grv_min  = NUMBER_FORMATTER(grv_min1,DECIMALS=3)

   ;; Create the new array's with REBIN
   X_ARRAY	=	FLTARR(NX,NY,/NO)
   Y_ARRAY	= 	FLTARR(NY,NY,/NO)
   GRAVITY  	=	FLTARR(NX,NY,/NO)
   xptl	= FLTARR(npoint,/NO)        ; LONGITUDE
   yptl	= FLTARR(npoint,/NO)        ; LATITUDE
  
   ;; Convert the x- and y-km gravity grid values to 
   ;; longitude values using CAL_LATLON.pro
   LATLON_CONV, xpt, ypt, npoint, xptl, yptl, REFLON, REFLAT

   ;; Assign values to the new arrays
   m = long(-1)
   FOR i=0L,NX-1 DO BEGIN
      FOR j = 0L, NY-1 DO BEGIN
         m = long(m + 1)
         X_ARRAY(i,j) = xptl(m)
         Y_ARRAY(i,j) = yptl(m)
         GRAVITY(i,j) = grv(m)
      ENDFOR
   ENDFOR

   ;; REBIN into a finer mesh
   rebinfac        =  15
   newx            =  NX * rebinfac
   newy            =  NY * rebinfac
   NEW_X_ARRAY     =  REBIN(X_ARRAY,newx,newy)
   NEW_Y_ARRAY     =  REBIN(Y_ARRAY,newx,newy)
   NEW_GRAVITY     =  REBIN(GRAVITY,newx,newy)

   ;; Create the colour scheme for the contoured gravity data


   ;; ======================================================================
   ;; Set up for screen presentation (Black background!)
   ;; USAGE: Uncomment the code below if you intend to use 
   ;; in a presentation
   ;; BACKGROUND ===> BLACK
   ;; FOREGROUND ===> YELLOW 
   ;; FOR POWERPOINT USE "CONVERT" AS BELOW:
   ;; convert -verbose -antialias -density 600 -quality 95 foo.eps
   ;;     foo.png
   ;; 
   ;; Or do a search on the Google group "comp.lang.idl-pvwave" for 
   ;;  "IDL eps outpt converted to png for MS PowerPoint"
   ;; 
   ;; My posting there has a script to convert to PNG for PowerPoint
   ;; ======================================================================
   resp_PLOT = ''
   PRINT,' '
   PRINT,'Please enter: '
   PRINT,' '
   PRINT,'[S] if you wish to create plots for a screen '
   PRINT,'    presentation (i.e., BLACK BACKGROUND). '
   PRINT,' '
   PRINT,'[P] if you wish to create COLOUR plots for print in '
   PRINT,'    a LaTeX document  (i.e., WHITE BACKGROUND). '
   PRINT,' '
   PRINT,'[B] if you wish to create B&W plots for print '
   PRINT,'    in a publication (i.e., WHITE BACKGROUND). '
   PRINT,' '
   PRINT,'Otherwise press <ENTER> for DEFAULT (i.e., [P])'
   PRINT,' '
   PRINT,'Nota bene: CASE SENSITIVE!'

   READ, FORMAT='(a1)', resp_PLOT

   PRINT,' '
   plotchoice = ['S','P', 'B']
   pc = where(strupcase(resp_PLOT) EQ plotchoice)
   CASE pc[0] OF
      0: BEGIN
         ctable = 26            ; Eos A <- I like it....
         LOADCT, ctable
         GAMMA_CT, 2.089, /CURRENT ; this stretches the colours to have 
                                ; a more "even" distribution
         !P.CHARTHICK=3.0       ; To plot annotations thick
         TVLCT, 95, 95, 95, 2   ; Redefine my own lighter grey
         TVLCT, 18, 18, 18, 3   ; Redefine my own darker grey
         TVLCT, 255, 255, 0, 158; Redefine my own yellow
         backColour  = 0        ; Black background
         !P.COLOR    = 158      ; Sets Fonts to yellow using CT = 26
         image = REPLICATE(backColour, 10, 10)
         TV, image, XSIZE=!D.X_SIZE, YSIZE=!D.Y_SIZE
         IF (!P.COLOR EQ 158) THEN BEGIN ; Reset for black background
            rval = 65           ; River colour
            cval = 2            ; Continent outline colour
            gval = 3            ; Gridline colour
         ENDIF
      END
      1: BEGIN
         ctable = 26            ; Eos A <- I like it....
         LOADCT, ctable
         GAMMA_CT, 2.089, /CURRENT ; this stretches the colours to have 
                                ; a more "even" distribution
         !P.CHARTHICK=3.0       ; To plot annotations thick
         PRINT,'Creating COLOUR plot suitable for print or LaTeX'
         TVLCT, 95, 95, 95, 26  ; Redefine my own lighter grey
         !P.COLOR    = 0        ; Sets Fonts & outlines to black
         rval = 65              ; River colour
         cval = 26              ; Continent outline colour
         gval = 0               ; Grid colour
      END
      2: BEGIN
         ctable = 0             ; B&W colour scheme....
         LOADCT, ctable
         !P.CHARTHICK=3.0       ; To plot annotations thick
         PRINT,'Creating B&W plot suitable for print or LaTeX'
         ;; Shift middle colours only
         FOR i = 1, 30 DO BEGIN
            ncol = (i*8) - 1
            ocol = (i-1)*8 - 1
            IF (ocol LT 0) THEN ocol = 0
            TVLCT, ocol, ocol, ocol, ncol-1
            TVLCT, ocol, ocol, ocol, ncol+1
            TVLCT, ocol, ocol, ocol, ncol
         ENDFOR
         !P.COLOR    = 0        ; Sets Fonts & outlines to black
         rval = 65              ; River colour
         cval = 26              ; Continent outline colour
         gval = 0               ; Grid colour
      END
      ELSE: BEGIN
         PRINT,'That is an invalid choice guy...'
         PRINT,'Creating plot suitable for print or LaTeX'
         TVLCT, 95, 95, 95, 26  ; Redefine my own lighter grey
         !P.COLOR    = 0        ; Sets Fonts & outlines to black
         rval = 65              ; River colour
         cval = 26              ; Continent outline colour
         gval = 0               ; Grid colour
      END
   ENDCASE
   ;; END OF SCREEN PRESENTATION CODE
   ;; ===================================================================

   ;; Choose the zoomed up region for the contour plot
   resp_EXMPL = ''
   PRINT,' '
   PRINT,'PLEASE CHOOSE THE EXAMPLE YOU WISH TO PLOT: '
   PRINT,' '
   PRINT,'[A] ALASKA'
   PRINT,'[B] JOSHUA TREE-LANDERS-HECTOR MINE'
   PRINT,' '
   PRINT,' OR '
   PRINT,' '
   PRINT,'[Z] CHOOSE YOUR OWN LATITUDE/LONGITUDE VALUES'   
   PRINT,' '
   READ, FORMAT='(A1)', resp_EXMPL
   PRINT,' '
   PRINT,' NOTE 1: You will want to edit this file for your own'
   PRINT,'         plotting purposes. These values are used for the'
   PRINT,'         example files.'
   PRINT,' '
   PRINT,' NOTE 2: You may have to adjust the latitude and longitude'
   PRINT,'         spacing variables LATDEL and LONDEL internally'
   PRINT,'         to make nicer plots.'
   PRINT,' '
   example_choice = ['A','B','Z']
   zc = WHERE(STRUPCASE(resp_EXMPL) EQ example_choice)

   CASE zc[0] OF
      0: BEGIN
         ;; ALASKA
         X0z  = -166.
         X1z  = -135.
         Y0z  =  51.
         Y1z  =  62.
         xgrv = [X0z, X1z] ; LONGITUDE
         ygrv = [Y0z, Y1z]   ;   LATITUDE
      END
      1: BEGIN
         ;; BIG BEAR
         X0z  = -118. 
         X1z  = -115.
         Y0z  =  33.
         Y1z  =  36.
         xgrv = [X0z, X1z]      ; LONGITUDE
         ygrv = [Y0z, Y1z]   ;   LATITUDE
      END
      2: BEGIN 
         ;; USER CHOICE
         PRINT, 'ENTER BOTTOM LEFT LONGITUDE'
         PRINT, ' '
         READ, X0z
         PRINT, 'ENTER BOTTOM LEFT LATITUDE'
         PRINT, ' '
         READ, Y0z
         PRINT, 'ENTER UPPER RIGHT LONGITUDE'
         PRINT, ' '
         READ, X1z
         PRINT, 'ENTER UPPER RIGHT LATITUDE'
         PRINT, ' '
         READ, Y1z
         xgrv = [X0z, X1z]
         ygrv = [Y0z, Y1z]
      END
      ELSE: BEGIN
         print,'DEFAULTING TO ALASKA'
         X0z  = -166.
         X1z  = -135.
         Y0z  =  51.
         Y1z  =  62.
         xgrv = [X0z, X1z] ; LONGITUDE
         ygrv = [Y0z, Y1z]   ;   LATITUDE
      END
   ENDCASE 

   ;; Centre of the plot region
   rot_map    =  0.0
   lon_centre = (xgrv[0] + xgrv[1])/2
   lat_centre = (ygrv[0] + ygrv[1])/2  
   ;; LIMIT and CLIP vectors
   plimit = [ygrv[0], xgrv[0], ygrv[1], xgrv[1]] ; For the PLOT

   ;; Get the main title for the map
   mtitle1 = 'GRAVITY CONTOUR PLOT!C'
   ;; PRINT, 'Please enter the MAIN title for the map'
   ;; READ, mtitle1
   ;; mtitle = mtitle1 + 'C!'      ; Include a space for the 
   ;; Initialise the PLOT parameters
   !P.POSITION = [0.1, 0.24, 0.9, 0.9] ; For use with Colour bar
   ;; !P.POSITION = [0.1, 0.1, 0.9, 0.9] ; For no colour bar on plot

   ;; Set INITIAL mapping plot & style to be used
   ;; Add "/ISOTROPIC" here if you desire even lat/lon ratio for plot
   MAP_SET, lat_centre, lon_centre, /MERCATOR, LIMIT = plimit, $
    TITLE = mtitle1, /NOERASE
   
   ;; NOW CONTOUR THE GRAVITY VALUES
   ;; Must define the levels used for contours
   nlev      = 10               ; Total number of contours
   lev       = nlev - 2         ; For the contours which are plotted
                                ; -2 b/c MAX and MIN are not plotted
   cont_lev  = lev + 1          ; For the colour bar "tickmarks"
                                ; cont_lev = the number of INTERVALS

   ;; Label all values except the max & min using vector of 1's
   label_tmp = REPLICATE(0,lev) ; change 0 to 1 for labels
   g_labels  = [0, label_tmp, 0] 

   ;; Create matrix for levels
   grv_tot = ABS(grv_max) + ABS(grv_min) ; B/c range is (+) to (-)
   grv_int = grv_tot/(cont_lev)
   grv_lev = FLTARR(nlev,/no)
   grv_lev = INDGEN(nlev)*grv_int + grv_min

   ;; Test that the contour levels are correct
   ;; Uncomment below if you wish to see 10 even spaced contours
   PRINT, '--------------------------------------------'
   PRINT, 'The minimum gravity value is:', grv_min
   PRINT, 'The maximum gravity value is:', grv_max
   ;; PRINT, ''
   ;; PRINT,'Ten equally spaced contour intervals are:'
   ;; FOR j = 0L,9L DO BEGIN
   ;;   PRINT, grv_lev(j)         ;, g_labels(j)
   ;; ENDFOR
   PRINT, '--------------------------------------------'
   PRINT, ' '
   PRINT, 'NOW CONTOURING GRAVITY DATA. MAY TAKE SOME TIME...'

   ;; Now must define num_colours for the shading between contours
   ;; bott        = 3    ; Starting colour value from 0 to 255 (no gamma)
   bott        = 31   ; Starting colour value from 0 to 255 (w/ gamma)
   nc_lev      = 256 - bott   ;  The number of gravity contours
   num_colours = nc_lev - 1   ;  The number of colours (the intervals)
   col_int     = (ABS(grv_max) + ABS(grv_min)) / (num_colours) 
   col_lev     = FINDGEN(nc_lev) * col_int + grv_min
   c_colours   = FINDGEN(nc_lev) * (nc_lev/num_colours) + bott

   ;; Plot the gravity colours using c_colours and col_lev
   CONTOUR, NEW_GRAVITY, NEW_X_ARRAY, NEW_Y_ARRAY,$
    LEVELS = col_lev, C_COLORS=c_colours, $
    /FILL, /OVERPLOT, /NOERASE

   ;; Overplot the the desired MAP datasets
   ;; MAP_CONTINENTS, /HIRES
   MAP_CONTINENTS, /HIRES, /COASTS, COLOR = cval 

   ;; Place gridlines on the plot in DARK colour
   MAP_GRID, /BOX_AXES, LATDEL = 0.5, LONDEL = 0.5, CHARSIZE = 0.75, $
    COLOR = gval
   ;; Overlay the box axes in the ANNOTATION colour
   MAP_GRID, /BOX_AXES, LATDEL = 0.5, LONDEL = 0.5, CHARSIZE = 0.75, $
    COLOR = !P.COLOR, /NO_GRID

   ;; Create LaTeX title for the colourbar
   ctitle = TEXTOIDL('\muGal (10^{-8}m/s^{2})')
   ;; Place colourbar below the figures. David Fanning is a dude.
   COLORBAR, NCOLORS=num_colours,TITLE=ctitle,$ 
    POSITION=[0.1, 0.1, 0.9, 0.14], DIVISIONS=n_lev,$
    RANGE=[MIN(grv),MAX(grv)], BOTTOM=bott,$
    FORMAT='(F7.1)'

   ;; Set unit string value with LaTeX style
   units = TEXTOIDL('\muGal')
   ;; Choose the title name from the the original gravity choice
   CASE choice[0] OF
      0:BEGIN
         astr0 = 'Total Gravity Over Years:'
         astr1 = strcompress(string(tbase)) + $
          ' -' + strcompress(string(tlater))
         astr2 = 'Max. =' + strcompress(string(grv_max)) + units
         astr3 = ' Min. =' + strcompress(string(grv_min)) + units
         astr01 = astr0 + astr1
         astr23 = astr2 + '  ' + astr3                 	
      END
      1:BEGIN
         astr0 = 'Dilatational Gravity Over Years:'
         astr1 = strcompress(string(tbase)) + $
          ' -' + strcompress(string(tlater))
         astr2 = 'Max. =' + strcompress(string(grv_max)) + units
         astr3 = ' Min. =' + strcompress(string(grv_min)) + units
         astr01 = astr0 + astr1
         astr23 = astr2 + '  ' +  astr3
      END
      2:BEGIN 
         astr0 = 'Free-Air & Bouguer Gravity Over Years:'
         astr1 = strcompress(string(tbase)) + $
          ' -' + strcompress(string(tlater))
         astr2 = 'Max. =' + strcompress(string(grv_max)) + units
         astr3 = ' Min. =' + strcompress(string(grv_min)) + units
         astr01 = astr0 + astr1
         astr23 = astr2 + '  ' + astr3                
      END
;;        3:BEGIN                ; Cavity case  
;;           astr0 = 'Gravity From Cavity Formation Over Years:'
;;           astr1 = strcompress(string(tbase)) + $
;;            ' - ' + strcompress(string(tlater))
;;           astr2 = 'Max. = ' + strcompress(string(grv_max))
;;           astr3 = 'Min. = ' + strcompress(string(grv_min)) + units
;;           astr23 = astr2 + '  ' + astr3
;;        END
      ELSE:BEGIN                ; Default case
         astr0 = 'Total Gravity Over Years:'
         astr1 = strcompress(string(tbase)) + $
          ' -' + strcompress(string(tlater))
         astr2 = 'Max. =' + strcompress(string(grv_max)) + units
         astr3 = ' Min. =' + strcompress(string(grv_min)) + units
         astr01 = astr0 + astr1
         astr23 = astr2 + '  ' + astr3      
      END
   ENDCASE

   ;; Append Subtitle, max/min, and date information
   ;; Below are the LAT/LON coords for the plot annotation
   xs0 = X0z - (X0z *.001)
   ys0 = Y0z + (Y0z *.005)
   XYOUTS, xs0, ys0, astr01, charsize  = 1.0, COLOR=0
   ;; Forego MAX/MIN annotation for this map b/c it becomes too busy
   ;; looking, not to mention that the colourbar has them anyways.
   ;; XYOUTS, xs2, ys2, astr23, charsize = 1.00, COLOR=!P.COLOR

   ;; This is the part which adds each fault patches over the data
   ;; First Change Colour Table back to Black and White for patches
   LOADCT, 0
   fval    = 0
   slpval  = 255
   l_thick = 14.0      ; line thickness for oplot <-- slipped patches
   l_thin  = 3.0       ; line thickness for oplot <-- non-slipped patches
   xle	= FLTARR(NF,/NO)        ; EAST segment LONGITUDE
   yle	= FLTARR(NF,/NO)        ; EAST segment LATITUDE
   xlw	= FLTARR(NF,/NO)        ; WEST segment LONGITUDE
   ylw	= FLTARR(NF,/NO)        ; WEST segment LATITUDE
   x1   = FLTARR(2,/NO)         ; Individual OPLOT value
   y1   = FLTARR(2,/NO)         ; Individual OPLOT value

   ;; Convert the fault segment data into the proper latitude and
   ;; longitude values using LATLON_CONV.PRO if KMTST = 1 but if they
   ;; are already in LAT/LON notation continue on
   IF (KMTST EQ 0) THEN BEGIN
      LATLON_CONV, xfe, yfe, nfault, xle, yle, REFLON, REFLAT
      LATLON_CONV, xfw, yfw, nfault, xlw, ylw, REFLON, REFLAT
      ;; Initialise the segment vectors output from above to use OPLOT
      XL = xlw
      XR = xle
      YL = ylw
      YR = yle
   ENDIF ELSE BEGIN
      ;;  Initialise the segment vectors already in LAT/LON for OPLOT
      XL = xfw
      XR = xfe
      YL = yfw
      YR = yfe
   ENDELSE

   ;; Call the OPLOT routine to overplot ALL faults
   PRINT, 'Plotting all the faults'
   PRINT, ' '
   FOR i=0L,NF-1 DO BEGIN
      x1(0) = XL(i)
      y1(0) = YL(i)
      x1(1) = XR(i)
      y1(1) = YR(i)
      OPLOT, x1, y1, LINESTYLE = 0, THICK = l_thin, COLOR = fval
   ENDFOR

   ;; Draw faults that slipped
   time_mat = FLTARR(NF,/no)
   init_vec = FLTARR(NF,/no)

   FOR fault_index = 0L,NF-1 DO BEGIN
      time_mat(fault_index) = 0
      init_vec(fault_index) = 0 
   ENDFOR

   ;; Find segments which slipped	
   FOR fault_index = 0L,NF-1 DO BEGIN
      FOR time_index = 0L,islip(fault_index) -1 DO BEGIN
         time_test_1 = tslip(fault_index,time_index)- tstep*.001
         time_test_2 = tslip(fault_index,time_index)- tstep*.001
         IF (time_test_1 GT tbase AND time_test_2 LE tlater) THEN BEGIN
            time_mat(fault_index) = 1
            tmp  = time_mat(fault_index)
            ;; Initiator code below.
            IF (tmp EQ 1 AND initeq(fault_index,time_index) EQ 1) THEN $
             init_vec(fault_index) = 1
         ENDIF		
      ENDFOR
   ENDFOR

   ;; Draw slipped segments in THICK BLACK
   PRINT, 'Plotting the slipped faults'
   PRINT, ' '
   FOR i=0L,NF-1 DO BEGIN
      IF (time_mat(i) EQ 1) THEN BEGIN
         x1(0) = XL(i)
         y1(0) = YL(i)
         x1(1) = XR(i)
         y1(1) = YR(i)
         OPLOT, x1, y1, LINESTYLE=0, THICK=l_thick, COLOR = fval
      ENDIF
   ENDFOR

   ;; Plot J-L-H earthquake epicentres in new colour
   ;; LOADCT, 12
   PRINT, 'Plotting the earthquake epicentre(s)'
   PRINT, ' '
   PRINT, 'EXAMINE CODE TO COMMENT OUT J-L-H EPICENTRES. '
   PRINT, 'See documentation in code lines 890 -- 937 '
   PRINT, 'of this file. '
   PRINT, ' '

   jlat =   33.9600
   jlon = -116.3200
   llat =   34.2000
   llon = -116.4317
   blat =   34.1833
   blon = -116.8167
   hlat =   34.5923
   hlon = -116.2750
   
   ;; Label locations for J-L-H names
   jlatt = jlat - 0.030
   jlont = jlon + 0.055
   llatt = llat
   llont = llon + 0.035
   blatt = blat + 0.030
   blont = blon - 0.12
   hlatt = hlat
   hlont = hlon + 0.025
   
   ;; Diamond shape
   XDIA = [-1, 0, 1, 0, -1]
   YDIA = [0, 1, 0, -1, 0]  
   USERSYM, XDIA, YDIA, /FILL

   ;; Star shape if you prefer
   ;; XSTAR = [0.0, 0.5, -0.8, 0.8, -0.5, 0.0]*2.
   ;; YSTAR = [1.0, -0.8, 0.3, 0.3, -0.8, 1.0]*2.   
   ;; USERSYM, XSTAR, YSTAR, /FILL

   ;; Use either fucshia or white if you want.
   ;; 128 = fuschia
   ;; 255 = white
   f = 128
   w = 255

   ;; Use either star or diamond, not both
   PLOTS, jlon, jlat, PSYM=8, COLOR=w
   PLOTS, llon, llat, PSYM=8, COLOR=w
   PLOTS, blon, blat, PSYM=8, COLOR=w
   PLOTS, hlon, hlat, PSYM=8, COLOR=w

   XYOUTS, jlont, jlatt, 'JT', charsize  = 1.0, COLOR=0
   XYOUTS, llont, llatt, 'L' , charsize  = 1.0, COLOR=0
   XYOUTS, blont, blatt, 'BB', charsize  = 1.0, COLOR=0
   XYOUTS, hlont, hlatt, 'HM', charsize  = 1.0, COLOR=0
   ;; END OF THE J-L-H STUFF --------------------------

   ;; Draw slipped segments in THICK BLACK and the 
   ;; INITATORS in THICK WHITE

   ;; Un-/Comment this part of the code if you want the 
   ;; initiators to show up as white fault patehes. 
   ;; For the J-L-H sequence, we use the white diamonds 
   ;; instead, so no need for this. But nice to have just 
   ;; in case.
;;    PRINT, 'Plotting the initiator faults'
;;    PRINT, ' '
;;    FOR i=0L,NF-1 DO BEGIN
;;       IF (init_vec(i) EQ 1) THEN BEGIN
;;          x1(0) = XL(i)
;;          y1(0) = YL(i)
;;          x1(1) = XR(i)
;;          y1(1) = YR(i)
;;          OPLOT, x1, y1, LINESTYLE=0, THICK=l_thick, COLOR = slpval
;;       ENDIF
;;    ENDFOR

   ;; Must close the PS device to get proper output
   DEVICE,/CLOSE_FILE	

   ;; Use below if you wish to have IDL visualize the EPS
   ;; file right away. I think if you switch to XTERM as the
   ;; device, it might default to creating a window, but I
   ;; prefer EPS figures for LaTeX.
   SPAWN, view 

   ;; Exit program
END
