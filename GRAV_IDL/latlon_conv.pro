;   PROCEDURE LATLON_CONV.PRO
;   
;   This procedure is designed to be called by other IDL procedures in
;   order to convert any xkm and ykm values used by VIRTUAL CALIFORNIA
;   to LATITUDE and LONGITUDE values.  This will allow the IDL mapping
;   plots to be called with ease.
;
;   ASSUMPTIONS:
;   (1) Spherical earth with R=6371.0[KM]
;   (2) All points given are in KM from 31N 121W
;   (3) Uses standard arclength formula: DIST = ANGDIST*R*pi/180
;       where ANGDIST is in units of RADIANS
;
;   Created by Tyler Hayes on May 9, 2006
;
;   INPUTS:
; 
;   xkms   - The xvalue to be converted to a LONGITUDE
;   ykms   - The yvalue to be converted to a LATITUDE
;   numpt  - The size of the vector being converted
;            NUMPT is typically:
;            1       = single value
;            nfault  = for fault segemnt conversions ==> 650
;            npoint  = for gravity/strain data observations
;            (10000,6400) points respectively
;
;   OUTPUT:
;  
;   xlon   - The LONGITUDE of the corresponding x-value XKMS
;   ylat   - The LATITUDE of the correspondin y-value YKMS
;   
;   NOTE:  It is assumed that the arrays xlon and ylat have been
;          initialized in the calling procedure
;           
;----------------------------------------------------------------------

PRO latlon_conv, xkms, ykms, numft, xlon, ylat, reflon, reflat

   ;; Set the variables
   NFT    = long(numft)
   radius = 6371.0 
   ;; reflat = 31.0
   ;; reflon = -121.0
   pi     = 3.141592654
   fac1   = 180.0/(radius*pi)
   fac2   = pi/180.0
   
; The main loop 
   FOR i = 0L,NFT-1 DO BEGIN
      ylat(i) = (ykms(i) * fac1) + reflat
      ;; NOTE:  I had to use -121 so as to account to the way IDL
      ;; refers to the WEST LONGITUDE.
      xlon(i) = reflon + ((xkms(i)*fac1)/(cos(ylat(i)*fac2)))
   ENDFOR
  
; Return to calling procedure
   RETURN
END
