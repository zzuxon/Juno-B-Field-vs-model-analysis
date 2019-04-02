;+
; NAME:
; 
;   FOOTPOINT
;
; PURPOSE:
; 
;   This function calculates the magnetic surface footpoint of a 
;   location around Jupiter in System III coordinates.  It can also
;   return the trace of a magnetic field line.
;
; CALLING SEQUENCE:
; 
;   fp = FOOTPOINT(RR, Theta, Phi, Model, InStepSize, [Time])
;
; INPUTS:
;   
;   RR:         Radial distance from Jupiter's center to location to map.
;                 Must be given in Jovian radii.
;   Theta:      The co-latitude of location to map (0 to 180 degrees).
;   Phi:        The east longitude of location to map (0 to 360 degrees).
;   Model:      The name of the Jovian magnetic field model to use to
;                 trace the magnetic field.  Currently supports 'VIP4', 
;                 'VIT4', 'VIPAL', and 'Dipole'.
;   InStepSize: Either a scalar or a two element array.  If the step size
;                 does not vary, provide a scalar.  If it does, provide a
;                 two element array of the form [minStepSize, maxStepSize].
;                 If two element array given, but neither keyword BVARY nor
;                 RVARY are set, then the constant value minStepSize will 
;                 be used.  All step sizes must be given in Jovian radii.
;                 
; OPTIONAL INPUTS:
;  
;   Time:  The 'Khurana' model requires a time in J2000 seconds.  Can use
;            the ctimer or doyToJ2000 functions to obtain this time.
; 
; KEYWORD PARAMETERS:
; 
;   FULL:        If set, this function will return all coordinates of the 
;                  traced magnetic field line instead of just the final
;                  surface footpoint.
;   REVDIR:      If set, the direction along the magnetic field line is 
;                  reversed.  The default is to initially trace towards 
;                  the planet.
;   BVARY:       If set, the step size will vary between minStepSize and
;                  maxStepSize relative to the magnetic field strength.
;   RVARY:       If set, the step size will vary between minStepSize and
;                  maxStepSize relative to the radial distance to Jupiter's
;                  center.
;   RMIN:        If keyword RVARY is set, then this gives the radial distance,
;                  in Jovian radii, where the step size is equal to minStepSize.  
;                  If the radial distance goes below RMIN, then the step size can
;                  go below minStepSize.  (Default: 1)
;   RMAX:        If keyword RVARY is set, then this gives the radial distance,
;                  in Jovian radii, where the step size is equal to maxStepSize.  
;                  If the radial distance goes above RMAX, then the step size can
;                  go above maxStepSize.  (Default: 40)  
;   BMIN:        If keyword BVARY is set, then this gives the magnetic field
;                  strength, in nanoTesla, where the step size is equal to 
;                  maxStepSize.  If the magnetic field strength goes below BMIN, 
;                  then the step size can go above maxStepSize.  
;                  (Default: Magnetic Field Strength at Initial Location)
;   BMAX:        If keyword BVARY is set, then this gives the magnetic field
;                  strength, in nanoTesla, where the step size is equal to 
;                  minStepSize.  If the magnetic field strength goes above BMAX, 
;                  then the step size can go below minStepSize. (Default: 100000)
;   LINEAR:      If set, then the step size will scale linearly between minStepSize
;                  and maxStepSize. (If BVARY is set, LINEAR is default)         
;   LOGARITHMIC: If set, then the step size will scale logarithmically between 
;                  minStepSize and maxStepSize. 
;                  (If RVARY is set, LOGARITHMIC is default) 
;   INPLANET:    If set, then the initial RR is less than the final radius, generally
;                  the planet's surface.
;   FINALR:      The distance, in Jovian radii, to the final footpoint location. 
;                  If set to 0, then the oblate spheroid surface is used.  
;                  (Default: 0)
;   OPENRDIST:   The distance, in Jovian radii, that the
;                  trace will assume it is on an open field line, break from the
;                  trace and return [999, 999, 999].  (Default: 100)               
;   EQUATOR:     If set, then the field trace will go away from the planet until Br 
;                  changes direction indicating the magnetic equator or current sheet.
;
; OUTPUTS:
; 
;   A three element array of the form [rr, theta, phi] with the System III
;   coordinate where the magnetic field line that passes through the initial
;   location crosses the Jovian surface (an oblate spheroid).  As with the 
;   input coordinates, rr is in Jovian radii, theta in colatitude degrees,
;   and phi in east longitude degrees.
;
; OPTIONAL OUTPUTS:
; 
;   If keyword FULL is set, then the output is a 2D array with all the locations
;   of the magnetic field line between the initial location and the final surface
;   footpoint.  The array is Nx3, where N is the number of steps of the form
;   [stepNum, coordinate].  The array [*,0] gives all rr.  The 
;   array [*,1] gives all theta.  And the array [*,2] gives all phi.  The array
;   [stepNum, *] gives the three element coordinate array for the given stepNum.
;
; PROCEDURE:
; 
;   This function uses Euler's method to trace along vectors of the magnetic field
;   with a variable step size.  The step size varies in one of three ways:
;     Constant - If inStepSize is a single scalar or neither BVARY nor RVARY are set,
;                  then the step size will be constant throughout.
;     Linear   - If linear scaling is chosen, the step size will vary according to:
;                  stepSize = slope*(currentValue-minValue) + minStepSize, where
;                  slope = (maxStepSize-minStepSize)/(maxValue - minValue)
;                  For BVAR, minValue and maxValue are swapped.
;     Logarithmic - If logarithmic scaling is chosen, the step size will vary
;                  according to:
;                  stepSize = minStepSize * 10^(kk*(currentValue-minValue)), where
;                  kk = ALOG10(maxStepSize/minStepSize)/(maxValue - minValue)
;                  For BVAR, minValue and maxValue are swapped.
;
;   The required functions and procedures are:
;     jovMagField
;     
;   If keyword FINALR not specified or set to 0, this function requires:
;     jovSurface
;     
;   The models 'VIP4', 'VIT4', and 'VIPAL' require:
;     getSchmidtCoefs
;     getLegendre
;     dPdtheta
;     
;   The 'Khurana' model requires:
;     kk_2009
;
; EXAMPLE:
; 
;   To find the footpoint of the location at 3 Jovian radii, -70 degrees latitude,
;   and 30 degrees longitude using the 'VIPAL' model and a constant step size of 0.001,
;   use:
;   
;     fp = footpoint(3, 160, 30, 'VIPAL', 0.001)
;     fpR     = fp[0]
;     fpTheta = fp[1]
;     fpPhi   = fp[2]
;     
;     
;   To find the footpoint from location [10, 20, 0] using the 'VIT4' model while 
;   varying the step size logarithmically according to distance, with maximum step 
;   size (0.1) happening at apojove of 38 Jovian radii and minimum (0.001) at 
;   1 Jovian radius, use:
;   
;     fp = footpoint(10, 20, 0, 'VIT4', [0.001, 0.1], /RVARY, RMAX=38, /LOG)
;     
;     
;   To find all coordinates of an entire magnetic field line that extends to an
;   equatorial radius of 5 Jovian Radii, while varying the step size linearly 
;   according to magnetic field strength, with the model 'VIP4', use:
;   
;     fp1 = footpoint(6, 90, 0, 'VIP4', [0.0001, 0.01], /BVARY, /LIN, /FULL)
;     fp2 = footpoint(6, 90, 0, 'VIP4', [0.0001, 0.01], /BVARY, /LIN, /FULL, /REV)
;     
;   To extend a point on the jovian surface (say 70 degrees latitude and 180
;   degrees longitude) out onto a spherical surface of one Jovian radius, 
;   with the 'Khurana' model at the day of year 2007/103 07:27:30 (perijove 17 of
;   Juno), use:
;   
;     surfaceR = jovSurface(20)
;     time = doyToJ2000(2007, 103, 7, 27, 30)
;     fp = footpoint(surfaceR, 20, 180, 'Khurana', 0.0001, time, /INPLANET, FINALR=1)
;
;   To convert to xyz coordinates, assuming full trace in spherical coordinates 
;   stored in variable fp, use:
;   
;     fpXYZ = [[fp[*,0] * SIN(fp[*,1]/!radeg) * COS(fp[*,2]/!radeg)], $
;              [fp[*,0] * SIN(fp[*,1]/!radeg) * SIN(fp[*,2]/!radeg)], $
;              [fp[*,0] * COS(fp[*,1]/!radeg)]]
;     fpX = fpXYZ[*,0]
;     fpY = fpXYZ[*,1]
;     fpZ = fpXYZ[*,2]
;   
;
; MODIFICATION HISTORY:
;   Written by: Drake A. Ranquist, 12 Jun 2013
;     Added EQUATOR keyword.  26 May 2015
;-


function footpoint, rr, theta, phi, model, inStepSize, time, $
  FULL=full, REVDIR=revDir, BVARY=bvary, RVARY=rvary, $
  RMIN=Rmin, RMAX=rmax, BMIN=Bmin, BMAX=Bmax, $
  LINEAR=lin, LOGARITHMIC=log, $
  INPLANET=inPlanet, FINALR=finalR, MAXR=maxR, EQUATOR=equator, $
  OPENRDIST=openRDist, $
  GNM=Gnm, HNM=Hnm

;-------------Set Defaults------------------
IF NOT KEYWORD_SET(full) THEN full=0
IF NOT KEYWORD_SET(revDir) THEN revDir=0
IF NOT KEYWORD_SET(bvary) THEN bvary=0
IF NOT KEYWORD_SET(rvary) THEN rvary=0
IF NOT KEYWORD_SET(Rmin) THEN Rmin=1
IF NOT KEYWORD_SET(Rmax) THEN Rmax=40
;Bmin is set to Binit as default later in code
IF NOT KEYWORD_SET(bmax) THEN Bmax=100000
IF NOT KEYWORD_SET(lin) THEN lin=0
IF NOT KEYWORD_SET(log) THEN log=0
IF NOT KEYWORD_SET(inPlanet) THEN inPlanet=0
IF NOT KEYWORD_SET(finalR) THEN finalR=0
IF NOT KEYWORD_SET(maxR) THEN maxR=100
IF NOT KEYWORD_SET(openRDist) THEN openRDist=100 
equator= KEYWORD_SET(equator)

IF inPlanet OR equator THEN revDir=1

;---Set Constants---
DEG2RAD = !Pi/180

;-------Prepare Full Trace Array-------
IF full THEN BEGIN
  ;Array will automatically increase in size if more room is needed
  traceSize = 10000l
  nn = 0
  tracePos = DBLARR(traceSize, 3)
  ;Store initial position
  tracePos[nn,*] = [rr, theta, phi]
  nn++
ENDIF

;-----------------Can't Vary By Both B and R-----------------
IF bvary AND rvary THEN BEGIN
  print, "Can't vary both parameters!  Will vary by radius."
  bvary = 0
ENDIF

;---------Can't Scale Both Linearly and Logarithmically--------
IF lin AND log THEN BEGIN
  print, "Can't scale both linearly and logarithmically!"
  print, "If vary by B, will scale linearly."
  print, "If vary by R, will scale logramithically."
  lin = 0
  log = 0
ENDIF

;------Set Minimum and Maximum Step Sizes------
IF N_ELEMENTS(inStepSize) gt 1 THEN BEGIN
  minStepSize = MIN(inStepSize)
  maxStepSize = MAX(inStepSize)
ENDIF ELSE BEGIN
  minStepSize = inStepSize
  maxStepSize = inStepSize
ENDELSE

iceGiant = 0
If model eq 'Neptune' Or model eq 'Uranus' THEN iceGiant=1

;---Find Direction to Trace Field Toward Planet---
IF iceGiant THEN Binitial = magField(rr, theta, phi, model, Gnm, Hnm) $
 ELSE Binitial = jovMagField(rr, theta, phi, model, time)
IF Binitial[0] gt 0 THEN dir = -1 ELSE dir = 1
IF revDir THEN dir *= -1

;--------------Find Initial Field Magnitude-----------------
IF NOT KEYWORD_SET(Bmin) THEN BEGIN
  Bmin = SQRT(Binitial[0]^2 + Binitial[1]^2 + Binitial[2]^2)
ENDIF

;--Find Slope or Rate for Varying Step Size by Magnetic Field Strength--
IF bvary THEN BEGIN
  ;----Default to Scale Linearly----
  IF NOT (lin AND log) THEN lin=1
  IF lin THEN BEGIN
    slope = (maxStepSize-minStepSize)/DOUBLE(Bmin - Bmax)
  ENDIF ELSE BEGIN
    kk = ALOG(maxStepSize/DOUBLE(minStepSize))/(Bmin - Bmax)
  ENDELSE
ENDIF

;----Find Slope or Rate for Varying Step Size by Radial Distance----
IF rvary THEN BEGIN
  ;----Default to Scale Logarithmically----
  IF NOT (lin AND log) THEN log=1
  IF lin THEN BEGIN
    slope = (maxStepSize-minStepSize)/DOUBLE(Rmax - Rmin)
  ENDIF ELSE BEGIN
    kk = ALOG10(maxStepSize/DOUBLE(minStepSize))/(Rmax - Rmin)
  ENDELSE
ENDIF

;---Prepare Temporary Variables of While Loop---
nextR = rr
nextTheta = theta
nextPhi = phi
cont = 1

;---------Keep Tracing Field Until Surface Is Hit-----------
WHILE cont DO BEGIN
  
  ;------Find Magnetic Field Vector and Strength-------
  IF iceGiant THEN Bfield = magField(nextR, nextTheta, nextPhi, model, Gnm, Hnm) $
    ELSE Bfield = jovMagField(nextR, nextTheta, nextPhi, model, time)
  Bmag = SQRT(Bfield[0]^2 + Bfield[1]^2 + Bfield[2]^2)
  
  ;---Exit early if change in Br detected at equator---
  IF equator THEN BEGIN
    IF Bfield[0]*dir lt 0 THEN BREAK
  ENDIF
  
  ;-------Vary the Step Size if Necessary-------
  stepFactor = minStepSize
  IF bvary THEN BEGIN
    IF lin THEN BEGIN
      stepFactor = slope*(Bmag-Bmax) + minStepSize
    ENDIF ELSE IF log THEN BEGIN
      stepFactor = minStepSize * 10^(kk*(Bmag-Bmax))
    ENDIF
  ENDIF ELSE IF rvary THEN BEGIN
    IF lin THEN BEGIN
      stepFactor = slope*(nextR-Rmin) + minStepSize
    ENDIF ELSE IF log THEN BEGIN
      IF nextR GT Rmax THEN BEGIN
        Rmax += 10
        kk = ALOG10(maxStepSize/DOUBLE(minStepSize))/(Rmax - Rmin)
      ENDIF
      stepFactor = minStepSize * 10^(kk*(nextR-Rmin))
    ENDIF
  ENDIF
  
  ;--Vector of Correct Direction and Step Size for Next Step--
  stepVector = dir * Bfield/Bmag * ABS(stepFactor)
  
  ;---------Calculate the Next Position in Spherical Coordinates-----------
  ;nextR must be last because used for calculations of nextPhi and nextTheta
  nextPhi += stepVector[2]/(DEG2RAD * nextR * SIN(nextTheta*DEG2RAD)) 
  nextTheta += stepVector[1]/(DEG2RAD * nextR)
  nextR += stepVector[0]
  
  ;-------Add Next Position to the Full Trace Array----------
  IF full THEN BEGIN
    tracePos[nn,*] = [nextR, nextTheta, nextPhi]
    nn++
    
    ;------Increase Array Size if Full-------
    IF nn ge traceSize THEN BEGIN
      newArr = DBLARR(9*traceSize, 3)
      tracePos = [tracePos, newArr]
      traceSize *= 10
    ENDIF
  ENDIF
  
  
  ;Determine whether to continue the while loop
  IF nextR gt openRDist THEN return, [999, 999, 999]
  
  IF inPlanet THEN BEGIN
    IF finalR eq 0 THEN cont = nextR lt jovSurface(nextTheta) $
      ELSE cont = nextR lt finalR
      
  ENDIF ELSE IF equator THEN BEGIN
    ;Determined earlier in loop
    
  ENDIF ELSE IF finalR eq 0 THEN BEGIN
    cont = nextR gt jovSurface(nextTheta)
  
  ENDIF ELSE cont = nextR gt finalR AND nextR lt maxR

ENDWHILE


;----Return the Full Magnetic Field Line Trace----
IF full THEN BEGIN
  return, tracePos[0:nn-1, *]
ENDIF

;----Return the Footpoint on the Surface----
return, [nextR, nextTheta, nextPhi]



END







