;+
; NAME:
; 
;   JOVMAGFIELD
;
; PURPOSE:
; 
;   This function returns the radial, theta, and phi components of the magnetic
;   field at a specific point given a Jovian magnetic field model.
;
; CALLING SEQUENCE:
; 
;   Bfield = JOVMAGFIELD(R, InTheta, InPhi, Model)
;   Br = Bfield[0]
;   Btheta = Bfield[1]
;   Bphi = Bfield[2]
; 
; INPUTS:
; 
;   R:       Radial coordinate of desired position (in Jupiter Radii)
;   INTHETA: Co-latitude coordinate of desired position (in degrees unless
;             keyword INRAD set)
;   INPHI:   East longitude coordinate of desired position (in degrees unless
;             keyword INRAD set)
;   INMODEL: Name of the internal field model to use.  Only models supported by 
;             getSchmidtCoefs.pro allowed.  At first writing, only 'VIP4', 'VIT4', 
;             'VIPAL', and 'Dipole'. Now supports the 'Khurana' and
;             'JM09' (New Juno Model) models. (Default: 'VIPAL')
; 
; OPTIONAL INPUTS:
;   
;   TIME:    The J2000 time in seconds.  Necessary with the 'Khurana' model.  
;              (Default: Use ctimeOfDPLT to make dipole point towards noon.)
; 
; KEYWORD PARAMETERS:
; 
;   INRAD:         If set, InTheta and InPhi are given in radians rather than degrees.
;   CURRENT_SHEET: Apply a current sheet to the selected model.  Currently only supports
;                    'CANSheet', which is the Connerney, Acuna, and Ness Current Sheet model.
;                    Do not apply a current sheet to 'Khurana' model because it already
;                    applies a current sheet.
;
; OUTPUTS:
; 
;   A three index array with the radial, theta, and phi components of the magnetic
;   field in nanotesla at the requested position and using the given model.  If above
;   one of the poles, the phi component of the magnetic field will be 0.
;
; PROCEDURE:
; 
;   Solves for the magnetic field using spherical harmonics with Schmidt quasi-
;   normalized Legendre functions.  
;   
;   The following functions are required to run:
;     getSchmidtCoefs
;     getLegendre
;     getLegendre2
;     dPdtheta
;     dPdtheta2
;     kk_2009
;     FBfield -> CAN_SHEET
;     ctimeOfDPLT
;
; EXAMPLE:
; 
;   To get the magnetic field at 2 Jupiter Radii, at the equator and 0 degrees
;   longitude, according to the 'VIP4' model, use:
;     Bfield = jovMagField(2, 90, 0, 'VIP4')
;     Br = Bfield[0]
;     Btheta = Bfield[1]
;     Bphi = Bfield[2]
;     
;   To add the Connerny current sheet to the 'VIP4' model use:
;     Bfield = jovMagField(2, 90, 0, 'VIP4', CURRENT_SHEET='CANSheet')
;     
;   To get the magnitude of the magnetic field at 5 Jupiter radii, -30 degrees latittude,
;   and 45 degrees west longitude, according to the 'VIPAL' model, use:
;     Bfield = jovMagField(5, 120, 315, 'VIPAL')
;     Bmag = SQRT(Bfield[0]^2 + Bfield[1]^2 + Bfield[2]^2)
;     
;
; MODIFICATION HISTORY:
;   Written by: Drake A. Ranquist, 7 Jun 2013
;     8 Jul 2013, Added compatibility with the Khurana Field Model
;     29 Sep 2016, Added Connerney, Acuna, and Ness Current Sheet Model
;     16 Feb 2018, Peter Kollmann: added getLegendre2 and dPdtheta2 that can handle models to infinite order
;-

function jovMagField, r, inTheta, inPhi, inModel, time, INRAD=inRad, CURRENT_SHEET=currentSheet

;----------------Set Defaults-------------------
IF NOT KEYWORD_SET(model) THEN model='VIPAL'
IF NOT KEYWORD_SET(inRad) THEN inRad=0
IF NOT KEYWORD_SET(currentSheet) THEN currentSheet=''

;Check that multiple current sheets aren't being applied
IF currentSheet EQ 'CANSheet' AND model EQ 'Khurana' THEN BEGIN
   MESSAGE, 'Khurana model already has a current sheet applied!'
ENDIF   

;Make it so that 'CANSheet' or 'JRM09+CAN' can be given as the model
IF inModel EQ 'CANSheet' AND currentSheet EQ '' THEN BEGIN
  model = 'VIP4'
  currentSheet = 'CANSheet'
ENDIF ELSE IF inModel EQ 'JRM09+CAN' AND currentSheet EQ '' THEN BEGIN
   model = 'JRM09'
   currentSheet = 'CANSheet'
ENDIF ELSE BEGIN
  model = inModel
ENDELSE

;-------Change inTheta and inPhi to Radians------
theta=inTheta
phi=inPhi
IF NOT inRad THEN BEGIN
  theta *= !Pi/180
  phi *= !Pi/180
ENDIF

;---------Model Specific Settings--------------
CASE model OF
  'VIPAL':   maxN = 5
  'VIP4':    maxN = 4
  'VIT4':    maxN = 4
  'JRM09':   maxN = 10
  'Khurana': BEGIN
               IF NOT KEYWORD_SET(time) THEN time = ctimeOfDPLT(12)
               kk_2009, time, r, theta, phi, Br, Btheta, Bphi
               return, [Br, Btheta, Bphi]
             END
  'Dipole':  return, [2*COS(theta)/r^3, SIN(theta)/r^3, 0]
  'VIP4Reduced': BEGIN
     maxN = 1
     model = 'VIP4'
  END
  ELSE:      MESSAGE, model + ' model not found.'
ENDCASE


;--Initialize Three B Field Components--
Br = 0
Btheta = 0
Bphi = 0

;Perform Summations for All Components of B
FOR nn=1, maxN DO BEGIN
  radialFactor = r^(-(nn+2.))
  FOR mm=0, nn DO BEGIN
    ;-------------Find Radial B Field Component---------------
    gh = getSchmidtCoefs(model, nn, mm)
    phiFactor = gh[0]*COS(mm*phi) + gh[1]*SIN(mm*phi)
    IF nn gt 5 THEN Pnm = getLegendre2(nn, mm, theta, /INRAD) ELSE $
      Pnm = getLegendre(nn, mm, theta, /INRAD)
    Br += (nn+1) * radialFactor * phiFactor * Pnm
    
    ;--------------Find Theta B Field Component---------------
    IF nn gt 5 THEN dPnmdtheta = dPdtheta2(nn, mm, theta, /INRAD) ELSE $
       dPnmdtheta = dPdtheta(nn, mm, theta, /INRAD)
    Btheta -= radialFactor * phiFactor * dPnmdtheta
    
    ;----------------Find Phi B Field Component-------------------
    Bphi += mm * radialFactor * (gh[0]*SIN(mm*phi) - gh[1]*COS(mm*phi)) * Pnm
  ENDFOR
ENDFOR

;--------Correct Bphi and Don't Divide by 0!------------
sinTheta = SIN(theta)
IF sinTheta ne 0 THEN Bphi /= sinTheta ELSE Bphi=0

;--------Apply Curent Sheet Model-------------
IF currentSheet EQ 'CANSheet' THEN BEGIN
  ;---------Put position in r, latitude, west longitude form
  IF inRad THEN BEGIN
    pos = [r, 90-inTheta*!radeg, 360-inPhi*!radeg]    
  ENDIF ELSE BEGIN
    pos = [r, 90-inTheta, 360-inPhi]
  ENDELSE
  
  ;Get Current Sheet perturbations and convert to nT
  CANB = CAN_SHEET(pos)*1.e5
  
  ;CAN_SHEET puts theta and phi vectors in opposite direction
  Br += CANB[0]
  Btheta -= CANB[1]
  Bphi -= CANB[2]
  
ENDIF

return, [Br, Btheta, Bphi]

END
