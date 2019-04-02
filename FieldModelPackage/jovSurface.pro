;+
; NAME:
; 
;   JOVSURFACE
;
; PURPOSE:
; 
;   Since Jupiter is an oblate spheroid, this function returns the radial 
;   distance from the center of Jupiter to the surface at a given
;   co-latitude.  Returned as a fraction of the equatorial radius (Rj).
;
; CALLING SEQUENCE:
; 
;   R_Surface = JOVSURFACE(Theta)
;
; INPUTS:
; 
;   Theta:  Co-latitude (0 to 180) in degrees (unless INRAD keyword set)
; 
; KEYWORD PARAMETERS:
; 
;   INRAD:  If set, Theta is assumed to be in radians.
;   OBLATE: Inverse of the oblateness of Jupiter.  (Default: 15.41)
;   RADIUS: Equatorial radius of Jupiter.  (Default: 1 Jovian Radius)
;
; OUTPUTS:
; 
;   Distance from the center of Jupiter to the surface at the given co-latitude.
;
; EXAMPLE:
; 
;   To get the ratio of the polar radius to the equatorial radius, use:
;     eqRadius = jovSurface(0)
;     
;   To get the distance to surface in km at -30 degrees latitude, use:
;     distance = jovSurface(120, RADIUS=71492)
;
; MODIFICATION HISTORY:
;   Written by: Drake A. Ranquist, 10 Jun 2013
;-

function jovSurface, theta, INRAD=inRad, OBLATE=oblate, RADIUS=radius

;-------------Set Defaults------------------
IF NOT KEYWORD_SET(inRad) THEN inRad=0
IF NOT KEYWORD_SET(oblate) THEN oblate=15.41
IF NOT KEYWORD_SET(radius) THEN radius=1

;----------Convert theta to Radians-----------
thisTheta = theta
IF NOT inRad THEN thisTheta = theta * !Pi/180

;-------Solve for the Polar Radius-----------
polRadius = radius - radius/DOUBLE(oblate)

;---Solve for X and Z Coordinates of Surface---
xx = radius*SIN(thisTheta)
zz = polRadius*COS(thisTheta)

;----Return Distance to Surface----
return, SQRT(xx^2 + zz^2)

END