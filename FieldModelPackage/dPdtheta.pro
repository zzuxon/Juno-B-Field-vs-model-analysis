;+
; NAME:
; 
;   DPDTHETA
;
; PURPOSE:
; 
;   Solves for the theta derivative of the Schmidt quasi-normalized 
;   Legendre functions related to Associated Legendre functions for n < 6.  
;   Uses a case look up.  These are the derivatives of functions in
;   getLegendre.pro.
;
; CALLING SEQUENCE:
; 
;   dPnmdTheta = dPdTheta(N, M, InTheta)
;
; INPUTS:
; 
;   N:       Degree of the Schmidt quasi-normalized Legendre function (1 to 5)
;   M:       Order of the Schmidt quasi-normalized Legendre function (0 to n)
;   InTheta: The angle (in degrees, unless specified by INRAD keyword)
; 
; KEYWORD PARAMETERS:
; 
;   INRAD: Set if InTheta is given in radians.
;
; OUTPUTS:
; 
;   Solution to the theta derivative of the Schmidt quasi-normalized Legendre 
;   function of degree N, order M, and angle InTheta.  If given an invalid N or M, 
;   returns 'NaN'
;
; PROCEDURE:
; 
;   The Schmidt quasi-normalized Legendre functions are related to 
;   Associated Legendre functions by:
;     Snm = Pnm, if m=0
;     Snm = (-1)^m * (2(n-m)!/(n+m)!)^(1/2) * Pnm, if m>0
;   The derivatives were solved using Mathematica
;
; EXAMPLE:
; 
;   dPnmdTheta = dPdtheta(3, 2, 30)
;   
;   dPnmdTheta = dPdtheta(5, 0, 0.5, /INRAD)
;
; MODIFICATION HISTORY:
;   Written by: Drake A. Ranquist, 5 Jun 2013
;-

;The derivatives, with respect to theta, of all functions in getLegendre.pro
function dPdtheta, n, m, inTheta, INRAD=inRad

IF NOT KEYWORD_SET(inRad) THEN inRad=0

theta = inTheta
IF NOT inRad THEN theta = inTheta*!Pi/180

CASE n OF
  1: BEGIN
      CASE m OF
        0: return, -SIN(theta)
        1: return, COS(theta)
      ENDCASE
    END
  2: BEGIN
      CASE m OF
        0: return, -1.5 * SIN(2*theta)
        1: return, SQRT(3) * COS(2*theta)
        2: return, SQRT(3) * COS(theta) * SIN(theta)
      ENDCASE
    END
  3: BEGIN
      CASE m OF
        0: return, -3/8. * (SIN(theta) + 5*SIN(3*theta))
        1: return, 1/8. * SQRT(1.5) * (COS(theta) + 15*COS(3*theta))
        2: return, SQRT(15)/4. * SIN(theta) * (3*COS(2*theta) + 1)
        3: return, 1.5 * SQRT(2.5) * COS(theta) * SIN(theta)^2
      ENDCASE
    END
  4: BEGIN
      CASE m OF
        0: return, -5/16. * (2*SIN(2*theta) + 7 * SIN(4*theta))
        1: return, SQRT(2.5)/4. *(COS(2*theta) + 7*COS(4*theta))
        2: return, SQRT(5)/4. * SIN(theta) * (5*COS(theta) +7*COS(3*theta))
        3: return, SQRT(17.5)/2. * SIN(theta) * SIN(3*theta)
        4: return, SQRT(35)/2. * COS(theta) * SIN(theta)^3
      ENDCASE
    END
  5: BEGIN
      CASE m OF
        0: return, -15/128. * (2*SIN(theta) + 7*SIN(3*theta) + 21*SIN(5*theta))
        1: return, SQRT(15)/128.*(2*COS(theta) + 21*COS(3*theta) + 105*COS(5*theta))
        2: return, -SQRT(105)/64. * (2*SIN(theta) + 3*SIN(3*theta) - 15*SIN(5*theta))
        3: return, 3/64. * SQRT(17.5) * SIN(theta) * (2*SIN(2*theta) + 15*SIN(4*theta))
        4: return, 3/16. * SQRT(35) * SIN(theta)^3 * (5*COS(2*theta) + 3)
        5: return, 15/8. * SQRT(3.5) * COS(theta) * SIN(theta)^4
      ENDCASE
    END
ENDCASE

return, 'NaN'

END