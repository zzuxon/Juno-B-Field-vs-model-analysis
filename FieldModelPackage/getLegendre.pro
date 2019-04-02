;+
; NAME:
; 
;   GETLEGENDRE
;
; PURPOSE:
; 
;   Solves for the Schmidt quasi-normalized Legendre functions
;   related to Associated Legendre functions for n < 6.  Uses a case
;   look up.
;
; CALLING SEQUENCE:
; 
;   Pnm = GETLEGENDRE(N, M, InTheta)
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
;   Solution to the Schmidt quasi-normalized Legendre function of degree N,
;   order M, and angle InTheta.  If given an invalid N or M, returns 'NaN'
;
; PROCEDURE:
; 
;   The Schmidt quasi-normalized Legendre functions are related to 
;   Associated Legendre functions by:
;     Snm = Pnm, if m=0
;     Snm = (-1)^m * (2(n-m)!/(n+m)!)^(1/2) * Pnm, if m>0
;
; EXAMPLE:
; 
;   Pnm = getLegendre(3, 2, 30)
;   
;   Pnm = getLegendre(5, 0, 0.5, /INRAD)
;
; MODIFICATION HISTORY:
;   Written by: Drake A. Ranquist, 5 Jun 2013
;-

function getLegendre, n, m, inTheta, INRAD=inRad

IF NOT KEYWORD_SET(inRad) THEN inRad=0

theta = inTheta
IF NOT inRad THEN theta = inTheta*!Pi/180

CASE n OF
  1: BEGIN
      CASE m OF
        0: return, COS(theta)
        1: return, SIN(theta)
      ENDCASE
    END
  2: BEGIN
      CASE m OF
        0: return, 1.5 * (COS(theta)^2 - 1/3.)
        1: return, SQRT(3) * COS(theta) * SIN(theta)
        2: return, SQRT(3)/2. * SIN(theta)^2
      ENDCASE
    END
  3: BEGIN
      CASE m OF
        0: return, 2.5 * COS(theta) * (COS(theta)^2 - 9/15.)
        1: return, 2.5 * SQRT(1.5) * SIN(theta) * (COS(theta)^2 - 3/15.)
        2: return, SQRT(15)/2. * COS(theta) * SIN(theta)^2
        3: return, SQRT(2.5)/2. * SIN(theta)^3
      ENDCASE
    END
  4: BEGIN
      CASE m OF
        0: return, 35/8. * (COS(theta)^4 - (6/7.)*COS(theta)^2 + 3/35.)
        1: return, 3.5*SQRT(2.5)*COS(theta)*SIN(theta)*(COS(theta)^2 - 3/7.)
        2: return, 7/4. * SQRT(5) * SIN(theta)^2 * (COS(theta)^2 - 1/7.)
        3: return, SQRT(17.5)/2. * COS(theta) * SIN(theta)^3
        4: return, SQRT(35)/8. * SIN(theta)^4
      ENDCASE
    END
  5: BEGIN
      CASE m OF
        0: return, 1/8. * (63*COS(theta)^5 - 70*COS(theta)^3 + 15*COS(theta))
        1: return, SQRT(15)/8.*SIN(theta)*(21*COS(theta)^4 - 14*COS(theta)^2 + 1)
        2: return, SQRT(105)/16. * SIN(theta)^2 * (5*COS(theta) + 3*COS(3*theta))
        3: return, SQRT(17.5)/16. * SIN(theta)^3 * (9*COS(2*theta) + 7)
        4: return, 3/8. * SQRT(35) * COS(theta) * SIN(theta)^4
        5: return, 3/8. * SQRT(3.5) * SIN(theta)^5
      ENDCASE
    END
ENDCASE

return, 'NaN'

END