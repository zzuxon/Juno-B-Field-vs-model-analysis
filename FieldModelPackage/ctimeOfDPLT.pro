;+
; NAME:
; 
;   CTIMEOFDPLT
;
; PURPOSE:
; 
;   This function gives a default time to use in the Khurana kk_2009
;   Jupiter magnetic field model for the dipole pointing towards
;   different local times.  The variable in kk_2009 is ctime, which is
;   the number of seconds since Jan 1, 1966.
;
; CALLING SEQUENCE:
; 
;   ctime = ctimeOfDPLT(localTime)
; 
; INPUTS:
; 
;   LOCALTIME: The local time you wish the dipole to point.  0 refers
;                to midnight, 6 to dawn, 12 to noon, and 18 to dusk.
;                Fractional inputs allowed.  (0 to 24).
;
; OUTPUTS:
; 
;   A default time (in seconds from Jan 1, 1966) to use in the Khurana
;   kk_2009 model where the dipole is pointing towards the given local time.
;
; PROCEDURE:
; 
;   No other procedures needed.
;
; EXAMPLE:
; 
;   To get a time with the dipole pointing towards noon, use:
;
;      ctime = ctimeOfDPLT(12)
;     
;
; MODIFICATION HISTORY:
;   Written by: Drake A. Ranquist, 26 Jun 2016
;-

function ctimeOfDPLT, localTime

;Length of one Jovian day in seconds used by Khurana model
jovDay = 35734

;A default ctime for when dipole points towards midnight (away from sun)
midTime = 31290  

;Calculate ctime of given local time
frac = localTime/24.
ctime = midTime + frac*jovDay

return, ctime




;Following section was used to try to identify correct local times
;
;;RH System III Dipole Longitude (Phi angle)
;dpPhi = 158
;
;Rotation rate of Jupiter used by Khurana model in degrees/sec
;omega = 870.536d/86400d
;jovDay = 360./omega
;
;nn = ceil(jovDay)
;times = findgen(nn) + jovDay
;thetas = fltarr(nn)
;phis = fltarr(nn)
;phases = fltarr(nn)
;
;FOR ii=0, nn-1 DO BEGIN
;  
;  JSun, times[ii], stheta, sphi, phase
;  thetas[ii] = stheta*!radeg
;  phis[ii] = sphi*!radeg
;  phases[ii] = phase*!radeg
;  
;ENDFOR
;
;void = min(abs(phis-dpPhi),dpLoc)
;void = min(abs(phis-202),midLoc)
;
;
;cgplot, times, phis
;
;STOP
;
;return, ''

end
