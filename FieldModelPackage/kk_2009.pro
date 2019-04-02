;--------------------------------------------------------------------------------------
;---------------Krishan Khurana's Jupiter magnetic field model for IDL-----------------
;--------------------------------------------------------------------------------------
; NAME:
;   kk_2009
;
; PURPOSE:
;   Given right-handed system III spherical coordinates, this program
;    outputs the strength of jupiter's magnetic field in spherical coordinates
;    at the specified coordinates. 
;
; DESCRIPTION:
;   Calculates the magnetic field of a tilted and warped shielded magnetosphere
;    using general deformation technique.
;
; CALLING SEQUENCE:
;   kk_2009, ctime, R, theta, phi, Brm, Btm, Bpm, mpCheck
;
; INPUTS:
;   ctime         - J2000 time in seconds (can be calculated by the
;                   ctimer subroutine located in this code)
;   R, theta, phi - right-handed system III coordinates (angles in radians)
;   
; KEYWORDS:
;   QUIET     - If set, then nothing will be printed to the screen.
;
; OUTPUTS:
;   Brm, Btm, Bpm    - r, theta, phi components of the magnetic field
;                      in nanoteslas [nT]
;   mpCheck          - equals 1 if  inside the magnetopause
;                    - equals 0 if outside the magnetopause
;  
;   note: as a procedure the output to kk_2009 can easily be changed
;         to include more than just the components to the magnetic field
;         - for example, the following can also be included:
;
;   XJSO, YJSO, ZJSO - x, y, z components of the magnetic field in the 
;
;                      juptier-sun-orbital coordinate system (see JROT subroutine),
;   ZNS3             - distance to the current sheet from the
;                      jovigraphic equator in system III coordinates
;
;
; MODIFICATION HISTORY:
;   2003, written by Krishan Khurana in Fortran 77/90
;   2008, ctimer function added by Mariel Desroche
;         (converts year, month, day, hour and second to J2000 seconds)
;   2009, converted to IDL by Adam Shinn
;   2009, leapSecond function added by Adam Shinn and Rob Wilson
;         (replaces eetime function to easily append for future leap seconds)
;   2013, modified ctimer to actually give J2000 seconds instead of 1966 seconds,
;         JSun to receive the correct ctime input, updated leapSecond and moved
;         above ctimer, and changed kk_2009, ctimer, leapSecond, and JRot's use
;         of system variables.  Added QUIET keyword.  Changes by Drake A. Ranquist.
;   2014, Fixed bug in JSUN for times when stheta was larger than variable three by
;         Drake A. Ranquist.
;
; WEBSITE:
;   This code is available on the website for the Magnetospheres of
;    Outer Planets group at the Univeristy of Colorado at Boulder,
;    associated with the Laboratory for Atmospheric and Space Physics:
;    http://lasp.colorado.edu/MOP/resources/
;--------------------------------------------------------------------------------------
pro kk_2009, ctime, R, theta, phi, Brm, Btm, Bpm, mpCheck, XJSO, YJSO, ZJSO, ZNS3, $
             QUIET=quiet

  quiet = KEYWORD_SET(quiet)

;STOP

  vecin = dblarr(3)
  vecou = dblarr(3)

  ; constants
  parmod = dblarr(10)
  M = 8
  M1 = 64
  M2 = 16
  Mode = 7
  pi = 3.14159265358979d
  radian = pi/180d

  rho = R*sin(theta)  
  xs3 = rho*cos(phi)
  ys3 = rho*sin(phi)
  zs3 = R*cos(theta)

  JSun, ctime, stheta, sphi, phase ; given ctime, then stheta, sphi, and phase are returned
  defsysv, "!ctime", ctime
  defsysv, "!stheta", stheta
  defsysv, "!sphi", sphi
  defsysv, "!phase", phase

  ; Now define the nine derivatives
  ; These are calculated at the original location
  ; First get the unit normals of the zp axis
  csheet_N, RNx, RNy, RNz, rho, phi, zs3, stheta, ctime


  ; The z axis of the mapped location is given by
  zpx = RNx
  zpy = RNy
  zpz = RNz
  zp = xs3*zpx + ys3*zpy + zs3*zpz

  ; The y axis of the mapped location is given by
  ypx = zpy
  ypy = -zpx
  ypz = 0d
  yp1 = sqrt(ypx*ypx + ypy*ypy + ypz*ypz)
  ypx = ypx/yp1
  ypy = ypy/yp1
  ypz = ypz/yp1
  yp = xs3*ypx + ys3*ypy + zs3*ypz

  ; The x axis of the mapped location is given by
  xpx = ypy*zpz - ypz*zpy
  xpy = ypz*zpx - ypx*zpz
  xpz = ypx*zpy - ypy*zpx
  xp = xs3*xpx + ys3*xpy + zs3*xpz
  
  ; Now define the nine derivatives
  dxpdx = xpx
  dypdx = ypx
  dzpdx = zpx
  dxpdy = xpy
  dypdy = ypy
  dzpdy = zpy
  dxpdz = xpz
  dypdz = ypz
  dzpdz = zpz

  R = sqrt(xs3*xs3 + ys3*ys3 + zs3*zs3)
  Rp = sqrt(xp*xp + yp*yp + zp*zp)

  ; Now calculate the T matrix
  Txx = (dypdy*dzpdz - dypdz*dzpdy)
  Txy = (dxpdz*dzpdy - dxpdy*dzpdz)
  Txz = (dxpdy*dypdz - dxpdz*dypdy)
  Tyx = (dypdz*dzpdx - dypdx*dzpdz)
  Tyy = (dxpdx*dzpdz - dxpdz*dzpdx)
  Tyz = (dxpdz*dypdx - dxpdx*dypdz)
  Tzx = (dypdx*dzpdy - dypdy*dzpdx)
  Tzy = (dxpdy*dzpdx - dxpdx*dzpdy)
  Tzz = (dxpdx*dypdy - dxpdy*dypdx)

  ; The mapped location is:
  vecin = [xs3, ys3, zs3]
  JROT, 's3c', 'jso', vecin, vecou, ctime
  XJSO = vecou[0]
  YJSO = vecou[1]
  ZJSO = vecou[2]

  csheet_struc, ZNS3, rho, phi, XJSO, YJSO, stheta, ctime

  RLT = atan(YJSO, XJSO)
  rmap = sqrt(xp*xp + yp*yp)
  pmap = atan(yp, xp)
  zmap = zs3 - ZNS3


  ; Now calculate the field at the mapped location in cylindrical coordinates
  scol = pi/2d - stheta

  get_mapped_sunangle, scol, sphi, scolOut, sphiOut, xpx, xpy, xpz, ypx, ypy, ypz, zpx, zpy, zpz

  dipole_shield_cyl_S3, PARMOD, ctime, rmap, pmap, zmap, Brds, Bpds, Bzds, sphiOut

  ;problem not here
  
  tail_mag_notilt_all_modes, rmap, pmap, zmap, RLT, Brcs, Bpcs, Bzcs, Mode

  ;problem is here

  tail_mag_shield_cyl_S3, ctime, rmap, pmap, zmap, Brcss, Bpcss, Bzcss, M, sphiOut, Mode
  
  

  Brmap = Brds + Brcs + Brcss 
  Bpmap = Bpds + Bpcs + Bpcss 
  Bzmap = Bzds + Bzcs + Bzcss 
  
  cyl2car_vect, Brmap, Bpmap, Bzmap, BXcarMap, BYcarMap, BZcarMap, pmap

  cyl2car_vect, Brds, Bpds, Bzds, BXdsMap, BYdsMap, BZdsMap, pmap

  cyl2car_vect, Brcs, Bpcs, Bzcs, BXcsMap, BYcsMap, BZcsMap, pmap

  cyl2car_vect, Brcss, Bpcss, Bzcss, BXcssMap, BYcssMap, BZcssMap, pmap
  
  Bxfinal =  Txx*BXcarMap + Txy*BYcarMap + Txz*BZcarMap 
  Byfinal =  Tyx*BXcarMap + Tyy*BYcarMap + Tyz*BZcarMap 
  Bzfinal =  Tzx*BXcarMap + Tzy*BYcarMap + Tzz*BZcarMap

  Bxdsfinal =  Txx*BXdsMap + Txy*BYdsMap + Txz*BZdsMap  
  Bydsfinal =  Tyx*BXdsMap + Tyy*BYdsMap + Tyz*BZdsMap  
  Bzdsfinal =  Tzx*BXdsMap + Tzy*BYdsMap + Tzz*BZdsMap

  Bxcsfinal =  Txx*BXcsMap + Txy*BYcsMap + Txz*BZcsMap  
  Bycsfinal =  Tyx*BXcsMap + Tyy*BYcsMap + Tyz*BZcsMap  
  Bzcsfinal =  Tzx*BXcsMap + Tzy*BYcsMap + Tzz*BZcsMap
  
  Bxcssfinal =  Txx*BXcssMap + Txy*BYcssMap + Txz*BZcssMap  
  Bycssfinal =  Tyx*BXcssMap + Tyy*BYcssMap + Tyz*BZcssMap  
  Bzcssfinal =  Tzx*BXcssMap + Tzy*BYcssMap + Tzz*BZcssMap

  CAR2SPH_MAG, Bxfinal, Byfinal, Bzfinal, Brm, Btm, Bpm, theta, phi

  CAR2SPH_MAG, Bxdsfinal, Bydsfinal, Bzdsfinal, Brdsm, Btdsm, Bpdsm, theta, phi

  CAR2SPH_MAG, Bxcsfinal, Bycsfinal, Bzcsfinal, Brcsm, Btcsm, Bpcsm, theta, phi

  CAR2SPH_MAG, Bxcssfinal, Bycssfinal, Bzcssfinal, Brcssm, Btcssm, Bpcssm, theta, phi 

  JOVIAN_VIP4_no_dipole, 4, R, theta, phi, BR, BT, BF

  Brm = Brm + BR
  Btm = Btm + BT
  Bpm = Bpm + BF

  Brdsm = Brdsm + BR
  Btdsm = Btdsm + BT
  Bpdsm = Bpdsm + BF

  ; now calculate the solar wind magnetic field
  getBimf, ctime, theta, phi, BrS3IMF, BtS3IMF, BpS3IMF

  CheckIfInsideMappedMp, ctime, XS3, YS3, ZS3, ZNS3, mpCheck
  ;checkIfInsideMagnetop, ctime, XS3, YS3, ZS3, mpCheck


  ;print, BrS3IMF, BtS3IMF, BpS3IMF
  if (mpCheck eq 1) then begin
     Brm = Brm + BrS3IMF
     Btm = Btm + BtS3IMF
     Bpm = Bpm + BpS3IMF
  endif else begin
     Brm = BrS3IMF
     Btm = BtS3IMF
     Bpm = BpS3IMF
  endelse
end   ; kkk_2008 main program

;--------------------------------------------------------------------------------------
; function doy
;
; day of year calculator
;--------------------------------------------------------------------------------------
function doy, iyr, imon, iday
  mon = [31,28,31,30,31,30,31,31,30,31,30,31]
  dy = 0L

  for i = 2, imon do begin
     dy = dy + mon[i - 2]
     ; Add an extra day for February
     if (i eq 3) then begin
        if ((iyr mod 100) ne 0) || ((iyr mod 400) eq 0) && ((iyr mod 4) eq 0) then dy = dy + 1
     endif
  endfor

  dy = dy + iday
  return, dy
end


;--------------------------------------------------------------------------------------
; function leapSecond(ctime)
;   allows ctime to account for leap seconds
;--------------------------------------------------------------------------------------
function leapSecond, ctime, year, month, QUIET=quiet

  quiet = KEYWORD_SET(quiet)

 time  = year*100L + month ; formats year and month as YYYYMM
 if (time lt 197207) then Tcor = 10d else $ ; leap second jun 30 1972
 if (time lt 197301) then Tcor = 11d else $ ; leap second dec 31 1972
 if (time lt 197401) then Tcor = 12d else $ ; leap second dec 31 1973
 if (time lt 197501) then Tcor = 13d else $ ; leap second dec 31 1974
 if (time lt 197601) then Tcor = 14d else $ ; leap second dec 31 1975
 if (time lt 197701) then Tcor = 15d else $ ; leap second dec 31 1976
 if (time lt 197801) then Tcor = 16d else $ ; leap second dec 31 1977
 if (time lt 197901) then Tcor = 17d else $ ; leap second dec 31 1978
 if (time lt 198001) then Tcor = 18d else $ ; leap second dec 31 1979
 if (time lt 198107) then Tcor = 19d else $ ; leap second jun 30 1981
 if (time lt 198207) then Tcor = 20d else $ ; leap second jun 30 1982
 if (time lt 198307) then Tcor = 21d else $ ; leap second jun 30 1983
 if (time lt 198507) then Tcor = 22d else $ ; leap second jun 30 1985
 if (time lt 198801) then Tcor = 23d else $ ; leap second dec 31 1987
 if (time lt 199001) then Tcor = 24d else $ ; leap second dec 31 1989
 if (time lt 199101) then Tcor = 25d else $ ; leap second dec 31 1990
 if (time lt 199207) then Tcor = 26d else $ ; leap second jun 30 1992
 if (time lt 199307) then Tcor = 27d else $ ; leap second jun 30 1993
 if (time lt 199407) then Tcor = 28d else $ ; leap second jun 30 1994
 if (time lt 199601) then Tcor = 29d else $ ; leap second dec 31 1995
 if (time lt 199707) then Tcor = 30d else $ ; leap second jun 30 1997
 if (time lt 199901) then Tcor = 31d else $ ; leap second dec 31 1998
 if (time lt 200601) then Tcor = 32d else $ ; leap second dec 31 2005  
 if (time lt 200901) then Tcor = 33d else $ ; leap second dec 31 2008
 if (time lt 201207) then Tcor = 34d else $ ; leap second jun 30 2012 
 if (time lt 201507) then Tcor = 35d else $ ; leap second jun 30 2015 
  begin
     Tcor = 35d 
     IF NOT quiet THEN $
       print, "Warning: For any time after 31 Dec 2015, leap second function needs updating!"
     ;STOP
  endelse
  thirty_four_years = 0.1072958367816d10
  leaptime = ctime + Tcor - thirty_four_years ; sets time to J2000 epoch
  return, leaptime
end 


;--------------------------------------------------------------------------------------
; function ctimer
; 
; The following lines (and others between comment bars below) were added
; by M. Desroche for use with data sets that use DOY instead of mon/day format
; The similiar lines above need to be commented out and these lines commented in
; for use with this data format
;**************************************************************
; function ctimer, iyr, iday, ihr, imin, sec
;**************************************************************
;--------------------------------------------------------------------------------------
function ctimer, iyr, imon, iday, ihr, imin, sec, QUIET=quiet
  ;dayofyear = julday(imon, iday, iyr, 0, 0, 0) - julday(12, 31, iyr-1, 0, 0, 0)
  ctime = 0d
  ndays = 0L
  ; First calculate the number of days from Jan 1, 1966
  if (iyr gt 1966) then begin
     for i = 1966, iyr - 1 do begin
        ndays = ndays + 365
        if ((i mod 100) ne 0) || ((i mod 400) eq 0) && ((i mod 4) eq 0) then ndays = ndays + 1
     endfor
     ; Now add the number of days of the current year
     ndays = ndays + doy(iyr, imon, iday) - 1 
     ;**************************************************************
     ; ndays = ndays + iday - 1
     ;**************************************************************

     ctime = ndays*86400d + ihr*3600d + imin*60d + sec
  endif else begin
     ; Calculate the seconds for the negative years
     for i = iyr, 1965 do begin
        ndays = ndays - 365
        if ((i mod 100) ne 0) || ((i mod 400) eq 0) && ((i mod 4) eq 0) then ndays = ndays - 1
     endfor
     ; Now subtract the number of days of the current year
     ndays = ndays + doy(iyr, imon, iday)
     ;***************************************************************
     ; ndays = ndays + iday
     ;***************************************************************

     ctime = (ndays - 1)*86400d + ihr*3600d + imin*60d + sec
  endelse

  ;jultime = julday(imon, iday, iyr, ihr, imin, sec)*24d*3600d
  ;print, "ctime", ctime
  ;print, "jultime", jultime

  ;defsysv, "!ctime", ctime
  ;defsysv, "!year", iyr
  ;defsysv, "!month", imon
  return, leapSecond(ctime, iyr, imon, QUIET=quiet)
end


;--------------------------------------------------------------------------------------
; procedure JSun
;
; INPUT:  ctime of the data point 
; OUTPUTS: stheta, sphi, latitude and longitude  (in radians) of the Sun in system III (RH).  
; OUTPUTS: phase, Orbital phase angle of Jupiter (in radians) 
; The equations are written in eetime the J2000 time convention followed by PDS. 
; We first convert ctime to eetime
;--------------------------------------------------------------------------------------
pro JSun, ctime, stheta, sphi, phase

  ; Last updated August 12, 2002
  ; The program first calculates the direction of the sun in non-rotating  
  ; coordinates from equations of the type:
  ;   theta=a1*cos(omegay*t)+a2*sin(omegay*t)+a3*cos(2.*omegay*t)+a4*sin(2.*omegay*t)+a5
  ;   fphi=b1*cos(omegay*t)+b2*sin(omegay*t)+b3*cos(2.*omegay*t)+b4*sin(2.*omegay*t)+b5
  ; Then we rotate into the System III coordinates

  ; ctime is in double precision. theta and phi are single precision variables  
  ; omega is jupiter's rotation rate
  ; omegay is jupiter's yearly orbital rate

  ; Initialize variables
  pi = 3.14159265358979d
  twopi = 2d*PI
  radian = PI/180d
  degree = 180d/PI
  yrjup = 11.85652502d*86400d*365.25d
  omega = 870.536d/86400d
  omegay = 2*PI/yrjup
  year = 86400d*0.36525d3
  D360 = 360d
  etime1 = -8.25767955817479d8
  three = 3.123d*radian
  tan3 = 0.054560676d
  sin3 = 0.054479647d
  cos3 = 0.99851488d

  aa = [0.14347029d, 3.1145815d, -0.12025561d, 0.093909436d, -0.39321884d-5, 0.10194945d-3, -0.12799464d]
  bb = [-4.5467523d, 3.1848875d, -0.16329986d, -0.09776818d, 0.17556527d-3, -0.01978317d, 44.55915d]

  ; First calculate the latitude and longitude in non-rotating Jupiter coordinates
  ; Calculate the best fit theta and fphi
  ;t = eetime(ctime) - etime1
  ;t = leapSecond(ctime) - etime1 ;modified ctimer to use leapSecond
  t = ctime - etime1
  
  cos_omeg = cos(omegay*t)
  sin_omeg = sin(omegay*t)
  
  ;x = dblarr(7)
  x5 = t/year
  x = [cos_omeg, sin_omeg, cos(2*omegay*t), $
          2*cos_omeg*sin_omeg, x5*x5, x5, 1d]

  ; fphi is phi in Jupiter fixed (non-rotating) coordinate
  fphi = total(bb*x)
  stheta = total(aa*x)
  
  ; Now rotate the longitude to Jupiter System III
  ; First Add the rotation of Jupiter around the Sun
  ; fphi is the phi of the Sun as unspinning Jupiter goes around the sun
  fphi = (fphi + t/yrjup*360d) mod D360

  ; Next add the rotation of Jupiter around its axis.
  sphi = (fphi - t*omega) mod D360
  if (sphi lt 0d) then sphi = sphi + 360d
  sphi = sphi*radian
  stheta = stheta*radian
  acostanstheta = acos(tan(stheta)/tan3)

  ; Now compute the orbital phase (called phi2 or phase here)
  ; There are two solutions to the problem. Only one that is close to phi2b is correct
  if (stheta ge three) then begin
    stheta = three
    acostanstheta = 0d
  endif
  if (-stheta ge three) then begin 
    stheta = -three
    acostanstheta = PI
  endif
  

  phi21 = (sphi + PI) + acostanstheta
  phi22 = (sphi + PI) - acostanstheta
  phi21 = phi21 mod twopi
  phi22 = phi22 mod twopi
  dphi = fphi + 48.23012d
  phi2b = sphi - dphi*radian
  phi2b = phi2b mod twopi
  phase = phi21
  dif2 = abs(phi22 - phi2b)
  if (dif2 gt 350d*radian) then dif2 = twopi - dif2
  dif1 = abs(phi21 - phi2b)
  if (dif1 gt 350d*radian) then dif1 = twopi - dif1
  if (dif2 lt dif1) then phase = phi22
  ;phase = phi2b
  phase = phase mod twopi
end

;--------------------------------------------------------------------------------------
; procedure JROT
;
; INPUTS: From:   a chracter*3 variable denoting the incoming coordinate system
;   To:     a chracter*3 variable denoting the outgoing coordinate system
;   Vecin   a variable of dimension 3 containing the incoming vector
;   Vecout: a variable of dimension 3 containing the outgoing vector
;   ctime:  a double precision variable denoting Cline time which can be 
;           calculated by using the function program ctimer.
;
; The "From" coordinate system is rotated to system III, then rotated
;   into the "To" coordinate system
;
; The supported coordinate systems are: S3C System III Cartesian (right-handed)
;                                       JSO Jupiter-Sun-Orbital
;                     JSM Jupite-Sun-Magnetic
;                                       DIP Dipole (cartesian)
;                                       JSS Jupiter-Sun-Spin rewritten 7/2/09 by M.D.
;--------------------------------------------------------------------------------------
pro JROT, From, To, Vecin, Vecout, ctime
  vector = dblarr(3)
  dipole = dblarr(3,3)
  first = (second = dblarr(3,3) )
  dummy = dblarr(3,3)
  ;Identity = dblarr(3,3)
  Identity = [ [1d, 0d, 0d], $
               [0d, 1d, 0d], $
               [0d, 0d, 1d]  ]
  PI = 3.14159265358979d
  twopi = 2d*PI
  radian = PI/180d
  degree = 180d/PI
  three = 3.123*radian
  tan3 = 0.054560676d
  sin3 = 0.054479647d
  cos3 = 0.99851488d
  TH_DIP = 0.16755161d
  PH_DIP = -3.5255651d
  dipole[2,*] = [-0.154625290d, 0.0624726744d,  0.98599604d]
  dipole[1,*] = [  -0.3746066d,  -0.92718385d,           0d]
  dipole[0,*] = [  -0.9141996d,   0.36936062d, -0.16676875d]
  ;dipole = transpose(dipole)
  
  defsysv, "!stheta", exists = exist
  if exist then begin
    if ctime eq !ctime then begin
      stheta = !stheta
      sphi = !sphi
      phase = !phase
    endif else begin
      JSun, ctime, stheta, sphi, phase
    endelse
  endif else begin
     JSun, ctime, stheta, sphi, phase
  endelse

  sin_stheta = sin(stheta)
  cos_stheta = cos(stheta)
  sin_sphi = sin(sphi)
  cos_sphi = cos(sphi)

  case from of 
     "s3c" : first = Identity
     "jso" : begin
                ; print, " Incoming coordinate system is JSO"
                ; Calculate the rotation matrix to go from JSO to S3R
                ; Define X component in system III of XJSO unit vector etc.
                dummy[0,*] = [cos_stheta*cos_sphi, cos_stheta*sin_sphi, sin_stheta]
                ; Calculate the Z axis of the JSO coordinate system from the fact that the Z axis is tilted
                ;   by 3.12 degrees from the Z axis of SIII and is normal to the X axis of the JSO system.
                ; Define X component in system III of ZJSO unit vector etc.
                dummy[2,*] = [sin3*cos(phase), sin3*sin(phase), cos3]
                ; Define X component in system III of YJSO unit vector etc.
                dummy[1,0] = dummy[2,1]*dummy[0,2] - dummy[2,2]*dummy[0,1]
                dummy[1,1] = dummy[2,2]*dummy[0,0] - dummy[2,0]*dummy[0,2]
                dummy[1,2] = dummy[2,0]*dummy[0,1] - dummy[2,1]*dummy[0,0]
                first = transpose(dummy)
             end ; jso case
     "jsm" : begin
                ; print, " Incoming coordinate system is JSM"
                ; Now define the JSM transpose vector
                ; Define X component in system III of XJSM unit vector etc.
                dummy[0,0] = cos_stheta*cos_sphi
                dummy[0,1] = cos_stheta*sin_sphi
                dummy[0,2] = sin_stheta
                ; Now define the Y vector so that it is perpendicular to the dipole vector and X
                dummy[1,0] = dummy[0,2]*dipole[2,1] - dummy[0,1]*dipole[2,2]
                dummy[1,1] = dummy[0,0]*dipole[2,2] - dummy[0,2]*dipole[2,0]
                dummy[1,2] = dummy[0,1]*dipole[2,0] - dummy[0,0]*dipole[2,1]
                denom = sqrt(dummy[1,0]*dummy[1,0] + dummy[1,1]*dummy[1,1] + dummy[1,2]*dummy[1,2])
                dummy[1,0] = dummy[1,0]/denom
                dummy[1,1] = dummy[1,1]/denom
                dummy[1,2] = dummy[1,2]/denom
                ; Now define the z vector
                dummy[2,0] = dummy[0,1]*dummy[1,2] - dummy[0,2]*dummy[1,1]
                dummy[2,1] = dummy[0,2]*dummy[1,0] - dummy[0,0]*dummy[1,2]
                dummy[2,2] = dummy[0,0]*dummy[1,1] - dummy[0,1]*dummy[1,0]
                first = transpose(dummy)
             end ; jsm case
     "dip" : first = transpose(dipole)
     "jss" : begin
                ; print, " Incoming coordinate system is JSS"
                ; Calculate the rotation matrix to go from JSS to S3R
                ; Define X component in system III of ZJSS unit vector etc.
                dummy[2,0] = 0d
                dummy[2,1] = 0d
                dummy[2,2] = 1d
                ; Define X component in system III of YJSS unit vector etc.
                ;   so that it is perpendicular to the system 3 Z spin axis and the 
                ;   Jupiter-Sun line
                dummy[1,0] = -cos_stheta*sin_sphi
                dummy[1,1] = cos_stheta*cos_sphi
                dummy[1,2] = 0d
                denom = sqrt(dummy[1,0]*dummy[1,0] + dummy[1,1]*dummy[1,1])
                dummy[1,0] = dummy[1,0]/denom
                dummy[1,1] = dummy[1,1]/denom
                dummy[1,2] = dummy[1,2]/denom
                ; Define X component in system III of XJSS unit vector etc.
                ;   so that it is perpendicular to ZJSS and YJSS
                dummy[0,0] = dummy[1,1]
                dummy[0,1] = -dummy[1,0]
                dummy[0,2] = 0d
                ; The transpose of the dummy matrix takes JSS to System 3.
                first = transpose(dummy)
             end ; jss case
  endcase
  
  case to of
     "s3c" : second = identity
     "jso" : begin
                ; print, " The outgoing coordinate system is JSO"
                ; Get the matrix that Rotates from System III cartesian into JSO. 
                ; Define X component in system III of XJSO unit vector etc.
                Second[0,0] = cos_stheta*cos_sphi
                Second[0,1] = cos_stheta*sin_sphi
                Second[0,2] = sin_stheta
                ; Calculate the Z axis of the JSO coordinate system from the fact that the Z axis is tilted
                ;   by 3.12 degrees from the Z axis of SIII and is normal to the X axis of the JSO system.
                ; Define X component in system III of ZJSO unit vector etc.
                Second[2,0] = sin3*cos(phase)
                Second[2,1] = sin3*sin(phase)
                Second[2,2] = cos3
                ; Define X component in system III of YJSO unit vector etc.
                Second[1,0] = Second[2,1]*Second[0,2] - Second[2,2]*Second[0,1]
                Second[1,1] = Second[2,2]*Second[0,0] - Second[2,0]*Second[0,2]
                Second[1,2] = Second[2,0]*Second[0,1] - Second[2,1]*Second[0,0]
             end ; jso case
     "jsm" : begin
                ; print, " The outgoing coordinate system is JSM"
                ; Get the matrix that Rotates from System III cartesian into JSM. 
                ; Define X component in system III of XJSM unit vector etc.
                Second[0,0] = cos_stheta*cos_sphi
                Second[0,1] = cos_stheta*sin_sphi
                Second[0,2] = sin_stheta
                ; Now define the Y vector so that it is perpendicular to  the dipole vector and X 
                Second[1,0] = second[0,2]*dipole[2,1] - second[0,1]*dipole[2,2]
                Second[1,1] = second[0,0]*dipole[2,2] - second[0,2]*dipole[2,0]
                Second[1,2] = second[0,1]*dipole[2,0] - second[0,0]*dipole[2,1]
                denom = sqrt(second[1,0]*second[1,0] + second[1,1]*second[1,1] + second[1,2]*second[1,2])
                Second[1,0] = second[1,0]/denom
                Second[1,1] = second[1,1]/denom
                Second[1,2] = second[1,2]/denom
                ; Now define the z vector
                Second[2,0] = second[0,1]*second[1,2] - second[0,2]*second[1,1]
                Second[2,1] = second[0,2]*second[1,0] - second[0,0]*second[1,2]
                Second[2,2] = second[0,0]*second[1,1] - second[0,1]*second[1,0]
             end ; jsm case
     "dip" : second = transpose(dipole)
     "jss" : begin
                ; print, " The outgoing coordinate system is JSS"
                ; Define X component in system III of ZJSS unit vector etc.
                second[2,0] = 0d
                second[2,1] = 0d
                second[2,2] = 1d
                ; Define X component in system III of YJSS unit vector etc.
                ;   so that it is perpendicular to the system 3 Z spin axis and the 
                ; Jupiter-Sun line
                second[1,0] = -cos_stheta*sin_sphi
                second[1,1] = cos_stheta*cos_sphi
                second[1,2] = 0d
                denom = sqrt(second[1,0]*second[1,0] + second[1,1]*second[1,1] + second[1,2]*second[1,2])
                second[1,0] = second[1,0]/denom
                second[1,1] = second[1,1]/denom
                second[1,2] = second[1,2]/denom
                ; Define X component in system III of XJSS unit vector etc.
                ;   so that it is perpendicular to ZJSS and YJSS
                second[0,0] = second[1,1]*second[2,2] - second[1,2]*second[2,1]
                second[0,1] = second[1,2]*second[2,0] - second[1,0]*second[2,2]
                second[0,2] = second[1,0]*second[2,1] - second[1,1]*second[2,0]
             end ; jss case
  endcase

  ; Now multimply vecin with first and second matrices to get the vecout
  vector = matrix_multiply(first, vecin)
  vecout = matrix_multiply(second, vector)
end   ; jrot procedure

;--------------------------------------------------------------------------------------
; procedure cy12car_pos
;--------------------------------------------------------------------------------------
pro cyl2car_pos, rho, phi, Z, Xcar, Ycar, Zcar
  Xcar = rho*cos(phi)
  Ycar = rho*sin(phi)
  Zcar = Z
end

;--------------------------------------------------------------------------------------
; procedure car2cyl_pos
;--------------------------------------------------------------------------------------
pro car2cyl_pos, rho, phi, Z, Xcar, Ycar, Zcar
  phi = atan(Ycar,Xcar)
  
  rho = sqrt(Xcar*Xcar + Ycar*Ycar)
  Z = Zcar
end

;--------------------------------------------------------------------------------------
; procedure car2cyl_vect
;--------------------------------------------------------------------------------------
pro car2cyl_vect, Brho, Bphi, BZ, BXcar, BYcar, BZcar, phi
  cos_phi = cos(phi)
  sin_phi = sin(phi)
  Brho = BXcar*cos_phi + BYcar*sin_phi
  Bphi = -BXcar*sin_phi + BYcar*cos_phi
  BZ = BZcar
end

;--------------------------------------------------------------------------------------
; procedure cy12car_vect
;--------------------------------------------------------------------------------------
pro cyl2car_vect, Brho, Bphi, BZ, BXcar, BYcar, BZcar, phi
  cos_phi = cos(phi)
  sin_phi = sin(phi)
  BXcar = Brho*cos_phi - Bphi*sin_phi
  BYcar = Brho*sin_phi + Bphi*cos_phi
  BZcar = BZ
end

;--------------------------------------------------------------------------------------
; procedure CAR2SPH_pos
;--------------------------------------------------------------------------------------
pro CAR2SPH_pos, X, Y, Z, R, TH, PHI
  ; BARTSCH
  R = sqrt(X*X + Y*Y + Z*Z)
  PHI = atan(Y, X)
  TH = asin(Z/R)
  
end

;--------------------------------------------------------------------------------------
; procedure SPH2CAR_pos
;--------------------------------------------------------------------------------------
pro SPH2CAR_pos, X, Y, Z, R, TH, PHI
  ; BARTSCH
  sin_th = sin(TH)
  X = R*sin_th*cos(PHI)
  Y = R*sin_th*sin(PHI)
  Z = R*cos(TH)
end

;--------------------------------------------------------------------------------------
; procedure SPH2CAR_MAG
;--------------------------------------------------------------------------------------
pro SPH2CAR_MAG, BX, BY, BZ, BR, BTH, BPHI, TH, PHI
  ; ARFKEN
  sin_th = sin(TH)
  cos_th = cos(TH)
  sin_phi = sin(PHI)
  cos_phi = cos(PHI)
  BX = BR*sin_th*cos_phi + BTH*cos_th*cos_phi - BPHI*sin_phi
  BY = BR*sin_th*sin_phi + BTH*cos_th*sin_phi + BPHI*cos_phi
  BZ = BR*cos_th - BTH*sin_th
end

;--------------------------------------------------------------------------------------
; procedure CAR2SPH_MAG
;--------------------------------------------------------------------------------------
pro CAR2SPH_MAG, BX, BY, BZ, BR, BTH, BPHI, TH, PHI
  ; ARFKEN
  sin_th = sin(TH)
  cos_th = cos(TH)
  sin_phi = sin(PHI)
  cos_phi = cos(PHI)
  BR = BX*sin_th*cos_phi + BY*sin_th*sin_phi + BZ*cos_th
  BTH = BX*cos_th*cos_phi + BY*cos_th*sin_phi - BZ*sin_th
  BPHI = -BX*sin_phi + BY*cos_phi
end

;--------------------------------------------------------------------------------------
; procedure PJSO2S3
;--------------------------------------------------------------------------------------
pro PJSO2S3, BXPJSO, BYPJSO, BZPJSO, BXS3, BYS3, BZS3, phi
  sin_phi = sin(phi)
  cos_phi = cos(phi)
  BXS3 = BXPJSO*cos_phi - BYPJSO*sin_phi
  BYS3 = BXPJSO*sin_phi + BYPJSO*cos_phi
  BZS3 = BZPJSO
end

;--------------------------------------------------------------------------------------
; procedure S32PJSO
;--------------------------------------------------------------------------------------
pro S32PJSO, BXPJSO, BYPJSO, BZPJSO, BXS3, BYS3, BZS3, phi
  sin_phi = sin(phi)
  cos_phi = cos(phi)
  BXPJSO = BXS3*cos_phi + BYS3*sin_phi
  BYPJSO = -BXS3*sin_phi + BYS3*cos_phi
  BZPJSO = BZS3
end

;--------------------------------------------------------------------------------------
; function bessj1(x)
; Taken from: http://audiolab.uwaterloo.ca/~jeffb/thesis/node50.html
;--------------------------------------------------------------------------------------
function bessj1, x
  a = abs([min(x), max(x)])
  if (a[0] gt a[1]) then b = [a[1],a[0]] else b = a
  ; everything here is between -8 and 8
  if (b[1] lt 8d) then begin
     y = x*x
     bessj1 = x*(72362614232d + y*(-7895059235d + y*(242396853.1d + y*(-2972611.439d + y*( 15704.48260d + y*(-30.160366606d) ))))) $
              /(144725228442d + y*(2300535178d + y*(18583304.74d + y*(99447.43394d + y*(376.9991397d + y*1.0d)))))
     return, bessj1
  endif
  ; everything is either < or = to -8 or > or = to +8
  if (b[0] ge 8d) then begin
     ax = abs(x)
     z = 8d/ax
     y = z*z
     xx = ax - 2.356194491d
     one = x
     one[*] = 1d
     negX = WHERE(x lt 0d)
     If negX[0] ne -1 THEN one[negX] = -1d
     ;if x ge 0d then one = 1d else one = -1d
     bessj1 = sqrt(0.636619772d/ax)*(cos(xx)*(1.0d + y*(0.183105d-2 + y*(-0.3516396496d-4 + y*( 0.2457520174d-5 + y*(-0.240337019d-6) )))) $
                                     -z*sin(xx)*(0.04687499995d0 + y*(-0.88228987d-6 + y*(0.105787412d-6 + y*(-0.88228987d-6 + y*0.105787412d-6)))))*one
     return, bessj1
  endif
  ; below is a mixture of the two cases above
  y = x
  for ii = 0, n_elements(x) - 1 do y[ii] = bessj1(x[ii])
  bessj1 = y
  return, bessj1
end

;--------------------------------------------------------------------------------------
; function bessj0(x)
; Taken from: http://audiolab.uwaterloo.ca/~jeffb/thesis/node50.html
;--------------------------------------------------------------------------------------
function bessj0, x
  a = abs([min(x), max(x)])
  if (a[0] gt a[1]) then b = [a[1],a[0]] else b = a
  ; everything here is between -8 and 8
  if (b[1] lt 8) then begin
     y = x*x
     bessj0 = (57568490574d + y*(-13362590354d + y*(651619640.7d + y*(-11214424.18d + y*( 77392.33017d - 184.9052456d*y))))) $
              /(57568490411d + y*(1029532985d + y*(9494680.718d + y*(59272.64853d + y*(267.8532712d + y*1.0d)))))
     return, bessj0
  endif
  ; everything is either < or = to -8 or > or = to +8
  if (b[0] ge 8) then begin
     ax = abs(x)
     z = 8d/ax
     y = z*z
     xx = ax - 0.785398164d
     bessj0 = sqrt(0.636619772d/ax)*(cos(xx)*(1.0d + y*(-.1098628627d-2 + y*(.2734510407d-4 + y*(-.2073370639d-5 + y*.2093887211d-6)))) $
                                     - z*sin(xx)*(-.1562499995d-1 + y*(.1430488765d-3 + y*(-.6911147651d-5 + y*(.7621095161d-6 - .934945152d-7*y)))))
     return, bessj0
  endif
  ; below is a mixture of the two cases above
  y = x
  for ii = 0, n_elements(x) - 1 do y[ii] = bessj0(x[ii])
  bessj0 = y
  return, bessj0
  
end

;--------------------------------------------------------------------------------------
; function DBSJ2(x)
;--------------------------------------------------------------------------------------
function DBSJ2, x
  DBSJ2 = (2d/x)*bessj1(x) - bessj0(x)
  return, DBSJ2
end

;--------------------------------------------------------------------------------------
; function DBSJ3(x)
;--------------------------------------------------------------------------------------
function DBSJ3, x
  DBSJ3 = (4d/x)*DBSJ2(x) - bessj1(x)
  ;print, (4d/x)*DBSJ2(x) - bessj1(x)
  ;print, (8d/(x*x) - 1)*bessj1(x) - (4d/x)*bessj0(x)
  ;print, ''
  ;DBSJ3 = (8d/(x*x) - 1)*bessj1(x) - (4d/x)*bessj0(x)  ; RJW
  return, DBSJ3
end

;--------------------------------------------------------------------------------------
; procedure B_mp_perp
;--------------------------------------------------------------------------------------
pro B_mp_perp, rho, phi, x, Bperpr, Bperpf, Bperpx, Nmodes
  ; A
  A = [14.8636286417325d, -20.6733061430101, -68.9799158183312d, -1.42668213471797D-006, $
       -2.10151726837109D-004, 7.79886406933501D-006, 1.68914984127525D-008, -1.07213724493682D-007, $
       5.22885344671751D-007, -4.46205598371400D-008, 2.83331673403534D-007, -1.38338264037040D-006, $
       7.94736850421990d, 23.4906454129448d, -17.3725909963207d, -3.95338621079697D-004, $
       -4.24194742905271D-003, 1.92089873646750D-003, -6.94191156816239D-002, -0.170169494909972d, $
       -0.171632981837208d, 4.57144670517197D-009,  6.53500870745233D-007, 6.56440399411187D-007, $
       -7.67358819268301D-002, -0.133394608915796d, -0.127893292811889d, 3.60644094744725D-006 ,  $
       2.04677410835679D-005, 3.26595444981956D-005]

  ; the b(J) are
  B = [32.9198722839323d, 64.1556701660173d, 128.854949951183d, 30.4416656494141d, $
       45.2417564392090d, 107.303962707520d, 28.7235050200276d, 63.9741210937499d, $
       128.002441406250d, 33.8896331786845d, 60.2880744933923d, 119.902839660605d, $
       32.1036300655818d, 62.8152580261215d, 127.073867797851d]
  
  Modes1 = Nmodes
  Modes2 = 2*Modes1
  Modes3 = 3*Modes1
  Modes4 = 4*Modes1
  Modes5 = 5*Modes1
  ;Modes6 = 6*Modes1 ; not used

  rhom = rho - 0.1d
  rhop = rho + 0.1d

  ;phim = phi - 0.001d
  ;phip = phi + 0.001d

  phi2 = 2d*phi
  ;phi2m = 2d*phi - 0.001d
  ;phi2p = 2d*phi + 0.001d

  phi3 = 3d*phi
  ;phi3m = 3d*phi - 0.001d
  ;phi3p = 3d*phi + 0.001d

  xm = x - 0.1d
  xp = x + 0.1d

  sinf = sin(phi)
  cosf = cos(phi)
  ;sin2f = sin(phi2)
  sin2f = 2*sinf*cosf ; double angle 
  ;cos2f = cos(phi2)
  cos2f = 2*cosf*cosf - 1 ; double angle
  sin3f = sin(phi3)
  cos3f = cos(phi3)

  expb = exp(x/b)
  exp_xm_minus_xp = exp(xp/b) - exp(xm/b)
  
  sin_001 = sin(0.001d)

  ;sinP_minus_sinM = sin(phip) - sin(phim)
  sinP_minus_sinM = 2d*cos(phi)*sin_001

  ;cosP_minus_cosM = cos(phip) - cos(phim)
  cosP_minus_cosM = -2d*sin(phi)*sin_001

  ;sin2P_minus_sin2M = sin(phi2p) - sin(phi2m)
  sin2P_minus_sin2M = 2d*cos(phi2)*sin_001

  ;cos2P_minus_cos2M = cos(phi2p) - cos(phi2m)
  cos2P_minus_cos2M = -2d*sin(phi2)*sin_001

  ;sin3P_minus_sin3M = sin(phi3p) - sin(phi3m)
  sin3P_minus_sin3M = 2d*cos(phi3)*sin_001

  ;cos3P_minus_cos3M = cos(phi3p) - cos(phi3m)
  cos3P_minus_cos3M = -2d*sin(phi3)*sin_001

  rhobb = rho/b
  xmb = (x - b)
  fiveHundredoverRho = (500d/rho)

  Bperpr = 0d
  Bperpf = 0d
  Bperpx = 0d
  KA = 0

  ;--------------------------------------------------------------------------------------
  ; Now compute the terms associated with phi
  Nmodes_indgen = indgen(Nmodes)
  Nmodes_x2 = 2*Nmodes
    
  ind     = [Nmodes_indgen]
  KA1     = KA + Nmodes_indgen  ; KA = 0 on entering here - just here for clarity with other bits below
  KA_Temp = ind + Nmodes
  bessj1PminusM = bessj1(rhop/b[ind]) - bessj1(rhom/b[ind])
  bessj1RHO     = bessj1(rhobb[ind])
  angle_diff    = sinf*a[KA1] + cosf*a[KA_Temp]

  Bperpr += expb[ind]*bessj1PminusM       *angle_diff
  Bperpf += expb[ind]*bessj1RHO           *(sinP_minus_sinM*a[KA1] + cosP_minus_cosM*a[KA_Temp])
  Bperpx += bessj1RHO*exp_xm_minus_xp[ind]*angle_diff
  KA += 2*Nmodes
    
  ; Now compute the terms associated with 2*phi
  ind     = [Nmodes_indgen + Nmodes]
  KA1     = KA + Nmodes_indgen
  KA_Temp = KA1 + Nmodes
  dbsj2PminusM = DBSJ2(rhop/b[ind]) - DBSJ2(rhom/b[ind])
  dbsj2Rhobb   = DBSJ2(rhobb[ind])
  angle_diff = sin2f*a[KA1] + cos2f*a[KA_Temp]

  Bperpr += expb[ind] *dbsj2PminusM        *angle_diff
  Bperpf += expb[ind] *dbsj2Rhobb          *(sin2P_minus_sin2M*a[KA1] + cos2P_minus_cos2M*a[KA_Temp])
  Bperpx += dbsj2Rhobb*exp_xm_minus_xp[ind]*angle_diff
  KA += Nmodes_x2
    
  ; Now compute the terms associated with 3*phi
  ind     = [Nmodes_indgen + Modes2]
  KA1     = KA + Nmodes_indgen
  KA_Temp = KA1 + Nmodes
  dbsj3PminusM = DBSJ3(rhop/b[ind]) - DBSJ3(rhom/b[ind])
  dbsj3Rhobb   = DBSJ3(rhobb[ind])
  angle_diff = sin3f*a[KA1] + cos3f*a[KA_Temp]

  Bperpr += expb[ind] *dbsj3PminusM        *angle_diff
  Bperpf += expb[ind] *dbsj3Rhobb          *(sin3P_minus_sin3M*a[KA1] + cos3P_minus_cos3M*a[KA_Temp])
  Bperpx += dbsj3Rhobb*exp_xm_minus_xp[ind]*angle_diff
  KA += Nmodes_x2
    
  ; Now include the terms associated with the derivative term
  ind     = [Nmodes_indgen + Modes3]
  KA1     = KA + Nmodes_indgen
  KA_Temp = KA1 + Nmodes
  besj0 = bessj0(rhobb[ind])
  besj1 = bessj1(rhobb[ind])
  BJ0minusBJ1       = rhop*bessj0(rhop/b[ind]) - rhom*bessj0(rhom/b[ind]) - xmb[ind]*(bessj1(rhom/b[ind]) - bessj1(rhop/b[ind]))
  expxp_minus_expxm = exp(xp/b[ind])*(rho*besj0 + (xp - b[ind])*besj1) - exp(xm/b[ind])*(rho*besj0 + (xm - b[ind])*besj1)
  BJ0plusBJ1        = rho*besj0 + xmb[ind]*besj1
  angle_diff = sinf*a[KA1] + cosf*a[KA_Temp]

  Bperpr += expb[ind]        *BJ0minusBJ1 *angle_diff
  Bperpf += expb[ind]        *BJ0plusBJ1  *(sinP_minus_sinM*a[KA1] + cosP_minus_cosM*a[KA_Temp])
  Bperpx += expxp_minus_expxm*             angle_diff
  KA += Nmodes_x2
    
  ind = [Nmodes_indgen + Modes4]
  KA1  = KA + Nmodes_indgen
  KA_Temp = KA1 + Nmodes
  besj2 = DBSJ2(rhobb[ind])
  besj3 = DBSJ3(rhobb[ind])
  dbsj2minusdbsj3   = rhop*DBSJ2(rhop/b[ind]) - rhom*DBSJ2(rhom/b[ind]) - (x - 3*b[ind])*(DBSJ3(rhom/b[ind]) - DBSJ3(rhop/b[ind]))
  expxp_minus_expxm = exp(xp/b[ind])*(rho*besj2 + (xp - 3*b[ind])*besj3) - exp(xm/b[ind])*(rho*besj2 + (xm - 3*b[ind])*besj3)
  BJ2minusBJ3       = rho*besj2 + (x - 3*b[ind])*besj3
  angle_diff = sin3f*a[KA1] + cos3f*a[KA_Temp]

  Bperpr += expb[ind]        *dbsj2minusdbsj3*angle_diff
  Bperpf += expb[ind]        *BJ2minusBJ3    *(sin3P_minus_sin3M*a[KA1] + cos3P_minus_cos3M*a[KA_Temp])
  Bperpx += expxp_minus_expxm                *angle_diff
  ;--------------------------------------------------------------------------------------

  Bperpr = total(Bperpr) * 5d
  Bperpf = total(Bperpf) * fiveHundredoverRho
  Bperpx = total(Bperpx) * 5d
end

;--------------------------------------------------------------------------------------
; procedure dipole_cyl
;   not used
;--------------------------------------------------------------------------------------
pro dipole_cyl, rmap, pmap, zmap, Brd, Bpd, Bzd, counter, ctime
  vecin = (vecou = dblarr(3))

  cyl2car_pos, rmap, pmap, zmap, Xmap, Ymap, Zmap

  ; convert to jso
  vecin = [Xmap, Ymap, Zmap]
  JROT, "s3c", "jso", vecin, vecou, ctime
  B0x = 0d
  B0y = 0d
  zz = 420543d*420543d + 65920d*65920d + 24992d*24992d
  B0z = sqrt(zz) ;Vip4 dipole magnitude

  dipole, B0x, B0y, B0z, vecou[0], vecou[1], vecou[2], Bx, By, Bz

  ; convert back to s3c
  vecin = [BX, BY, BZ]
  JROT, "jso", "s3c", vecin, vecou, ctime
  car2cyl_vect, Brd, Bpd, Bzd, vecou[0], vecou[1], vecou[2], pmap
end

;--------------------------------------------------------------------------------------
; procedure dipole
;
; Calculates the field of Jupiter's dipole for shielding in the magnetopause
;
;   (B0x, B0y, B0z) is the dipole moment
;   x, y, z is the position vector
;   Bx, By, Bz is the output field vector
;--------------------------------------------------------------------------------------
pro dipole, B0x, B0y, B0z, x, y ,z, Bx, By, Bz
  r = sqrt(x*x + y*y + z*z)
  ;a = dblarr(3,3)
  a = [  [3*x*x - r*r, 3*x*y, 3*x*z], $
         [3*x*y, 3*y*y - r*r, 3*y*z], $
         [3*x*z, 3*y*z, 3*z*z - r*r]  ]
  r5 = r*r*r*r*r
  a = a/r5
  Bx = a[0,0]*B0x + a[0,1]*B0y + a[0,2]*B0z
  By = a[1,0]*B0x + a[1,1]*B0y + a[1,2]*B0z
  Bz = a[2,0]*B0x + a[2,1]*B0y + a[2,2]*B0z
end

;--------------------------------------------------------------------------------------
; procedure dipole_shielded
;
; WRITTEN BY K. K. KHURANA   11/2002
; MODIFIED BY H. K. SCHWARZL 11/2003
; PARMOD is an input array (real*8) that contains the model parameters
;   Dimension of PARMOD is: DIMENSION PARMOD(10)
;   currently just PARMOD[0] is used for the dipole tilt angle
; x,y,z input position
; OUTPUT: mag .filed Bx, By, Bz at x, y, z
;--------------------------------------------------------------------------------------
pro dipole_shielded, PARMOD, x, y, z, Bx, By, Bz
  psir = PARMOD[0]
  ;B0 = sqrt(420543d*420543d + 65920d*65920d + 24992d*24992d) ; Vip4 dipole magnitude
  B0 = 426411.1411689427332021d
  B0x = B0*sin(psir)
  B0y = 0d
  B0z = B0*cos(psir)

  ; We will  first calculate the field of the dipole
  dipole, B0x, B0y, B0z, x, y, z, Bxd, Byd, Bzd

  ; Now calculate the parallel dipole shielding field
  if (z eq 0d) && (y eq 0d) then begin
     cos_phi = 1d
     sin_phi = 0d
  endif else begin
     phi = atan(z,y)
     cos_phi = cos(phi)
     sin_phi = sin(phi)
     
  endelse
  rho = y*cos_phi + z*sin_phi
  ; Number of dipole modes
  Nmodes = 3

  ; Call B_mp_par(rho,phi,x,Brho1,Bphi1,Bx1,Nmodes) !no par dipole anymore
  Brho1 = 0d
  Bphi1 = 0d
  Bx1 = 0d
  B_mp_perp, rho, phi, x, Brho2, Bphi2, Bx2, Nmodes

  By2 = Brho2*cos_phi - Bphi2*sin_phi
  Bz2 = Brho2*sin_phi + Bphi2*cos_phi

  Bx = Bxd + Bx2 
  By = Byd + By2
  Bz = Bzd + Bz2
end

;--------------------------------------------------------------------------------------
; procedure dipole_shield_cyl_s3
;--------------------------------------------------------------------------------------
pro dipole_shield_cyl_S3, PARMOD, ctime, rmap, pmap, zmap, Brds, Bpds, Bzds, sphiOut
  cyl2car_pos, rmap, pmap, zmap, Xcar, Ycar, Zcar

  S32PJSO, xpJSO, ypJSO, zpJSO, Xcar, Ycar, Zcar, sphiOut 

  dipole_shielded, PARMOD, xpJSO, ypJSO, zpJSO, BxpJSO, BypJSO, BzpJSO

  PJSO2S3, BxpJSO, BypJSO, BzpJSO, BxS3, ByS3, BzS3, sphiOut

  car2cyl_vect, Brds, Bpds, Bzds, BxS3, ByS3, BzS3, pmap
end

;--------------------------------------------------------------------------------------
; funciton U
;
; calculates the potential terms
;--------------------------------------------------------------------------------------
;function U, a, c, M, rpvecP, rpvecR, sin1, sin3, cos1, sin2, M2, M2_mM
;  U = 0d
;  for i = 0, M - 1 do begin
;     U += total(a[i,0:M2_mM]*(rpvecP[i,M:M2])*cos1[i]*sin1[M:M2] $
;                + c[i,0:M2_mM]*(rpvecR[i,M:M2])*sin2[i]*sin3[M:M2])
;  endfor
;  return, U
;end ; U function


function U_B, x1, y1, z1, a, c, p_1, p_2, r_1, r_2

  Up = [0d, 0d, 0d]
  Um = [0d, 0d, 0d]
  B = [0d, 0d, 0d]

  xp = [x1 + 0.001d, x1, x1]
  yp = [y1, y1 + 0.001d, y1]
  zp = [z1, z1, z1 + 0.001d]
  xm = [x1 - 0.001d, x1, x1]
  ym = [y1, y1 - 0.001d, y1]
  zm = [z1, z1, z1 - 0.001d]

  term_p = (term_r = dblarr(8,8))
  for i = 0, 7 do begin
        term_p[i,0:7] = 1d/(p_1[i]*p_1[i]) + 1d/(p_2[0:7]*p_2[0:7])
        term_r[i,0:7] = 1d/(r_1[i]*r_1[i]) + 1d/(r_2[0:7]*r_2[0:7])
  endfor
  term_p = sqrt(term_p)
  term_r = sqrt(term_r)

  for i = 0, 7 do begin
     cos_yp = cos(yp/p_1[i])
     cos_ym = cos(ym/p_1[i])
     sin_yp = sin(yp/r_1[i])
     sin_ym = sin(ym/r_1[i])
     for k = 0, 7 do begin
        Up = Up + a[i,k]*exp(term_p[i,k]*xp)*cos_yp*sin(zp/p_2[k]) + $
             c[i,k]*exp(term_r[i,k]*xp)*sin_yp*sin(zp/r_2[k])

        Um = Um + a[i,k]*exp(term_p[i,k]*xm)*cos_ym*sin(zm/p_2[k]) + $
             c[i,k]*exp(term_r[i,k]*xm)*sin_ym*sin(zm/r_2[k])
     endfor
  endfor

  B = Up - Um
  return, B
end

;--------------------------------------------------------------------------------------
; procedure B_tail_shield
;--------------------------------------------------------------------------------------
pro B_tail_shield, x, y, z, Bxout, Byout, Bzout, M, ModeIn
  ; declare constants a, p, c, r
  aMode = (cMode = dblarr(6, 8, 8))
  p = dblarr(8, 16)
  r = dblarr(8, 16)

  ; coeff a
  aMode[0,0,*] = [0.13217647646022754326d-20, 0.50643011679165539362d-16, -0.62311994734606299672d-13, 0.73000393989897114366d-10, $
                  -0.10249726242696675093d-07, 0.32914667786334130816d-06, -0.24800841788149776689d-05, 0.37244082145067776146d-05]
  aMode[0,1,*] = [0.22584041437827493403d-17, 0.20229608874904410065d-12, -0.13496084238172885161d-09, 0.19512264195191708182d-06, $
                  -0.28536292479920302156d-04, 0.92636894313957363067d-03, -0.69993880953881300044d-02, 0.10518459214365802889d-01]
  aMode[0,2,*] = [0.36523294131124144357d-14, -0.64420435736454777497d-09, -0.57096883444695674114d-07, -0.27055379124230958254d-04, $
                  0.44272407289722757184d-02, -0.14889477801070045259d+00, 0.11356290643917379412d+01, -0.17106683807405858033d+01]
  aMode[0,3,*] = [-0.18416420070040647516d-11, 0.27046079272174843310d-06, 0.44203411790099522704d-04, 0.12944858989375875779d-03, $ 
                  -0.51975317855123419619d-02, 0.20369979936316231494d+00, -0.16254799676715984801d+01, 0.24787158145629919481d+01]
  aMode[0,4,*] = [0.13636392185761578854d-09, -0.20156479561298437097d-04, -0.39633156929923650579d-02, -0.13365044136678470021d-01, $
                  -0.10442458756051709034d-01, 0.11380270847337552453d+00, -0.88061814759885734815d+00, 0.13114182967621741404d+01]
  aMode[0,5,*] = [-0.22305468438407642928d-08, 0.33146982079580284974d-03, 0.69006040956576439882d-01, 0.27321338618451562751d+00, $
                  0.15256440935574380191d+00, 0.82160642450780390078d-01, -0.38978676634608163453d+00, 0.47434195602775641731d+00]
  aMode[0,6,*] = [0.84447121731127730015d-08, -0.12569305425012970989d-02, -0.26569849211645015785d+00, -0.11069710175423441711d+01, $
                  -0.64840056779980832502d+00, -0.97562574536155111104d-02, 0.15631498731831953818d+00, -0.16203940584912954747d+00]
  aMode[0,7,*] = [-0.63486912141303140089d-08, 0.94534756045114320954d-03, 0.20061320014222592256d+00, 0.84751706801178663397d+00, $
                  0.50107973765980506897d+00, -0.14177955455321484379d-03, 0.56476656376811078530d-01, -0.52516865628955500255d-02]
  
  aMode[1,0,*] = [0.13619937133704040910d-23, 0.43621985467324879692d-18, -0.15098830327683961272d-14, 0.63123957003196471404d-11, $
                  -0.11978885144275235319d-08, 0.41934992058305047279d-07, -0.31951294540338555094d-06, 0.48135741458667444803d-06]
  aMode[1,1,*] = [0.24119037788698203250d-20, 0.40881122046673903369d-14, -0.21545235983615116381d-11, 0.11640656169328844615d-07, $
                  -0.21959806946493167778d-05, 0.76662481447507042631d-04, -0.58366072779967188566d-03, 0.87913032319456423380d-03]
  aMode[1,2,*] = [0.31144339871300497080d-16, -0.22620785858629064435d-10, -0.74029445112059493183d-08, -0.93859788335645291112d-07, $
                  0.82317630609657417295d-04, -0.36569184498810431982d-02, 0.29432171359264206245d-01, -0.44946922231727173269d-01]
  aMode[1,3,*] = [-0.58080626963236774429d-14, 0.39682408099390080735d-08, 0.16840660739515575627d-05, 0.38335076286793983157d-05, $
                  0.20836381731764990199d-03, -0.33920492396785180133d-02, 0.15820849562165701485d-01, -0.19644945238949279797d-01]
  aMode[1,4,*] = [0.41076898294619867968d-12, -0.27989015436962660920d-06, -0.13840368413368968703d-03, -0.70869170102476930495d-03, $
                  0.26033168308504186505d-03, -0.24526926326940952094d-02, -0.35013567157541127805d-01, 0.69275925576020052076d-01]
  aMode[1,5,*] = [-0.70799772156343383500d-11, 0.48274807820722331896d-05, 0.24974637159216332982d-02, 0.15464887864241698700d-01, $
                  0.11439264646130791192d-01, 0.38097770315818785036d-01, -0.13590064465091589163d+00, 0.17180997533186371128d+00]
  aMode[1,6,*] = [0.26783558161963370025d-10, -0.18266964017400193043d-04, -0.95586630597219102156d-02, -0.62918902435418315732d-01, $
                  -0.64039968185764291064d-01, -0.22222194307812768165d-03, 0.50813426126672718297d-01, -0.46991112128761134414d-01]
  aMode[1,7,*] = [-0.20108572992030198101d-10, 0.13715428126943252085d-04, 0.71979427421967230316d-02, 0.48199320155821636646d-01, $
                  0.49407199622570425745d-01, 0.92225425694026483824d-02, 0.53857653772015723347d-02, 0.91367943357959795491d-03]
  
  aMode[2,0,*] = [-0.13119415239474609968d-23, -0.60256122420913360571d-18, 0.37920759812307474057d-14, -0.18206762505356291370d-10, $
                  0.26774304351790334521d-08, -0.87684408227198336049d-07, 0.66457031757235851543d-06, -0.99948895864718618753d-06]
  aMode[2,1,*] = [-0.25960636282069127211d-19, -0.42000012814022671392d-16, 0.74796457397241065123d-11, -0.25318674869395683124d-07, $
                  0.34976874196721219334d-05, -0.11242351418125537954d-03, 0.84791358984972013956d-03, -0.12736554117114038398d-02]
  aMode[2,2,*] = [0.76343432773888704190d-16, 0.64378927275727946266d-12, 0.16598929643645348619d-08, -0.92837331289370776943d-05, $
                  0.14285390095735186477d-02, -0.47503134695174127344d-01, 0.36153292729534660665d+00, -0.54431063875700962384d+00]
  aMode[2,3,*] = [-0.25391252497671850107d-13, -0.31368629291300997863d-09, 0.14987071996970728449d-05, 0.17824032262389087222d-04, $
                  0.11355427076600370650d-02, -0.37610982091768150326d-01, 0.28513697308043282063d+00, -0.42869687938445020236d+00]
  aMode[2,4,*] = [0.17211923850683952252d-11, 0.31101556003157209140d-07, -0.12650545583107459801d-03, -0.30769692479640595728d-02, $
                  -0.45994232482610790668d-02, 0.15393025113972294448d-01, 0.86726269138905074385d-02, -0.41633915235090732664d-01]
  aMode[2,5,*] = [-0.27638896781422706006d-10, -0.55187436783116812222d-06, 0.21199958794900246594d-02, 0.61391435739846613728d-01, $
                  0.52329651593674135767d-01, -0.83734155027368615265d-01, 0.12671335310678004670d+00, -0.18213430508768142956d+00]
  aMode[2,6,*] = [0.10413809392985686752d-09, 0.21324280383467479893d-05, -0.80597222736636293660d-02, -0.24311900160452664110d+00, $
                  -0.13586934681330335994d+00, 0.10699875405203342904d-02, -0.67114222283495701404d-01, 0.56104987912055399590d-02]
  aMode[2,7,*] = [-0.78195074731598692707d-10, -0.16113432047677196834d-05, 0.60646416575606716392d-02, 0.18453527959746301334d+00, $
                  0.97748021711249979404d-01, -0.86931770826402736673d-01, -0.25859340706079221305d-01, 0.30535054676239381521d-01]
  
  aMode[3,0,*] = [-0.11168229337158983582d-20, -0.96888835482058599524d-16, 0.28484097342054432999d-12, -0.88991894939124396302d-09, $
                  0.13436311859750933450d-06, -0.43949492953948432472d-05, 0.33267991354851349505d-04, -0.50016760854385964307d-04]
  aMode[3,1,*] = [-0.20765635100251089717d-17, -0.36197386357992140659d-12, 0.35732517136988684036d-09, -0.11320141902699047964d-05, $
                  0.16737866108943116216d-03, -0.54399134640077706492d-02, 0.41108560676656589194d-01, -0.61778272173930721677d-01]
  aMode[3,2,*] = [-0.49227064641270947831d-14, 0.13773911494082902162d-08, 0.28000536764000152345d-06, -0.60982707102378110874d-05, $
                  -0.45545917365779216012d-03, 0.30474617144202409413d-01, -0.26240435518843665541d+00, 0.40666851745242809101d+00]
  aMode[3,3,*] = [0.24046710130055206633d-11, -0.56928231643345874601d-06, -0.14658425595003152785d-03, -0.54669311790752939117d-03, $
                  0.63289641147589756897d-01, -0.23416266019622873351d+01, 0.18376229679993183907d+02, -0.27890706711124808592d+02]
  aMode[3,4,*] = [-0.17790734879856566763d-09, 0.41757078162568088686d-04, 0.12661549470502249103d-01, 0.30791263485754920559d-01, $
                  0.52395704373591174274d-01, -0.15856752943903087427d+01, 0.11848914793373501730d+02, -0.17521753728113509396d+02]
  aMode[3,5,*] = [0.29086256117761690731d-08, -0.68275932316891383422d-03, -0.21784801486997875663d+00, -0.66623989758662736093d+00, $
                  -0.28050159504711622560d+00, -0.67774272696429793683d+00, 0.10271422822879689995d+01, -0.63392178592402839143d+00]
  aMode[3,6,*] = [-0.11009990780568283952d-07, 0.25849182195501789749d-02, 0.83614032826616710991d+00, 0.27685938326621157834d+01, $
                  0.16426046025686344975d+01, 0.78572426172170208857d+00, -0.18258375964853595263d+01, 0.12152196449017165225d+01]
  aMode[3,7,*] = [0.82768727273709110647d-08, -0.19433482187601436308d-02, -0.63081127374393464180d+00, -0.21366815638955616307d+01, $
                  -0.13446774521063316054d+01, -0.58592309940003861612d+00, 0.31827661237739279798d+00, -0.20240366730373042791d+00]
  
  aMode[4,0,*] = [-0.13947208148057608312d-20, -0.20210120268149411870d-15, 0.45539475964762736737d-12, -0.11170239869341296312d-08, $
                  0.16610385153669340319d-06, -0.54155353001624000341d-05, 0.40961114502502580236d-04, -0.61571027127369957199d-04]
  aMode[4,1,*] = [-0.15010983747727246750d-17, -0.56103448860629763217d-12, 0.43831990483680147718d-09, -0.10609291877688533656d-05, $
                  0.15387190085182278487d-03, -0.49799395553291745386d-02, 0.37594128280158018995d-01, -0.56482482158755402679d-01] 
  aMode[4,2,*] = [-0.84949700310157272298d-14, 0.19131227654780431635d-08, 0.23797873673175664599d-06, 0.67357172013098089990d-04, $
                  -0.11412327840147937774d-01, 0.38787752218711162299d+00, -0.29661484451688613361d+01, 0.44710237927891309794d+01]
  aMode[4,3,*] = [0.36903300990875096410d-11, -0.78266544732189755606d-06, -0.14356129913638142170d-03, -0.52616304038202654780d-03, $
                  0.50792590874790795041d-01, -0.18789041423837609556d+01, 0.14745664699188227864d+02, -0.22381820221529649117d+02]
  aMode[4,4,*] = [-0.26744771000491356360d-09, 0.57361959249991540943d-04, 0.12466277937525278574d-01, 0.35752121631968263315d-01, $
                  0.49443768069258187125d-01, -0.11592738471082895124d+01, 0.85575941704881230975d+01, -0.12621255590568500881d+02]
  aMode[4,5,*] = [0.43525625631567042006d-08, -0.93805114409436534117d-03, -0.21484367014734804257d+00, -0.75606598389073855770d+00, $
                  -0.38916312244605659742d+00, -0.45744351551847461934d+00, 0.91714655267490634571d+00, -0.70503052173487583687d+00]
  aMode[4,6,*] = [-0.16457579884296458239d-07, 0.35517012004150618764d-02, 0.82496531512910689087d+00, 0.31057726963483678339d+01, $
                  0.18813656203693941648d+01, 0.51092359576032500001d+00, -0.12428990388842335867d+01, 0.89617638629098603786d+00]
  aMode[4,7,*] = [0.12368783158269132105d-07, -0.26702315450896003667d-02, -0.62244835948407777337d+00, -0.23879836056356684714d+01, $
                  -0.14962346641346604414d+01, -0.32180960449888429408d+00, 0.14080714553385571541d+00, -0.12793599129792307955d+00]
  
  aMode[5,0,*] = [-0.71007031654422592126d-20, -0.82331351329343824829d-15, 0.14860700245019204501d-11, -0.18336776863420102046d-08, $
                  0.26381273327308489839d-06, -0.85305767261458509409d-05, 0.64390806859274407614d-04, -0.96740352003902945199d-04]
  aMode[5,1,*] = [-0.47759968145964668551d-17, -0.11110703618965760419d-11, 0.93331659872542349631d-09, -0.11361217227825839426d-05, $
                  0.16101720194965896126d-03, -0.51853068878838977084d-02, 0.39098877702660601585d-01, -0.58726418277052969685d-01]
  aMode[5,2,*] = [-0.20561632805425666958d-13, 0.30615328304455808883d-08, 0.24271749535760407390d-06, 0.15937085189299786236d-03, $
                  -0.25318551119527441528d-01, 0.84552595243637060917d+00, -0.64376874089793245659d+01, 0.96932895116823303283d+01]
  aMode[5,3,*] = [0.87503551825887093684d-11, -0.13163374744119336057d-05, -0.21315095274792228430d-03, -0.75679626156015320503d-03, $
                  0.43875888066817978483d-01, -0.16187067723465442981d+01, 0.12704254977048992092d+02, -0.19288846975099222191d+02]
  aMode[5,4,*] = [-0.62980025164513007140d-09, 0.99033147185834717873d-04, 0.19286483194829802556d-01, 0.67190391627066974322d-01, $
                  0.59018569821895905391d-01, -0.57675212643809139478d+00, 0.41759144712554761014d+01, -0.60616007357605896643d+01]
  aMode[5,5,*] = [0.10235134109675949609d-07, -0.16329320547735449054d-02, -0.33676325501000334838d+00, -0.13741116541224409619d+01, $
                  -0.82897929864310349046d+00, -0.29064544069947594095d+00, 0.99716723490406433683d+00, -0.11955304085306053352d+01]
  aMode[5,6,*] = [-0.38688354152446052580d-07, 0.61964020776766650655d-02, 0.12977027325262353585d+01, 0.55698639583852660450d+01, $
                  0.35055306218486181890d+01, 0.24501453852569592406d+00, -0.51229945875859312920d+00, 0.42642747170716877036d+00]
  aMode[5,7,*] = [0.29074290348871754119d-07, -0.46611910769336493132d-02, -0.98002188003461245813d+00, -0.42652335402388423801d+01, $
                  -0.27213526220765031915d+01, -0.21604814310942872523d+00, -0.40325674485435270000d+00, 0.12987860657487151350d+00]

  ; coeff c
  cMode[0,0,*] = [0.62608370709224905326d-18, -0.69296927042311144973d-15, 0.53793142672880289723d-22, 0.55953376048336647130d-13, $
                  -0.13644845214970875435d-11, 0.51560868494277940499d-10, -0.38841223167697958018d-09, 0.58169693048278663383d-09]
  cMode[0,1,*] = [0.79948744245788159190d-02, 0.28880230453715999061d-01, 0.18816017618889148366d-03, 0.30432220614412108794d-01, $
                  0.20335251312180755434d-01, 0.12253559136086547010d+01, -0.11162104253917091156d+02, 0.17710281832472645646d+02]
  cMode[0,2,*] = [-0.30562966462902028119d+03, -0.11990795349370895195d+04, -0.65802994518949660118d+01, -0.12273805927634358070d+04, $
                  -0.65349097657709229736d+03, -0.12824769943213845024d+03, -0.88590404530202224719d+02, -0.92014824378423334394d+02]
  cMode[0,3,*] = [-0.95152195796331646704d+02, -0.36956353969242154988d+03, -0.20695307247843501841d+01, -0.38362497103275590148d+03, $
                  -0.22776414558915512031d+03, -0.46624811555373160132d+02, 0.45685596890036768158d+02, -0.25467538299896208542d+03]
  cMode[0,4,*] = [0.32918009134216310584d+03, 0.12851583684685703445d+04, 0.71222223072393315845d+01, 0.13249069943127731452d+04, $
                  0.74436590714498462872d+03, 0.14749387222546454623d+03, 0.42968129090712094964d+02, 0.33666974003731491293d+03]
  cMode[0,5,*] = [0.20156315877857942098d+03, 0.80514059433962259504d+03, 0.42634222620010087112d+01, 0.79880139299917756190d+03, $
                  0.34272486825284445011d+03, 0.52606372363987690121d+02, 0.74826789775518696146d+01, -0.49206154796489274261d+01]
  cMode[0,6,*] = [-0.53421252140832686805d+03, -0.21573846013695305856d+04, -0.11180286753450849879d+02, -0.20919394990520148169d+04, $
                  -0.77083430972020012816d+03, -0.74269446238765679524d+02, 0.86081879160377337001d+01, -0.33189293006747679903d+01]
  cMode[0,7,*] = [0.64042667076955872573d+03, 0.25941583461638844099d+04, 0.13364255832056319839d+02, 0.24982768986941823463d+04, $
                  0.87945971890665362025d+03, 0.70070351276938351858d+02, -0.86139878258589916981d+01, 0.30011846979293261838d+01]

  cMode[1,0,*] = [-0.63013334128716040893d-13, 0.60865772431920603935d-05, -0.40016262502772708131d-09, -0.23069864565108759713d-04, $
                  0.23255146051530801720d-03,-0.51951144258714805346d-02, 0.36136455088219365805d-01, -0.53017682607232261560d-01]
  cMode[1,1,*] = [0.55436122531818750047d-03, -0.66992043381081387565d-02, 0.83803970952918511727d-02, 0.80754210845424374554d-02, $
                  -0.61427751801047527635d-02, 0.54122156026262873140d+00, -0.43010880044939874267d+01, 0.65777678112627633311d+01]
  cMode[1,2,*] = [-0.23084376754368056694d+02, 0.20899246010739616075d+03, -0.36836390749521057408d+03, -0.52914996781002994197d+03, $
                  -0.20007855689466294002d+03, -0.36056061195965001253d+02, -0.16132602175333610183d+02, 0.21643026623710106548d+01]
  cMode[1,3,*] = [-0.30489331973819799870d+01, 0.33304902492752548326d+02, -0.47908018506183607243d+02, -0.66968634395671120529d+02, $
                  -0.34606055921065168590d+02, -0.74942275618948936966d+01, 0.20818848022147209420d+02, -0.68532587218682348151d+02]
  cMode[1,4,*] = [0.44092752992499013586d+01, -0.47413797913825259655d+02, 0.69401964003076273002d+02, 0.97614272148050673649d+02, $
                  0.48704912842393142113d+02, 0.99941378554753921292d+01, -0.14570861805562320689d+02, 0.65068431858518227528d+02]
  cMode[1,5,*] = [0.27754841758044097588d+02, -0.24286715291266780525d+03, 0.44379959765035561503d+03, 0.63697225078825798760d+03, $
                  0.23079424266757904149d+03, 0.39146925311929221535d+02, 0.12991883575647105075d+02, -0.34552170740359109402d+01]
  cMode[1,6,*] = [-0.20272909693943255149d+02, 0.14757207634388509021d+03, -0.32708628601917539846d+03, -0.46225115463536488036d+03, $
                  -0.13753186496272293837d+03, -0.14358017591372269627d+02, 0.41733787941502891172d+00, -0.14829734668187728452d+00]
  cMode[1,7,*] = [0.22054195737479753702d+02, -0.15104824369566767217d+03, 0.35667767962821925742d+03, 0.50036616092734789162d+03, $
                  0.14078440554104824755d+03, 0.12025910667067751802d+02, -0.70328949883151059552d+00, 0.24434145951697048282d+00]

  cMode[2,0,*] = [-0.28862183650305684778d-01, 0.16202748718849406373d-02, 0.62452442623701900359d-08, -0.10292552636752272388d-04, $
                  0.14899979598450072693d-10, 0.54701873831866327790d+00, -0.33801096158365382393d+01, 0.46663975792709919687d+01]
  cMode[2,1,*] = [0.82086900572783427776d+00, 0.66781995371601139410d+00, 0.17689517609606757453d+01, 0.18458732991639696052d+01, $
                  -0.13408206456871818446d-02, -0.55045820392109447993d+01, 0.34057914707574474810d+02, -0.44172744618314503384d+02]
  cMode[2,2,*] = [-0.33440561402073085695d+04, -0.63006668870531044035d+04, -0.16148889754716360123d+05, -0.16734685574071548330d+05, $
                  0.98894432761575998824d+01, -0.10391741179165414621d+04, -0.17857836173455567951d+04, 0.29803489709704852117d+04]
  cMode[2,3,*] = [0.31915064954348162373d+04, 0.59270193788343030760d+04, 0.15182397380443129364d+05, 0.15734311019448650625d+05, $
                  -0.93434105903445914265d+01, 0.99260524095183413351d+03, 0.20341329769010583206d+04, -0.32962143597458424260d+04]
  cMode[2,4,*] = [-0.92169284592923155230d+02, -0.14814666530740792538d+03, -0.37839662711403616590d+03, -0.39255805350361141492d+03, $
                  0.24567231705288499199d+00, -0.11195431933835682247d+02, -0.27154689478860825069d+03, 0.37730991124653248114d+03]
  cMode[2,5,*] = [0.54075863200485683179d+03, 0.13068947201920562140d+04, 0.34015030753731947399d+04, 0.35222027482034166112d+04, $
                  -0.19372427035775936943d+01, 0.10115813738135293053d+03, -0.23832476149210184424d+02, -0.16890808415607296844d+02]
  cMode[2,6,*] = [-0.77613909679256050111d+03, -0.22744934351169181496d+04, -0.60171788754328519033d+04, -0.62294563978631467549d+04, $
                  0.32775011290202917813d+01, -0.72131812461142832404d+02, 0.27116603294455088324d+02, -0.71224576605140779150d+01]
  cMode[2,7,*] = [0.69369255218888747904d+03, 0.22078310233442990373d+04, 0.58858136079721381506d+04, 0.60932655114670426499d+04, $
                  -0.31519955877359868701d+01, 0.44446156981790601037d+02, -0.14233922505288183479d+02, 0.34387488557126792976d+01]
  
  cMode[3,0,*] = [-0.36868002035455376130d-21, 0.10541748946585167701d-17, 0.27755566547052099579d-13, -0.34780372211895671519d-15, $
                  -0.99328472275559516191d-12, 0.37378357885053432596d-10, -0.27315035342737905565d-09, 0.40965931239852784173d-09]
  cMode[3,1,*] = [-0.33737087677076624814d-04, -0.11396625664731094840d-02, -0.10789923154822660400d-01, -0.80796849963401182748d-02, $
                  0.32166227895723737972d-01, -0.25571424523548520468d+01, 0.21267119641439142796d+02, -0.33078925310689939465d+02]
  cMode[3,2,*] = [-0.17479871537929662750d+01, -0.72444600080193755076d+02, -0.55033093275560522883d+03, -0.50210145669063086515d+03, $
                  -0.25844521786360195036d+03, -0.55113408951058726614d+02, 0.20187135793924890769d+03, -0.81289982453232934034d+03]
  cMode[3,3,*] = [0.78603614447709357904d+00, 0.32329428441630612134d+02, 0.24762820448825531016d+03, 0.22413360241783610860d+03, $
                  0.11996249854833778059d+03, 0.30448969122459486058d+02, -0.23988687682537666034d+03, 0.67671098556149500424d+03]
  cMode[3,4,*] = [0.12159934205840530196d+01, 0.51016427583112839982d+02, 0.38232595903016939331d+03, 0.35346430190826185757d+03, $
                  0.17043627550273194870d+03, 0.36343814970370997841d+02, 0.25302217119464871508d+02, 0.16017628158741096910d+03]
  cMode[3,5,*] = [-0.13415948717761450037d+01, -0.59032589809690589888d+02, -0.41777999557453684431d+03, -0.40902457469604014406d+03, $
                  -0.14475183395103458749d+03, -0.29478536348031818548d+02, -0.13387094046024579085d+01, 0.23831368497822609242d+01]
  cMode[3,6,*] = [0.55621147775238517496d+01, 0.24890779339281028370d+03, 0.17228008107309571883d+04, 0.17256605312262474072d+04, $
                  0.52806415548478540245d+03, 0.64772764588200661961d+02, -0.12923255027738644784d+02, 0.58752559976531326668d+01]
  cMode[3,7,*] = [-0.73285857760141963623d+01, -0.32941251542733325230d+03, -0.22662989686584209536d+04, -0.22842752627361360140d+04, $
                  -0.66988125859857134969d+03, -0.67321294156948754405d+02, 0.12094082762406057618d+02, -0.50344380663030463551d+01]

  cMode[4,0,*] = [0.57951308146739251014d+01, -0.14999087019121093433d+03, -0.19559326802070165385d+01, -0.66720804039446477418d+03, $
                  -0.75756207610197607849d+02, -0.79695623597154439110d+02, 0.12741669724908954997d+04, -0.25779398744795645193d+04]
  cMode[4,1,*] = [0.25559373539070939784d+03, -0.26954324048080380293d+04, 0.42741834334120607508d+02, -0.13609107834915077361d+05, $
                  -0.18219129315199742435d+04, -0.49540503909243787106d+03, -0.86525145704320394202d+03, -0.19349682251463748983d+04]
  cMode[4,2,*] = [-0.87923543175638503299d+05, 0.30861003822923382955d+06, -0.25217316373834348652d+05, 0.21333098058144601694d+07, $
                  0.33404584481156045505d+06, 0.42354436002087458845d+05, 0.58637814837743640339d+04, -0.30291614862820930298d+04]
  cMode[4,3,*] = [0.21501887419077996277d+06, -0.98416860975998154970d+06, 0.50756266605703403982d+05, -0.58885368624678218196d+07, $
                  -0.88230632458050202160d+06, -0.24671074506661172520d+06, -0.14496629047056444505d+06, 0.68193611837692298394d+05]
  cMode[4,4,*] = [0.18016800297232990146d+02, 0.11854731605905291047d+00, 0.41242056558479642802d-02, 0.10788775915866732901d+00, $
                  0.76539029541206167195d+00, 0.67292188204585867694d+02, -0.31151519702859546967d+03, 0.44202304203357769551d+03]
  cMode[4,5,*] = [0.20068140169309098830d+06, -0.93099363016509286694d+06, 0.47215756436253171202d+05, -0.55473705494319487385d+07, $
                  -0.82963829864856943885d+06, -0.23595446376207669381d+06, -0.14387731979865256981d+06, 0.66726697732453104094d+05]
  cMode[4,6,*] = [-0.10659833147754620430d+06, 0.36525155976544767533d+06, -0.32972599682804299980d+05, 0.26280915665374022793d+07, $
                  0.41399542203693986408d+06, 0.40200322950430580348d+05, 0.33714384362519780324d+04, -0.16170695162643912823d+04]
  cMode[4,7,*] = [0.45385852289903034773d+05, -0.15367285853252436567d+06, 0.15938904084451054998d+05, -0.11759967065602889846d+07, $
                  -0.18617563074881875451d+06, -0.11629254312624910383d+05, -0.37356632983270783299d+03, 0.13368622059973152005d+03]
  
  cMode[5,0,*] = [0.37071935082976326114d-15, 0.46201169190899653571d-08, -0.62075147292222574435d-11, -0.36146708876274598054d-06, $
                  0.87046861427880060091d-05, -0.10076975633470059978d-03, 0.89298240723925399464d-03, -0.14060364310276510124d-02]
  cMode[5,1,*] = [-0.56188781974792609830d-01, -0.40650854381275376425d+01, -0.14988738284605299000d+01, -0.38632728626301734209d+01, $
                  -0.22261011246466999580d+01, -0.14212548821848476343d+02, 0.19383430976405556123d+03, -0.34806064213022160913d+03]
  cMode[5,2,*] = [0.27593841470781566016d+03, 0.22140385471263797079d+05, 0.78362796461845718454d+04, 0.19668830633326862766d+05, $
                  0.10181292569179314355d+05, 0.39096643822536609746d+04, 0.26876686999693699675d+04, 0.33310618692009836827d+04]
  cMode[5,3,*] = [0.44061551736067183782d+02, 0.34746817101944778016d+04, 0.12384774155446498511d+04, 0.31312056838986070950d+04, $
                  0.17260274467090328087d+04, 0.70922681375271823256d+03, -0.30968222241454737009d+03, 0.29976750308301927105d+04]
  cMode[5,4,*] = [-0.29826526924997253331d+03, -0.23828019858653441964d+05, -0.84485266842828430355d+04, -0.21247202725657969857d+05, $
                  -0.11181459854792739072d+05, -0.43708622712852429614d+04, -0.24560179552632419586d+04, -0.60645360750106398484d+04]
  cMode[5,5,*] = [-0.74893303991611572811d+02, -0.62909773314430736945d+04, -0.21848887114479706994d+04, -0.53368828632132494504d+04, $
                  -0.22162146435679854761d+04, -0.52579256909996399116d+03, -0.11309723650631768876d+03, 0.63863048387678587047d+02]
  cMode[5,6,*] = [0.19717911116782794067d+03, 0.16899744672548959734d+05, 0.58199397638247170050d+04, 0.13983841291452328015d+05, $
                  0.50676093636595025415d+04, 0.75432123378980122652d+03, 0.10996364977646411276d+02, -0.48990730626876937137d+01]
  cMode[5,7,*] = [-0.22796489432205970793d+03, -0.19665709629971770411d+05, -0.67538957248075544015d+04, -0.16130393271007925193d+05, $
                  -0.55557948391593710013d+04, -0.67829337377268315023d+03, 0.63691256079676055179d+01, -0.36820923261738309761d+01]

  ; coeff p
  p[0,*] = [0.32901733177624610249d-01, 0.54475681390900216882d-01, 0.11033542666255162778d+00, $
            0.21824017019698205288d+00, 0.43540342906595910221d+00, 0.87057139406808232706d+00, $
            0.17411041057694376377d+01, 0.34822027229232799250d+01, 0.29391078654639329670d-01, $
            0.53092570701901200536d-01, 0.10068022495463000431d+00, 0.21707901052684266396d+00, $
            0.43526892089758080217d+00, 0.87054714666225887498d+00, 0.17411005963420354447d+01, $
            0.34822023079824373503d+01]
  p[1,*] = [0.32535198971349430507d-01, 0.57792553225847100861d-01, 0.11708104989648597804d+00, $
            0.21135795367821805790d+00, 0.42489217791799518408d+00, 0.86257818103884389415d+00, $
            0.17181130605554559842d+01, 0.34339371403040641617d+01, 0.25660275820190880935d-01, $
            0.45569108067616728163d-01, 0.86363227093303382986d-01, 0.20636215029847653212d+00, $
            0.42568882948609800820d+00, 0.85912749746654117899d+00, 0.17163445496343952001d+01, $
            0.34350464850833977159d+01]
  p[2,*] = [0.32114531391568896800d-01, 0.54668167014354631660d-01, 0.10999415084820487464d+00, $
            0.21830735865316683863d+00, 0.43522952245694588313d+00, 0.87059254539023882557d+00, $
            0.17410412227916809868d+01, 0.34821651743153156921d+01, 0.25685870467860096866d-01, $
            0.45698157579046663201d-01, 0.90356977415432240263d-01, 0.21884901172766766386d+00, $
            0.43511983572804924236d+00, 0.87029548334738358050d+00, 0.17410538620391950992d+01, $
            0.34822549594214162738d+01]
  p[3,*] = [0.32356337372804593321d-01, 0.54152800183054496940d-01, 0.10961184933750356407d+00, $
            0.21794853515003347332d+00, 0.43533813536787366871d+00, 0.87056003928681793269d+00, $
            0.17411018194244403112d+01, 0.34822024269907854154d+01, 0.27588798509226557520d-01, $
            0.48812593148800704767d-01, 0.93213066985517798457d-01, 0.21740973446223286202d+00, $
            0.43526939652742839825d+00, 0.87054487704825334049d+00, 0.17411007942523434977d+01, $
            0.34822023168779896451d+01]
  p[4,*] = [0.32735344319694603676d-01, 0.54347222488222799441d-01, 0.10989084818706484902d+00, $
            0.21805120636819852464d+00, 0.43536325140487690532d+00, 0.87056578462259590622d+00, $
            0.17411028067887455605d+01, 0.34822028139735676788d+01, 0.27779830466381087994d-01, $
            0.50806515286805851161d-01, 0.95276608648916276678d-01, 0.21745769582596392588d+00, $
            0.43527860088967376128d+00, 0.87054808117661632849d+00, 0.17411000172672461694d+01, $
            0.34822021788641484008d+01]
  p[5,*] = [0.33568077258163553366d-01, 0.53905902840230535133d-01, 0.10958409128391679487d+00, $
            0.21675427898182215713d+00, 0.43240630217346920360d+00, 0.86456015794363612059d+00, $
            0.17290777196751532684d+01, 0.34581494323600296958d+01, 0.28810823097181398111d-01, $
            0.53751387895676616679d-01, 0.10088751776373805491d+00, 0.21567365023402529367d+00, $
            0.43226560363373209838d+00, 0.86453474685580449232d+00, 0.17290747052830268692d+01, $
            0.34581487532588752742d+01]

  ; coeff r
  r[0,*] = [0.23400582376918621640d-01, 0.29906910603763133593d+00, 0.59419217080359816307d+00, $
            0.51683309796288279258d+00, 0.55280210585007472090d+00, 0.89024092184836405294d+00, $
            0.17236578859735887547d+01, 0.34342717240262339295d+01, 0.89287729231317083389d-01, $
            0.15439963383324775136d+00, 0.50639243810157408276d-01, 0.24456407532792816539d+00, $
            0.39627924482769603997d+00, 0.86233612861692066076d+00, 0.17188649618875711411d+01, $
            0.34339261899212027984d+01]
  r[1,*] = [0.70036641254021372304d-01, 0.33297901220247387854d+00, 0.80631162662909527938d+00, $
            0.54483599484759341891d+00, 0.56290677049809243470d+00, 0.87422019588133075274d+00, $
            0.17295204193052738261d+01, 0.34651766994561397083d+01, 0.57115996571135090320d-01, $
            0.23559291883706992010d+00, 0.11690775021980281955d+00, 0.26929992715503474620d+00, $
            0.42836762405666322095d+00, 0.87049839288489376798d+00, 0.17350924596759142559d+01, $
            0.34557476257228634253d+01]
  r[2,*] = [0.11127103098086079668d+00, 0.42445608810651771491d+00, 0.71058432225687280237d+00, $
            0.69760620618057709307d+00, 0.58851822769682966551d+00, 0.10410712531090018373d+01, $
            0.18016710591015945297d+01, 0.35162043268667866335d+01, 0.50015594405195891170d+00, $
            0.32951760674056020938d+00, 0.11251791947458271714d+00, 0.19597524022502758711d+00, $
            -0.52856178098938624287d-01, 0.95266905942729742662d+00, 0.18279298809137124237d+01, $
            0.35144044692591625000d+01]
  r[3,*] = [0.22860180522694295568d-01, 0.25001972739726898709d+00, 0.46402414947596950511d+00, $
            0.44398927627952060603d+00, 0.50272011857581864191d+00, 0.88445628872250843244d+00, $
            0.17389166922733014786d+01, 0.34834653146383227628d+01, 0.49049989399817732760d-01, $
            0.85147557568237122183d-01, 0.23424821227315173466d+00, 0.14250169643289223309d+00, $
            0.40586458823694400166d+00, 0.87361464455904407344d+00, 0.17279679104606691097d+01, $
            0.34828619038557646625d+01]
  r[4,*] = [0.48558255088791124620d+00, -0.60422350819111194653d+00,-0.13683535645881987896d+01, $
            -0.85855064913451286656d+00, 0.28100294887963377377d+00, 0.84995521582071482669d+00, $
            0.17002422172314077819d+01, 0.34601959118742100507d+01, -0.60269701842058500673d+00, $
            -0.31233918670609019940d+00, -0.71139286069860601102d-01, 0.13800610282920311533d+00, $
            0.38218789451073229557d+00, 0.81980250822178710734d+00, 0.17183601833584307705d+01, $
            0.34710558209881763325d+01]
  r[5,*] = [0.38821262813796547419d-01, 0.35983028187700409894d+00, 0.57018233005798357737d+00, $
            0.51153414239116576922d+00, 0.55367851931772325002d+00, 0.92841885168264770555d+00, $
            0.17036409347489702703d+01, 0.34382652552776309384d+01, 0.54101616716924096905d-01, $
            0.17168542829213128797d+00, 0.10028211367824637623d+00, 0.27153034480893150082d+00, $
            0.43311645903980116045d+00, 0.75314774033930556029d+00, 0.17038798689307688150d+01, $
            0.34370447884477832722d+01]

  pvec = (rvec = dblarr(16))
  ain = (cin = dblarr(M,M))
  CS = [1.068627d, 1.572350d, 1.015095d, 1.714487d, -5.232800d, 6.528898d]/0.002d

  Bx = 0d
  By = 0d
  Bz = 0d

  Bxout = 0d
  Byout = 0d
  Bzout = 0d

  xp = [x + 0.001d, x, x]
  yp = [y, y + 0.001d, y]
  zp = [z, z, z + 0.001d]
  xm = [x - 0.001d, x, x]
  ym = [y, y - 0.001d, y]
  zm = [z, z, z - 0.001d]

  for Mode = 0, 5 do begin
     ain = reform(aMode[Mode,*,*])
     cin = reform(cMode[Mode,*,*])
     p_1 = reform(p[Mode,0:7])
     p_2 = reform(p[Mode,8:15])
     r_1 = reform(r[Mode,0:7])
     r_2 = reform(r[Mode,8:15])

     Up = dblarr(3)
     Um = Up
     
     term_p = (term_r = dblarr(8,8))
     for i = 0, 7 do begin
        term_p[i,0:7] = 1d/(p_1[i]*p_1[i]) + 1d/(p_2[0:7]*p_2[0:7])
        term_r[i,0:7] = 1d/(r_1[i]*r_1[i]) + 1d/(r_2[0:7]*r_2[0:7])
     endfor
     term_p = exp(sqrt(term_p))
     term_r = exp(sqrt(term_r))

     for i = 0, 7 do begin
        cos_yp = cos(yp/p_1[i])
        cos_ym = cos(ym/p_1[i])
        sin_yp = sin(yp/r_1[i])
        sin_ym = sin(ym/r_1[i])
        for k = 0, 7 do begin
           Up = Up + ain[i,k]*(term_p[i,k]^xp)*cos_yp*sin(zp/p_2[k]) + cin[i,k]*(term_r[i,k]^xp)*sin_yp*sin(zp/r_2[k])
           Um = Um + ain[i,k]*(term_p[i,k]^xm)*cos_ym*sin(zm/p_2[k]) + cin[i,k]*(term_r[i,k]^xm)*sin_ym*sin(zm/r_2[k])
        endfor
     endfor
     B = Up - Um
     Bxout += B[0]*CS[Mode]
     Byout += B[1]*CS[Mode]
     Bzout += B[2]*CS[Mode]
  endfor
end

;--------------------------------------------------------------------------------------
; procedure tail_mag_shield_cyl_S3
;--------------------------------------------------------------------------------------
pro tail_mag_shield_cyl_S3, ctime, rmap, pmap, zmap, Brcss, Bpcss, Bzcss, M, sphiOut, Mode
  ; transform to cartesian
  cyl2car_pos, rmap, pmap, zmap, Xcar, Ycar, Zcar

  S32PJSO, xJSO, yJSO, zJSO, Xcar, Ycar, Zcar, sphiOut

  ; Divide position by 100 to match the coefficients
  ; checkIfInsideMagnetop, xJSO, yJSO, zJSO, ianswer
  xJSO = xJSO/100d
  yJSO = yJSO/100d
  zJSO = zJSO/100d

  ; works in the original coordinate system S3 ***
  B_tail_shield, xJSO, yJSO, zJSO, BxJSO, ByJSO, BzJSO, M, Mode

  PJSO2S3, BxJSO, ByJSO, BzJSO, BxS3, ByS3, BzS3, sphiOut

  ; transform back to cylindrical
  car2cyl_vect, Brcss, Bpcss, Bzcss, BxS3, ByS3, BzS3, pmap
end   ; tail_mag_shield_cyl_S3

;--------------------------------------------------------------------------------------
; procedure csheet_struc
;
; NEW VERSION NOV 2003
; This program calculates the distance of current sheet from the 
; Jovigraphic equator in system III coordinates
;--------------------------------------------------------------------------------------
pro csheet_struc, ZNS3, rho, phi, XJSO, YJSO, stheta, ctime ; ctime?
  radian = 0.01745329252d
  period = 9.927953d
  OMEGAJ = 36.26125d
  X0 = -45d
  phip0 = 6.12611d
  ;ten = tan(0.167551606d)
  ten = 0.16913733749783d
  C1 = 0.005973d
  C2 = 5.114d-5
  c3 = 1.59d-5
  psi2 = -1.201201d
  C4 = 0.313244d
  C5 = -0.366166d
  psi4 = 2.2604522d

  if (abs(XJSO) lt 1d-6) then XJSO = 1d-6

  RLT = atan(YJSO, XJSO) ; calculate RLT here
  
  delay1 = C1*rho + 0.5d*C2*(rho*rho)*cos(RLT - psi2) + 0.5d*C3*rho*rho
  Alfven = -360d/( period*(C4*cos(RLT - psi4) + C5) )
  delay2 = rho/Alfven*Omegaj*radian

  phip = phip0 - delay1 - delay2

  hyptan = x0*tanh(XJSO/x0)
  rho1 = sqrt(hyptan*hyptan + YJSO*YJSO)
  
  ZNS3 = rho1*ten*cos(phi - phip) + XJSO*(1d - tanh( abs(x0/xjso) ))*tan(stheta)
end

;--------------------------------------------------------------------------------------
; procedure csheet_N
;
; This procedure uses csheet_structure to calculate the normal 
;   direction to the current sheet in system III coordinates
;--------------------------------------------------------------------------------------
pro csheet_N, RNx, RNy, RNz, rho, phi, zs3, stheta, ctime
  vecin = (vecou = dblarr(3))
  delta = 0.1d
  delta2 = 2d*delta
  xs3 = rho*cos(phi)
  ys3 = rho*sin(phi)
  xp = xs3 + delta
  xm = xs3 - delta
  yp = ys3 + delta
  ym = ys3 - delta

  ; First calculate the x derivatives
  ; Define the positive parameters
  rhop = sqrt(xp*xp + ys3*ys3)
  phip = atan(ys3, xp)
  
  vecin = [xp, ys3, zs3]
  JROT, 's3c', 'jso', vecin, vecou, ctime
  XJSOP = vecou[0]
  YJSOP = vecou[1]
  csheet_struc, ZNS3P, rhop, phip, XJSOP, YJSOP, stheta, ctime
  ; Now define the negative parameters
  rhom = sqrt(xm*xm + ys3*ys3)
  phim = atan(ys3,xm)
  
  vecin = [xm, ys3, zs3]
  JROT, 's3c', 'jso', vecin, vecou, ctime
  XJSOM = vecou[0]
  YJSOM = vecou[1]

  csheet_struc, ZNS3M, rhom, phim, XJSOM, YJSOM, stheta, ctime
  dzdx = (ZNS3P - ZNS3M)/delta2

  ; Next the y derivative
  ; Define the positive parameters
  rhop = sqrt(xs3*xs3 + yp*yp)
  phip = atan(yp, xs3)
  
  vecin = [xs3, yp, zs3]
  JROT, 's3c', 'jso', vecin, vecou, ctime
  XJSOP = vecou[0]
  YJSOP = vecou[1]

  csheet_struc, ZNS3P, rhop, phip, XJSOP, YJSOP, stheta, ctime
  ; Now define the negative parameters
  rhom = sqrt(xs3*xs3 + ym*ym)
  phim = atan(ym, xs3)
  
  vecin = [xs3, ym, zs3]
  JROT, 's3c', 'jso', vecin, vecou, ctime
  XJSOM = vecou[0]
  YJSOM = vecou[1]

  csheet_struc, ZNS3M, rhom, phim, XJSOM, YJSOM, stheta, ctime
  dzdy = (ZNS3P - ZNS3M)/delta2
  RNx = -dzdx
  RNy = -dzdy
  RNz = 1d
  RN = sqrt(RNx*RNx + RNy*RNy + RNz*RNz)
  RNx = RNx/RN
  RNy = RNy/RN
  RNz = RNz/RN
end
 
;--------------------------------------------------------------------------------------
; procedure get_mapped_sunangle
;--------------------------------------------------------------------------------------
pro get_mapped_sunangle, sthetaIN, sphiIN, sthetaOut, sphiOut, $
                         xpx, xpy, xpz, ypx, ypy, ypz, zpx, zpy, zpz

  SPH2CAR_pos, Xin, Yin, Zin, 1d, sthetaIN, sphiIN

  Xout = Xin*xpx + Yin*xpy + Zin*xpz
  Yout = Xin*ypx + Yin*ypy + Zin*ypz
  Zout = Xin*zpx + Yin*zpy + Zin*zpz

  CAR2SPH_pos, Xout, Yout, Zout, Rnotused, sthetaOut, sphiOut
end

;--------------------------------------------------------------------------------------
; procedure getBimf
; 
; Modified by Fran Bagenal Nov. 2008 - input background field
;   Br=Bx = +-0.2, Bp=By = -+1.0, Bt=Bz = 0.0 - in nT - OR Zero...
;--------------------------------------------------------------------------------------
pro getBimf, ctime, thetaS3, phiS3, BrS3IMF, BtS3IMF, BpS3IMF
  vecin = (vecou = dblarr(3))
  BxIMFgsm = 0d
  ByIMFgsm = 0d
  BzIMFgsm = 0d
      
  ; go from JSM to S3
  vecin = [BxIMFgsm, ByIMFgsm, BzIMFgsm]
  JROT, "jsm", "s3c", vecin, vecou, ctime
  BxIMFs3 = vecou[0]
  ByIMFs3 = vecou[1]
  BzIMFs3 = vecou[2]
      
  ; go to spherical
  CAR2SPH_MAG, BxIMFs3, ByIMFs3, BzIMFs3, BrS3IMF, BtS3IMF, BpS3IMF, thetaS3, phiS3
end

;--------------------------------------------------------------------------------------
; procedure CheckIfInsideMappedMp
;--------------------------------------------------------------------------------------
pro CheckIfInsideMappedMp, ctime, XS3, YS3, ZS3, ZNS3, ianswer
  vecin = (vecou = dblarr(3))

  vecin = [XS3, YS3, ZS3 + ZNS3]
  JROT, 's3c', 'jsm', vecin, vecou, ctime
  XJSOold = vecou[0]
  YJSOold = vecou[1]
  ZJSOold = vecou[2]
;print, ZJSOold
;
  checkIfInsideMagnetop, XJSOold, YJSOold, ZJSOold, ianswer
end

;--------------------------------------------------------------------------------------
; procedure checkIfInsideMagnetop
; 
; created by Mariel Desroche 11-17-2008
; checks whether the position x0, y0, z0 is in the magnetopause
;   defined by Joy et al 2002 polynomial
; Pd is the solar wind dynamic pressure (nPa)
; returns 0 = outside mp, 1 = inside mp
;--------------------------------------------------------------------------------------
;pro checkIfInsideMagnetop, x, y, z, ianswer
;  ; initialize answer 
;  ianswer = 0

;  Pd = 0.04d
;  Pd4 = 1d/sqrt(sqrt(Pd))
;  ; Magnetopause Parameters
;  A = -0.134d + 0.488d*Pd4
;  B = -0.581d - 0.225d*Pd4
;  C = -0.186d - 0.016d*Pd4
;  D = -0.014d + 0.096d*Pd
;  E = -0.814d - 0.811d*Pd
;  F = -0.050d + 0.186d*Pd

;  x0 = x/120d
;  y0 = y/120d
;  z0 = z/120d

;  aa = -0.186d - 0.016d*Pd4
;  bb = -0.581d - 0.225d*Pd4
;  cc = -0.134d + 0.488d*Pd4 - z0*z0 

;  x_max = (-bb - sqrt(bb*bb - 4*aa*cc))/(2*aa)
       
;  aa = -0.814d - 0.811d*Pd
;  bb = -0.014d + 0.096d*Pd + (-0.050d + 0.186d*Pd)*x0
;  cc = (-0.186d - 0.016d*Pd4)*x0*x0 + (B = -0.581d - 0.225d*Pd4)*x0 - 0.134d + 0.488d*Pd4 - z0*z0

;  quad = bb*bb - 4*aa*cc
;  y_pos = (-bb + sqrt(quad))/(2*aa)
;  y_neg = (-bb - sqrt(quad))/(2*aa)

;  if (x0 lt x_max) && (y0 gt y_pos) && (y0 lt y_neg) then ianswer = 1
;end

;--------------------------------------------------------------------------------------
; procedure checkIfInsideMagnetop
;
; original magnetopause check program... still being tested and under construction
;--------------------------------------------------------------------------------------
pro checkIfInsideMagnetop, x, y, z, answer
  ; init answer
  answer = 1

  ; Symmetric Magnetopause Parameter
  A = 13779.1546267827d
  B = -130.0469911647737d
  C = -0.2216972473764874d
  D = 0d
  E = -0.8453162863101464d
  F = 0d
    
  phi = atan(z,y)
  rho_in = sqrt(z*z + y*y)
  cos_phi = cos(phi)
  sin_phi = sin(phi)
  
  ; check if we are on the dayside
  if (x ge 0d) then begin
     aa = E*cos_phi*cos_phi - sin_phi*sin_phi + F*cos_phi*sin_phi
     bb = D*cos_phi
     cc = A + B*x + c*x*x
     term = bb*bb - 4*aa*cc
     sqrt_term = sqrt(term)
     rho1 = (-bb + sqrt_term)/(2d*aa)
     rho2 = (-bb - sqrt_term)/(2d*aa)
     ; The right solution has rho positive (i.e. the bigger of the two)
     rho = rho1
     if (rho2 gt rho1) then rho = rho2
     if (term lt 0d) || (rho lt 0d) || (rho_in gt rho) then answer = 0
  endif else begin
     ; Now the nightside 
     aa = E*cos_phi*cos_phi - sin_phi*sin_phi + F*cos_phi*sin_phi
     bb = D*cos_phi
     cc = A + B*x + c*x*x
     term = bb*bb - 4*aa*cc
     sqrt_term = sqrt(term)
     rho1 = (-bb + sqrt_term)/(2d*aa)
     rho2 = (-bb - sqrt_term)/(2d*aa)
     ; The right solution has rho positive (i.e. the bigger of the two)
     rho = rho1
     if (rho2 gt rho1) then rho = rho2
     if (term lt 0d) || (rho lt 0d) || (rho_in gt rho) then answer = 0
  endelse
end

;--------------------------------------------------------------------------------------
; procedure tail_mag_notilt_allmodes
;
; Modified from csheet_deform to calculate the field of only one point
; The subroutine first calculates the Brho and Bz components at the stretched location 
;   and then multiplies with the stretch matrices to calculate Bphi
; alat,colat RLT are in radians, RLT is measured from Noon
;--------------------------------------------------------------------------------------
pro tail_mag_notilt_all_modes, rho, phi, zz, RLT, Brcs, Bpcs, Bzcs, Mode
  f = (beta = dblarr(6,6))
  X = dblarr(6)
  isumno = 6

  ; begin Hannes Modes
  C = [-28.402612d, -18.160631d, 14.543719d, 12.873892d, -5.590983d, 3.461748d] ; short model

  ; ring current
  f[0,*] = [188.880535100774d, 170.216250249208d, -232.209387453753d, 17729.2564306361d, -17426.2836367108d, -176.740242190339d]
  beta[0,*] = [59.2554626464844d, 14.3016757965088d, 3.41575503349304d, 5.47725486755371d, 5.55856323242188d, 30.7563076019287d]
  ; L = 6
  f[1,*] = [34.5422110518d, 6.0536537495d, 414.2147150061d, -155.3647545840d, 316.3317439470d, -73.2913775361d]
  beta[1,*] = [59.2781028748d, 16.4449157715d, 4.1893663406d, 3.0382030010d, 9.1110754013d, 30.9898948669d]
  ; L = 12
  f[2,*] = [-165.1020233300d, 554.5217109532d, 315.5180924885d, -125.0952135958d, 449.1102649649d, 371.9172713735d]
  beta[2,*] = [55.9542961121d, 14.5898332596d, 3.9975905418d, 2.8588902950d, 8.6061649323d, 21.9151973724d]
  ; L = 24
  f[3,*] = [591.9706482108d, 1376.0166168630d, 298.8401265534d, -156.9445820815d, 521.9904943354d, 2387.3415724416d]
  beta[3,*] = [59.0788803101d, 16.2336521149d, 3.5862164497d, 2.7877552509d, 8.4433155060d, 30.6951904297d]
  ; L = 48
  f[4,*] = [9180.3492664995d, 1679.4592666407d, 274.9968669893d, -144.4824637207d, 532.5683191252d, 4540.5253277291d]
  beta[4,*] = [59.2788314819d, 16.4751930237d, 3.5809681416d, 2.7559804916d, 8.5185470581d, 31.0038948059d]
  ; L = 96
  f[5,*] = [28927.9037953983d, 2214.0596352902d, 242.5405186300d, -108.6915959578d, 609.8932996003d, 7437.5789048160d]
  beta[5,*] = [72.2351608276d, 18.0498085022d, 3.7080323696d, 2.6697711945d, 8.9703912735d, 35.2125587463d]
  ; end Hannes Modes

  ;! P = 3.0d-3
  ;! P = 0d
  P = 5.73d-3

  if (Mode gt 6) || (Mode lt 1) then begin
     istartLoop = 0
     istopLoop = 5
  endif else begin
     istartLoop = Mode - 1
     istopLoop = istartLoop
  endelse

  ; We work in the dipole magnetic coordinate system
  ; In this subroutine:
  ;   index I = 1,7 corresponds to each component of a mode
  ;   index L = 1,6 denotes the six L modes
  ;   index j sums over observations (Br and Btheta) and is twice the number of observations
  ;   index M is used to obtain least squares equations (same dimension as L)
  PI = 3.14159265358979d
  PI2 = PI/2d
  DEGREE = 180d/PI
  RADIAN = PI/180d
  drho = 0.05d
  dz = 0.05d
  phi = phi
  ; D is half thickness of current sheet
  ;! D = 2.0
  D0 = 11.0d
  D1 = 9.0d
  ;D = D0 + D1*cos(RLT - 1.5d*PI)
  D = 4.0d

  RMSBr = 0d
  RMSBt = 0d
  RMSBp = 0d
  ; Start calculating the field
  Brhodm = 0d
  Bzdm = 0d
  Bpdm = 0d

  rhomag = rho
  ZMAG = zz
  Z = ZMAG


  count = 0

  ;problem not present

  ; First calculate Brhod
  ;! Do 2 L = mode,mode
  for L = istartLoop, istopLoop do begin    
     
     ;THESE IF STATEMENT ARE THE PROBLEM!!!!
     ZM = abs(Z - dz)
     if (ZM lt D) then begin
      ZM = 0.5*(ZM*ZM/D + D)
      count++
     ENDIF
     
     ZP = abs(Z + dz)
     if (ZP lt D) then ZP = 0.5*(ZP*ZP/D + D)
     
     xlpp = 0d
     xlpm = 0d
     for i = 0, isumno - 1 do begin
        S1p = sqrt( (beta(L,i) + ZP)*(beta(L,i) + ZP) + (RHOMAG + Beta(L,i))*(RHOMAG + Beta(L,i)) )
        S2p = sqrt( (beta(L,i) + ZP)*(beta(L,i) + ZP) + (RHOMAG - Beta(L,i))*(RHOMAG - Beta(L,i)) )
        S1m = sqrt( (beta(L,i) + ZM)*(beta(L,i) + ZM) + (RHOMAG + Beta(L,i))*(RHOMAG + Beta(L,i)) )
        S2m = sqrt( (beta(L,i) + ZM)*(beta(L,i) + ZM) + (RHOMAG - Beta(L,i))*(RHOMAG - Beta(L,i)) )
        tp = 2d*Beta(L,i)/(S1p + S2p)
        tm = 2d*Beta(L,i)/(S1m + S2m)
        AAp = tp*sqrt(1d - tp*tp)/(S1p*S2p)
        AAm = tm*sqrt(1d - tm*tm)/(S1m*S2m)
        xlpp = xlpp + f(L,i)*AAp*rhomag
        xlpm = xlpm + f(L,i)*AAm*rhomag
     endfor
     dxpldz = (xlpp - xlpm)/(2d*dz)
     X[L] = -dxpldz
     Brhodm = Brhodm + X[L]*C[L]
     
     
  endfor
  
  ;STOP
  ;problem in Brhodm

  
  defsysv, "!Brm", count
  ;defsysv, "!Btm", RLT
  ;defsysv, "!Bpm", phi

  ; Now calculate BZD
  ; Now corrected for the current sheet derivative term. 11/2002.
  ;! Do 20 L = Mode,Mode
  for L = istartLoop, istopLoop do begin
     rhom = RHOMAG - drho
     rhop = RHOMAG + drho

     ; First calculate it for rhom
     rhos3m = rho - drho
     ZNM = 0d
     ZM = abs(ZMAG - ZNM)
     if (ZM lt D) then ZM = 0.5d*(ZM*ZM/D + D)

     ; Now calculate it for rhop
     rhos3p = rho + drho
     ZNP = 0d
     ZP = abs(ZMAG - ZNP)
     if (ZP lt D) then ZP = 0.5d*(ZP*ZP/D + D)
     xlpp = 0d
     xlpm = 0d
     for i = 0, isumno - 1 do begin
        S1p = sqrt( (beta(L,i) + ZP)*(beta(L,i) + ZP) + (rhop + Beta(L,i))*(rhop + Beta(L,i)) )
        S2p = sqrt( (beta(L,i) + ZP)*(beta(L,i) + ZP) + (rhop - Beta(L,i))*(rhop - Beta(L,i)) )
        S1m = sqrt( (beta(L,i) + ZM)*(beta(L,i) + ZM) + (rhom + Beta(L,i))*(rhom + Beta(L,i)) )
        S2m = sqrt( (beta(L,i) + ZM)*(beta(L,i) + ZM) + (rhom - Beta(L,i))*(rhom - Beta(L,i)) )
        tp = 2d*Beta(L,i)/(S1p + S2p)
        tm = 2d*Beta(L,i)/(S1m + S2m)
        AAp = tp*sqrt(1d - tp*tp)/(S1p*S2p)
        AAm = tm*sqrt(1d - tm*tm)/(S1m*S2m)
        xlpp = xlpp + f(L,i)*AAp*rhop
        xlpm = xlpm + f(L,i)*AAm*rhom
     endfor
     dxpldr = (rhop*xlpp - rhom*xlpm)/(2d*drho)
     X[L] = dxpldr/RHOMAG
     Bzdm = Bzdm + X[L]*C[L]
  endfor
  Bpdm = 0d
  ; Now add the stretch field to Bpdm 
  Bpdm = Bpdm - p*rhomag*Brhodm
  ;! PH_MAG=PH_MAG-p*rhomag

  ; Now rotate into system III coordinates.
  ;! Bphi_N=0  !changed Hannes
  Brcs = Brhodm
  Bpcs = Bpdm
  Bzcs = Bzdm
  ;! Bpdifm=BPHM
  ;! Btdifm=BTHM
end

;--------------------------------------------------------------------------------------
; procedure JOVIAN_VIP4_no_dipole
; 
; created by K.K.Khurana 1996
; BASED ON THE SUBROUTINE IGRF WRITTEN BY N.A. TSYGANENKO (1979)
; THE LEFT HANDED COODINATE SYSTEM SYSTEM III IS CHANGED TO A RIGHT
; HANDED SYSTEM BY FEEDING IN  -F.
;
; INPUT:  NM (INTEGER)- MAXIMUM ORDER OF HARMONICS TAKEN
;         INTO ACCOUNT (NM.LE.12)
;
;         R,T,F (REAL)- POSITION OF DESIRED FIELD VALUE IN
;         SPHERICAL GEOGRAPHIC COORDINATE SYSTEM
;         (R IN EARTH RADII, COLATITUDE T AND LONGITUDE F IN RADIANS)
;
; OUTPUT: BR,BT,BF (REAL)- COMPONENTS OF THE INTERNAL PORTION
;         OF THE MAIN MAGNETIC FIELD IN
;         SPHERICAL GEOGRAPHICAL COORD SYSTEM (VALUES GIVEN IN GAMMA)
;
; CALCULATES COMPONENTS OF MAIN JOVIAN FIELD IN SPHERICAL SYSTEM III (1965)
;   COORD SYSTEM BASING ON SPHERICAL HARMONIC COEFFICIENTS GIVEN BY ACUNA
;   ET.AL. (1983) IN DESSLER'S BOOK.
; Now changed to O6 model. April, 1996.
; MAXIMUM ORDER OF HARMONICS TAKEN INTO ACCOUNT (NOT MORE THAN 03)
; R,T,F ARE SPHERICAL COORDS OF THE POINT (R IN RJ,COLATITUDE T AND
;   LONGITUDE F IN RADIANS), FIELD COMPONENTS BR,BT,BF
;   ARE CALCULATED IN NANOTESLA(= 1 GAMMA)
;--------------------------------------------------------------------------------------
pro JOVIAN_VIP4_no_dipole, NM, r, theta, phi, BR, BT, BF  
  REC = replicate(1d, 91)
  A = (B = dblarr(13))
  d76 = dblarr(76)

  ; VIP4 no dipole component
  G = [0d, 0d, 0d, -5100d, -61900d, 49700d, -1600d, -52000d, 24400d, $
       -17600d, -16800d, 22200d, -6100d, -20200d, 6600d, d76]
  H = [0d, 0d, 0d, 0d, -36100d, 5300d, 0d, -8800d, 40800d, -31600d, $
       0d, 7600d, 40400d, -16600d, 3900d, d76]  

  ;! VIP4
  ;! G = [0.,420500.,-65900.,-5100.,-61900.,49700.,-1600.,-52000.,24400., $
  ;!      -17600.,-16800.,22200.,-6100.,-20200.,6600.,76*0.0]
  ;! H = [0.,0.,25000.,0.,-36100.,5300.,0.,-8800.,40800.,-31600.,0.,7600., $
  ;!      40400.,-16600.,3900,76*0.0]
  ;! O6
  ;! G = [0.,421800.,-66400.,-20300.,-73500.,51300.,-23300.,-7600.,16800.,-23100.,81*0.0]
  ;! H = [0.,0.,26400.,0.,-46900.,8800.,0.,-58000.,48700.,-29400.,81*0.0]

  G[0] = 0d
  H[0] = 0d
  KNM = 15
  for N = 1, 13 do begin
     N2 = (2d*N - 1d)*(2d*N - 3d)
     for M = 1, N do begin
        MN = N*(N - 1)/2 + M
        REC[MN-1] = (N - M)*(N + M - 2d)/N2
     endfor
  endfor
  S = 1d
  for N = 2, 13 do begin
     MN = N*(N - 1d)/2d
     S = S*(2d*N - 3d)/(N - 1d)
     G[MN] = G[MN]*S
     H[MN] = H[MN]*S
     P = S
     for M = 2, N do begin
        AA = 1d
        if (M eq 2) then AA = 2d
        P = P*sqrt(AA*(N - M + 1d)/(N + M - 2d))
        MNN = MN + M - 1
        G[MNN] = G[MNN]*P
        H[MNN] = H[MNN]*P
     endfor
  endfor
  if (KNM ne NM) then begin
     KNM = NM
     K = KNM + 1
  endif
  PP = 1d/r
  P = PP
  for N = 1, K do begin
     P = P*PP
     A[N-1] = P
     B[N-1] = P*N
  endfor
  P = 1d
  D = 0d
  BBR = 0d
  BBT = 0d
  BBF = 0d
  cos_phi = cos(phi)
  sin_phi = sin(phi)
  cos_theta = cos(theta)
  sin_theta = sin(theta)
  for M = 1, K do begin
     if (M eq 1) then begin
        X = 0d
        Y = 1d
     endif else begin
        MM = M - 1
        W = X
        X = W*cos_phi + Y*sin_phi
        Y = Y*cos_phi - W*sin_phi
     endelse
     Q = P
     Z = D
     BI = 0d
     P2 = 0d
     D2 = 0d
     for N = M, K do begin
        ;AN = A[N - 1]
        MN = N*(N - 1)/2 + M
        ;E = G[MN - 1]
        ;HH = H[MN - 1]
        W = G[MN - 1]*Y + H[MN - 1]*X
        ;if (ABS(P2) lt 1d-38) then P2 = 0d
        ;if (ABS(Q) lt 1d-38) then  Q = 0d
        BBR = BBR + B[N-1]*W*Q
        BBT = BBT - A[N - 1]*W*Z
        if (M ne 1) then begin
           QQ = Q
           if (sin_theta lt 1d-5) then QQ = Z
           BI = BI + A[N - 1]*(G[MN - 1]*X - H[MN - 1]*Y)*QQ
        endif
        XK = REC[MN-1]
        DP = cos_theta*Z - sin_theta*Q - XK*D2
        PM = cos_theta*Q - XK*P2
        D2 = Z
        P2 = Q
        Z = DP
        Q = PM 
     endfor
     D = sin_theta*D + cos_theta*P
     P = sin_theta*P
     if (M ne 1) then begin
        BI = BI*MM
        BBF = BBF + BI
     endif
  endfor
  BR = BBR
  BT = BBT
  if (sin_theta gt 1d-5) then begin
     BF = BBF/sin_theta
  endif else begin
     if (cos_theta lt 0d) then BBF = -BBF
     BF = BBF
  endelse
end


