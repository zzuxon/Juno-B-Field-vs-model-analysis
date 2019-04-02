FUNCTION fgm_B_spherical, file
  ON_ERROR,2

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ; temporary playing with file
  IF N_ELEMENTS(file) EQ 0 THEN MESSAGE,'ERROR: Must pass on a file as an argument to this function'



  IF FILE_TEST(file) EQ 0 THEN MESSAGE,'ERROR: input file does not exist: '+file

  IF STREGEX(file, '\.gz$', /BOOLEAN) THEN COMPRESS = 1 ELSE COMPRESS = 0

  ; pc = Planetocentric coordinates = IAU_JUPITER
  IF STREGEX(file,'fgm_jno_(ql|l3)_[0-9]{7}pc(_r1s|_r60s)?_v[0-9]{2}\.sts', /BOOLEAN) EQ 0 THEN MESSAGE, 'Code ONLY takes planetocentric (pc) files'

  line = ''
  nlines = FILE_LINES(file,COMPRESS=COMPRESS)
  OPENR, lun, file, /GET_LUN ,COMPRESS = COMPRESS
  array = DBLARR(nlines, 14)
  cnt = -1L
  WHILE NOT EOF(lun) DO BEGIN
    READF, lun, line
    IF STREGEX(line, '[A-Za-z]', /BOOLEAN) THEN CONTINUE
    ; Must be all numbers, so a record, carry on
    cnt ++
    l = STRSPLIT(line, /EXTRACT)
    array[cnt,*] = DOUBLE(l)
  ENDWHILE
  FREE_LUN, lun
  array = array[0L:cnt,*]
  n_rec_m1 = cnt
  n_recs = cnt + 1L

  ; Columns of file are
  ; Year, DOY, HH, MM, SS, milli-s, Decimal day , Bx, By, Bz, Range, PosX, PosY, PosZ
  ;    0    1   2   3   4        5            6    7   8   9     10    11    12    13

  ; So far these work for all tpes
  ;T_MAG = JULDAY(12, 31, array[*,0]-1, 0, 0, 0) + array[*,6]
  T_MAG = JULDAY(12, 31, array[*,0]-1, ROUND(array[*,2]), ROUND(array[*,3]), ROUND(array[*,4])) + array[*,1] + array[*,5]/86400000d
  MAG_JSSxyz  = array[*,7:9]
  range = array[*,10]
  Pos_JSSxyz = array[*,11:13]
  
  ;UTC_MAG = Julian_to_time_struct(T_MAG, /utc) & UTC = UTC_MAG.UTC1
  ; Originally used above, but below is faster
  UTC = STRING(ROUND(array[*,0]),FORMAT = '(I4)')+'-'+$
    STRING(ROUND(array[*,1]),FORMAT = '(I03)')+'T'+$
    STRING(ROUND(array[*,2]),FORMAT = '(I02)')+':'+$
    STRING(ROUND(array[*,3]),FORMAT = '(I02)')+':'+$
    STRING(ROUND(array[*,4]),FORMAT = '(I02)')+'.'+$
    STRING(ROUND(array[*,5]),FORMAT = '(I03)')


  Pos_JSS_r = SQRT(Pos_JSSxyz[*,0]*Pos_JSSxyz[*,0] + Pos_JSSxyz[*,1]*Pos_JSSxyz[*,1] + Pos_JSSxyz[*,2]*Pos_JSSxyz[*,2])

  CASE 4 OF
    4 : BEGIN
      ; this is my case where I don't bother with the phi_unit as a variable... Rob knows what's going on...
      r_unit = [$
        [Pos_JSSxyz[*,0] / Pos_JSS_r],$
        [Pos_JSSxyz[*,1] / Pos_JSS_r],$
        [Pos_JSSxyz[*,2] / Pos_JSS_r] ]
      theta_unit = [$
        [ r_unit[*,0]*r_unit[*,2]], $
        [ r_unit[*,2]*r_unit[*,1]], $
        [-r_unit[*,1]*r_unit[*,1]-r_unit[*,0]*r_unit[*,0]] ]
      MAG_JSSrtp = [$
        [MAG_JSSxyz[*,0]*r_unit[    *,0] + MAG_JSSxyz[*,1]*r_unit[    *,1] + MAG_JSSxyz[*,2]*r_unit[    *,2]],$
        [MAG_JSSxyz[*,0]*theta_unit[*,0] + MAG_JSSxyz[*,1]*theta_unit[*,1] + MAG_JSSxyz[*,2]*theta_unit[*,2]],$
        [MAG_JSSxyz[*,1]*r_unit[*,0] - MAG_JSSxyz[*,0]*r_unit[*,1]] ]
    END
  ENDCASE

  IF 0 THEN BEGIN
    DUMMY = LABEL_DATE(DATE_FORMAT=['%H:%I:%S'])
    h1 = plot(T_MAG,MAG_JSSrtp[*,0],'r',CURRENT=0,OVERPLOT=0,NAME = 'Br')
    h2 = plot(T_MAG,MAG_JSSrtp[*,1],'g',CURRENT=1,OVERPLOT=1,NAME = 'Bth')
    h3 = plot(T_MAG,MAG_JSSrtp[*,2],'b',CURRENT=1,OVERPLOT=1,NAME = 'Bphi')
    h4 = plot(T_MAG,SQRT(TOTAL(MAG_JSSrtp^2,2)),'k',/CURRENT,OVERPLOT=1,NAME = 'Bmag')
    h1.XTICKFORMAT = 'LABEL_DATE'
    h1.TITLE = 'JSS rtp'
    leg = LEGEND(TARGET=[h1,h2,h3,h4], /AUTO_TEXT_COLOR)
  ENDIF
  
  RETURN, CREATE_STRUCT('T',T_MAG,'UTC',UTC,'MAG_RANGE',range,'MAG_VECTOR_JSSRTP',MAG_JSSrtp)
END
