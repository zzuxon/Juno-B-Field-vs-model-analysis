

; This code is to be run after lines_at_Io.pro. It takes the outputs from that and generates the magentic field at each location.


CSV_Data = READ_CSV("/Users/cowo5730/Documents/juno orbit 11/039rtp.csv", HEADER = CSVHeader)
; Input the file you created with position extraction.
r_pH =csv_Data.(0)
theta_pH =csv_Data.(1)
phi_pH =csv_Data.(2)
b = fgm_B_spherical('/Users/cowo5730/Documents/juno orbit 11/fgm_jno_l3_2018039pc_r1s_v01.sts')
;Input the original file. This program handles the complexities of converting magnetic field components
;into rtp.


parkers_lines_magB = []
parkers_lines_rB = []
parkers_lines_thetaB = []
parkers_lines_phiB = []
nel = n_elements(r_pH)
for i = 0, nel-1 DO BEGIN

  Bfield = jovMagField(r_pH(i), theta_pH(i), phi_pH(i), 'JRM09', CURRENT_SHEET='CANSheet')
  ;;Bfieldnosheet = jovMagField(r_pH(i), theta_pH(i), phi_pH(i), 'JRM09')
  Bmag = SQRT(Bfield[0]^2 + Bfield[1]^2 + Bfield[2]^2)
  ;;Bmagnosheet = SQRT(Bfieldnosheet[0]^2 + Bfieldnosheet[1]^2 + Bfieldnosheet[2]^2)
  parkers_lines_magB = [parkers_lines_magB,Bmag]
  parkers_lines_rB = [parkers_lines_rB,Bfield[0]]
  parkers_lines_thetaB = [parkers_lines_thetaB,Bfield[1]]
  parkers_lines_phiB = [parkers_lines_phiB,Bfield[2]] 
  ;This gives not only the magnitude of the model's predicted b-field, but each component.
  print, i ;This allows you to see the progress of the program in the console.
  
 
ENDFOR


parkersArray = [[parkers_lines_magB],[parkers_lines_rB],[parkers_lines_thetaB],[parkers_lines_phiB]]
parkersArray = TRANSPOSE(parkersArray) 
;Due to the poor construction of the idl write_csv command, all components must be put in a single object.


write_csv, '/Users/cowo5730/Documents/juno orbit 11/039JRM09magB.csv', parkersArray;, header = ["magB, rB thetaB, phiB"]
write_csv, '/Users/cowo5730/Documents/juno orbit 11/039magrtp.csv', B.MAG_VECTOR_JSSRTP[*,0], B.MAG_VECTOR_JSSRTP[*,1], B.MAG_VECTOR_JSSRTP[*,2];, header = ["rB, thetaB, phiB"]
;Thus, a file giving the model values and actual values, respectively for the b-field are created

END