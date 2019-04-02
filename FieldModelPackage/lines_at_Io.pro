

;This code will use Drake's code to produce 360 field lines, 1 for each integer longitude of Io.
;it will make 360 field lines if FOR longitude_counter=0, 359 DO BEGIN is set
;for 3,600 lines use FOR longitude_counter=1, 3600 DO BEGIN (edit made 8/21/18)
;Made by Parker Hinton 4/24


parkers_lines = []

FOR longitude_counter=0.0d, 360d, .1d DO BEGIN
  fp1 = footpoint(6, 90, longitude_counter, 'JRM09+CAN',0.01,/FULL) ;change JRM09+CAN to CANSheet, this is VIP4, or vice versa as needed
  fp2 = footpoint(6, 90, longitude_counter, 'JRM09+CAN',0.01,/REVDIR,/FULL);change JRM09+CAN to CANSheet, this is VIP4, or vice versa as needed
parkers_lines = [parkers_lines,fp1,fp2]
print, longitude_counter
ENDFOR


parkers_lines = TRANSPOSE(parkers_lines)
write_csv, '/Users/pahi9557/Desktop/Field_lines_JRM09.csv', parkers_lines, header = ["Radius","Theta","Phi"]

fillervariable = 0
lines_at_Io_magB, fillervariable

END
