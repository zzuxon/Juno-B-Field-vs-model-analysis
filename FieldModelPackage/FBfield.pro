; FB's version - 
; B is in GAUSS (multiply by 1.e5 to get nT) - 
; B(7): Br, Bt, Bp, Bx, By, Bz, Bmag
FUNCTION CAN_SHEET, Pos
        ;Magentic field from a current sheet based on
        ;the Connerney, Acuna and Ness model
        ;Pos(0) is radius, in the same units as the parameter, Rj.
        ;Pos(1) is latitude, _not_ co-latitude in degrees.
        ;Pos(2) is West longitude (i.e. System III longitude)
        ;This coordinate system is right handed, but the
        ;       unit vectors in the latitude and longitude
        ;       directions both point anti-parallel to  the
        ;       traditional theta and phi of spherical coordinates

        r0 = 5.0
        r1 = 50.0
        d = 2.5
        c = 2.25E-3
        xt = 9.6*!PI/180.
        xp = 202.*!PI/180.

        B = DBLARR(7)
        r = Pos(0)
        theta = !PI*(90-Pos(1))/180.
        phi = !PI*(360-Pos(2))/180.

        xp = 2*!PI-xp
        ct=COS(xt)
        st=SIN(xt)
        cp=COS(xp)
        sp=SIN(xp)

        x1 = r*SIN(theta)*cos(phi)
        y1 = r*SIN(theta)*sin(phi)
        z1 = r*COS(theta)

        x2 = cp*ct*x1 + sp*ct*y1 - st*z1
        y2 = -sp*x1 + cp*y1
        z2 = st*cp*x1 + st*sp*y1 + ct*z1

        rr = SQRT(x2^2+y2^2)
        pp = ATAN(y2,x2)
        IF pp LT 0 THEN phi = phi + 2*!PI
        zz = z2
        IF rr LT r0 THEN BEGIN
                F1 = SQRT( (zz-d)^2 + r0^2 )
                F2 = SQRT( (zz+d)^2 + r0^2 )
                F3 = SQRT( zz^2 + r0^2 )
                Br = 0.5*rr*(1/F1 - 1/F2)
                Bz = 2*d/F3 - 0.25*rr^2*( (zz-d)/F1^3 - (zz+d)/F2^3 )
                END ELSE BEGIN
                F1 = SQRT( (zz-d)^2 + rr^2 )
                F2 = SQRT( (zz+d)^2 + rr^2 )
                Br = (F1 - F2 + 2*D)/rr
                IF (ABS(zz) GE D) AND (zz LT 0) THEN Br = (F1 - F2 - 2*D)/rr
                IF (ABS(zz) LT D) THEN Br = (F1 - F2 + 2*zz)/rr
                Br = Br - 0.25*r0^2*rr*(1/F1^3 - 1/F2^3)
                Bz = 2*d/SQRT(zz^2+rr^2) - 0.25*r0^2*( (zz-d)/F1^3 - $
                        (zz+d)/F2^3 )
                END

        F1 = SQRT( (zz-d)^2 + r1^2 )
        F2 = SQRT( (zz+d)^2 + r1^2 )
        F3 = SQRT( zz^2 + r1^2 )
        Br2 = 0.5*rr*(1/F1 - 1/F2)
        Bz2 = 2*d/F3 - 0.25*rr^2*( (zz-d)/F1^3 - (zz+d)/F2^3 )

        Bcy = [Br-Br2,0,Bz-Bz2]
        Bcy(1) = SIN(pp)*Bcy(0)
        Bcy(0) = COS(pp)*Bcy(0)

        B1 = cp*ct*Bcy(0) - sp*Bcy(1) + st*cp*Bcy(2)
        B2 = sp*ct*Bcy(0) + cp*Bcy(1) + st*sp*Bcy(2)
        B3 =   -st*Bcy(0)                + ct*Bcy(2)

        cth = COS(theta)
        sth = SIN(theta)
        cph = COS(phi)
        sph = SIN(phi)

        B(0) = sth*cph*B1 + sth*sph*B2 + cth*B3
        B(1) = cth*cph*B1 + cth*sph*B2 - sth*B3
        B(2) = -sph*B1 + cph*b2

        B(1) = -1.*B(1)
        B(2) = -1.*B(2)

	B(3)=B(0)*sth*cph + B(1)*cth*cph - B(2)*sph
	B(4)=B(0)*sth*sph + B(1)*cth*sph + B(2)*cph
	B(5)=B(0)*cth - B(1)*sth
        B(6)=sqrt(B(0)*B(0)+B(1)*B(1)+B(2)*B(2) )
        B = B*c
        bb=sqrt(B(3)*B(3)+B(4)*B(4)+B(5)*B(5) )
        db=abs((bb-B(6))/B(6))
if db gt 1e-2 then print,'can',pos(0),bb,B(6),db
        RETURN,B
        END

FUNCTION B_field, Pos,ff
        ;Multipole expansion Jovian magnetic field, translated
        ;from fortran O6RFT and O4RFT subroutines
        ;Pos(0) is radius, in the same units as the parameter, Rj.
        ;Pos(1) is latitude, _not_ co-latitude (i.e. 0 is the
        ;       north pole) in degrees.
        ;Pos(2) is West longitude (i.e. System III longitude)
        ;This coordinate system is right handed, but the
        ;       unit vectors in the latitude and longitude
        ;       directions both point anti-parallel to  the
        ;       traditional theta and phi of spherical coordinates


        IF ff EQ 0 OR ff EQ 2 THEN o4o6 = 1 ELSE o4o6 = 0
        IF ff EQ 0 OR ff EQ 1 THEN ish = 0 ELSE ish = 1
; ff  - 0=O4, 1=O6, 2=O4+sheet, 3=O6+sheet 

        B = DBLARR(7)
        sq3 = 1.73205
        sq5 = 2.23607
        sq7 = 2.64575
        sq8 = 2.82843
        sq15 = 3.87298
        dtr = !PI/180.
        theta = dtr*(90-Pos(1))
        phi = dtr*(360-Pos(2))

        ST = SIN(theta)
        CT = COS(theta)
        SP = SIN(phi)
        CP = COS(phi)

        RR = Pos(0)
        R3 = 1./RR^3
        R4 = R3/RR

        S2P = 2.*SP*CP
        C2P = 2.*CP*CP - 1.
        ST2 = ST*ST
        CT2 = CT*CT
        P20 = 1.5*CT2 - 0.5
        P21 = sq3*ST*CT
        P22 = 0.5*sq3*ST2
        DP20 = -3.*ST*CT
        DP21 = sq3*(CT2 - ST2)
        DP22 = sq3*ST*CT

        P21S = sq3*CT
        P22S = 0.5*sq3*ST

        R5 = R4/RR
        P30 = CT*(2.5*CT2 - 1.5)
        P31 = sq3*ST*(5.0*CT2 - 1.)/sq8
        P32 = 0.5*sq15*CT*ST2
        P33 = sq5*ST*ST2/sq8
        DP30 = ST*(-7.5*CT2 + 1.5)
        DP31 = sq3*CT*(5.*CT2 - 10.*ST2 - 1.)/sq8
        DP32 = sq15*ST*(CT2 - 0.5*ST2)
        DP33 = 3.0*sq5*ST2*CT/sq8

        P31S = sq3*(5.*CT2 - 1.)/sq8
        P32S = 0.5*sq15*ST*CT
        P33S = sq5*ST2/sq8
        S3P = SP*C2P + CP*S2P
        C3P = CP*C2P - SP*S2P

        IF O4O6 THEN GOTO,O4

        ;Parameters for the O6 model
        G10 = 4.24202
        G11 =-0.65929
        H11 = 0.24116
        G20 =-0.02181
        G21 =-0.71106
        H21 =-0.40304
        G22 = 0.48714
        H22 = 0.07179
        G30 = 0.07565
        G31 =-0.15493
        H31 =-0.38824
        G32 = 0.19775
        H32 = 0.34243
        G33 =-0.17958
        H33 =-0.22439

        GOTO, CALC

        O4: ;Parameters for the O4 model
        G10 = 4.218
        G11 =-0.664
        H11 = 0.264
        G20 =-0.203
        G21 =-0.735
        H21 =-0.469
        G22 = 0.513
        H22 = 0.088
        G30 =-0.233
        G31 =-0.076
        H31 =-0.580
        G32 = 0.168
        H32 = 0.487
        G33 =-0.231
        H33 =-0.294

        CALC: B(*) = 0.

        B(0) = B(0) + 4.*R5*(G30*P30 $
                + (G31*CP + H31*SP)*P31 $
                + (G32*C2P + H32*S2P)*P32 $
                + (G33*C3P + H33*S3P)*P33)
        B(1) = B(1) - R5*(G30*DP30 $
                + (G31*CP + H31*SP)*DP31 $
                + (G32*C2P + H32*S2P)*DP32 $
                + (G33*C3P + H33*S3P)*DP33)
        B(2) = B(2) + R5*((G31*SP - H31*CP)*P31S $
                + 2.*(G32*S2P - H32*C2P)*P32S $
                + 3.*(G33*S3P - H33*C3P)*P33S)

        B(0) = B(0) + 3.*R4*(G20*P20 $
                + (G21*CP + H21*SP)*P21 $
                + (G22*C2P + H22*S2P)*P22)
        B(1) = B(1) - R4*(G20*DP20 $
                + (G21*CP + H21*SP)*DP21 $
                + (G22*C2P + H22*S2P)*DP22)
        B(2) = B(2) + R4*((G21*SP - H21*CP)*P21S $
                + 2.*(G22*S2P - H22*C2P)*P22S)

        B(0) = B(0) + 2.*R3*(G10*CT $
                + (G11*CP + H11*SP)*ST)
        B(1) = B(1) - R3*(-1.*G10*ST $
                +(G11*CP + H11*SP)*CT)
        B(2) = B(2) + R3*(G11*SP - H11*CP)

        B(1) = -1.*B(1)
        B(2) = -1.*B(2)

	B(3)=B(0)*st*cp + B(1)*ct*cp - B(2)*sp
        B(4)=B(0)*st*sp + B(1)*ct*sp + B(2)*cp
        B(5)=B(0)*ct - B(1)*st
        B(6)=sqrt(B(0)*B(0)+B(1)*B(1)+B(2)*B(2) )

 	bb=sqrt(B(3)*B(3)+B(4)*B(4)+B(5)*B(5) )
        db=abs((bb-B(6))/B(6))
if db gt 1e-2 then print,'mod',pos(0),bb,B(6),db

        IF ish THEN B = B + CAN_SHEET(Pos)

        B(6)=sqrt(B(0)*B(0)+B(1)*B(1)+B(2)*B(2) )
 	bb=sqrt(B(3)*B(3)+B(4)*B(4)+B(5)*B(5) )
        db=abs((bb-B(6))/B(6))
if db gt 1e-2 then print,'tot',pos(0),bb,B(6),db

        RETURN, B
        END

