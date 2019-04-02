;+
; NAME:
; 
;   GETSCHMIDTCOEFS
;
; PURPOSE:
; 
;   This function retrieves the Schmidt quasi-normalized coefficients for
;   three different spherical harmonics models of the Jovian magnetosphere.
;   The available models are 'VIT4', 'VIP4', and 'VIPAL'.  
;
; CALLING SEQUENCE:
; 
;   gh = GETSCHMIDTCOEFS(MODEL, N, M)
;   g = gh[0]
;   h = gh[1]
;
; INPUTS:
;   MODEL: Name of the model
;   N:     Degree of the Schmidt quasi-normalized Legendre function.  N can take
;           on the values between 1 and 4 for the 'VIT4' and 'VIP4' models.  
;           For 'VIPAL', N can be between 1 and 5.
;   M:     Order of the Schmidt quasi-normalized Legendre function (0 to n)
;
; KEYWORD PARAMETERS:
;
;   INGAUSS: If set, the return value will be in gauss rather than nanotesla.
;             10^5 nanotesla = 1 Gauss
;
; OUTPUTS:
;
;   A two element array with the g and h coefficients for the Schmidt
;   quasi-normalized Legendre functions of degree n and order m.  The array is
;   arranged like so: [g, h].  The coefficients are given in nanotesla unless
;   the INGAUSS keyword is set.
;   
; PROCEDURE:
; 
;   Uses a series of case lookup.
;
; EXAMPLE:
;
;   To retrieve the coefficients of degree 3 and order 2 of the 'VIP4' model, use:
;     gh = getSchmidtCoefs('VIP4', 3, 2)
;     g = gh[0]
;     h = gh[1]
;     
;   To retrieve, in Gauss, the coefficients of degree 5 and order 0 
;   of the 'VIPAL' model, use:
;     gh = getSchmidtCoefs('VIPAL', 5, 0, /INGAUSS)
;     g = gh[0]
;     h = gh[1]
;
; MODIFICATION HISTORY:
;   Written by: Drake A. Ranquist, 7 Jun 2013
;   2016, Peter Kollmann: added VIP4dip that only includes the highest order moments
;   16 Feb 2018, Peter Kollmann: added JRM09 model
;-

function getSchmidtCoefs, model, n, m, INGAUSS=inGauss

IF NOT KEYWORD_SET(inGauss) THEN inGauss=0

gh = [0, 0]

CASE model OF
  'VIP4': CASE n OF ;; consistent with (but more digits) Connerney 98, Table 1, VIP4 column as well as with kk_2009.pro/JOVIAN_VIP4_no_dipole. Connerney uses the tradiational SIII(65) that is explained on Boulder page to be LEFT-handed, which agrees with the comments in kk_2009. 
    1: CASE m OF
      0: gh = [420543, 0]
      1: gh = [-65920, 24992]
    ENDCASE
    2: CASE m OF
      0: gh = [-5118, 0]
      1: gh = [-61904, -36052]
      2: gh = [49690, 5250]
    ENDCASE 
    3: CASE m OF
      0: gh = [-1576, 0]
      1: gh = [-52036, -8804]
      2: gh = [24386, 40829]
      3: gh = [-17597, -31586]
    ENDCASE   
    4: CASE m OF
      0: gh = [-16758, 0]
      1: gh = [22210, 7557]
      2: gh = [-6074, 40411]
      3: gh = [-20243, -16597]
      4: gh = [6643, 3866]
    ENDCASE
  ENDCASE
  'VIP4dip': CASE n OF
    1: CASE m OF
      0: gh = [420543, 0]
      1: gh = [-65920, 24992]
    ENDCASE
    2: CASE m OF
      0: gh = [-5118, 0]*0
      1: gh = [-61904, -36052]*0
      2: gh = [49690, 5250]*0
    ENDCASE 
    3: CASE m OF
      0: gh = [-1576, 0]*0
      1: gh = [-52036, -8804]*0
      2: gh = [24386, 40829]*0
      3: gh = [-17597, -31586]*0
    ENDCASE   
    4: CASE m OF
      0: gh = [-16758, 0]*0
      1: gh = [22210, 7557]*0
      2: gh = [-6074, 40411]*0
      3: gh = [-20243, -16597]*0
      4: gh = [6643, 3866]*0
    ENDCASE
  ENDCASE  
  'VIT4': CASE n OF
    1: CASE m OF
      0: gh = [428077, 0]
      1: gh = [-75306, 24616]
    ENDCASE
    2: CASE m OF
      0: gh = [-4283, 0]
      1: gh = [-59426, -50154]
      2: gh = [44386, 38452]
    ENDCASE 
    3: CASE m OF
      0: gh = [8906, 0]
      1: gh = [-21447, -17187]
      2: gh = [21130, 40667]
      3: gh = [-1190, -35263]
    ENDCASE   
    4: CASE m OF
      0: gh = [-22925, 0]
      1: gh = [18940, 16088]
      2: gh = [-3851, 11807]
      3: gh = [9926, 6195]
      4: gh = [1271, 12641]
    ENDCASE
  ENDCASE
  'VIPAL': CASE n OF
    1: CASE m OF
      0: gh = [420000, 0]
      1: gh = [-69750, 19730]
    ENDCASE
    2: CASE m OF
      0: gh = [64410, 0]
      1: gh = [-86720, -40410]
      2: gh = [95980, 60300]
    ENDCASE 
    3: CASE m OF
      0: gh = [-10580, 0]
      1: gh = [-59000, -23100]
      2: gh = [63220, 51600]
      3: gh = [46710, -11310]
    ENDCASE   
    4: CASE m OF
      0: gh = [-74660, 0]
      1: gh = [32820, 32830]
      2: gh = [-33800, -21310]
      3: gh = [18260, -6060]
      4: gh = [-14290, -4860]
    ENDCASE   
    5: CASE m OF
      0: gh = [-6600, 0]
      1: gh = [7370, 20650]
      2: gh = [-17110, -11670]
      3: gh = [-17930, -2880]
      4: gh = [-770, -500]
      5: gh = [-7400, -22790]
    ENDCASE
  ENDCASE 
  
  'JRM09': CASE n OF
  1: CASE m OF
0: gh=[410244.70,0]
1: gh=[-71498.300,21330.500]
ENDCASE
2: CASE m OF
0: gh=[11670.400,0]
1: gh=[-56835.800,-42027.300]
2: gh=[48689.500,19353.200]
ENDCASE
3: CASE m OF
0: gh=[4018.6000,0]
1: gh=[-37791.100,-32957.300]
2: gh=[15926.300,42084.500]
3: gh=[-2710.5000,-27544.200]
ENDCASE
4: CASE m OF
0: gh=[-34645.400,0]
1: gh=[-8247.6000,31994.500]
2: gh=[-2406.1000,27811.200]
3: gh=[-11083.800,-926.10000]
4: gh=[-17837.200,367.10000]
ENDCASE
5: CASE m OF
0: gh=[-18023.600,0]
1: gh=[4683.9000,45347.900]
2: gh=[16160.000,-749.00000]
3: gh=[-16402.000,6268.5000]
4: gh=[-2600.7000,10859.600]
5: gh=[-3660.7000,9608.4000]
ENDCASE
6: CASE m OF
0: gh=[-20819.600,0]
1: gh=[9992.9000,14533.100]
2: gh=[11791.800,-10592.900]
3: gh=[-12574.700,568.60000]
4: gh=[2669.7000,12871.700]
5: gh=[1113.2000,-4147.8000]
6: gh=[7584.9000,3604.4000]
ENDCASE
7: CASE m OF
0: gh=[598.40000,0]
1: gh=[4665.9000,-7626.3000]
2: gh=[-6495.7000,-10948.400]
3: gh=[-2516.5000,2633.3000]
4: gh=[-6448.5000,5394.2000]
5: gh=[1855.3000,-6050.8000]
6: gh=[-2892.9000,-1526.0000]
7: gh=[2968.0000,-5684.2000]
ENDCASE
8: CASE m OF
0: gh=[10059.200,0]
1: gh=[1934.4000,-2409.7000]
2: gh=[-6702.9000,-11614.600]
3: gh=[153.70000,9287.0000]
4: gh=[-4124.2000,-911.90000]
5: gh=[-867.20000,2754.5000]
6: gh=[-3740.6000,-2446.1000]
7: gh=[-732.40000,1207.3000]
8: gh=[-2433.2000,-2887.3000]
ENDCASE
9: CASE m OF
0: gh=[9671.8000,0]
1: gh=[-3046.2000,-8467.4000]
2: gh=[260.90000,-1383.8000]
3: gh=[2071.3000,5697.7000]
4: gh=[3329.6000,-2056.3000]
5: gh=[-2523.1000,3081.5000]
6: gh=[1787.1000,-721.20000]
7: gh=[-1148.2000,1352.5000]
8: gh=[1276.5000,-210.10000]
9: gh=[-1976.8000,1567.6000]
ENDCASE
10: CASE m OF
0: gh=[-2299.5000,0]
1: gh=[2009.7000,-4692.6000]
2: gh=[2127.8000,4445.8000]
3: gh=[3498.3000,-2378.6000]
4: gh=[2967.6000,-2204.3000]
5: gh=[16.300000,164.10000]
6: gh=[1806.5000,-1361.6000]
7: gh=[-46.500000,-2031.5000]
8: gh=[2897.8000,1411.8000]
9: gh=[574.50000,-714.30000]
10: gh=[1298.9000,1676.5000]
ENDCASE
  ENDCASE
ENDCASE

IF inGauss THEN gh /= 100000.

return, gh

END