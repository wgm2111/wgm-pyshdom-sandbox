      PARAMETER ( NLEVEL=15 )
      REAL*4 PGAS,PEDGE,DELP,C1,H,D,T,O,Q,S,OCM,WCM,WCMLO,HLO
      INTEGER*4 I,J,NPHD,NATM
      DIMENSION PGAS(NLEVEL),PEDGE(NLEVEL),DELP(NLEVEL)
      DATA PGAS/915.271,738.263,595.486,480.322,387.430,312.503,252.066,
     +203.318,163.997,132.281,106.699,86.0639,69.4195,55.9942,25.0000/
      DATA PEDGE/1013.25,817.293,659.232,531.740,428.904,345.957,
     +279.050,225.083,181.553,146.442,118.121,95.2769,76.8508,61.9883,
     +50.0001/
      DATA DELP/195.957,158.060,127.492,102.836,82.9479,66.9062,
     +53.9669,43.5300,35.1115,28.3211,22.8440,18.4261,14.8625,11.9882,
     +50.0000/
      C1=(8.314*273.16*100.0)/(18.*9.80665*1.01325)
      NPHD=2
      NATM=2
C
      X=29.0
      CALL ATMOS(P,X,D,T,O,Q,S,OCM,WCM,2,NATM)
      PRINT *,X,P,T,Q,S,WCM*0.8031/1000.,OCM
C
      X=27.0
      CALL ATMOS(P,X,D,T,O,Q,S,OCM,WCM,2,NATM)
      PRINT *,X,P,T,Q,S,WCM*0.8031/1000.,OCM
C
      X=25.0
      CALL ATMOS(P,X,D,T,O,Q,S,OCM,WCM,2,NATM)
      PRINT *,X,P,T,Q,S,WCM*0.8031/1000.,OCM
C
      X=23.0
      CALL ATMOS(P,X,D,T,O,Q,S,OCM,WCM,2,NATM)
      PRINT *,X,P,T,Q,S,WCM*0.8031/1000.,OCM
C
      X=21.0
      CALL ATMOS(P,X,D,T,O,Q,S,OCM,WCM,2,NATM)
      PRINT *,X,P,T,Q,S,WCM*0.8031/1000.,OCM
C
      X=19.0
      CALL ATMOS(P,X,D,T,O,Q,S,OCM,WCM,2,NATM)
      PRINT *,X,P,T,Q,S,WCM*0.8031/1000.,OCM
C
      X=17.0
      CALL ATMOS(P,X,D,T,O,Q,S,OCM,WCM,2,NATM)
      PRINT *,X,P,T,Q,S,WCM*0.8031/1000.,OCM
C
      X=15.0
      CALL ATMOS(P,X,D,T,O,Q,S,OCM,WCM,2,NATM)
      PRINT *,X,P,T,Q,S,WCM*0.8031/1000.,OCM
C
      X=13.0
      CALL ATMOS(P,X,D,T,O,Q,S,OCM,WCM,2,NATM)
      PRINT *,X,P,T,Q,S,WCM*0.8031/1000.,OCM

C
      X=11.0
      CALL ATMOS(P,X,D,T,O,Q,S,OCM,WCM,2,NATM)
      PRINT *,X,P,T,Q,S,WCM*0.8031/1000.,OCM
C
      X=11.0
      CALL ATMOS(P,X,D,T,O,Q,S,OCM,WCM,2,NATM)
      PRINT *,X,P,T,Q,S,WCM*0.8031/1000.,OCM
C
      X=9.0
      CALL ATMOS(P,X,D,T,O,Q,S,OCM,WCM,2,NATM)
      PRINT *,X,P,T,Q,S,WCM*0.8031/1000.,OCM
C
      X=7.0
      CALL ATMOS(P,X,D,T,O,Q,S,OCM,WCM,2,NATM)
      PRINT *,X,P,T,Q,S,WCM*0.8031/1000.,OCM
C
      X=5.0
      CALL ATMOS(P,X,D,T,O,Q,S,OCM,WCM,2,NATM)
      PRINT *,X,P,T,Q,S,WCM*0.8031/1000.,OCM
C
      X=3.0
      CALL ATMOS(P,X,D,T,O,Q,S,OCM,WCM,2,NATM)
      PRINT *,X,P,T,Q,S,WCM*0.8031/1000.,OCM
C
      X=1.0
      CALL ATMOS(P,X,D,T,O,Q,S,OCM,WCM,2,NATM)
      PRINT *,X,P,T,Q,S,WCM*0.8031/1000.,OCM
C
      X=0.5
      CALL ATMOS(P,X,D,T,O,Q,S,OCM,WCM,2,NATM)
      PRINT *,X,P,T,Q,S,WCM*0.8031/1000.,OCM
C
      STOP
      END
C                                                                               
      SUBROUTINE ATMOS(P,H,D,T,O,Q,S,OCM,WCM,NPHD,NATM)                         
C                                                                               
C     ------------------------------------------------------------------        
C     -------------     MCCLATCHY (1972) ATMOSPHERE DATA     -----------        
C     ------------------------------------------------------------------        
C                                                                               
C        INPUT DATA                                                             
C------------------                                                             
C                  NATM=0  GIVES ABREVIATED DATA FOR  STANDARD ATMOSPHER        
C                                (INPUT: P OR H) (RETURNS: H OR P & D,T)        
C                                                                               
C                  NATM=1  GIVES ATMOSPHERE DATA FOR  TROPICAL LATITUDES        
C                  NATM=2  GIVES ATMOSPHERE DATA FOR  MIDLATITUDE SUMMER        
C                  NATM=3  GIVES ATMOSPHERE DATA FOR  MIDLATITUDE WINTER        
C                  NATM=4  GIVES ATMOSPHERE DATA FOR  SUBARCTIC SUMMER          
C                  NATM=5  GIVES ATMOSPHERE DATA FOR  SUBARCTIC WINTER          
C                  NATM=6  GIVES ATMOSPHERE DATA FOR  STANDARD ATMOSPHER        
C                  NAMT=7      STANDARD ATMOSPHERE RH CONSTANT AT 0%            
C                  NAMT=8      STANDARD ATMOSPHERE RH CONSTANT AT 17%           
C                  NAMT=9      STANDARD ATMOSPHERE RH CONSTANT AT 87%           
C                                                                               
C                  NPHD=1  RETURNS H,D,T,O,Q,S DATA FOR GIVEN PRESSURE P        
C                  NPHD=2  RETURNS P,D,T,O,Q,S DATA FOR GIVEN   HEIGHT H        
C                  NPHD=3  RETURNS P,H,T,O,Q,S DATA FOR GIVEN  DENSITY D        
C                                                                               
C       OUTPUT DATA                                                             
C------------------                                                             
C                  P = PRESSURE IN MILLIBARS                                    
C                  H = HEIGHT IN KILOMETERS                                     
C                  D = DENSITY IN GRAMS/METER**3                                
C                  T = TEMPERATURE (ABSOLUTE)                                   
C                  O = OZONE MIXING RATIO (GRAMS OZONE)/(GRAMS AIR)             
C                  Q = SPECIFIC HUMIDITY (GRAMS WATER VAPOR)/(GRAMS AIR)        
C                  S = SATURATION RATIO (GRAMS WATER VAPOR)/(GRAMS AIR)         
C                  OCM = OZONE (CM-STP) ABOVE GIVEN HEIGHT                      
C                  WCM = WATER VAPOR (CM-STP) ABOVE GIVEN HEIGHT                
C                                                                               
C           REMARKS                                                             
C------------------                                                             
C                  INPUT P,H,D PARAMETERS ARE NOT ALTERED                       
C                  P,D INTERPOLATION IS EXPONENTIAL WITH HEIGHT                 
C                  NO EXTRAPOLATION IS MADE OUTSIDE 0-100 KM INTERVAL           
C                  S  IS NOT COMPUTED ABOVE 40 KM (FORMULA NOT ACCURATE)        
C                                                                               
C                  R = Q/S          GIVES RELATIVE HUMIDITY                     
C                  W = Q/(1-Q)      GIVES WATER VAPOR MIXING RATIO              
C                  N = D*2.079E 16  GIVES NUMBER DENSITY PER CM**3              
C                                                                               
C                                                                               
C                                                                               
C                                                                               
C                                                                               
      DIMENSION    PRS1(33),PRS2(33),PRS3(33),PRS4(33),PRS5(33),PRS6(33)        
     1            ,DNS1(33),DNS2(33),DNS3(33),DNS4(33),DNS5(33),DNS6(33)        
     2            ,TMP1(33),TMP2(33),TMP3(33),TMP4(33),TMP5(33),TMP6(33)        
     3            ,WVP1(33),WVP2(33),WVP3(33),WVP4(33),WVP5(33),WVP6(33)        
     4            ,OZO1(33),OZO2(33),OZO3(33),OZO4(33),OZO5(33),OZO6(33)        
      DIMENSION   PRES(33,6),DENS(33,6),TEMP(33,6),WVAP(33,6),OZON(33,6)        
C                                                                               
      EQUIVALENCE                                                               
     +     (PRES(1,1),PRS1(1)),(DENS(1,1),DNS1(1)),(TEMP(1,1),TMP1(1))          
     +    ,(PRES(1,2),PRS2(1)),(DENS(1,2),DNS2(1)),(TEMP(1,2),TMP2(1))          
     +    ,(PRES(1,3),PRS3(1)),(DENS(1,3),DNS3(1)),(TEMP(1,3),TMP3(1))          
     +    ,(PRES(1,4),PRS4(1)),(DENS(1,4),DNS4(1)),(TEMP(1,4),TMP4(1))          
     +    ,(PRES(1,5),PRS5(1)),(DENS(1,5),DNS5(1)),(TEMP(1,5),TMP5(1))          
     +    ,(PRES(1,6),PRS6(1)),(DENS(1,6),DNS6(1)),(TEMP(1,6),TMP6(1))          
      EQUIVALENCE (WVAP(1,1),WVP1(1)),(OZON(1,1),OZO1(1))                       
      EQUIVALENCE (WVAP(1,2),WVP2(1)),(OZON(1,2),OZO2(1))                       
      EQUIVALENCE (WVAP(1,3),WVP3(1)),(OZON(1,3),OZO3(1))                       
      EQUIVALENCE (WVAP(1,4),WVP4(1)),(OZON(1,4),OZO4(1))                       
      EQUIVALENCE (WVAP(1,5),WVP5(1)),(OZON(1,5),OZO5(1))                       
      EQUIVALENCE (WVAP(1,6),WVP6(1)),(OZON(1,6),OZO6(1))                       
C                                                                               
C                                                                               
      DIMENSION HTKM(33)                                                        
      DATA HTKM/1.0E-09, 1., 2., 3., 4., 5., 6., 7., 8., 9.,10.,11.             
     1         ,12.,13.,14.,15.,16.,17.,18.,19.,20.,21.,22.,23.,24.             
     2         ,25.,30.,35.,40.,45.,50.,70.,99.9/                               
C                                                                               
C                                                                               
C----------------------------------------------------------------------         
C0000 GLOBAL   U.S. (1976) STANDARD ATMOSPHERE   P, T, GEO H  PARAMETERS        
C----------------------------------------------------------------------         
C                                                                               
      DIMENSION SPLB(8),STLB(8),SHLB(8),SDLB(8)                                 
      DATA SPLB/1013.25,226.32,54.748,8.6801,1.109,.66938,.039564               
     +         ,3.7338E-03/                                                     
      DATA STLB/288.15,216.65,216.65,228.65,270.65,270.65,214.65,186.87/        
      DATA SHLB/0.0,11.0,20.0,32.0,47.0,51.0,71.0,84.852/                       
      DATA SDLB/-6.5,0.0,1.0,2.8,0.0,-2.8,-2.0,0.0/                             
      DATA HPCON/34.16319/                                                      
C                                                                               
C                                                                               
C-----------------------------------------------------------------------        
C1111 TROPICAL LATITUDES      MCCLATCHY (1972) ATMOSPHERE DATA VS HEIGHT        
C-----------------------------------------------------------------------        
C                                                                               
      DATA PRS1/      1.013E03,9.040E02,8.050E02,7.150E02,6.330E02,             
     1      5.590E02,4.920E02,4.320E02,3.780E02,3.290E02,2.860E02,              
     2      2.470E02,2.130E02,1.820E02,1.560E02,1.320E02,1.110E02,              
     3      9.370E01,7.890E01,6.660E01,5.650E01,4.800E01,4.090E01,              
     4      3.500E01,3.000E01,2.570E01,1.220E01,6.000E00,3.050E00,              
     5      1.590E00,8.540E-01,5.790E-02,3.000E-04/                             
      DATA DNS1/      1.167E03,1.064E03,9.689E02,8.756E02,7.951E02,             
     1      7.199E02,6.501E02,5.855E02,5.258E02,4.708E02,4.202E02,              
     2      3.740E02,3.316E02,2.929E02,2.578E02,2.260E02,1.972E02,              
     3      1.676E02,1.382E02,1.145E02,9.515E01,7.938E01,6.645E01,              
     4      5.618E01,4.763E01,4.045E01,1.831E01,8.600E00,4.181E00,              
     5      2.097E00,1.101E00,9.210E-02,5.000E-04/                              
      DATA TMP1/  300.0,294.0,288.0,284.0,277.0,270.0,264.0,257.0,250.0,        
     1244.0,237.0,230.0,224.0,217.0,210.0,204.0,197.0,195.0,199.0,203.0,        
     2207.0,211.0,215.0,217.0,219.0,221.0,232.0,243.0,254.0,265.0,270.0,        
     3  219.0,210.0/                                                            
      DATA WVP1/1.9E01,1.3E01,9.3E00,4.7E00,2.2E00,1.5E00,8.5E-01,              
     1  4.7E-01,2.5E-01,1.2E-01,5.0E-02,1.7E-02,6.0E-03,1.8E-03,1.0E-03,        
     2  7.6E-04,6.4E-04,5.6E-04,5.0E-04,4.9E-04,4.5E-04,5.1E-04,5.1E-04,        
     3  5.4E-04,6.0E-04,6.7E-04,3.6E-04,1.1E-04,4.3E-05,1.9E-05,6.3E-06,        
     4  1.4E-07,1.0E-09/                                                        
      DATA OZO1/5.6E-05,5.6E-05,5.4E-05,5.1E-05,4.7E-05,4.5E-05,4.3E-05,        
     1  4.1E-05,3.9E-05,3.9E-05,3.9E-05,4.1E-05,4.3E-05,4.5E-05,4.5E-05,        
     2  4.7E-05,4.7E-05,6.9E-05,9.0E-05,1.4E-04,1.9E-04,2.4E-04,2.8E-04,        
     3  3.2E-04,3.4E-04,3.4E-04,2.4E-04,9.2E-05,4.1E-05,1.3E-05,4.3E-06,        
     4  8.6E-08,4.3E-11/                                                        
C                                                                               
C-----------------------------------------------------------------------        
C2222 MIDLATITUDE SUMMER      MCCLATCHY (1972) ATMOSPHERE DATA VS HEIGHT        
C-----------------------------------------------------------------------        
C                                                                               
      DATA PRS2/      1.013E03,9.020E02,8.020E02,7.100E02,6.280E02,             
     1      5.540E02,4.870E02,4.260E02,3.720E02,3.240E02,2.810E02,              
     2      2.430E02,2.090E02,1.790E02,1.530E02,1.300E02,1.110E02,              
     3      9.500E01,8.120E01,6.950E01,5.950E01,5.100E01,4.370E01,              
     4      3.760E01,3.220E01,2.770E01,1.320E01,6.520E00,3.330E00,              
     5      1.760E00,9.510E-01,6.710E-02,3.000E-04/                             
      DATA DNS2/      1.191E03,1.080E03,9.757E02,8.846E02,7.998E02,             
     1      7.211E02,6.487E02,5.830E02,5.225E02,4.669E02,4.159E02,              
     2      3.693E02,3.269E02,2.882E02,2.464E02,2.104E02,1.797E02,              
     3      1.535E02,1.305E02,1.110E02,9.453E01,8.056E01,6.872E01,              
     4      5.867E01,5.014E01,4.288E01,1.322E01,6.519E00,3.330E00,              
     5      1.757E00,9.512E-01,6.706E-02,5.000E-04/                             
      DATA TMP2/  294.0,290.0,285.0,279.0,273.0,267.0,261.0,255.0,248.0,        
     1242.0,235.0,229.0,222.0,216.0,216.0,216.0,216.0,216.0,216.0,217.0,        
     2218.0,219.0,220.0,222.0,223.0,224.0,234.0,245.0,258.0,270.0,276.0,        
     3  218.0,210.0/                                                            
      DATA WVP2/1.4E01,9.3E00,5.9E00,3.3E00,1.9E00,1.0E00,6.1E-01,              
     1  3.7E-01,2.1E-01,1.2E-01,6.4E-02,2.2E-02,6.0E-03,1.8E-03,1.0E-03,        
     2  7.6E-04,6.4E-04,5.6E-04,5.0E-04,4.9E-04,4.5E-04,5.1E-04,5.1E-04,        
     3  5.4E-04,6.0E-04,6.7E-04,3.6E-04,1.1E-04,4.3E-05,1.9E-05,6.3E-06,        
     4  1.4E-07,1.0E-09/                                                        
      DATA OZO2/6.0E-05,6.0E-05,6.0E-05,6.2E-05,6.4E-05,6.6E-05,6.9E-05,        
     1  7.5E-05,7.9E-05,8.6E-05,9.0E-05,1.1E-04,1.2E-04,1.5E-04,1.8E-04,        
     2  1.9E-04,2.1E-04,2.4E-04,2.8E-04,3.2E-04,3.4E-04,3.6E-04,3.6E-04,        
     3  3.4E-04,3.2E-04,3.0E-04,2.0E-04,9.2E-05,4.1E-05,1.3E-05,4.3E-06,        
     4  8.6E-08,4.3E-11/                                                        
C                                                                               
C-----------------------------------------------------------------------        
C3333 MIDLATITUDE WINTER      MCCLATCHY (1972) ATMOSPHERE DATA VS HEIGHT        
C-----------------------------------------------------------------------        
C                                                                               
      DATA PRS3/      1.018E03,8.973E02,7.897E02,6.938E02,6.081E02,             
     1      5.313E02,4.627E02,4.016E02,3.473E02,2.992E02,2.568E02,              
     2      2.199E02,1.882E02,1.610E02,1.378E02,1.178E02,1.007E02,              
     3      8.610E01,7.350E01,6.280E01,5.370E01,4.580E01,3.910E01,              
     4      3.340E01,2.860E01,2.430E01,1.110E01,5.180E00,2.530E00,              
     5      1.290E00,6.820E-01,4.670E-02,3.000E-04/                             
      DATA DNS3/      1.301E03,1.162E03,1.037E03,9.230E02,8.282E02,             
     1      7.411E02,6.614E02,5.886E02,5.222E02,4.619E02,4.072E02,              
     2      3.496E02,2.999E02,2.572E02,2.206E02,1.890E02,1.620E02,              
     3      1.388E02,1.188E02,1.017E02,8.690E01,7.421E01,6.338E01,              
     4      5.415E01,4.624E01,3.950E01,1.783E01,7.924E00,3.625E00,              
     5      1.741E00,8.954E-01,7.051E-02,5.000E-04/                             
      DATA TMP3/  272.2,268.7,265.2,261.7,255.7,249.7,243.7,237.7,231.7,        
     1225.7,219.7,219.2,218.7,218.2,217.7,217.2,216.7,216.2,215.7,215.2,        
     2215.2,215.2,215.2,215.2,215.2,215.2,217.4,227.8,243.2,258.5,265.7,        
     3  230.7,210.2/                                                            
      DATA WVP3/3.5E00,2.5E00,1.8E00,1.2E00,6.6E-01,3.8E-01,2.1E-01,            
     1  8.5E-02,3.5E-02,1.6E-02,7.5E-03,6.9E-03,6.0E-03,1.8E-03,1.0E-03,        
     2  7.6E-04,6.4E-04,5.6E-04,5.0E-04,4.9E-04,4.5E-04,5.1E-04,5.1E-04,        
     3  5.4E-04,6.0E-04,6.7E-04,3.6E-04,1.1E-04,4.3E-05,1.9E-05,6.3E-06,        
     4  1.4E-07,1.0E-09/                                                        
      DATA OZO3/6.0E-05,5.4E-05,4.9E-05,4.9E-05,4.9E-05,5.8E-05,6.4E-05,        
     1  7.7E-05,9.0E-05,1.2E-04,1.6E-04,2.1E-04,2.6E-04,3.0E-04,3.2E-04,        
     2  3.4E-04,3.6E-04,3.9E-04,4.1E-04,4.3E-04,4.5E-04,4.3E-04,4.3E-04,        
     3  3.9E-04,3.6E-04,3.4E-04,1.9E-04,9.2E-05,4.1E-05,1.3E-05,4.3E-06,        
     4  8.6E-08,4.3E-11/                                                        
C                                                                               
C-----------------------------------------------------------------------        
C4444 SUBARCTIC SUMMER        MCCLATCHY (1972) ATMOSPHERE DATA VS HEIGH         
C----------------------------------------------------------------------         
C                                                                               
      DATA PRS4/      1.010E03,8.960E02,7.929E02,7.000E02,6.160E02,             
     1      5.410E02,4.730E02,4.130E02,3.590E02,3.107E02,2.677E02,              
     2      2.300E02,1.977E02,1.700E02,1.460E02,1.250E02,1.080E02,              
     3      9.280E01,7.980E01,6.860E01,5.890E01,5.070E01,4.360E01,              
     4      3.750E01,3.227E01,2.780E01,1.340E01,6.610E00,3.400E00,              
     5      1.810E00,9.870E-01,7.070E-02,3.000E-04/                             
      DATA DNS4/      1.220E03,1.110E03,9.971E02,8.985E02,8.077E02,             
     1      7.244E02,6.519E02,5.849E02,5.231E02,4.663E02,4.142E02,              
     2      3.559E02,3.059E02,2.630E02,2.260E02,1.943E02,1.671E02,              
     3      1.436E02,1.235E02,1.062E02,9.128E01,7.849E01,6.750E01,              
     4      5.805E01,4.963E01,4.247E01,1.338E01,6.614E00,3.404E00,              
     5      1.817E00,9.868E-01,7.071E-02,5.000E-04/                             
      DATA TMP4/  287.0,282.0,276.0,271.0,266.0,260.0,253.0,246.0,239.0,        
     1232.0,225.0,225.0,225.0,225.0,225.0,225.0,225.0,225.0,225.0,225.0,        
     2225.0,225.0,225.0,225.0,226.0,228.0,235.0,247.0,262.0,274.0,277.0,        
     3  216.0,210.0/                                                            
      DATA WVP4/9.1E00,6.0E00,4.2E00,2.7E00,1.7E00,1.0E00,5.4E-01,              
     1  2.9E-01,1.3E-02,4.2E-02,1.5E-02,9.4E-03,6.0E-03,1.8E-03,1.0E-03,        
     2  7.6E-04,6.4E-04,5.6E-04,5.0E-04,4.9E-04,4.5E-04,5.1E-04,5.1E-04,        
     3  5.4E-04,6.0E-04,6.7E-04,3.6E-04,1.1E-04,4.3E-05,1.9E-05,6.3E-06,        
     4  1.4E-07,1.0E-09/                                                        
      DATA OZO4/4.9E-05,5.4E-05,5.6E-05,5.8E-05,6.0E-05,6.4E-05,7.1E-05,        
     1  7.5E-05,7.9E-05,1.1E-04,1.3E-04,1.8E-04,2.1E-04,2.6E-04,2.8E-04,        
     2  3.2E-04,3.4E-04,3.9E-04,4.1E-04,4.1E-04,3.9E-04,3.6E-04,3.2E-04,        
     3  3.0E-04,2.8E-04,2.6E-04,1.4E-04,9.2E-05,4.1E-05,1.3E-05,4.3E-06,        
     4  8.6E-08,4.3E-11/                                                        
C                                                                               
C----------------------------------------------------------------------         
C5555 SUBARCTIC WINTER        MCCLATCHY (1972) ATMOSPHERE DATA VS HEIGH         
C----------------------------------------------------------------------         
C                                                                               
      DATA PRS5/      1.013E03,8.878E02,7.775E02,6.798E02,5.932E02,             
     1      5.158E02,4.467E02,3.853E02,3.308E02,2.829E02,2.418E02,              
     2      2.067E02,1.766E02,1.510E02,1.291E02,1.103E02,9.431E01,              
     3      8.058E01,6.882E01,5.875E01,5.014E01,4.277E01,3.647E01,              
     4      3.109E01,2.649E01,2.256E01,1.020E01,4.701E00,2.243E00,              
     5      1.113E00,5.719E-01,4.016E-02,3.000E-04/                             
      DATA DNS5/      1.372E03,1.193E03,1.058E03,9.366E02,8.339E02,             
     1      7.457E02,6.646E02,5.904E02,5.226E02,4.538E02,3.879E02,              
     2      3.315E02,2.834E02,2.422E02,2.071E02,1.770E02,1.517E02,              
     3      1.300E02,1.113E02,9.529E01,8.155E01,6.976E01,5.966E01,              
     4      5.100E01,4.358E01,3.722E01,1.645E01,7.368E00,3.330E00,              
     5      1.569E00,7.682E-01,5.695E-02,5.000E-04/                             
      DATA TMP5/  257.1,259.1,255.9,252.7,247.7,240.9,234.1,227.3,220.6,        
     1217.2,217.2,217.2,217.2,217.2,217.2,217.2,216.6,216.0,215.4,214.8,        
     2214.1,213.6,213.0,212.4,211.8,211.2,216.0,222.2,234.7,247.0,259.3,        
     3  245.7,210.0/                                                            
      DATA WVP5/1.2E00,1.2E00,9.4E-01,6.8E-01,4.1E-01,2.0E-01,9.8E-02,          
     1  5.4E-02,1.1E-02,8.4E-03,5.5E-03,3.8E-03,2.6E-03,1.8E-03,1.0E-03,        
     2  7.6E-04,6.4E-04,5.6E-04,5.0E-04,4.9E-04,4.5E-04,5.1E-04,5.1E-04,        
     3  5.4E-04,6.0E-04,6.7E-04,3.6E-04,1.1E-04,4.3E-05,1.9E-05,6.3E-06,        
     4  1.4E-07,1.0E-09/                                                        
      DATA OZO5/4.1E-05,4.1E-05,4.1E-05,4.3E-05,4.5E-05,4.7E-05,4.9E-05,        
     1  7.1E-05,9.0E-05,1.6E-04,2.4E-04,3.2E-04,4.3E-04,4.7E-04,4.9E-04,        
     2  5.6E-04,6.2E-04,6.2E-04,6.2E-04,6.0E-04,5.6E-04,5.1E-04,4.7E-04,        
     3  4.3E-04,3.6E-04,3.2E-04,1.5E-04,9.2E-05,4.1E-05,1.3E-05,4.3E-06,        
     4  8.6E-08,4.3E-11/                                                        
C                                                                               
C----------------------------------------------------------------------         
C6666 GLOBAL   U.S. (1976) STANDARD ATMOSPHERE   P, T, GEO H  PARAMETER         
C----------------------------------------------------------------------         
C                                                                               
      DATA PRS6/    1.01325E+03,8.987E+02,7.950E+02,7.011E+02,6.164E+02,        
     1      5.402E+02,4.718E+02,4.106E+02,3.560E+02,3.074E+02,2.644E+02,        
     2      2.263E+02,1.933E+02,1.651E+02,1.410E+02,1.204E+02,1.029E+02,        
     3      8.787E+01,7.505E+01,6.410E+01,5.475E+01,4.678E+01,4.000E+01,        
     4      3.422E+01,2.931E+01,2.511E+01,1.172E+01,5.589E+00,2.775E+00,        
     5      1.431E+00,7.594E-01,4.634E-02,2.384E-04/                            
      DATA DNS6/      1.225E+03,1.112E+03,1.006E+03,9.091E+02,8.191E+02,        
     1      7.361E+02,6.597E+02,5.895E+02,5.252E+02,4.663E+02,4.127E+02,        
     2      3.639E+02,3.108E+02,2.655E+02,2.268E+02,1.937E+02,1.654E+02,        
     3      1.413E+02,1.207E+02,1.031E+02,8.803E+01,7.487E+01,6.373E+01,        
     4      5.428E+01,4.627E+01,3.947E+01,1.801E+01,8.214E+00,3.851E+00,        
     5      1.881E+00,9.775E-01,7.424E-02,4.445E-04/                            
      DATA TMP6/                                                                
     1         288.150,281.650,275.150,268.650,262.150,255.650,249.150,         
     2         242.650,236.150,229.650,223.150,216.650,216.650,216.650,         
     3         216.650,216.650,216.650,216.650,216.650,216.650,216.650,         
     4         217.650,218.650,219.650,220.650,221.650,226.650,237.050,         
     5         251.050,265.050,270.650,217.450,186.870/                         
      DATA WVP6/      1.083E+01,6.323E+00,3.612E+00,2.015E+00,1.095E+00,        
     1      5.786E-01,2.965E-01,1.469E-01,7.021E-02,3.226E-02,1.419E-02,        
     2      5.956E-03,5.002E-03,4.186E-03,3.490E-03,2.896E-03,2.388E-03,        
     3      1.954E-03,1.583E-03,1.267E-03,9.967E-04,8.557E-04,7.104E-04,        
     4      5.600E-04,4.037E-04,2.406E-04,5.404E-05,2.464E-05,1.155E-05,        
     5      5.644E-06,2.932E-06,2.227E-07,1.334E-09/                            
      DATA OZO6/      7.526E-05,3.781E-05,6.203E-05,3.417E-05,5.694E-05,        
     1      3.759E-05,5.970E-05,4.841E-05,7.102E-05,6.784E-05,9.237E-05,        
     2      9.768E-05,1.251E-04,1.399E-04,1.715E-04,1.946E-04,2.300E-04,        
     3      2.585E-04,2.943E-04,3.224E-04,3.519E-04,3.714E-04,3.868E-04,        
     4      3.904E-04,3.872E-04,3.728E-04,2.344E-04,9.932E-05,3.677E-05,        
     5      1.227E-05,4.324E-06,5.294E-08,1.262E-10/                            
C                                                                               
C     SKELETON STANDARD ATMOSPHERE CASE                                         
C     PROFILES VS HEIGHT                                                        
C                                                                               
      IF(NATM.GT.0) GO TO 200                                                   
      O=1.E-10                                                                  
      Q=1.E-10                                                                  
      S=1.E-10                                                                  
      OCM=1.E-10                                                                
      WCM=1.E-10                                                                
      IF(NPHD.LT.2) GO TO 150                                                   
      DO 110 N=2,8                                                              
      IF(H.LT.SHLB(N)) GO TO 120                                                
 110  CONTINUE                                                                  
      N=9                                                                       
 120  N=N-1                                                                     
      IF(ABS(SDLB(N)).LT.1.E-04) GO TO 130                                      
      P=SPLB(N)*(1.E0+SDLB(N)/STLB(N)*(H-SHLB(N)))**(-HPCON/SDLB(N))            
      GO TO 140                                                                 
 130  P=SPLB(N)*EXP(-HPCON/STLB(N)*(H-SHLB(N)))                                 
 140  T=STLB(N)+SDLB(N)*(H-SHLB(N))                                             
      D=P/T*28.9644E05/8.31432E03                                               
      RETURN                                                                    
C                                                                               
C     PROFILES VS PRESSURE                                                      
C                                                                               
 150  CONTINUE                                                                  
      DO 160 N=2,8                                                              
 160  IF(P.GT.SPLB(N)) GO TO 170                                                
      N=9                                                                       
 170  N=N-1                                                                     
      IF(ABS(SDLB(N)).LT.1.E-04) GO TO 180                                      
      H=SHLB(N)+STLB(N)/SDLB(N)*((SPLB(N)/P)**(SDLB(N)/HPCON)-1.E0)             
      GO TO 190                                                                 
 180  H=SHLB(N)+STLB(N)/HPCON*ALOG(SPLB(N)/P)                                   
 190  T=STLB(N)+SDLB(N)*(H-SHLB(N))                                             
      D=P/T*28.9644E05/8.31432E03                                               
      RETURN                                                                    
C                                                                               
C     FULL ATMOSPHERIC CASES                                                    
C                                                                               
 200  CONTINUE                                                                  
      IF(NPHD.EQ.1) GO TO 240                                                   
      IF(NPHD.EQ.2) GO TO 220                                                   
      XX=D                                                                      
      XI=DENS(1,NATM)                                                           
      IF(D.GT.XI) XX=XI                                                         
      IF(D.LT.5.0E-04) GO TO 280                                                
      DO 210 J=2,33                                                             
      XJ=DENS(J,NATM)                                                           
      IF(XX.GT.XJ) GO TO 260                                                    
 210  XI=XJ                                                                     
C                                                                               
C    PROFILE VS HEIGHT                                                          
C                                                                               
 220  XX=H                                                                      
      XI=HTKM(1)                                                                
      IF(H.LT.XI) XX=XI                                                         
      IF(H.GT.99.9E0) GO TO 280                                                 
      DO 230 J=2,33                                                             
      XJ=HTKM(J)                                                                
      IF(XX.LT.XJ) GO TO 260                                                    
 230  XI=XJ                                                                     
C                                                                               
C    PROFILE VS PRESSURE                                                        
C                                                                               
 240  XX=P                                                                      
      XI=PRES(1,NATM)                                                           
      IF(P.GT.XI) XX=XI                                                         
      IF(P.LT.3.0E-04) GO TO 280                                                
      DO 250 J=2,33                                                             
      XJ=PRES(J,NATM)                                                           
      IF(XX.GT.XJ) GO TO 260                                                    
 250  XI=XJ                                                                     
 260  DELTA=(XX-XI)/(XJ-XI)                                                     
      I=J-1                                                                     
      IF(NPHD.NE.2) H=HTKM(I)+(HTKM(J)-HTKM(I))*ALOG(XX/XI)/ALOG(XJ/XI)         
      PI=PRES(I,NATM)                                                           
      PJ=PRES(J,NATM)                                                           
      DI=DENS(I,NATM)                                                           
      DJ=DENS(J,NATM)                                                           
C     IF(NPHD.NE.1) P=PI*(PJ/PI)**DELTA                                         
C     IF(NPHD.NE.3) D=DI*(DJ/DI)**DELTA                                         
C                                                                               
C     LINEAR INTERPOLATION OF PRESSURE AND DENSITY WITH HEIGHT                  
C                                                                               
      IF(NPHD.NE.1) P=PI+DELTA*(PJ-PI)                                          
      IF(NPHD.NE.3) D=DI+DELTA*(DJ-DI)                                          
C                                                                               
      T=TEMP(I,NATM)+DELTA*(TEMP(J,NATM)-TEMP(I,NATM))                          
      O=OZON(I,NATM)/DI+DELTA*(OZON(J,NATM)/DJ-OZON(I,NATM)/DI)                 
      Q=WVAP(I,NATM)/DI+DELTA*(WVAP(J,NATM)/DJ-WVAP(I,NATM)/DI)                 
      ES=10.E0**(9.4051E0-2353.0E0/T)                                           
      IF(P.LT.PI) PI=P                                                          
      S=1.E+06                                                                  
      RS=(PI-ES+0.622E0*ES)/(0.622E0*ES)                                        
      IF(RS.GT.1.E-06) S=1.E0/RS                                                
      OI=O                                                                      
      QI=Q                                                                      
      OCM=0.E0                                                                  
      WCM=0.E0                                                                  
      DO 270 K=J,33                                                             
      PJ=PRES(K,NATM)                                                           
      DJ=DENS(K,NATM)                                                           
      OJ=OZON(K,NATM)/DJ                                                        
      QJ=WVAP(K,NATM)/DJ                                                        
      DP=PI-PJ                                                                  
      OCM=OCM+0.5E0*(OI+OJ)*DP                                                  
      WCM=WCM+0.5E0*(QI+QJ)*DP                                                  
      OI=OJ                                                                     
      QI=QJ                                                                     
 270  PI=PJ                                                                     
      WCM=WCM/0.980E0*22420.7E0/18.0E0                                          
      OCM=OCM/0.980E0*22420.7E0/48.0E0                                          
      RETURN                                                                    
 280  T=210.0E0                                                                 
      IF(NATM.EQ.6) T=186.87                                                    
      O=1.E-10                                                                  
      Q=1.E-10                                                                  
      S=1.E-10                                                                  
      OCM=1.E-10                                                                
      WCM=1.E-10                                                                
      IF(NPHD.NE.1) P=1.E-05                                                    
      IF(NPHD.NE.2) H=99.99E0                                                   
      IF(NPHD.NE.3) D=2.E-05                                                    
      RETURN                                                                    
      END                                                                       
C                                                                               
