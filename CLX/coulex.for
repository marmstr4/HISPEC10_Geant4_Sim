      PROGRAM ECCUL                                                             
C     MULTIPLE COULOMB EXCITATION PROGRAM  COULEX 1980                       900
C         DIMENSION OF THE ZETA ARRAY (NUMBER OF COUPLINGS) = ZZZZ          1000
C         DIMENSION OF THE CATALOGUE(NUMBER OF MAGN.SUBSTATES)=CCCC         1100
C         DIMENSION OF THE LEVEL SCHEME(NUMBER OF EXCITED STATES)=NNNN      1200
C                                                                           1300
C     FOR INFORMATION WRITE TO:                                             1400
C                                                                           1500
C       J. DE BOER                                                          1600
C       AM COULOMBWALL 1                                                    1700
C       D-8046 GARCHING                                                     1800
C       WEST GERMANY                                                        1900
C                                                                           2000
C *********************************************************************     2100
C *                                                                   *     2200
C *  THE INPUT MAY BE GIVEN AS INTEGERS OR DECIMAL NUMBERS            *     2300
C *  WITH OR WITHOUT AN EXPONENT INDICATED BY AN E.                   *     2400
C *  THE NUMBERS MAY BEGIN WITH A SIGN AND THEY MUST BE SEPARATED     *     2500
C *  BY ONE COMMA OR ONE BLANK.                                       *     2600
C *  AFTER THE NUMBERS THE CARD MAY BE FILLED WITH COMMENTS.          *     2700
C *                                                                   *     2800
C *  CONTENT OF THE DATA CARDS (DEFAULT VALUES IN PARENTHESIS)        *     2900
C *                                                                   *     3000
C *   1    NMAX (2)                                                   *     3100
C *   2    NCM (2)                                                    *     3200
C *   3    NTIME (3600)                                               *     3300
C *   4    XIMAX (6.0)                                                *     3400
C *   5    EMMAX1 (50.0)                                              *     3500
C *   6    ACCUR (0.0001)                                             *     3600
C *   7    OUXI (1)                                                   *     3700
C *   8    OUPSI (1)                                                  *     3800
C *   9    OUAMP (0)                                                  *     3900
C *  10    OUPROW (0)                                                 *     4000
C *  11    OURHOB (0)                                                 *     4100
C *  12    OURHOC (0)                                                 *     4200
C *  14    ALPHA(0)BETA(0)    GAMMA(0)                                *     4300
C *  15    DIPOL (0.005)                                              *     4400
C *  16    INTERV (10)                                                *     4500
C *  17    ZP      A1                                                 *     4600
C *  18    ZT      A2                                                 *     4700
C *  19    EP                                                         *     4800
C *  20    TLBDG (0.0)                                                *     4900
C *  21    THETA (0.0)                                                *     5000
C *  22    N       SPIN(N)    EN(N)           IPAR(N)(+1) BANDK(N)    *     5100
C *  23    N       M          ME(N,M,LA)(0.0) LA (2)                  *     5200
C *  24    N       M          MM(N,M,LA)(0.0) LA (1)                  *     5300
C *  25    DTHETA (10.0)                                              *     5400
C *  26    NGMAX (1)                                                  *     5500
C *  27    I     NSTART(I)(0) NSTOP(I)(0)    EMMAX(I) (50.0)          *     5600
C *  28    I1      I2         MASTER(I1,I2) (1)                       *     5700
C *  29    I       IEXCIT(I) (2)                                      *     5800
C *  30    INTM1(0)INTE1(0)   INTE2 (1.0)    INTE3(0)                 *     5900
C *  31-40 FREE TO BE ASSIGNED                                        *     6100
C *                                                                   *     6200
C *********************************************************************     6300
C                                                                           6400
C                                                                           6500
C                                                                           6600
C                                                                           6700
      INTEGER OURHOB,OURHOC,OURHOX                                          6800
      INTEGER OUAMP,OUPROW,OUXI,OUPSI                                       6900
      LOGICAL ERR                                                           7000
      REAL MEM,MULT                                                         7100
      COMMON /BL1/LAMDA(12),LEAD(NNNN,NNNN,12),LDNUM(12,NNNN),LAMMAX        7200
      COMMON /BL2/NGMAX,NSTART(10),NSTOP(10),MASTER(10,10),EMMAX(10)        7300
      COMMON /BL3/EN(NNNN),SPIN(NNNN),ACCUR,DIPOL,ZPOL,BANDK(NNNN)          7400
      COMMON /BL4/MAGEXC                                                    7500
      COMMON /BL7/ZETA( ZZZZ),NZMAX                                         7600
      COMMON /BL8/ASQRT,B(200),IEX(200)                                     7700
      COMMON /BL9/IEXCIT(10),IEXNUM,LMAX(2),IRSTA(2),IRSTO(2),IP            7800
      COMMON/BL10/IPAR(NNNN),IFAC(NNNN)                                     7900
      COMMON/BL11/INTEND                                                    8000
      COMMON/BL18/CAT(CCCC,3),ICATMX, ISMAX                                 8100
      COMMON/BL19/ERR                                                       8200
      COMMON/BL20/MEM(NNNN,NNNN,12),MULT(NNNN,NNNN,12)                      8300
      COMMON/BL21/XIMAX                                                     8400
      COMMON/BL23/IDN17,IDN18,IDN19,IDN25                                   8500
      COMMON/BL31/OURHOB,OURHOC,OURHOX                                      8600
      COMMON/BL32/NMAX,NDIM                                                 8700
      COMMON/BL33/NTIME                                                     8800
      COMMON/BL34/SIMP                                                      8900
      COMMON/BL35/ZP,ZT,A1,A2,EP,TLBDG                                      9000
      COMMON/BL36/INTM1,INTE1,INTE2,INTE3,INTE4,INTE5,INTE6                 9100
      COMMON/BL38/ALPHA,BETA,GAMMA                                          9200
      COMMON/BL48/NANGLE,THI,TLI                                            9300
      COMMON/BL49/OUAMP,OUPROW,INTERV,INTIN                                 9400
      COMMON/BL50/OUXI,OUPSI,NCM,EMMAX1                                     9500
      COMMON/BL51/SIGTOT(NNNN)                                              9600
      COMMON/BL52/IPURG                                                     9700
      COMMON/BL55/DTHETA,DTLBDG                                             9800
                                                                                
      OPEN(UNIT=20,DEVICE='DSK',FILE='PARA.DAT',MODE='ASCII',                   
     1ACCESS='SEQIN')             !******************                           
C                                                                           9900
C     THE FOLLOWING THREE QUANTITIES MUST BE ADJUSTED IF THE DIMENSIONS    10000
C     IN THE PROGRAM ARE CHANGED                                           10100
      NZMAX =  ZZZZ                                                        10200
      ICATMX = CCCC                                                        10300
      NDIM = NNNN                                                          10400
C                                                                          10500
C     INITIALIZE FACTORIAL ARRAYS  B  AND  IEX                             10600
      A =  85.0                                                            10700
      ASQRT = SQRT(A)                                                      10800
      F = 1.0                                                              10900
      B(1) = 1.                                                            11000
      IEX(1) = 0                                                           11100
      DO 80 I=1,199                                                        11200
      F = F * FLOAT(I)/A                                                   11300
      IEX(I+1) = ALOG10(F)                                                 11400
      B(I+1) = F * 10.0**(-IEX(I+1))                                       11500
  80  CONTINUE                                                             11600
C                                                                          11700
C     DEFAULT VALUES                                                       11800
 100  NMAX = 2                                                             11900
      NGMAX = 1                                                            12000
      NCM = 2                                                              12100
      NTIME = 3600                                                         12200
      INTERV = 10                                                          12300
      INTIN = INTERV                                                       12400
      XIMAX = 6.                                                           12500
      EMMAX1 = 50.0                                                        12600
      ACCUR = .0001                                                        12700
      DIPOL = 0.005                                                        12800
      OUXI = 1                                                             12900
      OUPSI = 1                                                            13000
      OUAMP = 0                                                            13100
      OUPROW = 0                                                           13200
      OURHOB = 0                                                           13300
      OURHOC = 0                                                           13400
      OURHOX = 0                                                           13500
      DTHETA = 10.                                                         13600
      TLI=0.0                                                              13700
      THI=0.0                                                              13800
      MAGEXC = 0                                                           13900
      DO 140 I=1,10                                                        14000
      EMMAX(I) = 50.0                                                      14100
      IEXCIT(I) = 2                                                        14200
      NSTART(I) = 0                                                        14300
      NSTOP(I) = 0                                                         14400
      DO 120 J=1,10                                                        14500
      MASTER(I,J) = 1                                                      14600
 120  CONTINUE                                                             14700
 140  CONTINUE                                                             14800
      LAMMAX = 0                                                           14900
      DO 160 LAM = 1,12                                                    15000
 160  LAMDA(LAM) = 1000                                                    15100
      DO 220 I=1,NDIM                                                      15200
      EN(I) = 0.0                                                          15300
      SPIN(I) = 0.0                                                        15400
      IPAR(I) = 1                                                          15500
      DO 200 J=1,NDIM                                                      15600
      DO 180 LAM=1,12                                                      15700
      MEM(I,J,LAM)=0.0                                                     15800
      MULT(I,J,LAM)=1.0                                                    15900
 180  CONTINUE                                                             16000
 200  CONTINUE                                                             16100
 220  CONTINUE                                                             16200
      ZP=0.0                                                               16300
      ZT=0.0                                                               16400
      A1=0.0                                                               16500
      A2=0.0                                                               16600
      EP=0.0                                                               16700
      INTE2 = 1.0                                                          16800
      INTM1=0.0                                                            16900
      INTE1=0.0                                                            17000
      INTE3=0.0                                                            17100
      ALPHA = 0.                                                           17200
      BETA = 0.                                                            17300
      GAMMA = 0.                                                           17400
      IDN17 = 0                                                            17500
      IDN18 = 0                                                            17600
      IDN19 = 0                                                            17700
      IDN25 = 0                                                            17800
C                                                                          17900
C     EXECUTION OF PROGRAM STARTS HERE                                     18000
C                                                                          18100
 250  ERR = .FALSE.                                                        18200
      IPURG = 0                                                            18300
      INTEND = 0                                                           18400
      CALL ZCLOCK(TIME01)                                                  18500
      CALL READER                                                          18600
      IF(IPURG.EQ.1)GOTO 100                                               18700
      IF(ERR) GOTO 250                                                     18800
      ASSIGN 600 TO ITT                                                    18900
      IF(NANGLE.GT.1) ASSIGN 500 TO ITT                                    19000
      CALL PREP1                                                           19100
      IF(ERR) GOTO 250                                                     19200
      CALL ZCLOCK(TIME02)                                                  19300
      TPREP = TIME02 - TIME01                                              19400
      WRITE(5,908) TPREP                                                   19500
 300  CALL ZCLOCK(TIME1)                                                   19600
      CALL PREP2                                                           19700
      IF(ERR) GOTO 250                                                     19800
      CALL ZCLOCK(T1)                                                      19900
      CALL INTG                                                            20000
      IF(ERR) GOTO 250                                                     20100
      CALL ZCLOCK(T2)                                                      20200
      TINTG = T2 - T1                                                      20300
      WRITE(5,904) TINTG                                                   20400
      IF(OURHOB.EQ.0.AND.OURHOC.EQ.0.AND.OURHOX.EQ.0)                      20500
     1GOTO 420                                                             20600
      CALL TENS1                                                           20700
      IF(ERR) GOTO 250                                                     20800
      CALL ZCLOCK(T3)                                                      20900
      TTENS = T3 - T2                                                      21000
 420  CALL SIGM                                                            21100
      IF(ERR) GOTO 250                                                     21200
      CALL ZCLOCK(TIME2)                                                   21300
      IF(OURHOB.NE.0.OR.OURHOC.NE.0.OR.OURHOX.NE.0)                        21400
     1WRITE(5,906) TTENS                                                   21500
      GOTO ITT, (600,500)                                                  21600
C                                                                          21700
C     ESTIMATION OF TOTAL INTEGRATION TIME AFTER FIRST ANGLE               21800
 500  ESTTIM = (TIME2-TIME1)*FLOAT(NANGLE)                                 21900
      WRITE(5,900) ESTTIM                                                  22000
      ASSIGN 600 TO ITT                                                    22100
      IF(ESTTIM.LT.FLOAT(NTIME))GOTO 600                                   22200
      WRITE(5,902)                                                         22300
      GOTO 250                                                             22400
 600  IF(INTEND.EQ.1) GOTO 250                                             22500
      IF(SIMP.NE.0.0) GOTO 300                                             22600
      GOTO 250                                                             22700
C                                                                          22800
 900  FORMAT(21H0ESTIMATED JOB TIME = F10.3,4H SEC)                        22900
 902  FORMAT(44H0EXCEEDS ALLOWED TIME - EXECUTION TERMINATED)              23000
 904  FORMAT(23H0TIME FOR INTEGRATION =F10.3,4H SEC//)                     23100
 906  FORMAT(44H0TIME FOR CALCULATING STATISTICAL TENSORS = F10.3,         23200
     14H SEC)                                                              23300
 908  FORMAT(24H0TIME FOR PREPARATION = F10.3,4H SEC)                      23400
      END                                                                  23500
      SUBROUTINE AMPDER(AMP,AMPDOT,K)                                      23600
C                                                                          23700
C     THIS ROUTINE DETERMINES THE VALUE IR, FOR WHICH AMPDOT(IR,L) IS      23800
C     CALCULATED                                                           23900
C                                                                          24000
      COMPLEX AMP,AMPDOT                                                   24100
      INTEGER SSTART,SSTOP                                                 24200
      DIMENSION AMP(CCCC,4),AMPDOT(CCCC,4)                                 24300
      COMMON /BL3/EN(NNNN),SPIN(NNNN),ACCUR,DIPOL,ZPOL,BANDK(NNNN)         24400
      COMMON /BL5/NZ                                                       24500
      COMMON/BL12/LMX,IR1,IR2,N1,N2                                        24600
      COMMON/BL18/CAT(CCCC,3),ICATMX, ISMAX                                24700
      COMMON/BL40/SSTART(31),SSTOP(30)                                     24800
C                                                                          24900
      IF(K.EQ.2)GOTO 100                                                   25000
C     SETS INDEX COUNTER FOR ZETA-ARRAY EQUAL 0                            25100
      NZ = 0                                                               25200
 100  IF(SPIN(N1).EQ.0.)GOTO 300                                           25300
C                                                                          25400
C     IR - LOOP FOR GROUND STATE SPIN .NE. 0                               25500
      DO 200 IR=IR1,IR2                                                    25600
      N = CAT(IR,1)                                                        25700
      CALL LAISUM(AMP,AMPDOT,IR,N,LMX)                                     25800
 200  CONTINUE                                                             25900
      RETURN                                                               26000
C                                                                          26100
C     IR - LOOP FOR GROUND STATE SPIN = 0                                  26200
 300  DO 500 N=N1,N2                                                       26300
      IR = SSTART(N) - 1                                                   26400
 400  IR = IR + 1                                                          26500
      CALL LAISUM(AMP,AMPDOT,IR,N,LMX)                                     26600
      IF(CAT(IR,3).LT.-0.1)GOTO 400                                        26700
 500  CONTINUE                                                             26800
      RETURN                                                               26900
      END                                                                  27000
      SUBROUTINE CHECK                                                     27100
C                                                                          27200
C     THIS ROUTINE CHECKS THE INPUT DATA                                   27300
C                                                                          27400
      LOGICAL ERR                                                          27500
      INTEGER OURHOB,OURHOC,OURHOX                                         27600
      COMMON /BL2/NGMAX,NSTART(10),NSTOP(10),MASTER(10,10),EMMAX(10)       27700
      COMMON /BL3/EN(NNNN),SPIN(NNNN),ACCUR,DIPOL,ZPOL,BANDK(NNNN)         27800
      COMMON /BL9/IEXCIT(10),IEXNUM,LMAX(2),IRSTA(2),IRSTO(2),IP           27900
      COMMON/BL19/ERR                                                      28000
      COMMON/BL23/IDN17,IDN18,IDN19,IDN25                                  28100
      COMMON/BL31/OURHOB,OURHOC,OURHOX                                     28200
      COMMON/BL32/NMAX,NDIM                                                28300
      COMMON/BL50/OUXI,OUPSI,NCM,EMMAX1                                    28400
C                                                                          28500
C     CHECK IF PROJECTILE GROUPS AND TARGET GROUPS ARE SEPARATED           28600
C     THIS PART ALSO DETERMINES THE VARIABLES  IEXNUM  AND  IP             28700
      J = IEXCIT(1)                                                        28800
      IEXNUM = 1                                                           28900
      IP = 1                                                               29000
      DO 100 I=1,NGMAX                                                     29100
      IF(IEXCIT(I).EQ.J)GOTO 100                                           29200
      IP = I                                                               29300
      IEXNUM = IEXNUM + 1                                                  29400
      J = IEXCIT(I)                                                        29500
 100  CONTINUE                                                             29600
      IF(IEXNUM.LE.2)GOTO 120                                              29700
      ERR = .TRUE.                                                         29800
      WRITE(5,900)                                                         29900
      RETURN                                                               30000
C                                                                          30100
C     CHECK OF GROUND STATE SPIN, GROUND STATE ENERGY AND EMMAX            30200
 120  DO 140 K=1,IEXNUM                                                    30300
      J = 1                                                                30400
      IF(K.EQ.2) J=NSTART(IP)                                              30500
      IF(SPIN(J).LE.3.5)GOTO 125                                           30600
      ERR = .TRUE.                                                         30700
      WRITE(5,902) J,SPIN(J)                                               30800
 125  IF(EN(J).EQ.0.)GOTO 130                                              30900
      ERR = .TRUE.                                                         31000
      WRITE(5,904) J,EN(J)                                                 31100
 130  DO 139 I=1,NGMAX                                                     31200
      NGRD = 1                                                             31300
      IF(IEXNUM.EQ.2.AND.I.GE.IP) NGRD=NSTART(IP)                          31400
      IF(EMMAX(I).GE.SPIN(NGRD))GOTO 139                                   31500
      ERR = .TRUE.                                                         31600
      WRITE(5,906) I,EMMAX(I)                                              31700
 139  CONTINUE                                                             31800
 140  CONTINUE                                                             31900
C                                                                          32000
C     CHECK OF NSTART- AND NSTOP-ARRAY                                     32100
      IF(NGMAX.EQ.1)GOTO 165                                               32200
      DO 160 I=2,NGMAX                                                     32300
      J = I-1                                                              32400
      IF(NSTART(I)-NSTOP(J).EQ.1)GOTO 160                                  32500
      ERR = .TRUE.                                                         32600
      WRITE(5,908) J,NSTOP(J),I,NSTART(I)                                  32700
 160  CONTINUE                                                             32800
C                                                                          32900
C     CHECK OF NCM AND NMAX                                                33000
 165  IF(NGMAX.LE.NMAX)GOTO 170                                            33100
      ERR = .TRUE.                                                         33200
      WRITE(5,910) NMAX,NGMAX                                              33300
 170  IF(NCM.LE.NMAX)GOTO 180                                              33400
      WRITE(5,912) NMAX,NCM                                                33500
      NCM = 2                                                              33600
C                                                                          33700
C     CHECK IF INPUT DATA CONTAIN EP,ZP,A1,ZT,A2                           33800
 180  IF(IDN17.EQ.1.AND.IDN18.EQ.1.AND.IDN19.EQ.1)GOTO 190                 33900
      ERR = .TRUE.                                                         34000
      WRITE(5,914)                                                         34100
C                                                                          34200
C     CHECK IF ALL SPINS ARE EITHER INTEGER OR HALF INTEGER                34300
 190  DO 220 K=1,IEXNUM                                                    34400
      N1 = 1                                                               34500
      N2 = NMAX                                                            34600
      IF(IEXNUM.EQ.2) N2=NSTART(IP)-1                                      34700
      IF(K.EQ.1)GOTO 195                                                   34800
      N1 = NSTART(IP)                                                      34900
      N2 = NMAX                                                            35000
 195  DO 210 N=N1,N2                                                       35100
      S = SPIN(N)-SPIN(N1)                                                 35200
      J = S                                                                35300
      DIFF = S - FLOAT(J)                                                  35400
      IF(DIFF.LT.0.01)GOTO 210                                             35500
      ERR = .TRUE.                                                         35600
      WRITE(5,916) N1,SPIN(N1),N,SPIN(N)                                   35700
 210  CONTINUE                                                             35800
 220  CONTINUE                                                             35900
      RETURN                                                               36000
C                                                                          36100
 900  FORMAT(64H0ERROR  -  TARGET GROUPS AND PROJECTILE GROUPS ARE NOT S   36200
     1EPARATED//21H0EXECUTION TERMINATED)                                  36300
 902  FORMAT(6H0SPIN(I2,4H) = F4.1/                                        36400
     147H0ERROR  -  GROUND STATE SPIN MUST BE  .LT.  3.5)                  36500
 904  FORMAT(4H0EN(I2,4H) = F8.4,4H MEV/                                   36600
     149H0ERROR  -  GROUND STATE ENERGY MUST BE  .EQ.  0.0)                36700
 906  FORMAT(7H0EMMAX(I2,4H) = F4.1/                                       36800
     149H0ERROR  -  EMMAX MUST BE  .GE.  GROUND STATE SPIN)                36900
 908  FORMAT(7H0NSTOP(I2,4H) = I2,3X,7HNSTART(I2,4H) = I2/                 37000
     151H0ERROR  -  NSTART(I+1) - NSTOP(I)  MUST BE  .EQ.  1)              37100
 910  FORMAT(8H0NMAX = I2,4X8HNGMAX = I2/                                  37200
     136H0ERROR  -  NGMAX MUST BE  .LE.  NMAX)                             37300
 912  FORMAT(8H0NMAX = I2,4X,6HNCM = I2/25H0PROGRAM ASSUMES  NCM = 2)      37400
 914  FORMAT('   ERROR  -  SOME OF THE FOLLOWING INPUT QUANTITIES ARE      37500
     1 MISSING  EP, ZP, A1, ZT, A2 ')                                      37600
 916  FORMAT(6H0SPIN(I2,4H) = F4.1,4X,5HSPIN(I2,4H) = F4.1/9H0ERROR  -     37700
     159H  FOR EACH NUCLEUS ALL SPINS MUST BE EITHER INTEGER OR HALF       37800
     28H INTEGER)                                                          37900
      END                                                                  38000
      SUBROUTINE CMLAB(A1,A2,EP,EN,NMAX,TLBDG,TCMDG,ZLBDG,R3,R4,EPP,ETP)   38100
C                                                                          38200
C     THIS ROUTINE TRANSFORMS FROM LAB-SYSTEM TO CM-SYSTEM                 38300
C     THE FORMULAS ARE FROM ALDER/WINTHER, "ELECTROMAGNETIC EXCITATION"    38400
C     NOTH-HOLLAND PUBL.CO., AMSTERDAM (1975), P. 263 - 268                38500
C                                                                          38600
      LOGICAL ERR                                                          38700
      DIMENSION EN(NNNN),R3(NNNN),R4(NNNN),EPP(NNNN),ETP(NNNN)             38800
      DIMENSION TCMDG(NNNN),TCMRAD(NNNN),ZLBDG(NNNN),ZLBRAD(NNNN)          38900
      COMMON/BL19/ERR                                                      39000
C                                                                          39100
      TLBRAD = TLBDG/57.2957795                                            39200
      ARED = 1.0 + A1/A2                                                   39300
      EMAX = EP/ARED                                                       39400
      DO 200 N=1,NMAX                                                      39500
      IF (EN(N).LT.EMAX) GOTO 100                                          39600
C     STOPS CALCULATION IF EXCITATION ENERGY IS KINEMATICALLY FORBIDDEN    39700
      WRITE(5,900) EMAX                                                    39800
      GOTO 850                                                             39900
 100  EPMIN = EP - EN(N)*ARED                                              40000
      TAUP = SQRT(EP/EPMIN)                                                40100
      TAU = TAUP * A1/A2                                                   40200
      IF (TAU.LE.1.0) GOTO 120                                             40300
      TMXDG = ASIN(1.0/TAU) * 57.2957795                                   40400
      IF (TMXDG.GE.TLBDG) GOTO 120                                         40500
C     STOPS CALCULATION IF LABORATORY SCATTERING ANGLE IS KINEMATICALLY    40600
C     FORBIDDEN                                                            40700
      WRITE(5,902) TMXDG                                                   40800
      GOTO 850                                                             40900
C                                                                          41000
C     COMPUTATION OF CM SCATTERING ANGLE TCMDG(N)                          41100
 120  TCMRAD(N) = TLBRAD + ASIN(TAU*SIN(TLBRAD))                           41200
      TCMDG(N) = TCMRAD(N)*57.2957795                                      41300
      IF (TAU.LE.1.0) GOTO 140                                             41400
      WRITE(5,904) TCMDG(N),N                                              41500
C     TCMDG(N) = 180.+2.*TLBDG-TCMDG(N)                                    41600
      TCMRAD(N) = TCMDG(N)/57.2957795                                      41700
C                                                                          41800
C     COMPUTATION OF LAB RECOIL ANGLE ZLBDG(N)                             41900
 140  ZCMDG = 180.0 - TCMDG(N)                                             42000
      ZCMRAD = ZCMDG/57.2957795                                            42100
      ZLBRAD(N) = ATAN(SIN(ZCMRAD)/(COS(ZCMRAD)+TAUP))                     42200
      ZLBDG(N) = ZLBRAD(N)*57.2957795                                      42300
      IF(ABS(TCMDG(N)-180.).LT.1.E-10)GOTO 160                             42400
C                                                                          42500
C     COMPUTATIION OF SOLID ANGLE RATIOS  R3(N)  AND  R4(N)                42600
      R3(N)=(SIN(TLBRAD)/SIN(TCMRAD(N)))**2*ABS(COS(TCMRAD(N)-TLBRAD))     42700
      R4(N)=(SIN(ZLBRAD(N))/SIN(ZCMRAD))**2*ABS(COS(ZCMRAD-ZLBRAD(N)))     42800
      GOTO 180                                                             42900
 160  R3(N) = 1.0/(1.0-TAU)**2                                             43000
      R4(N) = 1.0/(1.0+TAUP)**2                                            43100
 180  R3(N) = 1.0/R3(N)                                                    43200
      R4(N) = 1.0/R4(N)                                                    43300
C                                                                          43400
C     COMPUTATION OF LAB ENERGIES OF SCATTERED PARTICLES                   43500
      EPP(N)=(A2/(A1+A2))**2*(1.0+TAU*TAU+2*TAU*COS(TCMRAD(N)))*EPMIN      43600
      ETP(N)=A1*A2/(A1+A2)**2*(1.+TAUP*TAUP+2.*TAUP*COS(ZCMRAD))*EPMIN     43700
 200  CONTINUE                                                             43800
      RETURN                                                               43900
C                                                                          44000
C     ERROR EXIT                                                           44100
C                                                                          44200
 850  ERR = .TRUE.                                                         44300
      RETURN                                                               44400
C                                                                          44500
 900  FORMAT(45H0ERROR- MAXIMUM ALLOWED EXCITATION ENERGY IS F8.4,         44600
     14H MEV)                                                              44700
 902  FORMAT(44H0ERROR- MAXIMUM ALLOWED SCATTERING ANGLE IS F7.2,          44800
     18H DEGREES)                                                          44900
 904  FORMAT(39H0SECOND POSSIBLE CM SCATTERING ANGLE IS F7.2,              45000
     117H DEGREES FOR N = I3)                                              45100
      END                                                                  45200
      FUNCTION DJMM(BETA,RJ,RM,RMP)                                        45300
C                                                                          45400
C     THIS ROUTINE COMPUTES THE D-FUNCTION ACCORDING TO EQ. (I.4.6)        45500
C                                                                          45600
      COMMON /BL8/ASQRT,B(200),IEX(200)                                    45700
C                                                                          45800
      J=2.001*RJ                                                           45900
      M=2.001*RM                                                           46000
      MP=2.001*RMP                                                         46100
      IF(J.LT.IABS(M).OR.J.LT.IABS(MP).OR.J.LT.0) GOTO 200                 46200
      BE=BETA/57.295779                                                    46300
      CB=COS(BE/2.0)                                                       46400
      SB=SIN(BE/2.0)                                                       46500
      JA=(J+MP)/2+1                                                        46600
      JB=(J-MP)/2+1                                                        46700
      JC=(J+M)/2+1                                                         46800
      JD=(J-M)/2+1                                                         46900
      B1 = B(JA)*B(JB)*B(JC)*B(JD)                                         47000
      IEX1 = IEX(JA)+IEX(JB)+IEX(JC)+IEX(JD)                               47100
      MINSIG=0                                                             47200
      IF(M+MP.LT.0) MINSIG=-M-MP                                           47300
      MAXSIG=J-M                                                           47400
      IF(M.LT.MP) MAXSIG=J-MP                                              47500
      ISIG=MINSIG                                                          47600
      DJMM = 0.0                                                           47700
C                                                                          47800
C     SUMMATION OVER SIGMA                                                 47900
 100  JA=ISIG/2+1                                                          48000
      JB=(J-M-ISIG)/2+1                                                    48100
      JC=(J-MP-ISIG)/2+1                                                   48200
      JD=(M+MP+ISIG)/2+1                                                   48300
      FASE=(-1.0)**((J-MP-ISIG)/2)                                         48400
      IC=ISIG+(M+MP)/2                                                     48500
      IS=J-ISIG-(M+MP)/2                                                   48600
      B2 = B(JA)*B(JB)*B(JC)*B(JD)                                         48700
      IEX2 = IEX(JA)+IEX(JB)+IEX(JC)+IEX(JD)                               48800
      DJMM=DJMM+FASE*CB**IC*SB**IS*SQRT(B1/(B2*B2)*10.0**(IEX1-2*IEX2))    48900
      ISIG=ISIG+2                                                          49000
      IF(ISIG.LE.MAXSIG) GOTO 100                                          49100
      RETURN                                                               49200
C                                                                          49300
C     EXIT IF ERROR IN DJMM - ARGUMENTS                                    49400
C                                                                          49500
 200  WRITE(5,900)                                                         49600
      STOP                                                                 49700
C                                                                          49800
 900  FORMAT(25H0 ERROR IN DJMM-ARGUMENTS)                                 49900
      END                                                                  50000
      SUBROUTINE ERROR(I)                                                  50100
C                                                                          50200
C     THIS ROUTINE PRINTS ERROR MESSAGES AND SETS LOGICAL VARIABLE         50300
C     ERR TO .TRUE.                                                        50400
C                                                                          50500
      LOGICAL ERR                                                          50600
      COMMON/BL19/ERR                                                      50700
      COMMON/BL32/NMAX,NDIM                                                50800
C                                                                          50900
      ERR = .TRUE.                                                         51000
C                                                                          51100
C     PRINTS ERROR MESSAGE ACCORDING TO VALUE OF I                         51200
C                                                                          51300
      GOTO(101,102,103,104,105,106,107,108,109,110,111,112,113,114,        51400
     1 115,116,117,118,119,120,121,122),I                                  51500
 101  WRITE(5,901) NDIM                                                    51600
      RETURN                                                               51700
 102  WRITE(5,902)                                                         51800
      RETURN                                                               51900
 103  WRITE(5,903)                                                         52000
      RETURN                                                               52100
 104  WRITE(5,904)                                                         52200
      RETURN                                                               52300
 105  WRITE(5,905)                                                         52400
      RETURN                                                               52500
 106  WRITE(5,906)                                                         52600
      RETURN                                                               52700
 107  WRITE(5,907)                                                         52800
      RETURN                                                               52900
 108  WRITE(5,908)                                                         53000
      RETURN                                                               53100
 109  WRITE(5,909) NDIM                                                    53200
      RETURN                                                               53300
 110  WRITE(5,910)                                                         53400
      RETURN                                                               53500
 111  WRITE(5,911)                                                         53600
      RETURN                                                               53700
 112  WRITE(5,912)                                                         53800
      RETURN                                                               53900
 113  WRITE(5,913)                                                         54000
      RETURN                                                               54100
 114  WRITE(5,914)                                                         54200
      RETURN                                                               54300
 115  WRITE(5,915)                                                         54400
      RETURN                                                               54500
 116  WRITE(5,916) NDIM                                                    54600
      RETURN                                                               54700
 117  WRITE(5,917)                                                         54800
      RETURN                                                               54900
 118  WRITE(5,918)                                                         55000
      RETURN                                                               55100
 119  WRITE(5,919)                                                         55200
      RETURN                                                               55300
 120  WRITE(5,920)                                                         55400
      RETURN                                                               55500
 121  WRITE(5,921)                                                         55600
      RETURN                                                               55700
 122  WRITE(5,922)                                                         55800
      RETURN                                                               55900
C                                                                          56000
 901  FORMAT(46H INPUT ERROR  --  NMAX MUST BE .GE. 2 AND .LE. I2/)        56100
 902  FORMAT(54H INPUT ERROR  --  XIMAX MUST BE .GT. 0.0 AND .LT. 99.0/)   56200
 903  FORMAT(41H INPUT ERROR  --  EMMAX1 MUST BE .GE. 0.0/)                56300
 904  FORMAT(40H INPUT ERROR  --  ACCUR MUST BE .GT. 0.0/)                 56400
 905  FORMAT(53H INPUT ERROR  --  CHARGE NUMBER AND MASS NUMBER MUST       56500
     19HBE .GT. 0/)                                                        56600
 906  FORMAT(53H INPUT ERROR  --  CHARGE NUMBER MUST BE .LE. MASS NUM      56700
     13HBER/)                                                              56800
 907  FORMAT(53H INPUT ERROR  --  ENERGY OF PROJECTILE MUST BE .GT. 0/)    56900
 908  FORMAT(53H INPUT ERROR  --  SCATTERING ANGLE MUST BE BETWEEN 0       57000
     115HAND 180 DEGREES/)                                                 57100
 909  FORMAT(53H INPUT ERROR  --  LEVEL INDEX MUST BE .GE. 1 AND .LE.      57200
     11H I2/)                                                              57300
 910  FORMAT(54H INPUT ERROR  --  SPIN QUANTUM NUMBER MUST BE .GE. 0.0/)   57400
 911  FORMAT(41H INPUT ERROR  --  ENERGY MUST BE .GE. 0.0/)                57500
 912  FORMAT(51H INPUT ERROR  --  IPAR(N) MUST BE EQUAL 1 FOR EVEN         57600
     128HPARITY AND -1 FOR ODD PARITY/)                                    57700
 913  FORMAT(53H INPUT ERROR  --  MULTIPOLARITY MUST BE .GE. 1 AND .L      57800
     14HE. 6/)                                                             57900
 914  FORMAT(50H INPUT ERROR  --  NGMAX MUST BE .GE. 1 AND .LE. 10/)       58000
 915  FORMAT(53H INPUT ERROR  --  GROUP INDEX MUST BE .GE. 1 AND .LE.      58100
     13H 10/)                                                              58200
 916  FORMAT(53H INPUT ERROR  --  NSTART(I) AND NSTOP(I) MUST BE .GE.      58300
     112H 1 AND .LE. I2/)                                                  58400
 917  FORMAT(53H INPUT ERROR  --  NSTOP(I) MUST BE GREATER THAN NSTAR      58500
     14HT(I)/)                                                             58600
 918  FORMAT(43H INPUT ERROR  --  EMMAX(I) MUST BE .GE. 0.0/)              58700
 919  FORMAT(48H INPUT ERROR  --  IEXCIT(I) MUST BE EQUAL 1 (FOR           58800
     158H PROJECTILE EXCITATION) OR EQUAL 2 (FOR TARGET EXCITATION)/)      58900
 920  FORMAT(42H INPUT ERROR  --  MAXIMUM NUMBER OF LEVELS                 59000
     120H IN EACH GROUP IS 11)                                             59100
 921  FORMAT(38H INPUT ERROR  --  NTIME MUST BE .GT. 0)                    59200
 922  FORMAT(54H INPUT ERROR  --  DTHETA MUST BE POSITIVE AND .LT. 180)    59300
      END                                                                  59400
      SUBROUTINE INTG                                                      59500
C                                                                          59600
C     THIS ROUTINE INTEGRATES THE COUPLED DIFFERENTIAL EQUATIONS           59700
C                                                                          59800
      INTEGER SSTART,SSTOP,OUAMP,OUPROW                                    59900
      COMPLEX AMP,AMPDOT,AMPP,F,Q1,RK1,RK2,RK3                             60000
      LOGICAL ERR                                                          60100
      DIMENSION PROB(CCCC)                                                 60200
      DIMENSION Q1(CCCC,4),AMPP(CCCC,4),F(CCCC,4,4),AMPDOT(CCCC,4)         60300
      COMMON /BL2/NGMAX,NSTART(10),NSTOP(10),MASTER(10,10),EMMAX(10)       60400
      COMMON /BL3/EN(NNNN),SPIN(NNNN),ACCUR,DIPOL,ZPOL,BANDK(NNNN)         60500
      COMMON /BL9/IEXCIT(10),IEXNUM,LMAX(2),IRSTA(2),IRSTO(2),IP           60600
      COMMON/BL10/IPAR(NNNN),IFAC(NNNN)                                    60700
      COMMON/BL12/LMX,IR1,IR2,N1,N2                                        60800
      COMMON/BL18/CAT(CCCC,3),ICATMX, ISMAX                                60900
      COMMON/BL19/ERR                                                      61000
      COMMON/BL24/AMP(CCCC,4)                                              61100
      COMMON/BL32/NMAX,NDIM                                                61200
      COMMON/BL33/NTIME                                                    61300
      COMMON/BL39/UP,DW,ISTEP,D2W                                          61400
      COMMON/BL40/SSTART(31),SSTOP(30)                                     61500
      COMMON/BL49/OUAMP,OUPROW,INTERV,INTIN                                61600
      COMMON/BL53/P(NNNN)                                                  61700
      DOUBLEPRECISION FZR,FZI,FZ        !*************************              
C                                                                          61800
C     INITIAL CONDITIONS FOR INTEGRATION                                   61900
      KAST = 0                                                             62000
      ISTEPS = 0                                                           62100
C     ISTEPS COUNTS ACTUAL NUMBER OF STEPS                                 62200
      INTERV = INTIN                                                       62300
      W = -UP                                                              62400
C     INITIAL VALUES OF AMPLITUDES  AMP(W=-UP)                             62500
      DO 110 K=1,IEXNUM                                                    62600
      LMX = LMAX(K)                                                        62700
      IR1 = IRSTA(K)                                                       62800
      IR2 = IRSTO(K)                                                       62900
      DO 105 L=1,LMX                                                       63000
      DO 100 IR=IR1,IR2                                                    63100
      AMP(IR,L) = (0.0,0.0)                                                63200
      IF(L.EQ.IR-IR1+1) AMP(IR,L) = (1.0,0.0)                              63300
 100  CONTINUE                                                             63400
 105  CONTINUE                                                             63500
 110  CONTINUE                                                             63600
C     CONSTANTS USED IN RUNGE-KUTTA-GILL EQUATIONS                         63700
      B1=0.5857864                                                         63800
      C1=0.1213204                                                         63900
      B2=3.4142136                                                         64000
      C2=-4.1213204                                                        64100
C                                                                          64200
C     INITIAL TIMING OF INTEGRATION                                        64300
      CALL ZCLOCK(TIME1)                                                   64400
      ASSIGN 360 TO ITI                                                    64500
      IF(OUPROW.EQ.0) GOTO 120                                             64600
C                                                                          64700
C     HEADING FOR PRINT-OUT OF PROBABILITIES                               64800
      WRITE(5,900)(N,N=1,NMAX)                                             64900
      P(1)=1.0                                                             65000
      DO 115 N=2,NMAX                                                      65100
 115  P(N) = 0.0                                                           65200
      N1 = NSTART(IP)                                                      65300
      IF(IEXNUM.EQ.2) P(N1)=1.0                                            65400
      PTOT = P(1) + P(N1)                                                  65500
      WRITE(5,902) W,(P(N),N=1,NMAX)                                       65600
C                                                                          65700
C     ******************************************************************   65800
C     THE RUNGE-KUTTA-GILL INTEGRATION PROCEDURE                           65900
C     ******************************************************************   66000
C                                                                          66100
 120  CONTINUE                                                             66200
C     COMPUTATION OF START VALUES OF THE DERIVATIVES                       66300
      CALL Q(W)                                                            66400
      DO 145  K=1,IEXNUM                                                   66500
      CALL PREP3(K)                                                        66600
      CALL AMPDER(AMP,AMPDOT,K)                                            66700
      IF(SPIN(N1).EQ.0.)GOTO 130                                           66800
      DO 125 L=1,LMX                                                       66900
      DO 125 IR=IR1,IR2                                                    67000
      F(IR,L,1) = AMPDOT(IR,L)                                             67100
 125  CONTINUE                                                             67200
      GOTO 145                                                             67300
 130  DO 140 L=1,LMX                                                       67400
      DO 140 N=N1,N2                                                       67500
      IR = SSTART(N) - 1                                                   67600
 135  IR = IR + 1                                                          67700
      F(IR,L,1) = AMPDOT(IR,L)                                             67800
      IF(CAT(IR,3).LT.-0.1)GOTO 135                                        67900
 140  CONTINUE                                                             68000
 145  CONTINUE                                                             68100
C     ******************************************************************   68200
C     COMPUTATION OF AMPLITUDES ACCORDING TO EQUATION (II.7.3A)            68300
      DO 275 NAM=2,4                                                       68400
      DO 170 K=1,IEXNUM                                                    68500
      CALL PREP3(K)                                                        68600
      IF(SPIN(N1).EQ.0.)GOTO 155                                           68700
      DO 150 L=1,LMX                                                       68800
      DO 150 IR=IR1,IR2                                                    68900
      Q1(IR,L) = DW*AMPDOT(IR,L)                                           69000
      AMP(IR,L) = AMP(IR,L) + Q1(IR,L)                                     69100
 150  CONTINUE                                                             69200
      GOTO 170                                                             69300
 155  DO 165 L=1,LMX                                                       69400
      DO 165 N=N1,N2                                                       69500
      IR = SSTART(N) - 1                                                   69600
 160  IR = IR + 1                                                          69700
      Q1(IR,L) = DW*AMPDOT(IR,L)                                           69800
      AMP(IR,L) = AMP(IR,L) + Q1(IR,L)                                     69900
C     MAKE USE OF SYMMETRY RELATION (EQUATION (I.2.30))                    70000
      MIR = CAT(IR,3)                                                      70100
      IR1 = IR - 2*MIR                                                     70200
      AMP(IR1,L) = IFAC(N)*AMP(IR,L)                                       70300
      IF(CAT(IR,3).LT.-0.1)GOTO 160                                        70400
 165  CONTINUE                                                             70500
 170  CONTINUE                                                             70600
C     COMPUTATION OF AMPLITUDES ACCORDING TO EQUATION (II.7.3B)            70700
      W = W + DW                                                           70800
      CALL Q(W)                                                            70900
      DO 195 K=1,IEXNUM                                                    71000
      CALL PREP3(K)                                                        71100
      CALL AMPDER(AMP,AMPDOT,K)                                            71200
      IF(SPIN(N1).EQ.0.)GOTO 180                                           71300
      DO 175 L=1,LMX                                                       71400
      DO 175 IR=IR1,IR2                                                    71500
      RK1 = DW*AMPDOT(IR,L)                                                71600
      AMP(IR,L) = AMP(IR,L) + B1*(RK1-Q1(IR,L))                            71700
      Q1(IR,L) = B1*RK1 + C1*Q1(IR,L)                                      71800
 175  CONTINUE                                                             71900
      GOTO 195                                                             72000
 180  DO 190 L=1,LMX                                                       72100
      DO 190 N=N1,N2                                                       72200
      IR = SSTART(N) - 1                                                   72300
 185  IR = IR + 1                                                          72400
      RK1 = DW*AMPDOT(IR,L)                                                72500
      AMP(IR,L) = AMP(IR,L) + B1*(RK1-Q1(IR,L))                            72600
      Q1(IR,L) = B1*RK1 + C1*Q1(IR,L)                                      72700
C     MAKE USE OF SYMMETRIE RELATION (EQUATION (I.2.30))                   72800
      MIR = CAT(IR,3)                                                      72900
      IR1 = IR - 2*MIR                                                     73000
      AMP(IR1,L) = IFAC(N)*AMP(IR,L)                                       73100
      IF(FLOAT(MIR).LT.-0.1)GOTO 185                                       73200
 190  CONTINUE                                                             73300
 195  CONTINUE                                                             73400
C     ******************************************************************   73500
C     COMPUTATION OF AMPLITUDES ACCORDING TO EQUATION (II.7.3C)            73600
      DO 220 K=1,IEXNUM                                                    73700
      CALL PREP3(K)                                                        73800
      CALL AMPDER(AMP,AMPDOT,K)                                            73900
      IF(SPIN(N1).EQ.0.)GOTO 205                                           74000
      DO 200 L=1,LMX                                                       74100
      DO 200 IR=IR1,IR2                                                    74200
      RK2 = DW*AMPDOT(IR,L)                                                74300
      AMP(IR,L) = AMP(IR,L) + B2*(RK2-Q1(IR,L))                            74400
      Q1(IR,L) = B2*RK2 + C2*Q1(IR,L)                                      74500
 200  CONTINUE                                                             74600
      GOTO 220                                                             74700
 205  DO 215 L=1,LMX                                                       74800
      DO 215 N=N1,N2                                                       74900
      IR = SSTART(N) - 1                                                   75000
 210  IR = IR + 1                                                          75100
      RK2 = DW*AMPDOT(IR,L)                                                75200
      AMP(IR,L) = AMP(IR,L) + B2*(RK2-Q1(IR,L))                            75300
      Q1(IR,L) = B2*RK2 + C2*Q1(IR,L)                                      75400
C     MAKE USE OF SYMMETRIE RELATION (EQUATION (I.2.30))                   75500
      MIR = CAT(IR,3)                                                      75600
      IR1 = IR - 2*MIR                                                     75700
      AMP(IR1,L) = IFAC(N)*AMP(IR,L)                                       75800
      IF(FLOAT(MIR).LT.-0.1)GOTO 210                                       75900
 215  CONTINUE                                                             76000
 220  CONTINUE                                                             76100
C     ******************************************************************   76200
C     COMPUTATION OF AMPLITUDES ACCORDING TO EQUATION (II.7.3D)            76300
      W = W + DW                                                           76400
      CALL Q(W)                                                            76500
      DO 245 K=1,IEXNUM                                                    76600
      CALL PREP3(K)                                                        76700
      CALL AMPDER(AMP,AMPDOT,K)                                            76800
      IF(SPIN(N1).EQ.0.)GOTO 230                                           76900
      DO 225 L=1,LMX                                                       77000
      DO 225 IR=IR1,IR2                                                    77100
      RK3 = DW*AMPDOT(IR,L)                                                77200
      AMP(IR,L) = AMP(IR,L) + RK3/3.0 - 2.0*Q1(IR,L)/3.0                   77300
 225  CONTINUE                                                             77400
      GOTO 245                                                             77500
 230  DO 240 L=1,LMX                                                       77600
      DO 240 N=N1,N2                                                       77700
      IR = SSTART(N) - 1                                                   77800
 235  IR = IR + 1                                                          77900
      RK3 = DW*AMPDOT(IR,L)                                                78000
      AMP(IR,L) = AMP(IR,L) + RK3/3.0 - 2.0*Q1(IR,L)/3.0                   78100
C     MAKE USE OF SYMMETRIE RELATION (EQUATION (I.2.30))                   78200
      MIR = CAT(IR,3)                                                      78300
      IR1 = IR - 2*MIR                                                     78400
      AMP(IR1,L) = IFAC(N)*AMP(IR,L)                                       78500
      IF(FLOAT(MIR).LT.-0.1)GOTO 235                                       78600
 240  CONTINUE                                                             78700
 245  CONTINUE                                                             78800
C     ******************************************************************   78900
C     STORING OF DERIVATIVES IN ARRAY F                                    79000
      DO 270 K=1,IEXNUM                                                    79100
      CALL PREP3(K)                                                        79200
      CALL AMPDER(AMP,AMPDOT,K)                                            79300
      IF(SPIN(N1).EQ.0.)GOTO 255                                           79400
      DO 250 L=1,LMX                                                       79500
      DO 250 IR=IR1,IR2                                                    79600
      F(IR,L,NAM) = AMPDOT(IR,L)                                           79700
 250  CONTINUE                                                             79800
      GOTO 270                                                             79900
 255  DO 265 L=1,LMX                                                       80000
      DO 265 N=N1,N2                                                       80100
      IR = SSTART(N) - 1                                                   80200
 260  IR = IR + 1                                                          80300
      F(IR,L,NAM) = AMPDOT(IR,L)                                           80400
      IF(CAT(IR,3).LT.-0.1)GOTO 260                                        80500
 265  CONTINUE                                                             80600
 270  CONTINUE                                                             80700
      ISTEPS = ISTEPS + 1                                                  80800
      KAST = KAST + 1                                                      80900
 275  CONTINUE                                                             81000
C     WE NOW HAVE THE FOUR START VALUES FOR THE DERIVATIVES AND CAN        81100
C     PROCEED BY THE FASTER ADAMS-MOULTON-PROCEDURE                        81200
C                                                                          81300
C     ******************************************************************   81400
C     THE ADAMS-MOULTON-PROCEDURE                                          81500
C     ******************************************************************   81600
C                                                                          81700
 280  CONTINUE                                                             81800
C     COMPUTATION OF AMPLITUDES AMPP ACCORDING TO EQUATION (II.7.6)        81900
      DO 305 K=1,IEXNUM                                                    82000
      CALL PREP3(K)                                                        82100
      IF(SPIN(N1).EQ.0.)GOTO 290                                           82200
      DO 285 L=1,LMX                                                       82300
      DO 285 IR=IR1,IR2                                                    82400
      AMPP(IR,L) = AMP(IR,L) + DW/12.0*(55.0*F(IR,L,4)-59.0*F(IR,L,3)      82500
     1+37.0*F(IR,L,2)-9.0*F(IR,L,1))                                       82600
 285  CONTINUE                                                             82700
      GOTO 305                                                             82800
 290  DO 300 L=1,LMX                                                       82900
      DO 300 N=N1,N2                                                       83000
      IR = SSTART(N) - 1                                                   83100
 295  IR = IR + 1                                                          83200
      AMPP(IR,L) = AMP(IR,L) +DW/12.0*(55.0*F(IR,L,4)-59.0*F(IR,L,3)       83300
     1+37.0*F(IR,L,2)-9.0*F(IR,L,1))                                       83400
C     USING OF THE SYMMETRY RELATION EQUATION (I.2.30)                    83500 
      MIR = CAT(IR,3)                                                      83600
      IR1 = IR - 2*MIR                                                     83700
      AMPP(IR1,L) = IFAC(N)*AMPP(IR,L)                                     83800
      IF(FLOAT(MIR).LT.-0.1)GOTO 295                                       83900
 300  CONTINUE                                                             84000
 305  CONTINUE                                                             84100
C     ******************************************************************   84200
C     COMPUTATION OF AMPLITUDES AMP ACCORDING TO EQUATION (II.7.7)         84300
      W = W + DW + DW                                                      84400
      KAST = KAST + 1                                                      84500
      ISTEPS = ISTEPS + 1                                                  84600
      CALL Q(W)                                                            84700
      DO 330 K=1,IEXNUM                                                    84800
      CALL PREP3(K)                                                        84900
      CALL AMPDER(AMPP,AMPDOT,K)                                           85000
      IF(SPIN(N1).EQ.0.)GOTO 315                                           85100
      DO 310 L=1,LMX                                                       85200
      DO 310 IR=IR1,IR2                                                    85300
      AMP(IR,L) = AMP(IR,L) + DW/12.0*(9.0*AMPDOT(IR,L)+19.0*F(IR,L,4)     85400
     1-5.0*F(IR,L,3)+F(IR,L,2))                                            85500
 310  CONTINUE                                                             85600
      GOTO 330                                                             85700
 315  DO 325 L=1,LMX                                                       85800
      DO 325 N=N1,N2                                                       85900
      IR = SSTART(N) - 1                                                   86000
 320  IR = IR + 1                                                          86100
      AMP(IR,L) = AMP(IR,L) + DW/12.0*(9.0*AMPDOT(IR,L)+19.0*F(IR,L,4)     86200
     1-5.0*F(IR,L,3)+F(IR,L,2))                                            86300
C     MAKE USE OF SYMMETRIE RELATION (EQUATION (I.2.30))                   86400
      MIR = CAT(IR,3)                                                      86500
      IR1 = IR - 2*MIR                                                     86600
      AMP(IR1,L) = IFAC(N)*AMP(IR,L)                                       86700
      IF(FLOAT(MIR).LT.-0.1)GOTO 320                                       86800
 325  CONTINUE                                                             86900
 330  CONTINUE                                                             87000
C     ******************************************************************   87100
      DO 355 K=1,IEXNUM                                                    87200
      CALL PREP3(K)                                                        87300
      CALL AMPDER(AMP,AMPDOT,K)                                            87400
      IF(SPIN(N1).EQ.0.)GOTO 340                                           87500
      DO 335 L=1,LMX                                                       87600
      DO 335 IR=IR1,IR2                                                    87700
      F(IR,L,1) = F(IR,L,2)                                                87800
      F(IR,L,2) = F(IR,L,3)                                                87900
      F(IR,L,3) = F(IR,L,4)                                                88000
      F(IR,L,4) = AMPDOT(IR,L)                                             88100
 335  CONTINUE                                                             88200
      GOTO 355                                                             88300
 340  DO 350 L=1,LMX                                                       88400
      DO 350 N=N1,N2                                                       88500
      IR = SSTART(N) - 1                                                   88600
 345  IR = IR + 1                                                          88700
      F(IR,L,1) = F(IR,L,2)                                                88800
      F(IR,L,2) = F(IR,L,3)                                                88900
      F(IR,L,3) = F(IR,L,4)                                                89000
      F(IR,L,4) = AMPDOT(IR,L)                                             89100
      IF(CAT(IR,3).LT.-0.1)GOTO 345                                        89200
 350  CONTINUE                                                             89300
 355  CONTINUE                                                             89400
C     ******************************************************************   89500
C     ******************************************************************   89600
      IF(ISTEPS.LT.10)GOTO 370                                             89700
C                                                                          89800
C     TIMING OF INTEGRATION AFTER FIRST TEN STEPS                          89900
      GOTO ITI,(360,370)                                                   90000
 360  CALL ZCLOCK(TIME2)                                                   90100
      TLAPSE = TIME2-TIME1                                                 90200
      ESTTIM = TLAPSE*FLOAT(ISTEP)/10.                                     90300
      WRITE(5,926) ESTTIM                                                  90400
      ASSIGN 370 TO ITI                                                    90500
      IF(ESTTIM.LT.FLOAT(NTIME))GOTO 370                                   90600
      WRITE(5,928)                                                         90700
      GOTO 850                                                             90800
 370  IF(W+DW.GT.UP) GOTO 405                                              90900
C                                                                          91000
C     ACCURACY CONTROL                                                     91100
      IF(KAST.LT.INTERV) GOTO 280                                          91200
C     FIND LARGEST  AMPP - AMP                                             91300
      FF=0.0                                                               91400
      DO 390 K=1,IEXNUM                                                    91500
      LMX = LMAX(K)                                                        91600
      IR1 = IRSTA(K)                                                       91700
      IR2 = IRSTO(K)                                                       91800
      DO 385 L=1,LMX                                                       91900
      DO 380 IR=IR1,IR2                                                    92000
      FZR = REAL(AMPP(IR,L)) - REAL(AMP(IR,L))                             92100
      FZI = AIMAG(AMPP(IR,L)) - AIMAG(AMP(IR,L))                           92200
      FZ=(SQRT(FZR*FZR+FZI*FZI))/14.0                                      92300
      IF(FZ.GT.FF) FF=FZ                                                   92400
 380  CONTINUE                                                             92500
 385  CONTINUE                                                             92600
 390  CONTINUE                                                             92700
      ACC050 = ACCUR/50.0                                                  92800
      IF (FF .LT. ACC050) GOTO 395                                         92900
      IF (FF .GT. ACCUR) GOTO 400                                          93000
      GOTO 405                                                             93100
 395  DW=2.0*DW                                                            93200
      INTERV=INTERV/2                                                      93300
      IF(INTERV.EQ.0) INTERV = 1                                           93400
      D2W=2.*DW                                                            93500
      WRITE(5,904) W,D2W                                                   93600
      GOTO 405                                                             93700
 400  CONTINUE                                                             93800
      DW=DW/2.0                                                            93900
      INTERV=2*INTERV                                                      94000
      D2W=2.*DW                                                            94100
      WRITE(5,906) W,D2W                                                   94200
      KAST = 0                                                             94300
      GOTO 120                                                             94400
C                                                                          94500
C     THE EXCITATION PROBABILITIES DURING INTEGRATION                      94600
 405  DO 435 K=1,IEXNUM                                                    94700
      N1 = 1                                                               94800
      IF(K.EQ.2) N1=NSTART(IP)                                             94900
      LLMAX = 2*(SPIN(N1)+1.1)                                             95000
      IR1 = IRSTA(K)                                                       95100
      IR2 = IRSTO(K)                                                       95200
      DO 430 IR=IR1,IR2                                                    95300
      PROB(IR) = 0.0                                                       95400
C     SUMMATION OVER GROUND STATE POLARIZATIONS                            95500
      L = 1                                                                95600
      DO 425 LL=2,LLMAX,2                                                  95700
      IF(LL-LLMAX)410,415,425                                              95800
 410  FAC = 2.0/(2.0*SPIN(N1)+1.0)                                         95900
      GOTO420                                                              96000
 415  FAC = 1.0/(2.0*SPIN(N1)+1.0)                                         96100
 420  FZR = REAL(AMP(IR,L))                                                96200
      FZI = AIMAG(AMP(IR,L))                                               96300
      PROB(IR) = PROB(IR) + FAC * (FZR*FZR + FZI*FZI)                      96400
      L = L + 1                                                            96500
 425  CONTINUE                                                             96600
 430  CONTINUE                                                             96700
 435  CONTINUE                                                             96800
      DO 445 N=1,NMAX                                                      96900
      P(N) = 0.0                                                           97000
      IR1 = SSTART(N)                                                      97100
      IR2 = SSTOP(N)                                                       97200
      DO 440 IR=IR1,IR2                                                    97300
 440  P(N) = P(N) + PROB(IR)                                               97400
 445  CONTINUE                                                             97500
C                                                                          97600
C     TOTAL PROBABILITY                                                    97700
      PTOTP = 0.0                                                          97800
      PTOTT = 0.0                                                          97900
      DO 475 K=1,IEXNUM                                                    98000
      IF(K.EQ.2)GOTO 450                                                   98100
      N1 = 1                                                               98200
      N2 = NMAX                                                            98300
      IF(IEXNUM.EQ.2) N2=NSTART(IP)-1                                      98400
      IGROUP = 1                                                           98500
      GOTO 455                                                             98600
 450  N1 = NSTART(IP)                                                      98700
      N2 = NMAX                                                            98800
      IGROUP = IP                                                          98900
 455  IEXC = IEXCIT(IGROUP)                                                99000
      IF(IEXC.EQ.2)GOTO 465                                                99100
      DO 460 N=N1,N2                                                       99200
 460  PTOTP = PTOTP + P(N)                                                 99300
      GOTO 475                                                             99400
 465  DO 470 N=N1,N2                                                       99500
 470  PTOTT = PTOTT + P(N)                                                 99600
 475  CONTINUE                                                             99700
C                                                                          99800
C     UNITARITY CHECK OF THE TOTAL PROBABILITY                             99900
      DO 485 K=1,IEXNUM                                                   100000
      IGROUP = 1                                                          100100
      IF(K.EQ.2) IGROUP=IP                                                100200
      IEXC = IEXCIT(IGROUP)                                               100300
      IF(IEXC.EQ.2)GOTO 480                                               100400
      ABW = ABS(PTOTP-1.0)/20.0                                           100500
      IF(ABW.LT.ACCUR)GOTO 485                                            100600
      WRITE(5,908) W,PTOTP                                                100700
      GOTO 485                                                            100800
 480  ABW = ABS(PTOTT-1.0)/20.0                                           100900
      IF(ABW.LT.ACCUR)GOTO 485                                            101000
      WRITE(5,910) W,PTOTT                                                101100
 485  CONTINUE                                                            101200
      IF(OUPROW.EQ.0)GOTO 490                                             101300
C                                                                         101400
C     PRINT-OUT OF EXCITATION PROBABILITIES FOR CURRENT VALUE OF W        101500
      PTOT = PTOTP + PTOTT                                                101600
      WRITE(5,902) W,(P(N),N=1,NMAX)                                      101700
      WRITE(5,912) PTOT                                                   101800
 490  KAST=0                                                              101900
      IF(W+DW.GT.UP) GOTO 500                                             102000
      IF(FF.LT.ACC050) GOTO 120                                           102100
      GOTO 280                                                            102200
 500  CONTINUE                                                            102300
C                                                                         102400
C     INTEGRATION COMPLETED                                               102500
C                                                                         102600
      WRITE(5,914) ISTEPS                                                 102700
      IF(OUAMP.EQ.0)GOTO 550                                              102800
C                                                                         102900
C     PRINT-OUT OF THE FINAL AMPLITUDES AMP(W=+UP)                        103000
      DO 545 K=1,IEXNUM                                                   103100
      LMX = LMAX(K)                                                       103200
      IF(K.EQ.2)GOTO 510                                                  103300
      N1 = 1                                                              103400
      IGROUP = 1                                                          103500
      GOTO 515                                                            103600
 510  N1 = NSTART(IP)                                                     103700
      IGROUP = IP                                                         103800
 515  IEXC = IEXCIT(IGROUP)                                               103900
      DO 540 L=1,LMX                                                      104000
      IR1 = SSTART(N1)                                                    104100
      IRL = IR1 + L - 1                                                   104200
      IF(IEXC.EQ.2)GOTO 520                                               104300
      WRITE(5,916) CAT(IRL,3)                                             104400
      GOTO 525                                                            104500
 520  WRITE(5,918) CAT(IRL,3)                                             104600
 525  WRITE(5,920)                                                        104700
      VAL = N1                                                            104800
      IR1 = IRSTA(K)                                                      104900
      IR2 = IRSTO(K)                                                      105000
      DO 540 IZR=IR1,IR2                                                  105100
      POP=(REAL(AMP(IZR,L)))**2+(AIMAG(AMP(IZR,L)))**2                    105200
      IF(CAT(IZR,1).LE.VAL)GOTO 530                                       105300
      VAL = CAT(IZR,1)                                                    105400
      WRITE(5,922)                                                        105500
 530  WRITE (5,924) (CAT(IZR,LC),LC=2,3), AMP(IZR,L), POP                 105600
 540  CONTINUE                                                            105700
 545  CONTINUE                                                            105800
 550  RETURN                                                              105900
C                                                                         106000
C     ERROR EXIT                                                          106100
C                                                                         106200
 850  ERR = .TRUE.                                                        106300
      RETURN                                                              106400
C                                                                         106500
 900  FORMAT(1H0,5X,1HW4X,10(3X,2HP(I2,1H)4X),/3X,                        106600
     111(3X,2HP(I2,1H)4X),/3X,11(3X,2HP(I2,1H)4X))                        106700
 902  FORMAT(1H0,F10.3,10E12.4,/3X,11E12.4,/3X,11E12.4)                   106800
 904  FORMAT(8H0AT W = F7.3,36H, STEP WIDTH WAS DOUBLED TO BE D2W= F8.5)  106900
 906  FORMAT(8H0AT W = F7.3,35H, STEP WIDTH WAS HALVED TO BE D2W= F8.5)   107000
 908  FORMAT(1H0,26H  - -  WARNING  - - AT W= ,F6.3,10X,5HPTOT            107100
     116HFOR PROJECTILE =F10.6,3X,32HERROR IN PTOT EXCEEDS 20 X ACCUR)    107200
 910  FORMAT(1H0,26H  - -  WARNING  - - AT W= ,F6.3,10X,                  107300
     117HPTOT FOR TARGET =F10.7,3X,32HERROR IN PTOT EXCEEDS 20 X ACCUR)   107400
 912  FORMAT(7H PTOT =E11.4,//)                                           107500
 914  FORMAT(34H0ACTUAL NUMBER OF STEPS, ISTEPS = I4)                     107600
 916  FORMAT(13H1INITIAL M = F4.1,20X,21HPROJECTILE AMPLITUDES)           107700
 918  FORMAT(13H1INITIAL M = F4.1,20X,17HTARGET AMPLITUDES)               107800
 920  FORMAT (1H0,1X4HSPIN2X12HMAG.QUAN.NO.1X14HREAL AMPLITUDE1X          107900
     114HIMAG AMPLITUDE3X10HPOPULATION)                                   108000
 922  FORMAT(1H0)                                                         108100
 924  FORMAT(1XF5.1,F9.1,3X3E15.4)                                        108200
 926  FORMAT(18H0ESTIMATED TIME = F10.4,4H SEC)                           108300
 928  FORMAT(45H0 EXCEEDS ALLOWED TIME - EXECUTION TERMINATED)            108400
      END                                                                 108500
      SUBROUTINE LAISUM(AMP,AMPDOT,IR,N,LMX)                              108600
C                                                                         108700
C     THIS ROUTINE COMPUTES AMPDOT(IR,L)                                  108800
C                                                                         108900
      COMPLEX AMP,AMPDOT,RC,CI                                            109000
      COMPLEX CQU1,CQU2,CQU3,CQU4,CQU5,CQU6,CQU7,CQU8,CQU9,CQU10,CQU11,   109100
     1CQU12                                                               109200
      INTEGER SSTART,SSTOP                                                109300
      DIMENSION AMP(CCCC,4),AMPDOT(CCCC,4)                                109400
      COMMON /BL1/LAMDA(12),LEAD(NNNN,NNNN,12),LDNUM(12,NNNN),LAMMAX      109500
      COMMON /BL5/NZ                                                      109600
      COMMON /BL7/ZETA( ZZZZ),NZMAX                                       109700
      COMMON/BL13/CQU1(7,NNNN,NNNN),CQU2(7,NNNN,NNNN),CQU3(7,NNNN,NNNN)   109800
      COMMON/BL14/CQU4(7,NNNN,NNNN),CQU5(7,NNNN,NNNN),CQU6(7,NNNN,NNNN)   109900
      COMMON/BL15/CQU7(7,NNNN,NNNN),CQU8(7,NNNN,NNNN),CQU9(7,NNNN,NNNN)   110000
      COMMON/BL16/CQU10(7,NNNN,NNNN),CQU11(7,NNNN,NNNN),                  110100
     1CQU12(7,NNNN,NNNN)                                                        
      COMMON/BL18/CAT(CCCC,3),ICATMX, ISMAX                               110200
      COMMON/BL40/SSTART(31),SSTOP(30)                                    110300
      DATA CI/(0.0,1.0)/                                                  110400
C                                                                         110500
      RMIR = CAT(IR,3)                                                    110600
      DO 100 L=1,LMX                                                      110700
 100  AMPDOT(IR,L) = (0.0,0.0)                                            110800
      DO 580 I1=1,LAMMAX                                                  110900
      LAM = LAMDA(I1)                                                     111000
      LA = LAM                                                            111100
      IF(LAM.GT.6) LAM=LAM-6                                              111200
      LD = LDNUM(LA,N)                                                    111300
      IF(LD.EQ.0)GOTO 580                                                 111400
      DO 560 I2=1,LD                                                      111500
      M = LEAD(N,I2,LA)                                                   111600
      ISMIN = 0                                                           111700
      IS1 = SSTART(M)                                                     111800
      ISPLUS = IFIX(RMIR-CAT(IS1,3))-LAM                                  111900
      IF(ISPLUS.GE.0)GOTO 150                                             112000
      ISMIN = ISPLUS                                                      112100
      ISPLUS = 0                                                          112200
 150  IS2 = IS1+ISPLUS-1                                                  112300
      MRANGE = 2*LAM+1+ISMIN                                              112400
      IF(IS2+MRANGE.GT.SSTOP(M)) MRANGE=SSTOP(M)-IS2                      112500
      IF(MRANGE.LE.0)GOTO 560                                             112600
      DO 540 I3=1,MRANGE                                                  112700
      IS = IS2+I3                                                         112800
      NZ = NZ + 1                                                         112900
      Z = ZETA(NZ)                                                        113000
      RMU = CAT(IS,3) - RMIR                                              113100
      MUA = ABS(RMU) + 1.1                                                113200
      DO 520 L=1,LMX                                                      113300
C     COMPUTATION OF THE PARTIAL SUM FOR MULTIPOLARITY LA                 113400
      GOTO(210,220,230,240,250,260,310,320,330,340,350,360),LA            113500
 210  RC = CQU1(MUA,N,M)*AMP(IS,L)*Z                                      113600
      GOTO 400                                                            113700
 220  RC = CQU2(MUA,N,M)*AMP(IS,L)*Z                                      113800
      GOTO 400                                                            113900
 230  RC = CQU3(MUA,N,M)*AMP(IS,L)*Z                                      114000
      GOTO 400                                                            114100
 240  RC = CQU4(MUA,N,M)*AMP(IS,L)*Z                                      114200
      GOTO 400                                                            114300
 250  RC = CQU5(MUA,N,M)*AMP(IS,L)*Z                                      114400
      GOTO 400                                                            114500
 260  RC = CQU6(MUA,N,M)*AMP(IS,L)*Z                                      114600
      GOTO 400                                                            114700
 310  RC = CQU7(MUA,N,M)*AMP(IS,L)*Z                                      114800
      IF(RMU.LT.-0.1) RC=-RC                                              114900
      GOTO 400                                                            115000
 320  RC = CQU8(MUA,N,M)*AMP(IS,L)*Z                                      115100
      IF(RMU.LT.-0.1) RC=-RC                                              115200
      GOTO 400                                                            115300
 330  RC = CQU9(MUA,N,M)*AMP(IS,L)*Z                                      115400
      IF(RMU.LT.-0.1) RC=-RC                                              115500
      GOTO 400                                                            115600
 340  RC =CQU10(MUA,N,M)*AMP(IS,L)*Z                                      115700
      IF(RMU.LT.-0.1) RC=-RC                                              115800
      GOTO 400                                                            115900
 350  RC =CQU11(MUA,N,M)*AMP(IS,L)*Z                                      116000
      IF(RMU.LT.-0.1) RC=-RC                                              116100
      GOTO 400                                                            116200
 360  RC =CQU12(MUA,N,M)*AMP(IS,L)*Z                                      116300
      IF(RMU.LT.-0.1) RC=-RC                                              116400
 400  AMPDOT(IR,L) = AMPDOT(IR,L) - CI*RC                                 116500
 520  CONTINUE                                                            116600
 540  CONTINUE                                                            116700
 560  CONTINUE                                                            116800
 580  CONTINUE                                                            116900
      RETURN                                                              117000
      END                                                                 117100
      SUBROUTINE LEXEM(A,X,N,MODE,IERR)                                   117200
                                                                          117300
C     SUBROUTINE LEXEM SCANS THE 72 CHARACTERS STORED IN THE FIELD        117400
C     A(72). IF THE CHARACTERS ARE INTERPRETED AS DIGITS OF ONE OR        117500
C     MORE NUMBERS, THESE NUMBERS ARE STORED AS AN ELEMENT OF THE         117600
C     ARRAY X. THE SCANNING IS STOPPED BY A CHARACTER INTERPRETED         117700
C     AS START OF A TEXT-STRING.                                          117800
C         A IS THE INPUT FIELD FOR 72 CHARACTERS.                         117900
C         X IS THE OUTPUT FIELD FOR THE DETECTED NUMBERS.                 118000
C         N TELLS,  HOW MANY NUMBERS HAVE BEEN DETECTED IN THE ARRAY A.   118100
C         MODE = 0 ARRAY A CONTAINS ONLY NUMBERS AND NO TEXT.            118200 
C            .GE.1 ARRAY A CONTAINS TEXT STARTING AT THE LOCATION        118300 
C                   A(MODE).                                              118400
C     BY R. ZIMMERMANN, MUNICH,1979                                       118500
                                                                          118600
      DIMENSION X(72),A(72),B(16)                                         118700
      DATA B/1H0,1H1,1H2,1H3,1H4,1H5,1H6,1H7,1H8,1H9,                     118800
     11H ,1H,,1H+,1H-,1H.,1HE/                                            118900
      IERR=0                                                              119000
      MODE=0                                                              119100
      I=0                                                                 119200
      N=0                                                                 119300
      IZ=1                                                                119400
      SIGM=1.0                                                            119500
      XM=0.0                                                              119600
      SIGE=1.0                                                            119700
      XE=0.0                                                              119800
   90 I=I+1                                                               119900
      IF (IZ.EQ.1) NM=I                                                   120000
      IF (I.GE.73) GOTO 300                                               120100
      DO 92 K=1,16                                                        120200
      IF (A(I).EQ.B(K)) GOTO 94                                           120300
   92 CONTINUE                                                            120400
      GOTO 100                                                            120500
   94 CONTINUE                                                            120600
      GOTO(1,1,1,1,1,1,1,1,1,1,2,3,4,4,5,6),K                             120700
                                                                          120800
CHARACTER IS ONE DIGIT OF A NUMBER                                        120900
1     XX=K-1                                                              121000
      GOTO(21,22,23,24,25,26),IZ                                          121100
                                                                          121200
C     START OF A NUMBER, FIRST DIGIT DETECTED                             121300
   21 NM=I                                                                121400
      SIGM=+1.0                                                           121500
      XM=XX                                                               121600
      IZ=3                                                                121700
      GOTO 90                                                             121800
                                                                          121900
C     A SIGN WAS PRECEEDING THE FIRST DIGIT                               122000
   22 XM=XX                                                               122100
      IZ=3                                                                122200
      GOTO 90                                                             122300
                                                                          122400
C     ONE OF THE OTHER DIGITS BEFORE THE DECIMAL POINT                    122500
   23 XM=10.0*XM+XX                                                       122600
      GOTO 90                                                             122700
                                                                          122800
C     ONE OF THE DIGITS AFTER THE DECIMAL POINT                           122900
   24 XME=10.0*XME                                                        123000
      XM=XM+XX/XME                                                        123100
      GOTO 90                                                             123200
                                                                          123300
C     FIRST DIGIT OF THE EXPONENT                                         123400
   25 SIGE=+1.0                                                           123500
      XE=XX                                                               123600
      IZ=6                                                                123700
      GOTO 90                                                             123800
                                                                          123900
C     ONE OF THE OTHER DIGITS OF THE EXPONENT                             124000
   26 XE=10.0*XE+XX                                                       124100
      GOTO 90                                                             124200
                                                                          124300
CHARACTER IS A BLANK                                                      124400
2     CONTINUE                                                            124500
      GOTO(90,100,200,200,100,200),IZ                                     124600
                                                                          124700
CHARACTER IS A COMMA                                                      124800
3     CONTINUE                                                            124900
      GOTO(100,100,200,200,100,200),IZ                                    125000
                                                                          125100
CHARACTER IS THE SIGN OF A NUMBER                                         125200
4     CONTINUE                                                            125300
      GOTO(40,100,100,100,45,100),IZ                                      125400
   40 IZ=2                                                                125500
      NM=I                                                                125600
      SIGM=+1.0                                                           125700
      IF (A(I).EQ.B(14)) SIGM=-1.0                                        125800
      GOTO 90                                                             125900
   45 SIGE=+1.0                                                           126000
      IF (A(I).EQ.B(14)) SIGE=-1.0                                        126100
      IZ=6                                                                126200
      GOTO 90                                                             126300
                                                                          126400
CHARACTER IS A DECIMAL POINT                                              126500
5     CONTINUE                                                            126600
      GOTO(50,51,52,100,100,100),IZ                                       126700
   50 SIGM=+1.0                                                           126800
      NM=I                                                                126900
   51 XM=0.0                                                              127000
   52 IZ=4                                                                127100
      XME=1.0                                                             127200
      GOTO 90                                                             127300
                                                                          127400
CHARACTER IS AN E, THE FOLLOWING DIGITS ARE THE EXPONENT                  127500
6     CONTINUE                                                            127600
      GOTO(100,100,60,60,100,100),IZ                                      127700
   60 IZ=5                                                                127800
      XE=0.0                                                              127900
      GOTO 90                                                             128000
                                                                          128100
C     TEXT,STOP SCANNING                                                  128200
  100 MODE=NM                                                             128300
      RETURN                                                              128400
                                                                          128500
C     DELIMITER, STORE THE DETECTED NUMBER AS THE N-TH ELEMENT            128600
C     OF THE ARRAY X                                                      128700
  200 N=N+1                                                               128800
      X(N)=SIGM*XM*10.0**(SIGE*XE)                                        128900
      NM=I                                                                129000
      IZ=1                                                                129100
      SIGM=1.0                                                            129200
      XM=0.0                                                              129300
      SIGE=1.0                                                            129400
      XE=0.0                                                              129500
      GOTO 90                                                             129600
                                                                          129700
C     ERROR EXIT                                                          129800
  300 IF (IZ.EQ.1) RETURN                                                 129900
      IERR=1                                                              130000
      GOTO 100                                                            130100
      RETURN                                                              130200
      END                                                                 130300
      SUBROUTINE LSLOOP(IR,N,NZ)                                          130400
C                                                                         130500
C     THIS ROUTINE COMPUTES THE ZETA-ARRAY                                130600
C                                                                         130700
      INTEGER SSTART,SSTOP                                                130800
      COMMON /BL1/LAMDA(12),LEAD(NNNN,NNNN,12),LDNUM(12,NNNN),LAMMAX      130900
      COMMON /BL3/EN(NNNN),SPIN(NNNN),ACCUR,DIPOL,ZPOL,BANDK(NNNN)        131000
      COMMON /BL6/PSI(NNNN,NNNN,12)                                       131100
      COMMON /BL7/ZETA( ZZZZ),NZMAX                                       131200
      COMMON/BL18/CAT(CCCC,3),ICATMX, ISMAX                               131300
      COMMON/BL40/SSTART(31),SSTOP(30)                                    131400
C                                                                         131500
C     LAMDA - LOOP  AND  IS - LOOP  FOR THE INPUT-VALUE IR                131600
      SPINR = CAT(IR,2)                                                   131700
      RMIR = CAT(IR,3)                                                    131800
      DO 400 I1=1,LAMMAX                                                  131900
      LAM = LAMDA(I1)                                                     132000
      LA = LAM                                                            132100
      IF(LAM.GT.6) LAM=LAM-6                                              132200
      RLAM = FLOAT(LAM)                                                   132300
      SSQRT = SQRT(2.0*RLAM+1.0)                                          132400
      LD = LDNUM(LA,N)                                                    132500
      IF(LD.EQ.0)GOTO 400                                                 132600
      DO 300 I2=1,LD                                                      132700
      M = LEAD(N,I2,LA)                                                   132800
      ISMIN = 0                                                           132900
      SPINS = SPIN(M)                                                     133000
      IS1 = SSTART(M)                                                     133100
      ISPLUS = IFIX(RMIR-CAT(IS1,3))-LAM                                  133200
      IF(ISPLUS.GE.0)GOTO 100                                             133300
      ISMIN = ISPLUS                                                      133400
      ISPLUS = 0                                                          133500
 100  IS2 = IS1+ISPLUS-1                                                  133600
      MRANGE = 2*LAM+1+ISMIN                                              133700
      IF(IS2+MRANGE.GT.SSTOP(M)) MRANGE = SSTOP(M)-IS2                    133800
      IF(MRANGE.LE.0)GOTO 300                                             133900
      DO 200 I3=1,MRANGE                                                  134000
      IS = IS2+I3                                                         134100
      RMIS = CAT(IS,3)                                                    134200
      G1 = -RMIS                                                          134300
      G2 = RMIS - RMIR                                                    134400
      NZ = NZ + 1                                                         134500
      IF(NZ.GT.NZMAX) GOTO 200                                            134600
      IIEX = SPINS - RMIS                                                 134700
      PHZ = (-1.0)**IIEX                                                  134800
C     WRITE(5,910)SPINS,RLAM,SPINR                                        134900
C     WRITE(5,910)G1,G2,RMIR                                              135000
      TTT=THREEJ(SPINS,G1,RLAM,G2,SPINR,RMIR)                             135100
C     WRITE(5,920)TTT                                                     135200
      ZETA(NZ)=PHZ*PSI(N,M,LA)*SSQRT*TTT                                  135300
C     WRITE(5,900) NZ,IR,IS,LA,ZETA(NZ)                                   135400
 200  CONTINUE                                                            135500
 300  CONTINUE                                                            135600
 400  CONTINUE                                                            135700
      RETURN                                                              135800
C                                                                         135900
 900  FORMAT(1H ,4HNZ =I4,6X,4HIR =I4,6X,4HIS =I4,6X,5HLAM =I2,6X,        136000
     16HZETA =F12.4)                                                      136100
  910 FORMAT(3F10.1)                                                      136200
  920 FORMAT(35X,F10.6)                                                   136300
      END                                                                 136400
      SUBROUTINE MAT                                                      136500
C                                                                         136600
C     THIS ROUTINE CHECKS THE INPUT OF MATRIX ELEMENTS MEM(N,M,LAMDA)     136700
C     THE MEM-MATRIX IS SYMMETRIZED AS NEEDED                             136800
C                                                                         136900
      REAL MEM,MTR                                                        137000
      LOGICAL ERR                                                         137100
      COMMON /BL1/LAMDA(12),LEAD(NNNN,NNNN,12),LDNUM(12,NNNN),LAMMAX      137200
      COMMON /BL2/NGMAX,NSTART(10),NSTOP(10),MASTER(10,10),EMMAX(10)      137300
      COMMON /BL3/EN(NNNN),SPIN(NNNN),ACCUR,DIPOL,ZPOL,BANDK(NNNN)        137400
      COMMON/BL10/IPAR(NNNN),IFAC(NNNN)                                   137500
      COMMON/BL19/ERR                                                     137600
      COMMON/BL20/MEM(NNNN,NNNN,12)                                       137700
      COMMON/BL32/NMAX,NDIM                                               137800
C                                                                         137900
      DO 500 N=1,NMAX                                                     138000
      DO 500 M=N,NMAX                                                     138100
      DO 400 I=1,LAMMAX                                                   138200
      LAM = LAMDA(I)                                                      138300
      LA = LAM                                                            138400
      IF(LAM.GT.6) LA=LAM-6                                               138500
      MTR=MEM(N,M,LAM)                                                    138600
      IF(MTR.EQ.0.0) GOTO 300                                             138700
C                                                                         138800
C     TRIANGLE RELATION                                                   138900
      IF(ABS(SPIN(N)-LA).LE.SPIN(M).AND.SPIN(M).LE.(SPIN(N)+LA))GOTO 140  139000
      MEM(N,M,LAM)=0.0                                                    139100
C                                                                         139200
C     PARITY CONSERVATION                                                 139300
 140  IDP = 1                                                             139400
      IF(IPAR(N).NE.IPAR(M)) IDP=-1                                       139500
      L = LA                                                              139600
      IF(LAM.GT.6) L=LA+1                                                 139700
      IF(IDP.EQ.(-1)**L)GOTO 150                                          139800
      MEM(N,M,LAM)=0.0                                                    139900
 150  IF(M.EQ.N)GOTO 400                                                  140000
C                                                                         140100
C     SYMMETRISATION OF THE MEM-MATRIX                                    140200
 300  FAC = (-1.)**IFIX(ABS(SPIN(N)-SPIN(M))+0.001)                       140300
      MEM(M,N,LAM) = FAC*MTR                                              140400
 400  CONTINUE                                                            140500
 500  CONTINUE                                                            140600
      RETURN                                                              140700
      END                                                                 140800
      SUBROUTINE MATROT                                                   140900
C                                                                         141000
C     THIS SUBROUTINE COMPUTES THE MATRIX ELEMENTS                        141100
C     ACCORDING TO A ROTATIONAL MODEL WITH INTRINSIC PARAMETERS           141200
C     INTM1,INTE1,INTE2,INTE3,INTE4,INTE5,INTE6                           141300
C     DEFINED IN BOHR AND MOTTELSON, VOL II, PG.155                       141400
C                                                                         141500
      REAL INTM1                                                          141600
      REAL INTE1,INTE2,INTE3,INTE4,INTE5,INTE6                            141700
      REAL MEM,MULT                                                       141800
      COMMON /BL1/LAMDA(12),LEAD(NNNN,NNNN,12),LDNUM(12,NNNN),LAMMAX      141900
      COMMON /BL2/NGMAX,NSTART(10),NSTOP(10),MASTER(10,10),EMMAX(10)      142000
      COMMON /BL3/EN(NNNN),SPIN(NNNN),ACCUR,DIPOL,ZPOL,BANDK(NNNN)        142100
      COMMON/BL10/IPAR(NNNN),IFAC(NNNN)                                   142200
      COMMON/BL19/ERR                                                     142300
      COMMON/BL20/MEM(NNNN,NNNN,12),MULT(NNNN,NNNN,12)                    142400
      COMMON/BL32/NMAX,NDIM                                               142500
      COMMON/BL36/INTM1,INTE1,INTE2,INTE3,INTE4,INTE5,INTE6               142600
C                                                                         142700
      K1=0                                                                142800
      K2=0                                                                142900
      K3=0                                                                143000
      K4=0                                                                143100
      K5=0                                                                143200
      K6=0                                                                143300
      K7=0                                                                143400
      PI=3.1415926                                                        143500
      FPY1=1.0/SQRT(4.0*PI/3.0)                                           143600
      FPY2=1.0/SQRT(16.0*PI/5.0)                                          143700
      FPY1=1.0                                                            143800
      FPY2=1.0                                                            143900
      DO 600 N = 1,NMAX                                                   144000
      AJ1 = SPIN(N)                                                       144100
      AM1 = BANDK(N)                                                      144200
      ROOTN = SQRT(2.0*SPIN(N) + 1.0)                                     144300
      DO 600 M = N,NMAX                                                   144400
      AJ3 = SPIN(M)                                                       144500
      ROOTM = SQRT(SPIN(M)*2.0 + 1.0)                                     144600
      IEX=ABS(SPIN(N)+BANDK(M))+0.001                                     144700
      DELK=BANDK(M)-BANDK(N)                                              144800
      AM3=-BANDK(M)                                                       144900
      PHZ = (-1.0)**IEX                                                   145000
      IF(INTE2.EQ.0.0)GOTO 200                                            145100
C                                                                         145200
C     COMPUTATION OF E2-MATRIX ELEMENT MEM(N,M,2)                         145300
      LAM = 2                                                             145400
      IF(ABS(DELK).GT.2.0)GOTO 200                                        145500
      IF(IPAR(N).NE.IPAR(M))GOTO 200                                      145600
      TRJ2=THREEJ(AJ1,AM1,2.00,DELK,AJ3,AM3)                              145700
      MEM(N,M,2)=PHZ*ROOTN*ROOTM*INTE2*MULT(N,M,LAM)*TRJ2*FPY2            145800
      IF(MEM(N,M,2).EQ.0.0)GOTO 200                                       145900
      I = 1                                                               146000
      K2 = K2 + 1                                                         146100
      IF(K2.EQ.1)GOTO 400                                                 146200
 200  IF(INTE3.EQ.0.0)GOTO 300                                            146300
C                                                                         146400
C     COMPUTATION OF E3-MATRIX ELEMENT MEM(N,M,3)                         146500
      LAM = 3                                                             146600
      IF(ABS(DELK).GT.3.0)GOTO 300                                        146700
      IF(IPAR(N).EQ.IPAR(M))GO TO 300                                     146800
      TRJ3=THREEJ(AJ1,AM1,3.00,DELK,AJ3,AM3)                              146900
      MEM(N,M,3)=-PHZ*ROOTN*ROOTM*INTE3*TRJ3                              147000
      IF(MEM(N,M,3).EQ.0.0)GOTO 300                                       147100
      I = 2                                                               147200
      K3 = K3 + 1                                                         147300
      IF(K3.EQ.1)GOTO 400                                                 147400
 300  IF(INTM1.EQ.0.0)GOTO 600                                            147500
C                                                                         147600
C     COMPUTATION OF M1 MATRIX-ELEMENT MEM(N,M,7)                         147700
      LAM=7                                                               147800
      IF(ABS(DELK).GT.1.0) GO TO 600                                      147900
      IF(IPAR(N).NE.IPAR(M))GO TO 600                                     148000
      TRJ7=THREEJ(AJ1,AM1,1.00,DELK,AJ3,AM3)                              148100
      MEM(N,M,7)=-PHZ*ROOTN*ROOTM*INTM1*TRJ7*FPY1                         148200
      IF(MEM(N,M,7).EQ.0.0)GO TO 600                                      148300
      I = 3                                                               148400
      K7 = K7+1                                                           148500
      IF(K7.NE.1) GO TO 600                                               148600
C                                                                         148700
C     DETERMINATION OF  LAMMAX  AND  LAMDA(I)                             148800
 400  IF(LAM.EQ.LAMDA(1).OR.LAM.EQ.LAMDA(2).OR.LAM.EQ.LAMDA(3).OR.LAM.EQ  148900
     1.LAMDA(4).OR.LAM.EQ.LAMDA(5).OR.LAM.EQ.LAMDA(6).OR.LAM.EQ.LAMDA(7)  149000
     2.OR.LAM.EQ.LAMDA(8).OR.LAM.EQ.LAMDA(9).OR.LAM.EQ.LAMDA(10).OR.LAM.  149100
     3EQ.LAMDA(11).OR.LAM.EQ.LAMDA(12))GOTO 500                           149200
      LAMMAX = LAMMAX + 1                                                 149300
      LAMDA(LAMMAX) = LAM                                                 149400
 500  IF(I.EQ.1)GOTO 200                                                  149500
      IF(I.EQ.2)GOTO 300                                                  149600
 600  CONTINUE                                                            149700
C                                                                         149800
C     PRINT-OUT OF ROTATIONAL MODEL PARAMETERS AND OF MULT-MATRIX         149900
      WRITE(5,900)                                                        150000
      WRITE(6,902) INTM1,INTE1,INTE2,INTE3,INTE4,INTE5,INTE6              150100
      RETURN                                                              150200
C                                                                         150300
 900  FORMAT(47H0MATRIX ELEMENTS ARE CALCULATED FROM ROTATIONAL,          150400
     122H MODEL WITH INTRINSIC /12H PARAMETERS,9X5HINTM19X               150500 
     26HINTE1 ,8X6HINTE2 ,8X6HINTE3 8X6HINTE4 8X6HINTE5 8X6HINTE6 )       150600
 902  FORMAT(13X7(4XF10.4))                                               150700
      END                                                                 150800
      SUBROUTINE PREP1                                                    150900
C                                                                         151000
C     THIS ROUTINE PREPARES THE INTEGRATION                               151100
C     IT COMPUTES XI-MATRIX,PSI-MATRIX,LEAD-MATRIX,SSTART- AND SSTOP-ARR  151200
C                                                                         151300
      INTEGER SSTART,SSTOP,OUXI,OUPSI                                     151400
      LOGICAL ERR                                                         151500
      REAL MEM                                                            151600
      COMMON /BL1/LAMDA(12),LEAD(NNNN,NNNN,12),LDNUM(12,NNNN),LAMMAX      151700
      COMMON /BL2/NGMAX,NSTART(10),NSTOP(10),MASTER(10,10),EMMAX(10)      151800
      COMMON /BL3/EN(NNNN),SPIN(NNNN),ACCUR,DIPOL,ZPOL,BANDK(NNNN)        151900
      COMMON /BL6/PSI(NNNN,NNNN,12)                                       152000
      COMMON /BL7/ZETA( ZZZZ),NZMAX                                       152100
      COMMON /BL9/IEXCIT(10),IEXNUM,LMAX(2),IRSTA(2),IRSTO(2),IP          152200
      COMMON/BL18/CAT(CCCC,3),ICATMX, ISMAX                               152300
      COMMON/BL19/ERR                                                     152400
      COMMON/BL20/MEM(NNNN,NNNN,12),MULT(NNNN,NNNN,12)                    152500
      COMMON/BL32/NMAX,NDIM                                               152600
      COMMON/BL35/ZP,ZT,A1,A2,EP,TLBDG                                    152700
      COMMON/BL40/SSTART(31),SSTOP(30)                                    152800
      COMMON/BL41/XI(NNNN,NNNN)                                           152900
      COMMON/BL50/OUXI,OUPSI,NCM,EMMAX1                                   153000
      COMMON/BL54/DISTA                                                   153100
      DIMENSION ETAN(NNNN),CPSI(12)                                       153200
C                                                                         153300
C     CRITICAL SCATTERING ANGLE IN CM SYSTEM                              153400
      DISTR = 1.20*(A1**0.33333 + A2**0.33333) + 4.0                      153500
      DISTA=0.0719949*(1.0+A1/A2)*ZP*ZT/EP                                153600
      D2A = 20.0*DISTA                                                    153700
      VINF = 0.0463365*SQRT(EP/A1)                                        153800
      WRITE(5,942) VINF,D2A                                               153900
      IF(D2A.GT.DISTR) GOTO 100                                           154000
      TCRIT=2.0* ASIN(DISTA*10./(DISTR-DISTA))*57.295779                  154100
      WRITE(5,900) DISTR,TCRIT                                            154200
C                                                                         154300
C     APPROXIMATION OF THE VIRTUAL EXCITATION OF THE GIANT DIPOLE         154400
C     RESONANCE BY DIPOLE POLARIZATION EFFECTS                            154500
C     SEE ALDER/WINTHER, "ELECTROMAGNETIC EXCITATION", APPENDIX J         154600
C     NORTH-HOLLAND PUBL.CO., AMSTERDAM (1975)                            154700
 100  ZPOL = DIPOL*EP*A2/(ZT*ZT*(1.+A1/A2))                               154800
C                                                                         154900
C     XI- MATRIX                                                          155000
      ETA=ZP*ZT*SQRT(A1/EP)/6.349770                                      155100
      DO 110 M=1,NMAX                                                     155200
      DEP=(1.0+A1/A2)*EN(M)                                               155300
      ZET=DEP/EP                                                          155400
      SZET=SQRT(1.0-ZET)                                                  155500
      ETAN(M)=ETA/SZET                                                    155600
 110  CONTINUE                                                            155700
      DO 120 N=1,NMAX                                                     155800
      DO 120 M=1,NMAX                                                     155900
      XI(N,M) = 99.9                                                      156000
 120  CONTINUE                                                            156100
      DO 150 I1=1,LAMMAX                                                  156200
      LAM = LAMDA(I1)                                                     156300
      DO 140 N=1,NMAX                                                     156400
      DO 130 M=1,NMAX                                                     156500
      IF( MEM(N,M,LAM).EQ.0.0 ) GOTO 130                                  156600
      XI(N,M) = ETAN(N) - ETAN(M)                                         156700
 130  CONTINUE                                                            156800
 140  CONTINUE                                                            156900
 150  CONTINUE                                                            157000
      IF(OUXI.EQ.0)GOTO 200                                               157100
C     PRINT-OUT OF XI-MATRIX                                              157200
      WRITE(5,902)                                                        157300
      NN = 1                                                              157400
      IF(NMAX.GT.11) NN=NGMAX                                             157500
      DO 190 I=1,NN                                                       157600
      M1 = NSTART(I)                                                      157700
      M2 = NSTOP(I)                                                       157800
      IF(NMAX.GT.11)GOTO 165                                              157900
      M1 = 1                                                              158000
      M2 = NMAX                                                           158100
 165  WRITE(5,906)(M,M=M1,M2)                                             158200
      DO 170 N=1,NMAX                                                     158300
      WRITE(5,908) N,(XI(N,M),M=M1,M2)                                    158400
 170  CONTINUE                                                            158500
      WRITE(5,904)                                                        158600
 190  CONTINUE                                                            158700
C                                                                         158800
C     COMPUTATION OF LEAD-MATRIX AND LDNUM-MATRIX                         158900
 200  DO 250 I5=1,LAMMAX                                                  159000
      LAM = LAMDA(I5)                                                     159100
      DO 240 I1=1,NGMAX                                                   159200
      N1 = NSTART(I1)                                                     159300
      N2 = NSTOP(I1)                                                      159400
      DO 230 N=N1,N2                                                      159500
      J = 0                                                               159600
      DO 220 I2=1,NGMAX                                                   159700
      IF ( MASTER(I1,I2).EQ.0 ) GOTO 220                                  159800
      M1 = NSTART(I2)                                                     159900
      M2 = NSTOP(I2)                                                      160000
      DO 210 M=M1,M2                                                      160100
      IF ( MEM(N,M,LAM).EQ.0.0 ) GOTO 210                                 160200
      J = J + 1                                                           160300
      LEAD(N,J,LAM) = M                                                   160400
 210  CONTINUE                                                            160500
 220  CONTINUE                                                            160600
      LDNUM(LAM,N) = J                                                    160700
 230  CONTINUE                                                            160800
 240  CONTINUE                                                            160900
 250  CONTINUE                                                            161000
C                                                                         161100
C     PSI- MATRIX                                                         161200
      AAZZ=1./(1.+A1/A2)/ZP/ZT                                            161300
      CPSI(1) = 5.169286*AAZZ                                             161400
      CPSI(2) = 14.359366*AAZZ*AAZZ                                       161500
      CPSI(3) = 56.982577 *AAZZ**3                                        161600
      CPSI(4) = 263.812653*AAZZ**4                                        161700
      CPSI(5) =1332.409500*AAZZ**5                                        161800
      CPSI(6) =7117.691577*AAZZ**6                                        161900
      CPSI(7) = VINF*CPSI(1)/95.0981942                                   162000
      CPSI(8) = VINF*CPSI(2)/95.0981942                                   162100
      CPSI(9) = VINF*CPSI(3)/95.0981942                                   162200
      CPSI(10)= VINF*CPSI(4)/95.0981942                                   162300
      CPSI(11)= VINF*CPSI(5)/95.0981942                                   162400
      CPSI(12)= VINF*CPSI(6)/95.0981942                                   162500
      DO 290 I1=1,LAMMAX                                                  162600
      LAM = LAMDA(I1)                                                     162700
      LAM1 = LAM                                                          162800
      IF(LAM.GT.6) LAM1=LAM-6                                             162900
      DO 280 I2=1,NGMAX                                                   163000
      Z = ZP                                                              163100
C     FOR PROJECTILE EXCITATION   Z = ZT                                  163200
      IF(IEXCIT(I2).EQ.1) Z=ZT                                            163300
      ZSQA = Z*SQRT(A1)                                                   163400
      N1 = NSTART(I2)                                                     163500
      N2 = NSTOP(I2)                                                      163600
      DO 270 N=N1,N2                                                      163700
      PP1 = (EP-(1.0+A1/A2)*EN(N))**0.25                                  163800
      DO 270 M=1,NMAX                                                     163900
      IF(MEM(N,M,LAM).EQ.0.0)GOTO 260                                     164000
      PP2 = (EP-(1.0+A1/A2)*EN(M))**0.25                                  164100
      PSI(N,M,LAM)=CPSI(LAM)*ZSQA*(PP1*PP2)**(2*LAM1-1)*MEM(N,M,LAM)      164200
      GOTO 270                                                            164300
 260  PSI(N,M,LAM) = 0.0                                                  164400
 270  CONTINUE                                                            164500
 280  CONTINUE                                                            164600
 290  CONTINUE                                                            164700
      IF(OUPSI.EQ.0)GOTO 350                                              164800
C     PRINT-OUT OF PSI-MATRIX                                             164900
      WRITE(5,910)                                                        165000
      DO 340 I1=1,LAMMAX                                                  165100
      LAM = LAMDA(I1)                                                     165200
      IF(LAM.GT.6)GOTO 302                                                165300
      WRITE(5,912) LAM                                                    165400
      GOTO 303                                                            165500
 302  LA = LAM - 6                                                        165600
      WRITE(5,938) LA                                                     165700
 303  NN = 1                                                              165800
      IF(NMAX.GT.11) NN=NGMAX                                             165900
      DO 330 I=1,NN                                                       166000
      M1 = NSTART(I)                                                      166100
      M2 = NSTOP(I)                                                       166200
      IF(NMAX.GT.11)GOTO 305                                              166300
      M1 = 1                                                              166400
      M2 = NMAX                                                           166500
 305  WRITE(5,906) (M,M=M1,M2)                                            166600
      DO 310 N=1,NMAX                                                     166700
      K = 0                                                               166800
      DO 307 M=M1,M2                                                      166900
      P = PSI(N,M,LAM)                                                    167000
      IF(P.NE.0..AND.ABS(P).LT.0.001) K=1                                 167100
 307  CONTINUE                                                            167200
      IF(K.EQ.1)GOTO 308                                                  167300
      WRITE(5,908) N,(PSI(N,M,LAM),M=M1,M2)                               167400
      GOTO 310                                                            167500
 308  WRITE(5,940) N,(PSI(N,M,LAM),M=M1,M2)                               167600
 310  CONTINUE                                                            167700
      WRITE(5,904)                                                        167800
 330  CONTINUE                                                            167900
 340  CONTINUE                                                            168000
C                                                                         168100
C     CATALOG OF MAGNETIC SUBSTATES; SSTART- AND SSTOP-ARRAY              168200
 350  IS = 1                                                              168300
      SSTART(1) = 1                                                       168400
      DO 390 I1=1,NGMAX                                                   168500
      N1 = NSTART(I1)                                                     168600
      N2 = NSTOP(I1)                                                      168700
      DO 380 N= N1,N2                                                     168800
      QUAN = SPIN(N)                                                      168900
      IF(QUAN.GT.EMMAX(I1))QUAN = EMMAX(I1)                               169000
      MSTOP = 2.0*QUAN + 1.0001                                           169100
      QUAN = -QUAN                                                        169200
      DO 370 I = 1,MSTOP                                                  169300
      IF ( IS.GT.ICATMX ) GOTO 360                                        169400
      CAT(IS,1) = N                                                       169500
      CAT(IS,2)= SPIN(N)                                                  169600
      CAT(IS,3) = QUAN                                                    169700
 360  QUAN = QUAN +1.0                                                    169800
      IS = IS + 1                                                         169900
 370  CONTINUE                                                            170000
      SSTART(N+1) = IS                                                    170100
      SSTOP(N) = IS - 1                                                   170200
 380  CONTINUE                                                            170300
 390  CONTINUE                                                            170400
      ISMAX = IS -1                                                       170500
      WRITE(5,914)ISMAX                                                   170600
      IF(ISMAX.LE.ICATMX)GOTO 400                                         170700
C     ERROR EXIT IF NUMBER OF MAGNETIC SUBSTATES IS TOO LARGE             170800
      WRITE(5,916) ICATMX                                                 170900
      GOTO 850                                                            171000
C     PRINTOUT OF SSTART ARRAY AND SSTOP ARRAY                            171100
 400  WRITE(5,918) (SSTART(I),I=1,NMAX)                                   171200
      WRITE(5,920) (SSTOP(I),I=1,NMAX)                                    171300
C                                                                         171400
C     LMAX-, IRSTA- AND IRSTO-ARRAY                                       171500
      LMAX(1) = SPIN(1) + 1.1                                             171600
      IRSTA(1) = 1                                                        171700
      IF(IEXNUM.EQ.2)GOTO 410                                             171800
      IRSTO(1) = ISMAX                                                    171900
      GOTO 420                                                            172000
 410  N1 = NSTART(IP)                                                     172100
      IRSTO(1) = SSTART(N1) - 1                                           172200
      LMAX(2) = SPIN(N1) + 1.1                                            172300
      IRSTA(2) = IRSTO(1) + 1                                             172400
      IRSTO(2) = ISMAX                                                    172500
C                                                                         172600
C     ZETA-ARRAY                                                          172700
 420  NZ = 0                                                              172800
      DO 500 K=1,IEXNUM                                                   172900
      N1 = 1                                                              173000
      IF(K.EQ.2) N1=NSTART(IP)                                            173100
      IF(SPIN(N1).EQ.0.)GOTO 450                                          173200
      IR1 = IRSTA(K)                                                      173300
      IR2 = IRSTO(K)                                                      173400
C     IR - LOOP FOR GROUND STATE SPIN .NE. 0                              173500
      DO 440 IR=IR1,IR2                                                   173600
      N = CAT(IR,1)                                                       173700
      CALL LSLOOP(IR,N,NZ)                                                173800
 440  CONTINUE                                                            173900
      GOTO 500                                                            174000
 450  IF(K.EQ.2)GOTO 460                                                  174100
      N1 = 1                                                              174200
      N2 = NMAX                                                           174300
      IF(IEXNUM.EQ.2) N2=NSTART(IP)-1                                     174400
      GOTO 470                                                            174500
 460  N1 = NSTART(IP)                                                     174600
      N2 = NMAX                                                           174700
C     IR - LOOP FOR GROUND STATE SPIN = 0                                 174800
 470  DO 490 N=N1,N2                                                      174900
      IR = SSTART(N) - 1                                                  175000
 480  IR = IR + 1                                                         175100
      CALL LSLOOP(IR,N,NZ)                                                175200
C     ZETA ARRAY IS PRODUCED IN THIS SUBROUTINE                           175300
      IF(CAT(IR,3).LT.-0.1)GOTO 480                                       175400
 490  CONTINUE                                                            175500
 500  CONTINUE                                                            175600
      WRITE(5,922)NZ                                                      175700
C     ERROR EXIT IF NUMBER OF ELEMENTS IN ZETA-ARRAY IS TOO LARGE         175800
      IF ( NZ.GT.NZMAX ) GOTO 540                                         175900
      WRITE(5,924) ETA                                                    176000
C                                                                         176100
      RETURN                                                              176200
 540  WRITE(5,936) NZMAX                                                  176300
 850  ERR = .TRUE.                                                        176400
      RETURN                                                              176500
C                                                                         176600
 900  FORMAT(56H0SAFE DISTANCE (1.20*(AP**(1/3)+AT**(1/3)) + 4.0)FERMI=  176700 
     1F6.2,6H FERMI,//44H WARNING  --  BOMBARDING ENERGY TOO HIGH FOR     176800
     218H HEAD-ON COLLISION,/34H CRITICAL SCATTERING ANGLE THETA =        176900
     3F6.2,21H DEGREES IN CM SYSTEM/)                                     177000
 902  FORMAT(10H0XI MATRIX)                                               177100
 904  FORMAT(///)                                                         177200
 906  FORMAT(6X,I14,10I11)                                                177300
 908  FORMAT(1H05X3HN =I2,11F11.4)                                        177400
 910  FORMAT(1H1)                                                         177500
 912  FORMAT(/18H0PSI MATRIX FOR  EI1)                                    177600
 914  FORMAT(45H0TOTAL NUMBER OF MAGNETIC SUBSTATES, ISMAX = I3)          177700
 916  FORMAT(28H0 ERROR ISMAX EXCEEDS ICATMX,5X,9H(ICATMX =I4,1H))        177800
 918  FORMAT(17H0SSTART ARRAY    ,20I5/17X20I5)                           177900
 920  FORMAT(17H0SSTOP  ARRAY    ,20I5/17X20I5,//)                        178000
 922  FORMAT(36H0NUMBER OF ELEMENTS IN ZETA-ARRAY = I5)                   178100
 924  FORMAT(7H0ETA = F6.2)                                               178200
 936  FORMAT(57H0ERROR  -  NUMBER OF ELEMENTS IN ZETA ARRAY EXCEEDS NZMA  178300
     1X,5X,8H(NZMAX =I6,1H))                                              178400
 938  FORMAT(/18H0PSI MATRIX FOR  MI1)                                    178500
 940  FORMAT(1H05X3HN =I2,11E11.3)                                        178600
 942  FORMAT(//37H0INITIAL VELOCITY OF PROJECTILE V/C =F6.3,             178700 
     1/56H DISTANCE OF CLOSED APPROACH IN HEAD-ON COLLISION 2*A =,       178800 
     2F8.2,6H FERMI)                                                      178900
      END                                                                 179000
      SUBROUTINE PREP2                                                    179100
C                                                                         179200
C     THIS ROUTINE COMPUTES RANGE AND STEP WIDTH OF INTEGRATION           179300
C                                                                         179400
      INTEGER OUXI,OUPSI                                                  179500
      LOGICAL ERR                                                         179600
      COMMON /BL3/EN(NNNN),SPIN(NNNN),ACCUR,DIPOL,ZPOL,BANDK(NNNN)        179700
      COMMON/BL19/ERR                                                     179800
      COMMON/BL21/XIMAX                                                   179900
      COMMON/BL22/THETA                                                   180000
      COMMON/BL25/EPS,EROOT                                               180100
      COMMON/BL32/NMAX,NDIM                                               180200
      COMMON/BL35/ZP,ZT,A1,A2,EP,TLBDG                                    180300
      COMMON/BL39/UP,DW,ISTEP,D2W                                         180400
      COMMON/BL41/XI(NNNN,NNNN)                                           180500
      COMMON/BL42/TCMDG(NNNN),ZLBDG(NNNN),R3(NNNN),R4(NNNN),EPP(NNNN),    180600
     1ETP(NNNN)                                                                 
      COMMON/BL50/OUXI,OUPSI,NCM,EMMAX1                                   180700
      COMMON/BL51/SIGTOT(NNNN)                                            180800
      COMMON/BL54/DISTA                                                   180900
C                                                                         181000
C     TRANSITION TO CM COORDINATE SYSTEM IF SCATTERING ANGLE IS GIVEN     181100
C     IN LABORATORY SYSTEM                                                181200
      IF(TLBDG.EQ.0.0.AND.THETA.NE.0.0) GOTO 100                          181300
      CALL CMLAB(A1,A2,EP,EN,NMAX,TLBDG,TCMDG,ZLBDG,R3,R4,EPP,ETP)        181400
      IF(ERR) GOTO 850                                                    181500
      THETA=TCMDG(NCM)                                                    181600
      WRITE(5,902) THETA                                                  181700
C                                                                         181800
C     RANGE AND STEP WIDTH OF INTEGRATION                                 181900
 100  TRAD=THETA/57.295779                                                182000
      STR = SIN(TRAD/2.0)                                                 182100
      EPS = 1.0/STR                                                       182200
      EROOT = SQRT(EPS*EPS-1.)                                            182300
      CLOSE = DISTA*(1.+EPS)*10.0                                         182400
      UP=ALOG(1.0/(EPS*SQRT(ACCUR)))                                      182500
C     DETERMINATION OF THE LARGEST XI-VALUE USED IN THE CALCULATION       182600
      XIM = 0.0                                                           182700
      DO 200 N = 1,NMAX                                                   182800
      DO 150 M = 1,NMAX                                                   182900
      IF(XI(N,M).GT.XIMAX)GOTO 150                                        183000
      IF(XI(N,M).LE.XIM)GOTO 150                                          183100
      XIM = XI(N,M)                                                       183200
 150  CONTINUE                                                            183300
 200  CONTINUE                                                            183400
      WRITE(5,900)XIM                                                     183500
      DW=40.0*(ACCUR**0.2)/(10.0+48.0*XIM+16.0*XIM*EPS)                   183600
      ISTEP=UP/(DW*8.)+1.                                                 183700
      ISTEP=ISTEP*8                                                       183800
      DW=UP/(FLOAT(ISTEP)-0.25)                                           183900
      UP=DW*FLOAT(ISTEP)                                                  184000
      WRITE(5,904)EPS,CLOSE,UP,ISTEP                                      184100
      D2W = DW + DW                                                       184200
      WRITE(5,906)D2W                                                     184300
 850  RETURN                                                              184400
C                                                                         184500
 900  FORMAT(7H0XIM = F10.4)                                              184600
 902  FORMAT(29H0CM SCATTERING ANGLE USED IS F7.2,8H DEGREES)             184700
 904  FORMAT(7H0EPS = F7.3,/                                              184800
     133H0DISTANCE OF CLOSEST APPROACH IS F8.2,6H FERMI/                  184900
     228H0RANGE OF INTEGRATION, UP = F6.2 /                               185000
     336H0ESTIMATED NUMBER OF STEPS, ISTEP = I4)                          185100
 906  FORMAT(27H0INITIAL STEP WIDTH, D2W = F8.5)                          185200
      END                                                                 185300
      SUBROUTINE PREP3(K)                                                 185400
C                                                                         185500
C     THIS ROUTINE IS CALLED DURING THE INTEGRATION OF THE AMPLITUDES     185600
C     IT DETERMINES THE RANGE OF DO-LOOPS FOR TARGET AND PROJECTILE       185700
C     EXCITATION (K=1,2)                                                  185800
C                                                                         185900
      COMMON /BL2/NGMAX,NSTART(10),NSTOP(10),MASTER(10,10),EMMAX(10)      186000
      COMMON /BL9/IEXCIT(10),IEXNUM,LMAX(2),IRSTA(2),IRSTO(2),IP          186100
      COMMON/BL12/LMX,IR1,IR2,N1,N2                                       186200
      COMMON/BL32/NMAX,NDIM                                               186300
C                                                                         186400
      LMX = LMAX(K)                                                       186500
      IR1 = IRSTA(K)                                                      186600
      IR2 = IRSTO(K)                                                      186700
      IF(K.EQ.2)GOTO 100                                                  186800
      N1 = 1                                                              186900
      N2 = NMAX                                                           187000
      IF(IEXNUM.EQ.2) N2=NSTART(IP)-1                                     187100
      RETURN                                                              187200
 100  N1 = NSTART(IP)                                                     187300
      N2 = NMAX                                                           187400
      RETURN                                                              187500
      END                                                                 187600
      SUBROUTINE Q(W)                                                     187700
C                                                                         187800
C     THIS ROUTINE COMPUTES Q-FUNCTIONS                                   187900
C                                                                         188000
      IMPLICIT COMPLEX (A)                                                188100
      REAL ACCUR                                                          188200
      COMPLEX CI,EX,QLM(12,7),DEN7,DEN8,DEN9,DEN10,DEN11,DEN12            188300
      COMPLEX CQU1,CQU2,CQU3,CQU4,CQU5,CQU6,CQU7,CQU8,CQU9,CQU10,CQU11,   188400
     1CQU12                                                               188500
      COMMON /BL1/LAMDA(12),LEAD(NNNN,NNNN,12),LDNUM(12,NNNN),LAMMAX      188600
      COMMON /BL3/EN(NNNN),SPIN(NNNN),ACCUR,DIPOL,ZPOL,BANDK(NNNN)        188700
      COMMON /BL4/MAGEXC                                                  188800
      COMMON/BL13/CQU1(7,NNNN,NNNN),CQU2(7,NNNN,NNNN),CQU3(7,NNNN,NNNN)   188900
      COMMON/BL14/CQU4(7,NNNN,NNNN),CQU5(7,NNNN,NNNN),CQU6(7,NNNN,NNNN)   189000
      COMMON/BL15/CQU7(7,NNNN,NNNN),CQU8(7,NNNN,NNNN),CQU9(7,NNNN,NNNN)   189100
      COMMON/BL16/CQU10(7,NNNN,NNNN),CQU11(7,NNNN,NNNN),                  189200
     1CQU12(7,NNNN,NNNN)                                                        
      COMMON/BL25/EPS,EROOT                                               189300
      COMMON/BL32/NMAX,NDIM                                               189400
      COMMON/BL41/XI(NNNN,NNNN)                                           189500
      DATA CI/(0.0,1.0)/                                                  189600
C                                                                         189700
C      VARIABLES NEEDED FOR THE CALCULATION OF QLM(LAM,MY)                189800
      EW=EXP(W)                                                           189900
      COSHY=0.5*(EW+1.0/EW)                                               190000
      SINHY=0.5*(EW-1.0/EW)                                               190100
      DEN=EPS*COSHY+1.0                                                   190200
      POL = 1.0-ZPOL/DEN                                                  190300
C     POL ACCOUNTS FOR THE VIRTUAL EXCITATION OF THE DIPOLE GIANT RESONA  190400
      DEN1=DEN*DEN                                                        190500
      DEN2=DEN1*DEN1                                                      190600
      DEN3=DEN2*DEN1                                                      190700
      DEN4=DEN2*DEN2                                                      190800
      DEN5=DEN4*DEN1                                                      190900
      DEN6=DEN3*DEN3                                                      191000
      SH1=EROOT*SINHY                                                     191100
      SH2=SH1*SH1                                                         191200
      SH3=SH2*SH1                                                         191300
      SH4=SH2*SH2                                                         191400
      SH5=SH2*SH3                                                         191500
      SH6=SH2*SH4                                                         191600
      CH1=COSHY+EPS                                                       191700
      CH2=CH1*CH1                                                         191800
      CH4=CH2*CH2                                                         191900
      IF(MAGEXC.EQ.0) GOTO 100                                            192000
C     QUANTITIES NEEDED ONLY FOR MAGNETIC EXCITATIONS                     192100
      A = COSHY + EPS + CI*EROOT*SINHY                                    192200
      B = EPS*COSHY + 1.                                                  192300
      B2 = B*B                                                            192400
      A2 = A*A                                                            192500
      B4 = B2*B2                                                          192600
      A4 = A2*A2                                                          192700
      A2B2 = A2*B2                                                        192800
      B6 = B4*B2                                                          192900
      A6 = A4*A2                                                          193000
      A2B4 = A2B2 * B2                                                    193100
      A4B2 = A4 * B2                                                      193200
      B8 = B6 * B2                                                        193300
      A8 = A6 * A2                                                        193400
      A2B6 = A2B4 * B2                                                    193500
      A4B4 = A4B2 * B2                                                    193600
      A6B2 = A6 * B2                                                      193700
      B10 = B8 * B2                                                       193800
      A10 = A8 * A2                                                       193900
      A2B8 = A2B6 * B2                                                    194000
      A4B6 = A4B4 * B2                                                    194100
      A6B4 = A6B2 * B2                                                    194200
      A8B2 = A8 * B2                                                      194300
      A1B2 = A * B2                                                       194400
      DEN7 = B2                                                           194500
      DEN8 = DEN7 * A1B2                                                  194600
      DEN9 = DEN8 * A1B2                                                  194700
      DEN10 = DEN9 * A1B2                                                 194800
      DEN11 = DEN10 * A1B2                                                194900
      DEN12 = DEN11 * A1B2                                                195000
C                                                                         195100
C     COMPUTE FUNCTIONS QLM(LAM,MY)                                       195200
 100  DO 300 I1=1,LAMMAX                                                  195300
      LAM = LAMDA(I1)                                                     195400
      GOTO(110,120,130,140,150,160,210,220,230,240,250,260),LAM           195500
 110  QLM(1,1)= 0.5*CH1/DEN1                                              195600
      QLM(1,2)=-0.35355339*CI*SH1/DEN1                                    195700
      GOTO 300                                                            195800
 120  QLM(2,1)= 0.75*(2.0*CH2-SH2)/DEN2 * POL                             195900
      QLM(2,2)=-1.83711730*CI*CH1*SH1/DEN2 * POL                          196000
      QLM(2,3)=-0.91855865*SH2/DEN2 * POL                                 196100
      GOTO 300                                                            196200
 130  QLM(3,1)= 1.875*CH1*(2.0*CH2-3.0*SH2)/DEN3                          196300
      QLM(3,2)=-1.62379763*CI*(4.0*CH2-SH2)*SH1/DEN3                      196400
      QLM(3,3)=-5.13489890*CH1*SH2/DEN3                                   196500
      QLM(3,4)= 2.09631373*CI*SH3/DEN3                                    196600
      GOTO 300                                                            196700
 140  QLM(4,1)= 1.09375000*(8.0*CH4-24.0*CH2*SH2+3.0*SH4)/DEN4            196800
      QLM(4,2)=-4.89139867*CI*CH1*(4.0*CH2-3.0*SH2)*SH1/DEN4              196900
      QLM(4,3)=-3.45874113*(6.0*CH2-SH2)*SH2/DEN4                         197000
      QLM(4,4)=12.9414244 *CI*CH1*SH3/DEN4                                197100
      QLM(4,5)= 4.57548440*SH4/DEN4                                       197200
      GOTO 300                                                            197300
 150  QLM(5,1)=1.230468*CH1*(-14.*CH2*(9.*SH2+DEN1)+30.*DEN2)/DEN5        197400
      QLM(5,2)=-1.347911*CI*SH1*(35.*CH2*(-3.*SH2+DEN1)+5.*DEN2)/DEN5     197500
      QLM(5,3)=-35.662372*SH2*CH1*(-3.*SH2+2.*DEN1)/DEN5                  197600
      QLM(5,4)=7.279552*CI*SH3*(9.*CH2-DEN1)/DEN5                         197700
      QLM(5,5)=30.884521*SH4*CH1/DEN5                                     197800
      QLM(5,6)=-9.766543*CI*SH5/DEN5                                      197900
      GOTO 300                                                            198000
 160  QLM(6,1)=2.707031*(21.*CH2*(-CH2*(11.*SH2+4.*DEN1)+5.*DEN2)-        198100
     1 5.*DEN3)/DEN6                                                      198200
      QLM(6,2)=-17.543567*CI*SH1*CH1*(3.*CH2*(-11.*SH2+DEN1)+5.*DEN2)/    198300
     1 DEN6                                                               198400
      QLM(6,3)=-13.869408*SH2*(3.*CH2*(-11.*SH2+5.*DEN1)+DEN2)/DEN6       198500
      QLM(6,4)=-27.738815*CI*SH3*CH1*(-11.*SH2+8.*DEN1)/DEN6              198600
      QLM(6,5)=15.193177*SH4*(11.*CH2-DEN1)/DEN6                          198700
      QLM(6,6)=-71.262308*CI*SH5*CH1/DEN6                                 198800
      QLM(6,7)=-20.571656*SH6/DEN6                                        198900
      GOTO 300                                                            199000
 210  QLM(7,1) = (0.,0.)                                                  199100
      QLM(7,2) = -0.35355339*EROOT/DEN7                                   199200
      GOTO 300                                                            199300
 220  AFAC = EROOT/DEN8                                                   199400
      QLM(8,1)=(0.,0.)                                                    199500
      QLM(8,2)=-0.45927933*AFAC*(B2+A2)                                   199600
      QLM(8,3)=-0.45927933*AFAC*(B2-A2)                                   199700
      GOTO 300                                                            199800
 230  AFAC = EROOT/DEN9                                                   199900
      QLM(9,1)=(0.,0.)                                                    200000
      QLM(9,2)=-0.13531647*AFAC*(5.*B4+6.*A2B2+5.*A4)                     200100
      QLM(9,3)=-0.85581650*AFAC*(B4-A4)                                   200200
      QLM(9,4)=-0.52407843*AFAC*(B2-A2)**2                                200300
      GOTO 300                                                            200400
 240  AFAC = EROOT/DEN10                                                  200500
      QLM(10,1)=(0.,0.)                                                   200600
      QLM(10,2)=-.15285621*AFAC*(7.*B6+9.*A2B4+9.*A4B2+7.*A6)             200700
      QLM(10,3)=-.21617132*AFAC*(7.*B6+3.*A2B4-3.*A4B2-7.*A6)             200800
      QLM(10,4)=-1.21325855*AFAC*(B6-A2B4-A4B2+A6)                        200900
      QLM(10,5)=-.57193557*AFAC*(B2-A2)**3                                201000
      GOTO 300                                                            201100
 250  AFAC = EROOT/DEN11                                                  201200
      QLM(11,1)=(0.,0.)                                                   201300
      QLM(11,2)=-0.08424444*AFAC*(21.*B8+28.*A2B6+30.*A4B4+28.*A6B2       201400
     1+21.*A8)                                                            201500
      QLM(11,3)=-0.89155931*AFAC*(3.*B8+2.*A2B6-2.*A6B2-3.*A8)            201600
      QLM(11,4)=-0.27298317*AFAC*(9.*B8-4.*A2B6-10.*A4B4-4.*A6B2+9.*A8)   201700
      QLM(11,5)=-1.54422603*AFAC*(B8-2.*A2B6+2.*A6B2-A8)                  201800
      QLM(11,6)=-0.61040893*AFAC*(B2-A2)**4                               201900
      GOTO 300                                                            202000
 260  AFAC = EROOT/DEN12                                                  202100
      QLM(12,1)=(0.,0.)                                                   202200
      QLM(12,2)=-0.09137275*AFAC*(33.*B10+45.*A2B8+50.*A4B6+50.*A6B4      202300
     1+45.*A8B2+33.*A10)                                                  202400
      QLM(12,3)=-0.14447300*AFAC*(33.*B10+27.*A2B8+10.*A4B6-10.*A6B4      202500
     1-27.*A8B2-33.*A10)                                                  202600
      QLM(12,4)=-0.43341900*AFAC*(11.*B10-A2B8-10.*A4B6-10.*A6B4-A8B2     202700
     1+11.*A10)                                                           202800
      QLM(12,5)=-0.31652448*AFAC*(11.*B10-15.*A2B8-10.*A4B6+10.*A6B4      202900
     1+15.*A8B2-11.*A10)                                                  203000
      QLM(12,6)=-1.85578928*AFAC*(B10-3.*A2B8+2.*A4B6+2.*A6B4-3.*A8B2     203100
     1+A10)                                                               203200
      QLM(12,7)=-0.64286427*AFAC*(B2-A2)**5                               203300
 300  CONTINUE                                                            203400
C                                                                         203500
C     MULTIPLY QLM(LAM,MY) WITH THE KINETIC EXPONENTIAL-FACTOR            203600
      RALFA=EPS*SINHY+W                                                   203700
      DO 930 I1=1,LAMMAX                                                  203800
      LAM = LAMDA(I1)                                                     203900
      LA = LAM                                                            204000
      IF(LAM.GT.6) LAM=LAM-6                                              204100
      MYMAX = LAM + 1                                                     204200
      DO 920 N=1,NMAX                                                     204300
      LD = LDNUM(LA,N)                                                    204400
      IF(LD.EQ.0) GOTO 920                                                204500
      DO 910 I2=1,LD                                                      204600
      M = LEAD(N,I2,LA)                                                   204700
      IF(N.EQ.M)GOTO 610                                                  204800
      X = XI(N,M)                                                         204900
      EX = COS(X*RALFA) + CI*SIN(X*RALFA)                                 205000
      DO 600 MY=1,MYMAX                                                   205100
      GOTO(410,420,430,440,450,460,510,520,530,540,550,560),LA            205200
 410  CQU1(MY,N,M) = QLM(LA,MY)*EX                                        205300
      CQU1(MY,M,N) = QLM(LA,MY)/EX                                        205400
      GOTO 600                                                            205500
 420  CQU2(MY,N,M) = QLM(LA,MY)*EX                                        205600
      CQU2(MY,M,N) = QLM(LA,MY)/EX                                        205700
      GOTO 600                                                            205800
 430  CQU3(MY,N,M) = QLM(LA,MY)*EX                                        205900
      CQU3(MY,M,N) = QLM(LA,MY)/EX                                        206000
      GOTO 600                                                            206100
 440  CQU4(MY,N,M) = QLM(LA,MY)*EX                                        206200
      CQU4(MY,M,N) = QLM(LA,MY)/EX                                        206300
      GOTO 600                                                            206400
 450  CQU5(MY,N,M) = QLM(LA,MY)*EX                                        206500
      CQU5(MY,M,N) = QLM(LA,MY)/EX                                        206600
      GOTO 600                                                            206700
 460  CQU6(MY,N,M) = QLM(LA,MY)*EX                                        206800
      CQU6(MY,M,N) = QLM(LA,MY)/EX                                        206900
      GOTO 600                                                            207000
 510  CQU7(MY,N,M) = QLM(LA,MY)*EX                                        207100
      CQU7(MY,M,N) = QLM(LA,MY)/EX                                        207200
      GOTO 600                                                            207300
 520  CQU8(MY,N,M) = QLM(LA,MY)*EX                                        207400
      CQU8(MY,M,N) = QLM(LA,MY)/EX                                        207500
      GOTO 600                                                            207600
 530  CQU9(MY,N,M) = QLM(LA,MY)*EX                                        207700
      CQU9(MY,M,N) = QLM(LA,MY)/EX                                        207800
      GOTO 600                                                            207900
 540  CQU10(MY,N,M)= QLM(LA,MY)*EX                                        208000
      CQU10(MY,M,N)= QLM(LA,MY)/EX                                        208100
      GOTO 600                                                            208200
 550  CQU11(MY,N,M)= QLM(LA,MY)*EX                                        208300
      CQU11(MY,M,N)= QLM(LA,MY)/EX                                        208400
      GOTO 600                                                            208500
 560  CQU12(MY,N,M)= QLM(LA,MY)*EX                                        208600
      CQU12(MY,M,N)= QLM(LA,MY)/EX                                        208700
 600  CONTINUE                                                            208800
      GOTO 910                                                            208900
 610  DO 900 MY=1,MYMAX                                                   209000
      GOTO(710,720,730,740,750,760,810,820,830,840,850,860),LA            209100
 710  CQU1(MY,N,M) = QLM(LA,MY)                                           209200
      GOTO 900                                                            209300
 720  CQU2(MY,N,M) = QLM(LA,MY)                                           209400
      GOTO 900                                                            209500
 730  CQU3(MY,N,M) = QLM(LA,MY)                                           209600
      GOTO 900                                                            209700
 740  CQU4(MY,N,M) = QLM(LA,MY)                                           209800
      GOTO 900                                                            209900
 750  CQU5(MY,N,M) = QLM(LA,MY)                                           210000
      GOTO 900                                                            210100
 760  CQU6(MY,N,M) = QLM(LA,MY)                                           210200
      GOTO 900                                                            210300
 810  CQU7(MY,N,M) = QLM(LA,MY)                                           210400
      GOTO 900                                                            210500
 820  CQU8(MY,N,M) = QLM(LA,MY)                                           210600
      GOTO 900                                                            210700
 830  CQU9(MY,N,M) = QLM(LA,MY)                                           210800
      GOTO 900                                                            210900
 840  CQU10(MY,N,M)= QLM(LA,MY)                                           211000
      GOTO 900                                                            211100
 850  CQU11(MY,N,M)= QLM(LA,MY)                                           211200
      GOTO 900                                                            211300
 860  CQU12(MY,N,M)= QLM(LA,MY)                                           211400
 900  CONTINUE                                                            211500
 910  CONTINUE                                                            211600
 920  CONTINUE                                                            211700
 930  CONTINUE                                                            211800
      RETURN                                                              211900
      END                                                                 212000
      SUBROUTINE READER                                                   212100
C                                                                         212200
C     THIS ROUTINE READS AND CHECKS THE INPUT DATA                        212300
C                                                                         212400
      LOGICAL ERR                                                         212500
      REAL INTM1                                                          212600
      REAL INTE1,INTE2,INTE3,INTE4,INTE5,INTE6                            212700
      INTEGER OUAMP,OUPROW,OUXI,OUPSI                                     212800
      INTEGER OURHOB,OURHOC,OURHOX                                        212900
      REAL MEM                                                            213000
      COMMON /BL1/LAMDA(12),LEAD(NNNN,NNNN,12),LDNUM(12,NNNN),LAMMAX      213100
      COMMON /BL2/NGMAX,NSTART(10),NSTOP(10),MASTER(10,10),EMMAX(10)      213200
      COMMON /BL3/EN(NNNN),SPIN(NNNN),ACCUR,DIPOL,ZPOL,BANDK(NNNN)        213300
      COMMON /BL4/MAGEXC                                                  213400
      COMMON /BL9/IEXCIT(10),IEXNUM,LMAX(2),IRSTA(2),IRSTO(2),IP          213500
      COMMON/BL10/IPAR(NNNN),IFAC(NNNN)                                   213600
      COMMON/BL19/ERR                                                     213700
      COMMON/BL20/MEM(NNNN,NNNN,12),MULT(NNNN,NNNN,12)                    213800
      COMMON/BL21/XIMAX                                                   213900
      COMMON/BL22/THETA                                                   214000
      COMMON/BL23/IDN17,IDN18,IDN19,IDN25                                 214100
      COMMON/BL31/OURHOB,OURHOC,OURHOX                                    214200
      COMMON/BL32/NMAX,NDIM                                               214300
      COMMON/BL33/NTIME                                                   214400
      COMMON/BL35/ZP,ZT,A1,A2,EP,TLBDG                                    214500
      COMMON/BL36/INTM1,INTE1,INTE2,INTE3,INTE4,INTE5,INTE6               214600
      COMMON/BL38/ALPHA,BETA,GAMMA                                        214700
      COMMON/BL48/NANGLE,THI,TLI                                          214800
      COMMON/BL49/OUAMP,OUPROW,INTERV,INTIN                               214900
      COMMON/BL50/OUXI,OUPSI,NCM,EMMAX1                                   215000
      COMMON/BL51/SIGTOT(NNNN)                                            215100
      COMMON/BL52/IPURG                                                   215200
      COMMON/BL55/DTHETA,DTLBDG                                           215300
      DIMENSION FINPUT(72),A(72)                                          215400
C                                                                         215500
C     READ IN STARTS                                                      215600
C                                                                         215700
      WRITE(5,900)                                                        215800
 100  CONTINUE                                                            215900
      DO 101 NI=1,8                                                       216000
 101  FINPUT(NI)=0.0                                                      216100
      READ(20,901,END=105),(A(I), I=1,72)        !********                216200
      CALL LEXEM(A,FINPUT,N,MODE,IERR)                                    216300
      IF (N.LT.8 .AND. IERR.EQ.0) GOTO 102                                216400
      WRITE(5,912)                                                        216500
      ERR = .TRUE.                                                        216600
 102  IF( FINPUT(1) .EQ. 0. ) GOTO 600                                    216700
C     STATEMENT 600 MEANS  START EXECUTION OF THE JOB                     216800
      IF(FINPUT(1).LT.0.)GOTO 700                                         216900
C     STATEMENT 700 MEANS  PURGE ALL PREVIOUS INPUT                       217000
      WRITE(5,902) (FINPUT(I), I=1,7)                                     217100
      IF( FINPUT(1) .LE. 40. ) GOTO 110                                   217200
C                                                                         217300
C     PROGRAM STOPS HERE IF NO MORE DATA                                  217400
 105  STOP                                                                217500
C                                                                         217600
C     ASSIGN INPUT DATA TO PROPER VARIABLES                               217700
C                                                                         217800
 110  IJUMP = IFIX( FINPUT(1) + .00001 )                                  217900
      GOTO (501,502,503,504,505,506,507,508,509,510,511,512,513,514,      218000
     1 515,516,517,518,519,520,521,522,523,524,525,526,527,528,529,530,   218100
     2 531,532,533,534,535,536,537,538,539,540),IJUMP                     218200
 501  CONTINUE                                                            218300
      NMAX = FINPUT(2)                                                    218400
      IF(NMAX.LT.2.OR.NMAX.GT.NDIM) CALL ERROR(1)                         218500
      GOTO 100                                                            218600
 502  CONTINUE                                                            218700
      NCM = FINPUT(2)                                                     218800
      IF(NCM.GE.1.AND.NCM.LE.NDIM)GOTO 100                                218900
      NCM = 2                                                             219000
      WRITE(5,904)                                                        219100
      GOTO 100                                                            219200
 503  CONTINUE                                                            219300
      NTIME = FINPUT(2)                                                   219400
      IF(NTIME.LE.0) CALL ERROR(21)                                       219500
      GOTO 100                                                            219600
 504  CONTINUE                                                            219700
      XIMAX = FINPUT(2)                                                   219800
      IF(XIMAX.LE.0..OR.XIMAX.GE.99.)CALL ERROR(2)                        219900
      GOTO 100                                                            220000
 505  CONTINUE                                                            220100
      EMMAX1 = FINPUT(2)                                                  220200
      IF(EMMAX1.LT.0.)CALL ERROR(3)                                       220300
      GOTO 100                                                            220400
 506  CONTINUE                                                            220500
      ACCUR=FINPUT(2)                                                     220600
      IF(ACCUR.LE.0.)CALL ERROR(4)                                        220700
      IF(ACCUR.LE.0.0001)GOTO 100                                         220800
      ACCUR = 0.0001                                                      220900
      WRITE(5,906)                                                        221000
      GOTO 100                                                            221100
 507  CONTINUE                                                            221200
      OUXI = FINPUT(2)                                                    221300
      GOTO 100                                                            221400
 508  CONTINUE                                                            221500
      OUPSI = FINPUT(2)                                                   221600
      GOTO 100                                                            221700
 509  CONTINUE                                                            221800
      OUAMP = FINPUT(2)                                                   221900
      GOTO 100                                                            222000
 510  CONTINUE                                                            222100
      OUPROW = FINPUT(2)                                                  222200
      GOTO 100                                                            222300
 511  CONTINUE                                                            222400
      OURHOB = FINPUT(2)                                                  222500
      GOTO 100                                                            222600
 512  CONTINUE                                                            222700
      OURHOC = FINPUT(2)                                                  222800
      GOTO 100                                                            222900
 513  CONTINUE                                                            223000
      GOTO 100                                                            223100
 514  CONTINUE                                                            223200
      ALPHA = FINPUT(2)                                                   223300
      BETA = FINPUT(3)                                                    223400
      GAMMA = FINPUT(4)                                                   223500
      OURHOX = 1                                                          223600
      OURHOC = 1                                                          223700
      GOTO 100                                                            223800
 515  CONTINUE                                                            223900
      DIPOL = FINPUT(2)                                                   224000
      GOTO 100                                                            224100
 516  CONTINUE                                                            224200
      INTIN=FINPUT(2)                                                     224300
      IF(INTIN.GE.1)GOTO 100                                              224400
      INTIN = 1                                                           224500
      WRITE(5,908)                                                        224600
      GOTO 100                                                            224700
 517  CONTINUE                                                            224800
      IDN17 = 1                                                           224900
      ZP = FINPUT(2)                                                      225000
      A1 = FINPUT(3)                                                      225100
      IF(ZP.LE.0..OR.A1.LE.0.)CALL ERROR(5)                               225200
      IF(ZP.GT.A1)CALL ERROR(6)                                           225300
      GOTO 100                                                            225400
 518  CONTINUE                                                            225500
      IDN18 = 1                                                           225600
      ZT = FINPUT(2)                                                      225700
      A2 = FINPUT(3)                                                      225800
      IF(ZT.LE.0..OR.A2.LE.0.)CALL ERROR(5)                               225900
      IF(ZT.GT.A2)CALL ERROR(6)                                           226000
      GOTO 100                                                            226100
 519  CONTINUE                                                            226200
      IDN19 = 1                                                           226300
      EP = FINPUT(2)                                                      226400
      IF(EP.LE.0.)CALL ERROR(7)                                           226500
      GOTO 100                                                            226600
 520  CONTINUE                                                            226700
      TLI = FINPUT(2)                                                     226800
      IF(TLI.LT.0..OR.TLI.GT.180.)CALL ERROR(8)                           226900
      TLBDG = TLI                                                         227000
      GOTO 100                                                            227100
 521  CONTINUE                                                            227200
      THI = FINPUT(2)                                                     227300
      IF(THI.LT.0..OR.THI.GT.180.)CALL ERROR(8)                           227400
      THETA = THI                                                         227500
      GOTO 100                                                            227600
 522  CONTINUE                                                            227700
      N = FINPUT(2)                                                       227800
      IF(N.LT.1.OR.N.GT.NDIM) CALL ERROR(9)                               227900
      SPIN(N) = FINPUT(3)                                                 228000
      EN(N) = FINPUT(4)                                                   228100
      IPAR(N) = FINPUT(5)                                                 228200
      BANDK(N)=FINPUT(6)                                                  228300
      IF(SPIN(N).LT.0.)CALL ERROR(10)                                     228400
      IF(EN(N).LT.0.)CALL ERROR(11)                                       228500
      IF(IPAR(N).EQ.0.0)IPAR(N)=1.0                                       228600
      IF(IPAR(N).NE.1..AND.IPAR(N).NE.-1.)CALL ERROR(12)                  228700
      GOTO 100                                                            228800
 523  CONTINUE                                                            228900
      N = FINPUT(2)                                                       229000
      M = FINPUT(3)                                                       229100
      LAM=FINPUT(5)                                                       229200
      IF(LAM.EQ.0) LAM=2                                                  229300
      IF(N.LT.1.OR.N.GT.NDIM)CALL ERROR(9)                                229400
      IF(M.LT.1.OR.M.GT.NDIM)CALL ERROR(9)                                229500
      IF(LAM.LT.1.OR.LAM.GT.6)CALL ERROR(13)                              229600
      MEM(N,M,LAM)=FINPUT(4)                                              229700
      IF(LAM.EQ.LAMDA(1).OR.LAM.EQ.LAMDA(2).OR.LAM.EQ.LAMDA(3).OR.LAM.EQ  229800
     1.LAMDA(4).OR.LAM.EQ.LAMDA(5).OR.LAM.EQ.LAMDA(6).OR.LAM.EQ.LAMDA(7)  229900
     2.OR.LAM.EQ.LAMDA(8).OR.LAM.EQ.LAMDA(9).OR.LAM.EQ.LAMDA(10).OR.      230000
     3LAM.EQ.LAMDA(11).OR.LAM.EQ.LAMDA(12))GOTO 100                       230100
      LAMMAX = LAMMAX + 1                                                 230200
      LAMDA(LAMMAX) = LAM                                                 230300
      GOTO 100                                                            230400
 524  CONTINUE                                                            230500
      MAGEXC = 1                                                          230600
      N = FINPUT(2)                                                       230700
      M = FINPUT(3)                                                       230800
      LA= FINPUT(5)                                                       230900
      IF(LA.EQ.0) LA = 1                                                  231000
      LAM = LA + 6                                                        231100
      IF(N.LT.1.OR.N.GT.NDIM) CALL ERROR(9)                               231200
      IF(M.LT.1.OR.M.GT.NDIM) CALL ERROR(9)                               231300
      IF(LA.LT.1.OR.LA.GT.6) CALL ERROR(13)                               231400
      MEM(N,M,LAM)=FINPUT(4)                                              231500
      IF(LAM.EQ.LAMDA(1).OR.LAM.EQ.LAMDA(2).OR.LAM.EQ.LAMDA(3).OR.LAM.EQ  231600
     1.LAMDA(4).OR.LAM.EQ.LAMDA(5).OR.LAM.EQ.LAMDA(6).OR.LAM.EQ.LAMDA(7)  231700
     2.OR.LAM.EQ.LAMDA(8).OR.LAM.EQ.LAMDA(9).OR.LAM.EQ.LAMDA(10).OR.      231800
     3LAM.EQ.LAMDA(11).OR.LAM.EQ.LAMDA(12))GOTO 100                       231900
      LAMMAX = LAMMAX + 1                                                 232000
      LAMDA(LAMMAX) = LAM                                                 232100
      GOTO 100                                                            232200
 525  CONTINUE                                                            232300
      IDN25 = 1                                                           232400
      DTHETA = FINPUT(2)                                                  232500
      IF(DTHETA.LE.0..OR.DTHETA.GE.180.)CALL ERROR(22)                    232600
      GOTO 100                                                            232700
 526  CONTINUE                                                            232800
      NGMAX = FINPUT(2)                                                   232900
      IF(NGMAX.LT.1.OR.NGMAX.GT.10)CALL ERROR(14)                         233000
      GOTO 100                                                            233100
 527  CONTINUE                                                            233200
      I = FINPUT(2)                                                       233300
      IF(I.LT.1.OR.I.GT.10)CALL ERROR(15)                                 233400
      NSTART(I) = FINPUT(3)                                               233500
      NSTOP(I) = FINPUT(4)                                                233600
      EMMAX(I) = FINPUT(5)                                                233700
      IF(NSTART(I).LT.1.OR.NSTART(I).GT.NDIM)CALL ERROR(16)               233800
      IF(NSTOP(I).LT.1.OR.NSTOP(I).GT.NDIM)CALL ERROR(16)                 233900
      IF(NSTOP(I).LT.NSTART(I))CALL ERROR(17)                             234000
      IF(EMMAX(I).LT.0.)CALL ERROR(18)                                    234100
      IF(NSTOP(I)-NSTART(I).GT.10) CALL ERROR(20)                         234200
      GOTO 100                                                            234300
 528  CONTINUE                                                            234400
      I1 = FINPUT(2)                                                      234500
      I2 = FINPUT(3)                                                      234600
      IF(I1.LT.1.OR.I1.GT.10)CALL ERROR(15)                               234700
      IF(I2.LT.1.OR.I2.GT.10)CALL ERROR(15)                               234800
      MASTER(I1,I2) = FINPUT(4)                                           234900
      GOTO 100                                                            235000
 529  CONTINUE                                                            235100
      I = FINPUT(2)                                                       235200
      IF(I.LT.1.OR.I.GT.10)CALL ERROR(15)                                 235300
      IEXCIT(I) = FINPUT(3)                                               235400
      IF(IEXCIT(I).NE.1.AND.IEXCIT(I).NE.2)CALL ERROR(19)                 235500
      GOTO 100                                                            235600
 530  CONTINUE                                                            235700
      INTM1 = FINPUT(2)                                                   235800
      INTE1 = FINPUT(3)                                                   235900
      INTE2 = FINPUT(4)                                                   236000
      INTE3 = FINPUT(5)                                                   236100
      INTE4 = FINPUT(6)                                                   236200
      INTE5 = FINPUT(7)                                                   236300
      INTE6=FINPUT(8)                                                     236400
      CALL MATROT                                                         236500
      GOTO 100                                                            236600
 531  CONTINUE                                                            236700
 532  CONTINUE                                                            236800
 533  CONTINUE                                                            236900
 534  CONTINUE                                                            237000
 535  CONTINUE                                                            237100
 536  CONTINUE                                                            237200
 537  CONTINUE                                                            237300
 538  CONTINUE                                                            237400
 539  CONTINUE                                                            237500
 540  CONTINUE                                                            237600
      GOTO 100                                                            237700
C                                                                         237800
C     START EXECUTION OF THE JOB                                          237900
 600  IF(ERR)GOTO 710                                                     238000
      CALL CHECK                                                          238100
      IF(ERR)GOTO 710                                                     238200
      CALL START                                                          238300
      IF(ERR)GOTO 710                                                     238400
      RETURN                                                              238500
C                                                                         238600
C     RETURN TO MAIN PROGRAM AND PURGE ALL PREVIOUS INPUT                 238700
 700  IPURG = 1                                                           238800
      RETURN                                                              238900
C                                                                         239000
C     PRINT MESSAGE "EXECUTION TERMINATED"                                239100
 710  WRITE(5,910)                                                        239200
      RETURN                                                              239300
C                                                                         239400
 900  FORMAT( //, 32H1LIST OF NEW CARDS FOR THIS RUN. )                   239500
 901  FORMAT(72A1)                                                        239600
 902  FORMAT( 8F10.5 )                                                    239700
 904  FORMAT(24H PROGRAM ASSUMES NCM = 2)                                 239800
 906  FORMAT(32HPROGRAM ASSUMES  ACCUR = 0.0001,/)                        239900
 908  FORMAT(28H PROGRAM ASSUMES  INTERV = 1)                             240000
 910  FORMAT(21H0EXECUTION TERMINATED)                                    240100
 912  FORMAT(50H INPUT ERROR -- DATA CARD MUST CONTAIN LESS THAN 8        240200
     141H NUMBERS AND LESS THAN 72 DATA CHARACTERS)                       240300
      END                                                                 240400
      SUBROUTINE ROTATE(RHONEW,RHOOLD,ALPHA,BETA,GAMMA)                   240500
C                                                                         240600
C     THIS ROUTINE TRANSFORMS THE STATISTICAL TENSORS FROM                240700
C     ONE COORDINATE SYSTEM TO ANONTHER ONE GENERATED BY A                240800
C     ROTATION WITH THE EULERIAN ANGLES ALPHA, BETA, GAMMA                240900
C                                                                         241000
      COMPLEX CI,RHOOLD,RHONEW,TE                                         241100
      COMMON /BL2/NGMAX,NSTART(10),NSTOP(10),MASTER(10,10),EMMAX(10)      241200
      COMMON /BL3/EN(NNNN),SPIN(NNNN),ACCUR,DIPOL,ZPOL,BANDK(NNNN)        241300
      COMMON /BL9/IEXCIT(10),IEXNUM,LMAX(2),IRSTA(2),IRSTO(2),IP          241400
      COMMON/BL32/NMAX,NDIM                                               241500
      DIMENSION RHOOLD(NNNN,3,5),RHONEW(NNNN,3,5)                         241600
      DATA CI/(0.0,1.0)/                                                  241700
C                                                                         241800
      ALPHAR = ALPHA/57.295779                                            241900
      GAMMAR = GAMMA/57.295779                                            242000
      DO 590 K=1,IEXNUM                                                   242100
      N1 = 2                                                              242200
      N2 = NMAX                                                           242300
      IF(IEXNUM.EQ.2) N2=NSTART(IP)-1                                     242400
      IF(K.EQ.1)GOTO 510                                                  242500
      N1 = NSTART(IP) + 1                                                 242600
      N2 = NMAX                                                           242700
 510  DO 580 N=N1,N2                                                      242800
      KA = 0                                                              242900
      KAMAX = 2.02*SPIN(N)                                                243000
      IF(KAMAX.GT.4) KAMAX=4                                              243100
      KASTOP = KAMAX + 1                                                  243200
      DO 570 KAI=1,KASTOP,2                                               243300
C     KAI IS INDEX COUNTER ONLY                                           243400
      DJ = KA                                                             243500
      KAINDX = KA/2 + 1                                                   243600
      KAPPA = KA                                                          243700
      KAPSTP = KA + 1                                                     243800
      DO 560 KAPPI=1,KAPSTP                                               243900
C     KAPPI IS INDEX COUNTER ONLY                                         244000
      DMP = KAPPA                                                         244100
      KAPPIN = KAPPA + 1                                                  244200
      TE = (0.,0.)                                                        244300
C     SUMMATION OVER INDICES .GE. 0 IN EQUATION (I.4.4)                   244400
      KPR = 0                                                             244500
 520  DM = KPR                                                            244600
      KPRI = KPR + 1                                                      244700
      TE = TE+RHOOLD(N,KAINDX,KPRI)*(COS(KPR*ALPHAR) -                    244800
     1CI*SIN(KPR*ALPHAR))*DJMM(BETA,DJ,DM,DMP)                            244900
      KPR = KPR + 1                                                       245000
      IF(KPR.LE.KA)GOTO 520                                               245100
      IF(KA.EQ.0)GOTO 550                                                 245200
C     SUMMATION OVER NEGATIVE INDICES IN EQUATION (I.4.4)                 245300
      KPR = 1                                                             245400
 540  DM = -KPR                                                           245500
      KPRI = KPR + 1                                                      245600
      TE = TE+(-1)**KPR*CONJG(RHOOLD(N,KAINDX,KPRI))*(COS(-KPR*           245700
     1ALPHAR)-CI*SIN(-KPR*ALPHAR))*DJMM(BETA,DJ,DM,DMP)                   245800
      KPR = KPR + 1                                                       245900
      IF(KPR.LE.KA)GOTO 540                                               246000
 550  RHONEW(N,KAINDX,KAPPIN)=TE*(COS(KAPPA*GAMMAR)-CI*                   246100
     1SIN(KAPPA*GAMMAR))                                                  246200
      KAPPA = KAPPA - 1                                                   246300
 560  CONTINUE                                                            246400
      KA = KA + 2                                                         246500
 570  CONTINUE                                                            246600
 580  CONTINUE                                                            246700
 590  CONTINUE                                                            246800
      RETURN                                                              246900
      END                                                                 247000
      SUBROUTINE SIGM                                                     247100
C                                                                         247200
C     THIS ROUTINE COMPUTES THE DIFFERENTIAL CROSS SECTIONS AND PERFORMS  247300
C     THE INTEGRATION OVER SCATTERING ANGLES                              247400
C                                                                         247500
      LOGICAL ERR                                                         247600
      INTEGER OURHOB,OURHOC,OURHOX                                        247700
      COMMON /BL2/NGMAX,NSTART(10),NSTOP(10),MASTER(10,10),EMMAX(10)      247800
      COMMON /BL3/EN(NNNN),SPIN(NNNN),ACCUR,DIPOL,ZPOL,BANDK(NNNN)        247900
      COMMON/BL11/INTEND                                                  248000
      COMMON/BL19/ERR                                                     248100
      COMMON/BL21/XIMAX                                                   248200
      COMMON/BL22/THETA                                                   248300
      COMMON/BL25/EPS,EROOT                                               248400
      COMMON/BL31/OURHOB,OURHOC,OURHOX                                    248500
      COMMON/BL32/NMAX,NDIM                                               248600
      COMMON/BL34/SIMP                                                    248700
      COMMON/BL35/ZP,ZT,A1,A2,EP,TLBDG                                    248800
      COMMON/BL38/ALPHA,BETA,GAMMA                                        248900
      COMMON/BL41/XI(NNNN,NNNN)                                           249000
      COMMON/BL42/TCMDG(NNNN),ZLBDG(NNNN),R3(NNNN),R4(NNNN),EPP(NNNN),    249100
     1ETP(NNNN)                                                                 
      COMMON/BL45/DTHETR,DTLBR                                            249200
      COMMON/BL46/DSIG(NNNN),DSIGLB(NNNN),RU(NNNN)                        249300
      COMMON/BL48/NANGLE,THI,TLI                                          249400
      COMMON/BL51/SIGTOT(NNNN)                                            249500
      COMMON/BL53/P(NNNN)                                                 249600
      COMMON/BL54/DISTA                                                   249700
      COMMON/BL55/DTHETA,DTLBDG                                           249800
C                                                                         249900
C     COMPUTATION AND PRINT-OUT OF THE DIFFERENTIAL CROSS-SECTIONS        250000
      IF(TLBDG.NE.0.)GOTO 110                                             250100
C     HEADING FOR CM CROSS SECTIONS                                       250200
      WRITE(5,902) THETA                                                  250300
      WRITE(5,904)                                                        250400
      GOTO 120                                                            250500
C     HEADING FOR LAB CROSS SECTIONS                                      250600
 110  WRITE(5,906) TLBDG                                                  250700
      WRITE(5,908)                                                        250800
 120  DO 190 I=1,NGMAX                                                    250900
      WRITE(5,910) I                                                      251000
      N1 = NSTART(I)                                                      251100
      N2 = NSTOP(I)                                                       251200
      DO 180 N=N1,N2                                                      251300
      RU(N)=0.25*SQRT(EP/(EP-(1.0+A1/A2)*EN(N)))*DISTA**2*EPS**4          251400
      DSIG(N) = RU(N)*P(N)                                                251500
      IF(TLBDG.NE.0.)GOTO 160                                             251600
      WRITE(5,912)N,P(N),DSIG(N)                                          251700
      GOTO 180                                                            251800
 160  DSIGLB(N) = DSIG(N)*R3(N)                                           251900
      WRITE(5,916)N,TCMDG(N),P(N),DSIGLB(N),EPP(N),ETP(N),ZLBDG(N),R4(N)  252000
 180  CONTINUE                                                            252100
 190  CONTINUE                                                            252200
C                                                                         252300
C     PRINT-OUT OF THE STATISTICAL TENSORS                                252400
      IF(OURHOB.EQ.0)GOTO 250                                             252500
      WRITE(5,920)                                                        252600
      CALL TENS2(1)                                                       252700
 250  IF(OURHOC.EQ.0)GOTO 270                                             252800
      WRITE(5,926)                                                        252900
      CALL TENS2(2)                                                       253000
 270  IF(OURHOX.EQ.0)GOTO 310                                             253100
      WRITE(5,928)                                                        253200
      CALL TENS2(5)                                                       253300
      WRITE(5,930) ALPHA,BETA,GAMMA                                       253400
C                                                                         253500
C     RETURNS TO MAIN PROGRAM IF NO INTEGRATION OVER SCATTERING ANGLE     253600
C     IS PERFORMED                                                        253700
 310  IF(SIMP.EQ.0.0) RETURN                                              253800
C     INTEGRATION OVER SCATTERING ANGLE                                   253900
C                                                                         254000
      IF( TLBDG.NE.0.0) GOTO 380                                          254100
C     TOTAL CROSS SECTIONS FROM CM DIFFERENTIAL CROSS SECTIONS            254200
      DO 330 N=2,NMAX                                                     254300
      SIGTOT(N)=2.0*3.1415926*DTHETR/3.0*SIMP*DSIG(N)*SIN(THETA/57.29577  254400
     19) + SIGTOT(N)                                                      254500
 330  CONTINUE                                                            254600
C     TOTAL STATISTICAL TENSORS                                           254700
      IF(OURHOC.EQ.0)GOTO 350                                             254800
      CALL TENS2(3)                                                       254900
 350  THETA = THETA - DTHETA                                              255000
      IF(THETA.LT.0.5*DTHETA)GOTO 500                                     255100
      GOTO 440                                                            255200
C     TOTAL CROSS SECTIONS FROM LABORATORY DIFF CROSS SECTIONS            255300
 380  DO 390 N=1,NMAX                                                     255400
      SIGTOT(N)=2.*3.1415926*DTLBR/3.*SIMP*DSIGLB(N)*                     255500
     1SIN(TLBDG/57.295779) + SIGTOT(N)                                    255600
 390  CONTINUE                                                            255700
C     TOTAL STATISTICAL TENSORS                                           255800
      IF(OURHOC.EQ.0)GOTO 410                                             255900
      CALL TENS2(6)                                                       256000
 410  TLBDG = TLBDG-DTLBDG                                                256100
      IF(TLBDG.LT.THI-0.01)GOTO 500                                       256200
      IF(TLBDG.GT.THI+0.0001)GOTO 440                                     256300
      SIMP = 1.0                                                          256400
      GOTO 570                                                            256500
C     ADJUST THE VARIABLE SIMP                                            256600
 440  IF (SIMP.LE.2.001) GOTO 450                                         256700
      IF(SIMP.LE.3.001) GOTO 570                                          256800
      SIMP =2.0                                                           256900
      GOTO 570                                                            257000
 450  SIMP =4.0                                                           257100
      GOTO 570                                                            257200
C                                                                         257300
C     PRINT-OUT OF TOTAL CROSS SECTIONS AND TOTAL STATISTICAL TENSORS     257400
 500  CONTINUE                                                            257500
      IF(TLBDG.NE.0.0) GOTO 510                                           257600
      WRITE(5,938)                                                        257700
      GOTO 520                                                            257800
 510  WRITE(5,940)TLI,THI                                                 257900
 520  WRITE(5,942)                                                        258000
      NMIN=2                                                              258100
      IF(TLBDG.NE.0.0.AND.THI.GT.0.01) NMIN=1                             258200
      DO 530 N=NMIN,NMAX                                                  258300
 530  WRITE(5,944) N,SIGTOT(N)                                            258400
      INTEND = 1                                                          258500
      IF(OURHOC.EQ.0)GOTO 560                                             258600
      WRITE(5,946)                                                        258700
      CALL TENS2(4)                                                       258800
 560  IF(OURHOX.EQ.0)GOTO 570                                             258900
      WRITE(5,947)                                                        259000
      CALL TENS2(12)                                                      259100
 570  RETURN                                                              259200
C                                                                         259300
 902  FORMAT(52H0SCATTERING ANGLE IN CENTER OF MASS SYSTEM  THETA = F7.2  259400
     1,8H DEGREES)                                                        259500
 904  FORMAT(1H0,12H LEVEL INDEX10X10HEXCITATION11X16HCM CROSS SECTION/   259600
     122X13HPROBABILITIES9X15HBARNS/STERADIAN//7X1HN18X4HP(N)18X          259700
     27HDSIG(N))                                                          259800
 906  FORMAT(48H0SCATTERING ANGLE IN LABORATORY SYSTEM  TLBDG = F7.2,     259900
     1 8H DEGREES)                                                        260000
 908  FORMAT(6H0LEVEL4X10HSCATTERING8X10HEXCITATION8X9HLAB CROSS9X        260100
     19HENERGY OF9X9HERERGY OF9X10HLAB RECOIL8X11HSOLID ANGLE/6H INDEX    260200
     24X8HANGLE-CM10X13HPROBABILITIES5X7HSECTION11X9HSCATTERED9X6HRECOIL  260300
     33HING9X5HANGLE13X9HRATIO FOR/10X7HDEGREES29X15HBARNS/STERADIAN      260400
     43X14HPROJECTILE,MEV4X10HTARGET,MEV8X,7HDEGREES11X6HRECOIL//3X1HN6X  260500
     58HTCMDG(N)10X4HP(N)14X9HDSIGLB(N)9X6HEPP(N)12X6HETP(N)12X           260600
     68HZLBDG(N)10X5HR4(N))                                               260700
 910  FORMAT(6H0GROUPI3,/)                                                260800
 912  FORMAT (I8,6X,E19.4,5X,E19.4)                                       260900
 916  FORMAT(1XI3,5XF7.2,7XE14.4,4XE14.4,8XF9.4,9XF9.4,9XF9.4,7XF9.4)     261000
 920  FORMAT(1H0,9X40HTHE STATISTICAL TENSORS RHOB(N,KA,KAPPA))           261100
 926  FORMAT(1H0,9X40HTHE STATISTICAL TENSORS RHOC(N,KA,KAPPA))           261200
 928  FORMAT(1H0,9X40HTHE STATISTICAL TENSORS RHOX(N,KA,KAPPA))           261300
 930  FORMAT(31H0EULERIAN ANGLES USED   ALPHA = F7.2,3X,7HBETA = F7.2,    261400
     13X,8HGAMMA = F7.2)                                                  261500
 938  FORMAT(1H1,16X,29HTOTAL CROSS SECTIONS IN BARNS)                    261600
 940  FORMAT(56H1TOTAL LABORATORY CROSS SECTIONS IN RING COUNTER BETWEEN  261700
     1 F7.2,13H DEGREES AND F7.2, 8H DEGREES,11H (IN BARNS))              261800
 942  FORMAT(1H0,6X1HN18X9HSIGTOT(N) )                                    261900
 944  FORMAT(5XI3,11XE19.4)                                               262000
 946  FORMAT(44H1TOTAL STATISTICAL TENSORS RHOCT(N,KA,KAPPA)//)           262100
 947  FORMAT(44H1TOTAL STATISTICAL TENSORS RHOXT(N,KA,KAPPA)//)           262200
      END                                                                 262300
      SUBROUTINE START                                                    262400
C                                                                         262500
C     THIS ROUTINE STARTS EXECUTION OF PROGRAM                            262600
C                                                                         262700
      LOGICAL ERR                                                         262800
      INTEGER OUAMP,OUPROW,OUXI,OUPSI                                     262900
      INTEGER OURHOB,OURHOC,OURHOX                                        263000
      COMPLEX RHOCT                                                       263100
      REAL MEM                                                            263200
      COMMON /BL1/LAMDA(12),LEAD(NNNN,NNNN,12),LDNUM(12,NNNN),LAMMAX      263300
      COMMON /BL2/NGMAX,NSTART(10),NSTOP(10),MASTER(10,10),EMMAX(10)      263400
      COMMON /BL3/EN(NNNN),SPIN(NNNN),ACCUR,DIPOL,ZPOL,BANDK(NNNN)        263500
      COMMON /BL9/IEXCIT(10),IEXNUM,LMAX(2),IRSTA(2),IRSTO(2),IP          263600
      COMMON/BL10/IPAR(NNNN),IFAC(NNNN)                                   263700
      COMMON/BL19/ERR                                                     263800
      COMMON/BL20/MEM(NNNN,NNNN,12)                                       263900
      COMMON/BL21/XIMAX                                                   264000
      COMMON/BL22/THETA                                                   264100
      COMMON/BL23/IDN17,IDN18,IDN19,IDN25                                 264200
      COMMON/BL30/RHOCT(NNNN,3,5)                                         264300
      COMMON/BL31/OURHOB,OURHOC,OURHOX                                    264400
      COMMON/BL32/NMAX,NDIM                                               264500
      COMMON/BL33/NTIME                                                   264600
      COMMON/BL34/SIMP                                                    264700
      COMMON/BL35/ZP,ZT,A1,A2,EP,TLBDG                                    264800
      COMMON/BL45/DTHETR,DTLBR                                            264900
      COMMON/BL48/NANGLE,THI,TLI                                          265000
      COMMON/BL49/OUAMP,OUPROW,INTERV,INTIN                               265100
      COMMON/BL50/OUXI,OUPSI,NCM,EMMAX1                                   265200
      COMMON/BL51/SIGTOT(NNNN)                                            265300
      COMMON/BL55/DTHETA,DTLBDG                                           265400
C                                                                         265500
      INTERV = INTIN                                                      265600
      THETA = THI                                                         265700
      TLBDG = TLI                                                         265800
C                                                                         265900
C     CHECK AND SYMMETRIZATION OF MEM-MATRIX BY SUBROUTINE MAT            266000
      CALL MAT                                                            266100
      IF(ERR) RETURN                                                      266200
      IF (NGMAX.NE.1) GOTO 110                                            266300
C                                                                         266400
C     CONSTRUCTION OF GROUPS FOR NMAX .GT. 11 BY THE PROGRAM              266500
      NGMAX = NMAX/11                                                     266600
      NGMAX = NGMAX + 1                                                   266700
      IF(MOD(NMAX,11).EQ.0) NGMAX=NMAX/11                                 266800
      NSTA = 1                                                            266900
      NSTO = 11                                                           267000
      DO 100 I=1,NGMAX                                                    267100
      IEXCIT(I) = IEXCIT(1)                                               267200
      EMMAX(I) = EMMAX1                                                   267300
      NSTART(I) = NSTA                                                    267400
      NSTOP(I) = NSTO                                                     267500
      IF(I.EQ.NGMAX) NSTOP(I)=NMAX                                        267600
      NSTA = NSTA+11                                                      267700
      NSTO = NSTO+11                                                      267800
 100  CONTINUE                                                            267900
      GOTO 180                                                            268000
C                                                                         268100
C     SYMMETRIZATION OF MASTER - MATRIX                                   268200
 110  DO 170 I1 = 1,NGMAX                                                 268300
      DO 170 I2 = I1,NGMAX                                                268400
      MASTER(I2,I1) = MASTER(I1,I2)                                       268500
 170  CONTINUE                                                            268600
C                                                                         268700
C     COMPUTATION OF IFAC-ARRAY                                           268800
 180  DO 220 K=1,IEXNUM                                                   268900
      IF(K.EQ.2)GOTO 190                                                  269000
      N1 = 1                                                              269100
      N2 = NMAX                                                           269200
      IF(IEXNUM.EQ.2) N2 = NSTART(IP) - 1                                 269300
      GOTO 200                                                            269400
 190  N1 = NSTART(IP)                                                     269500
      N2 = NMAX                                                           269600
 200  DO 210 N=N1,N2                                                      269700
      ISPIN = SPIN(N) + 0.0001                                            269800
      IDPAR = 0                                                           269900
      IF(IPAR(N).NE.IPAR(N1)) IDPAR=1                                     270000
      IFAC(N) = (-1)**(IDPAR-ISPIN)                                       270100
 210  CONTINUE                                                            270200
 220  CONTINUE                                                            270300
      IF (NCM .GT. NMAX) NCM = NMAX                                       270400
C                                                                         270500
C     PRINT-OUT OF INPUT DATA                                             270600
      WRITE(5,900)                                                        270700
      WRITE(5,902) ZP,A1,EP                                               270800
      WRITE(5,904) DIPOL                                                  270900
      IF(TLI.EQ.0.0.AND.THI.EQ.0.0) GOTO 330                              271000
C     STATEMENT 330 PREPARES THE COMPUTATION OF TOTAL CROSS SECTIONS      271100
      IF(TLI.NE.0.0.AND.THI.NE.0.0) GOTO 350                              271200
C     STATEMENT 350 PREPARES THE COMPUTATION OF CROSS-SECTIONS FOR A      271300
C     RING-COUNTER                                                        271400
      SIMP = 0.0                                                          271500
      NANGLE=0                                                            271600
      THETA = THI                                                         271700
      TLBDG = TLI                                                         271800
      IF(TLBDG.NE.0.)GOTO 240                                             271900
      WRITE(5,906) THETA                                                  272000
      GOTO 250                                                            272100
 240  WRITE(5,908) TLBDG                                                  272200
 250  WRITE(5,910) ZT,A2                                                  272300
      WRITE(5,912) NGMAX                                                  272400
      WRITE(5,914)                                                        272500
      DO 260 I=1,NGMAX                                                    272600
      N1 = NSTART(I)                                                      272700
      N2 = NSTOP(I)                                                       272800
      WRITE(5,916) I, (N,N=N1,N2)                                         272900
      WRITE(5,918)(EN(N),N=N1,N2)                                         273000
      WRITE(5,920)(SPIN(N),N=N1,N2)                                       273100
      WRITE(5,922) (IPAR(N),N=N1,N2)                                      273200
 260  CONTINUE                                                            273300
      IF(NGMAX.EQ.1) GOTO 280                                             273400
C     PRINTOUT OF MASTER - MATRIX                                         273500
      WRITE(5,924)                                                        273600
      WRITE(5,926) ( I, I=1,NGMAX)                                        273700
      DO 270 N=1,NGMAX                                                    273800
      WRITE(5,928) N, ( MASTER(N,I), I = 1,NGMAX)                         273900
 270  CONTINUE                                                            274000
 280  CONTINUE                                                            274100
C     PRINT-OUT OF THE MATRIX ELEMENTS                                    274200
      WRITE(5,930)                                                        274300
      WRITE(5,960)                                                        274400
      DO 320 I5=1,LAMMAX                                                  274500
      LAM = LAMDA(I5)                                                     274600
      IF(LAM.GT.6)GOTO 282                                                274700
      WRITE(5,932) LAM                                                    274800
      GOTO 283                                                            274900
 282  LA = LAM - 6                                                        275000
      WRITE(5,962) LA                                                     275100
 283  NN = 1                                                              275200
      IF(NMAX.GT.11) NN=NGMAX                                             275300
      DO 310 I=1,NN                                                       275400
      M1 = NSTART(I)                                                      275500
      M2 = NSTOP(I)                                                       275600
      IF(NMAX.GT.11)GOTO 285                                              275700
      M1 = 1                                                              275800
      M2 = NMAX                                                           275900
 285  WRITE(5,936) (M,M=M1,M2)                                            276000
      DO 290 N=1,NMAX                                                     276100
      K = 0                                                               276200
      DO 286 M=M1,M2                                                      276300
      RMEM = MEM(N,M,LAM)                                                 276400
      IF(RMEM.NE.0..AND.ABS(RMEM).LT.0.001) K=1                           276500
 286  CONTINUE                                                            276600
      IF(K.EQ.1)GOTO 288                                                  276700
      WRITE(5,938) N,(MEM(N,M,LAM),M=M1,M2)                               276800
      GOTO 290                                                            276900
 288  WRITE(5,964) N,(MEM(N,M,LAM),M=M1,M2)                               277000
 290  CONTINUE                                                            277100
      WRITE(5,934)                                                        277200
 310  CONTINUE                                                            277300
 320  CONTINUE                                                            277400
      WRITE(5,940)                                                        277500
      WRITE(5,942)NMAX,INTERV,NCM,XIMAX,ACCUR,NTIME                       277600
      WRITE(5,944) (I,I=1,NGMAX)                                          277700
      WRITE(5,946) (EMMAX(I),I=1,NGMAX)                                   277800
      WRITE(5,948) (IEXCIT(I),I=1,NGMAX)                                  277900
      WRITE(5,950)                                                        278000
      WRITE(5,952) OUXI,OUPSI,OUAMP,OUPROW,OURHOB,OURHOC,OURHOX           278100
      GOTO 400                                                            278200
C                                                                         278300
C     PREPARE THE COMPUTATION OF TOTAL CROSS SECTIONS                     278400
 330  WRITE(5,954)                                                        278500
      RANGLE = 180.0/DTHETA                                               278600
      NANGLE = RANGLE                                                     278700
      IF(ABS(RANGLE-NANGLE).GT.1.E-6) NANGLE=NANGLE+1                     278800
      IF(MOD(NANGLE,2).EQ.0)GOTO 340                                      278900
      NANGLE = NANGLE + 1                                                 279000
 340  DTHETA = 180.0/FLOAT(NANGLE)                                        279100
      WRITE(5,966) DTHETA                                                 279200
      DTHETR = DTHETA/57.295779                                           279300
      NANGLE = NANGLE - 1                                                 279400
      THETA = 180.0 - DTHETA                                              279500
      SIMP = 4.0                                                          279600
      GOTO 380                                                            279700
C                                                                         279800
C     PREPARE THE COMPUTATION OF CROSS-SECTIONS FOR A                     279900
C     RING-COUNTER                                                        280000
 350  IF(TLI.EQ.THI) GOTO 410                                             280100
      IF(TLI.GT.THI)GOTO 355                                              280200
      T = TLBDG                                                           280300
      TLBDG = THETA                                                       280400
      THETA = T                                                           280500
      T = TLI                                                             280600
      TLI = THI                                                           280700
      THI = T                                                             280800
 355  WRITE(5,956)TLI,THI                                                 280900
      DELTHE = TLI - THI                                                  281000
      DTLBDG = DTHETA                                                     281100
      IF(IDN25.NE.0)GOTO 356                                              281200
      IF(DELTHE.LE.20.)GOTO 360                                           281300
 356  RANGLE = DELTHE/DTLBDG                                              281400
      NANGLE = RANGLE                                                     281500
      IF(ABS(RANGLE-NANGLE).GT.1.E-6) NANGLE=NANGLE+1                     281600
      IF(MOD(NANGLE,2).NE.0) NANGLE=NANGLE+1                              281700
      DTLBDG = DELTHE/FLOAT(NANGLE)                                       281800
      NANGLE = NANGLE + 1                                                 281900
      WRITE(5,966) DTLBDG                                                 282000
      DTLBR = DTLBDG/57.295779                                            282100
      IF(TLI.EQ.180.)GOTO 358                                             282200
      SIMP =1.0                                                           282300
      GOTO 380                                                            282400
 358  TLBDG = TLI - DTLBDG                                                282500
      SIMP = 4.0                                                          282600
      GOTO 380                                                            282700
 360  IF(DELTHE.LE.10.)GOTO 370                                           282800
      DTLBDG = 0.6*DELTHE                                                 282900
      DTLBR = DTLBDG/57.295779                                            283000
      TLBDG = TLI-0.2*DELTHE                                              283100
      NANGLE=2                                                            283200
      SIMP=2.5                                                            283300
      GOTO 380                                                            283400
 370  TLBDG = TLI-0.5*DELTHE                                              283500
      DTLBDG = DELTHE                                                     283600
      DTLBR = DTLBDG/57.295779                                            283700
      NANGLE=1                                                            283800
      SIMP=3.0                                                            283900
C                                                                         284000
C     INITIALIZE SIGTOT(N), RHOCT(N,K,KAPPA)                              284100
 380  DO 390 N=1,NMAX                                                     284200
      SIGTOT(N) = 0.0                                                     284300
      DO 390 L=1,5                                                        284400
      DO 390 K=1,3                                                        284500
      RHOCT(N,K,L) = 0.0                                                  284600
 390  CONTINUE                                                            284700
      GOTO 250                                                            284800
 400  RETURN                                                              284900
C                                                                         285000
C     ERROR EXIT                                                          285100
 410  WRITE(5,958)                                                        285200
      ERR = .TRUE.                                                        285300
      RETURN                                                              285400
C                                                                         285500
 900  FORMAT(12H1COULEX 1978,/28H MULTIPLE COULOMB EXCITATION,            285600
     18H PROGRAM,2X,26HFOR  E1 - E6  AND  M1 - M6//)                      285700
 902  FORMAT(32H0PROJECTILE CHARGE NUMBER  ZP = F6.2,14H,  MASS  A1 =     285800
     1F8.3,23HAMU,  LAB ENERGY  EP = F8.3,3HMEV)                          285900
 904  FORMAT('   POLARIZATION-PARAMETER FOR VIRTUAL EXCITATION OF THE,    286000
     1    GIANT DIPOLE RESONANCE DIPOLE',F7.4)                            286100
 906  FORMAT(52H0SCATTERING ANGLE IN CENTER OF MASS SYSTEM  THETA = F7.2  286200
     1,8H DEGREES)                                                        286300
 908  FORMAT(48H0SCATTERING ANGLE IN LABORATORY SYSTEM  TLBDG = F7.2,     286400
     1 8H DEGREES)                                                        286500
 910  FORMAT(28H0TARGET CHARGE NUMBER  ZT = F6.2,14H,  MASS  A2 = F8.3,   286600
     13HAMU/)                                                             286700
 912  FORMAT(27H0NUMBER OF GROUPS   NGMAX =I3)                            286800
 914  FORMAT(///16H0ENERGY SPECTRUM)                                      286900
 916  FORMAT(//3X,5HGROUPI3,9X,13HLEVEL INDEX N,11I8,/40X,11I8)           287000
 918  FORMAT(1H0,20X,13HENERGY IN MEV,11F8.4,/(40X,11F8.4))               287100
 920  FORMAT(1H0,20X,4HSPIN,11X,11(F5.1,3X),/(40X,11(F5.1,3X)))           287200
 922  FORMAT(1H0,20X,6HPARITY,5X,11I8)                                    287300
 924  FORMAT(////16H0MASTER - MATRIX)                                     287400
 926  FORMAT(7X,10I10)                                                    287500
 928  FORMAT(4H0N =I3,10I10)                                              287600
 930  FORMAT(///26H0MULTIPOLE MATRIX ELEMENTS)                            287700
 932  FORMAT(/11H0MATRIX  MEI1,5H(N,M))                                   287800
 934  FORMAT(///)                                                         287900
 936  FORMAT(6X,I14,10I11)                                                288000
 938  FORMAT(1H05X3HN =I2,11F11.4)                                        288100
 940  FORMAT(75H0PERFORMANCE CONTROLS  NMAX,   INTERV,   NCM,   XIMAX,    288200
     1   ACCUR,     NTIME)                                                288300
 942  FORMAT(1H0,19X,I6,2X,2I8,F9.2,F13.7,I9)                             288400
 944  FORMAT(12H0GROUP NO. I,10(6X,I2,2X))                                288500
 946  FORMAT(9H EMMAX(I),3X,10(5X,F5.1))                                  288600
 948  FORMAT(10H IEXCIT(I),2X,10(6X,I2,2X))                               288700
 950  FORMAT(73H0OUTPUT CONTROLS       OUXI, OUPSI, OUAMP, OUPROW, OURHO  288800
     1B, OURHOC, OURHOX)                                                  288900
 952  FORMAT(1H0,17X3I7,4I8)                                              289000
 954  FORMAT(69H0TOTAL CROSS SECTIONS FROM CENTER OF MASS DIFFERENTIAL C  289100
     1ROSS SECTIONS)                                                      289200
 956  FORMAT(49H0TOTAL LABORATORY CROSS SECTIONS BETWEEN TLBDG = F7.2,    289300
     120H DEGREES AND TLBDG = F7.2,8H DEGREES )                           289400
 958  FORMAT(69H0ERROR - ILLEGAL DEFINITION OF TLBDG AND THETA - EXECUTI  289500
     1ON TERMINATED)                                                      289600
 960  FORMAT(36H   MEM(N,M,LAM) IN  E*BARN**(LAM/2),                      289700
     144H  MM(N,M,LAM) IN  MAGNETON*BARN**((LAM-1)/2))                    289800
 962  FORMAT(/11H0MATRIX  MMI1,5H(N,M))                                   289900
 964  FORMAT(1H05X3HN =I2,11E11.3)                                        290000
 966  FORMAT(38H STEP WIDTH FOR INTEGRATION  DTHETA = F6.2,8H DEGREES)    290100
      END                                                                 290200
      SUBROUTINE TENS1                                                    290300
C                                                                         290400
C     THIS ROUTINE COMPUTES THE STATISTICAL TENSORS CHARACTERIZING        290500
C     THE ANGULAR DISTRIBUTION OF THE DECAY GAMMA-RAYS                    290600
C                                                                         290700
      COMPLEX AMP,TE,CI,R                                                 290800
      COMPLEX RHOB,RHOC                                                   290900
      INTEGER OURHOB,OURHOC,OURHOX                                        291000
      INTEGER SSTART,SSTOP                                                291100
      COMMON /BL2/NGMAX,NSTART(10),NSTOP(10),MASTER(10,10),EMMAX(10)      291200
      COMMON /BL3/EN(NNNN),SPIN(NNNN),ACCUR,DIPOL,ZPOL,BANDK(NNNN)        291300
      COMMON /BL9/IEXCIT(10),IEXNUM,LMAX(2),IRSTA(2),IRSTO(2),IP          291400
      COMMON/BL18/CAT(CCCC,3),ICATMX, ISMAX                               291500
      COMMON/BL22/THETA                                                   291600
      COMMON/BL23/IDN17,IDN18,IDN19,IDN25                                 291700
      COMMON/BL24/AMP(CCCC,4)                                             291800
      COMMON/BL29/RHOB(NNNN,3,5),RHOC(NNNN,3,5)                           291900
      COMMON/BL31/OURHOB,OURHOC,OURHOX                                    292000
      COMMON/BL32/NMAX,NDIM                                               292100
      COMMON/BL35/ZP,ZT,A1,A2,EP,TLBDG                                    292200
      COMMON/BL37/TGAMMA,PGAMMA                                           292300
      COMMON/BL38/ALPHA,BETA,GAMMA                                        292400
      COMMON/BL40/SSTART(31),SSTOP(30)                                    292500
      COMMON/BL42/TCMDG(NNNN),ZLBDG(NNNN),R3(NNNN),R4(NNNN),EPP(NNNN),    292600
     1ETP(NNNN)                                                                 
      COMMON/BL47/WDCP,WDCT                                               292700
      COMMON/BL50/OUXI,OUPSI,NCM,EMMAX1                                   292800
      DATA CI/(0.0,1.0)/                                                  292900
C                                                                         293000
C     CALCULATION OF THE STATISTICAL TENSORS RHOB(N,KA,KAPPA) IN COORDIN  293100
C     SYSTEM B                                                            293200
      DO 190 K=1,IEXNUM                                                   293300
      N1 = 2                                                              293400
      N2 = NMAX                                                           293500
      IF(IEXNUM.EQ.2) N2=NSTART(IP)-1                                     293600
      IF(K.EQ.1)GOTO 100                                                  293700
      N1 = NSTART(IP) + 1                                                 293800
      N2 = NMAX                                                           293900
 100  LLMAX = 2.0*(SPIN(N1-1)+1.0)                                        294000
      CE3 = 1.0/(2.0*SPIN(N1-1)+1.0)                                      294100
      DO 180 N=N1,N2                                                      294200
      CE2 = SQRT(2.0*SPIN(N)+1.0)                                         294300
      KAMAX = 2.02*SPIN(N)                                                294400
      IF(KAMAX.GT.4) KAMAX=4                                              294500
      KASTOP = KAMAX + 1                                                  294600
      KA = 0                                                              294700
      DO 170 KAI=1,KASTOP,2                                               294800
C     KAI IS INDEX COUNTER ONLY                                           294900
      KAPSTP = KA + 1                                                     295000
C     CALCULATION OF RHOB FOR KAPPA = K,K-1,......,0                      295100
      KAPPA = KA                                                          295200
      DO 160 KAPPI=1,KAPSTP                                               295300
C     KAPPI IS INDEX COUNTER ONLY                                         295400
      TE = (0.0,0.0)                                                      295500
      IR1 = SSTART(N)                                                     295600
      IR2 = SSTOP(N)                                                      295700
      DO 150 IR=IR1,IR2                                                   295800
      IRP = IR - KAPPA                                                    295900
      IF(IRP.LT.SSTART(N))GOTO 150                                        296000
C     DEFINITION OF THE ARGUMENTS OF THE THREE-J SYMBOL                   296100
      AA1 = SPIN(N)                                                       296200
      BB1 = -CAT(IR,3)                                                    296300
      BB2 = CAT(IRP,3)                                                    296400
      AA3 = KA                                                            296500
      BB3 = KAPPA                                                         296600
      IEX = SPIN(N) + CAT(IR,3) + 0.01                                    296700
      FAC = (-1)**IEX * THREEJ(AA1,BB1,AA1,BB2,AA3,BB3)                   296800
C     SUMMATION OVER THE MAGNETIC SUBSTATES OF THE GROUND STATE           296900
      L = 1                                                               297000
 110  LL = L + L                                                          297100
      LLL = LLMAX - LL                                                    297200
      IF (LLL) 150,130,120                                                297300
C     STATEMENT 120 STARTS SUMMATION FOR NEGATIVE QUANTUM NUMBERS         297400
 120  TE = TE + FAC * CONJG(AMP(IR,L)) * AMP(IRP,L)                       297500
C     STATEMENT 130 COMPLETES SUMMATION FOR MAGNETIC QUANTUM NUMBERS .GE  297600
 130  JR = 2.02*CAT(IR,3)                                                 297700
      JRP = 2.02*CAT(IRP,3)                                               297800
      IRPOS = IR - JR                                                     297900
      IRPPOS = IRP - JRP                                                  298000
      TE = TE + FAC * CONJG(AMP(IRPOS,L)) * AMP(IRPPOS,L)                 298100
      L = L + 1                                                           298200
      GOTO 110                                                            298300
 150  CONTINUE                                                            298400
      KAINDX = KA/2 + 1                                                   298500
      KAPPIN = KAPPA + 1                                                  298600
      RHOB(N,KAINDX,KAPPIN) = CE2*CE3*TE                                  298700
      KAPPA = KAPPA - 1                                                   298800
 160  CONTINUE                                                            298900
      KA = KA + 2                                                         299000
 170  CONTINUE                                                            299100
 180  CONTINUE                                                            299200
 190  CONTINUE                                                            299300
      IF(OURHOC.EQ.0)GOTO 490                                             299400
C                                                                         299500
C     CALCULATION OF THE STATISTICAL TENSORS RHOC(N,KA,KAPPA) IN COORDIN  299600
C     SYSTEM C                                                            299700
      DFARG = (180.0 + THETA)/2.0                                         299800
      DO 480 K=1,IEXNUM                                                   299900
      N1 = 2                                                              300000
      N2 = NMAX                                                           300100
      IF(IEXNUM.EQ.2) N2=NSTART(IP)-1                                     300200
      IEXC=IEXCIT(1)                                                      300300
      IF(K.EQ.1)GOTO 300                                                  300400
      N1 = NSTART(IP) + 1                                                 300500
      N2 = NMAX                                                           300600
      IEXC=IEXCIT(IP)                                                     300700
 300  DO 470 N=N1,N2                                                      300800
      KA = 0                                                              300900
      KAMAX = 2.02*SPIN(N)                                                301000
      IF(KAMAX.GT.4) KAMAX=4                                              301100
      KASTOP = KAMAX + 1                                                  301200
      DO 460 KAI=1,KASTOP,2                                               301300
C     KAI IS INDEX COUNTER ONLY                                           301400
      DJ = KA                                                             301500
      KAINDX = KA/2 + 1                                                   301600
      KAPPA = KA                                                          301700
      KAPSTP = KA + 1                                                     301800
      DO 450 KAPPI=1,KAPSTP                                               301900
C     KAPPI IS INDEX COUNTER ONLY                                         302000
      DMP = KAPPA                                                         302100
      KAPPIN = KAPPA + 1                                                  302200
      SUM = (0.,0.)                                                       302300
C     SUMMATION OVER INDICES .GE. 0  IN EQUATION (I.4.9)                  302400
      KPR = 0                                                             302500
 390  DM = KPR                                                            302600
      KPRI = KPR + 1                                                      302700
      SUM=SUM+CI**(-KPR)*RHOB(N,KAINDX,KPRI)*DJMM(DFARG,DJ,DM,DMP)        302800
      KPR = KPR + 1                                                       302900
      IF(KPR.LE.KA)GOTO 390                                               303000
      IF(KA.EQ.0)GOTO 420                                                 303100
C     SUMMATION OVER NEGATIVE INDICES IN EQUATION (I.4.9)                 303200
      KPR = 1                                                             303300
 410  DM = -KPR                                                           303400
      KPRI = KPR + 1                                                      303500
      SUM=SUM+(-1)**KPR*CI**KPR*CONJG(RHOB(N,KAINDX,KPRI))*               303600
     1DJMM(DFARG,DJ,DM,DMP)                                               303700
      KPR = KPR + 1                                                       303800
      IF(KPR.LE.KA)GOTO 410                                               303900
 420  RHOC(N,KAINDX,KAPPIN) = (-1)**KAPPA * SUM                           304000
      KAPPA = KAPPA - 1                                                   304100
 450  CONTINUE                                                            304200
      KA = KA + 2                                                         304300
 460  CONTINUE                                                            304400
 470  CONTINUE                                                            304500
 480  CONTINUE                                                            304600
 490  RETURN                                                              304700
      END                                                                 304800
      SUBROUTINE TENS2(I)                                                 304900
C                                                                         305000
C     THIS ROUTINE PRINTS THE STATISTICAL TENSORS AND COMPUTES TOTAL      305100
C     STATISTICAL TENSORS                                                 305200
C                                                                         305300
      COMPLEX RHOB,RHOC,RHOX,RHOXT,RHOCT                                  305400
      COMMON /BL2/NGMAX,NSTART(10),NSTOP(10),MASTER(10,10),EMMAX(10)      305500
      COMMON /BL3/EN(NNNN),SPIN(NNNN),ACCUR,DIPOL,ZPOL,BANDK(NNNN)        305600
      COMMON /BL9/IEXCIT(10),IEXNUM,LMAX(2),IRSTA(2),IRSTO(2),IP          305700
      COMMON/BL22/THETA                                                   305800
      COMMON/BL29/RHOB(NNNN,3,5),RHOC(NNNN,3,5)                           305900
      COMMON/BL30/RHOCT(NNNN,3,5)                                         306000
      COMMON/BL32/NMAX,NDIM                                               306100
      COMMON/BL34/SIMP                                                    306200
      COMMON/BL35/ZP,ZT,A1,A2,EP,TLBDG                                    306300
      COMMON/BL38/ALPHA,BETA,GAMMA                                        306400
      COMMON/BL42/TCMDG(NNNN),ZLBDG(NNNN),R3(NNNN),R4(NNNN),EPP(NNNN),    306500
     1ETP(NNNN)                                                                 
      COMMON/BL45/DTHETR,DTLBR                                            306600
      COMMON/BL46/DSIG(NNNN),DSIGLB(NNNN),RU(NNNN)                        306700
      COMMON/BL47/WDCP,WDCT                                               306800
      COMMON/BL51/SIGTOT(NNNN)                                            306900
      DIMENSION RHOX(NNNN,3,5),RHOXT(NNNN,3,5)                            307000
C                                                                         307100
      IF(I.EQ.5) CALL ROTATE(RHOX,RHOC,ALPHA,BETA,GAMMA)                  307200
      IF(I.EQ.12) CALL ROTATE(RHOXT,RHOCT,ALPHA,BETA,GAMMA)               307300
      DO 540 K=1,IEXNUM                                                   307400
      N1 = 2                                                              307500
      N2 = NMAX                                                           307600
      IF(IEXNUM.EQ.2) N2=NSTART(IP)-1                                     307700
      IEXC = IEXCIT(1)                                                    307800
      IF(K.EQ.1)GOTO 80                                                   307900
      N1 = NSTART(IP) + 1                                                 308000
      N2 = NMAX                                                           308100
      IEXC = IEXCIT(IP)                                                   308200
  80  IF(I.EQ.3.OR.I.EQ.6.OR.I.EQ.8.OR.I.EQ.9)GOTO 150                    308300
      IF(IEXC.EQ.1)GOTO 85                                                308400
      WRITE(5,900)                                                        308500
      GOTO 90                                                             308600
  85  WRITE(5,902)                                                        308700
  90  GOTO(100,103,150,104,105,150,150,150,150,150,150,109,150),I         308800
 100  WRITE(5,904)                                                        308900
      GOTO 150                                                            309000
 103  WRITE(5,910)                                                        309100
      GOTO 150                                                            309200
 104  WRITE(5,912)                                                        309300
      GOTO 150                                                            309400
 105  WRITE(5,914)                                                        309500
      GOTO 150                                                            309600
 109  WRITE(5,913)                                                        309700
 150  CONTINUE                                                            309800
      DO 530 N=N1,N2                                                      309900
      KA = 0                                                              310000
      KAMAX = 2.02*SPIN(N)                                                310100
      IF(KAMAX.GT.4) KAMAX=4                                              310200
      KASTOP = KAMAX + 1                                                  310300
      DO 520 KAI=1,KASTOP,2                                               310400
C     KAI IS INDEX COUNTER ONLY                                           310500
      KAINDX = KA/2 + 1                                                   310600
      KAPSTP = KA + 1                                                     310700
      KAPPA = KA                                                          310800
      DO 510 KAPPI=1,KAPSTP                                               310900
C     KAPPI IS INDEX COUNTER ONLY                                         311000
      KAPPIN = KAPPA + 1                                                  311100
      GOTO (200,203,204,205,206,207,208,208,208,208,208,213,208),I        311200
 200  WRITE(5,920) N,KA,KAPPA,RHOB(N,KAINDX,KAPPIN)                       311300
      GOTO 500                                                            311400
 203  RHO = REAL(RHOC(N,KAINDX,KAPPIN))                                   311500
      WRITE(5,920) N,KA,KAPPA,RHO                                         311600
      GOTO 500                                                            311700
 204  CONTINUE                                                            311800
      IF(KAPPA.NE.0) GO TO 500                                            311900
      RHOCT(N,KAINDX,KAPPIN) = RHOCT(N,KAINDX,KAPPIN) + DTHETR/3.0*       312000
     1SIMP*RHOC(N,KAINDX,KAPPIN)*RU(N)*SIN(THETA/57.295779)*2.*3.1415926  312100
      GOTO 500                                                            312200
 205  RHO = REAL(RHOCT(N,KAINDX,KAPPIN))                                  312300
      WRITE(5,920) N,KA,KAPPA,RHO                                         312400
      GOTO 500                                                            312500
 206  WRITE(5,920) N,KA,KAPPA,RHOX(N,KAINDX,KAPPIN)                       312600
      GOTO 500                                                            312700
 207  CONTINUE                                                            312800
      IF(KAPPA.NE.0) GO TO 500                                            312900
      RHOCT(N,KAINDX,KAPPIN)=RHOCT(N,KAINDX,KAPPIN)+DTLBR/3.*SIMP*RHOC    313000
     1(N,KAINDX,KAPPIN)*RU(N)*R3(N)*SIN(TLBDG/57.295779)*2.*3.1415926     313100
      GOTO 500                                                            313200
 208  CONTINUE                                                            313300
      GOTO 500                                                            313400
 213  WRITE(5,920) N,KA,KAPPA,RHOXT(N,KAINDX,KAPPIN)                      313500
      GOTO 500                                                            313600
 500  KAPPA = KAPPA - 1                                                   313700
 510  CONTINUE                                                            313800
      KA = KA + 2                                                         313900
 520  CONTINUE                                                            314000
      IF(I.NE.3.AND.I.NE.6.AND.I.NE.8.AND.I.NE.9)WRITE(5,916)             314100
 530  CONTINUE                                                            314200
 540  CONTINUE                                                            314300
      RETURN                                                              314400
 900  FORMAT(27H0TENSORS FOR TARGET NUCLEUS)                              314500
 902  FORMAT(31H0TENSORS FOR PROJECTILE NUCLEUS)                          314600
 904  FORMAT(7X1HN5X2HKA4X5HKAPPA9X9HREAL RHOB13X9HIMAG RHOB)             314700
 910  FORMAT(7X,1HN5X2HKA4X5HKAPPA9X9HREAL RHOC)                          314800
 912  FORMAT(7X1HN5X2HKA4X5HKAPPA9X10HREAL RHOCT/)                        314900
 913  FORMAT(7X1HN5X2HKA4X5HKAPPA9X10HREAL RHOXT11X10HIMAG RHOXT/)        315000
 914  FORMAT(7X,1HN5X2HKA4X5HKAPPA9X9HREAL RHOX13X9HIMAG RHOX)            315100
 916  FORMAT(1H0)                                                         315200
 920  FORMAT(1H 3I7,2E22.4)                                               315300
      END                                                                 315400
      FUNCTION THREEJ(RJ1,RM1,RJ2,RM2,RJ3,RM3)                            315500
C                                                                         315600
C     THIS ROUTINE COMPUTES THE THREE-J-SYMBOL ACCORDING TO EQ. (II.23.1  315700
C                                                                         315800
      COMMON /BL8/ASQRT,B(200),IEX(200)                                   315900
      J1=2.001*RJ1                                                        316000
      J2=2.001*RJ2                                                        316100
      J =2.001*RJ3                                                        316200
      M1=2.001*RM1                                                        316300
      M2=2.001*RM2                                                        316400
      M =2.001*RM3                                                        316500
      SUM=0.0                                                             316600
      IF((M1+M2+M.NE.0).OR.(J.GT.J1+J2).OR.(J.LT.IABS(J1-J2))) GOTO 400   316700
      IF((IABS(M1).GT.J1).OR.(IABS(M2).GT.J2).OR.(IABS(M).GT.J))GOTO 400  316800
      K2MIN=0                                                             316900
      IF(J-J2+M1.LT.0) K2MIN=-J+J2-M1                                     317000
      IF(J-J1-M2+K2MIN.LT.0) K2MIN=-J+J1+M2                               317100
      K2MAX=J1+J2-J                                                       317200
      IF(J2+M2-K2MAX.LT.0) K2MAX=J2+M2                                    317300
      IF(J1-M1-K2MAX.LT.0) K2MAX=J1-M1                                    317400
      JA=(J1+M1)/2+1                                                      317500
      JB=JA-M1                                                            317600
      JC=(J2+M2)/2+1                                                      317700
      JD=JC-M2                                                            317800
      JE=(J +M )/2+1                                                      317900
      JF=JE-M                                                             318000
      JG=(J1+J2-J)/2+1                                                    318100
      JH=JA+JB-JG                                                         318200
      JI=JC+JD-JG                                                         318300
      JJ=JE+JF+JG-1                                                       318400
      IF (JJ .GT. 200) GOTO 600                                           318500
C     STATEMENT 600 STOPS THE COMPUTATION IF THE RANGE OF FACTORIALS IS   318600
      B1=B(JA)*B(JB)*B(JC)*B(JD)*B(JE)*B(JF)*B(JG)*B(JH)*B(JI)/B(JJ)      318700
      IEX1=IEX(JA)+IEX(JB)+IEX(JC)+IEX(JD)+IEX(JE)+IEX(JF)+IEX(JG)+       318800
     1IEX(JH)+IEX(JI)-IEX(JJ)                                             318900
      IA=K2MIN/2                                                          319000
      IB=JG-IA+1                                                          319100
      IC=JB-IA+1                                                          319200
      ID=JC-IA+1                                                          319300
      IE=JA-JG+IA                                                         319400
      IF=JD-JG+IA                                                         319500
      FASE=1.0                                                            319600
      IF(MOD(IA,2).EQ.0) FASE=-FASE                                       319700
      K2 =K2MIN                                                           319800
C                                                                         319900
C     SUMMATION OVER K STARTS HERE                                        320000
 100  IA=IA+1                                                             320100
      IB=IB-1                                                             320200
      IC=IC-1                                                             320300
      ID=ID-1                                                             320400
      IE=IE+1                                                             320500
      IF=IF+1                                                             320600
      FASE=-FASE                                                          320700
      B2=B(IA)*B(IB)*B(IC)*B(ID)*B(IE)*B(IF)                              320800
      IEX2=IEX(IA)+IEX(IB)+IEX(IC)+IEX(ID)+IEX(IE)+IEX(IF)                320900
      STOR=SQRT(B1/(B2*B2) * 10.0**(IEX1-IEX2-IEX2))                      321000
      SUM = SUM + FASE*STOR                                               321100
      K2=K2+2                                                             321200
      IF(K2.LE.K2MAX) GOTO 100                                            321300
      N = ABS(RJ1-RJ2-RM3)+0.01                                           321400
      SUM = (-1)**N * SUM                                                 321500
 400  THREEJ=SUM/ASQRT                                                    321600
      RETURN                                                              321700
C                                                                         321800
C     ERROR EXIT                                                          321900
 600  WRITE(5,902)                                                        322000
      WRITE(5,904)RJ1,RJ2,RJ3,RM1,RM2,RM3                                 322100
      STOP                                                                322200
C                                                                         322300
 902  FORMAT(29H0 ERROR - FACTORIALS EXCEEDED)                            322400
 904  FORMAT(6H0RJ1 =F6.2,6H RJ2 =F6.2,6H RJ3 =F6.2,6H RM1 =F6.2,         322500
     16H RM2 =F6.2,6H RM3 =F6.2)                                          322600
      END                                                                 322700
      SUBROUTINE ZCLOCK(T)                                                322800
C                                                                         322900
C     THIS ROUTINE CALLS AN INTERNAL TIME ROUTINE SPECIFIC TO             323000
C     THE USED COMPUTER. FOR OTHER COMPUTERS CHANGE OR REPLACE            323100
C     THE TIME ROUTINE                                                    323200
C                                                                         323300
C     CALL SECOND(T)                                                      323400
      RETURN                                                              323500
      END                                                                 323600
