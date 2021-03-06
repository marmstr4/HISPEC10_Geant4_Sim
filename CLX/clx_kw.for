C
C	MULTIPLE COULOMB EXCITATION PROGRAM  C L X
C	
C	FORTRAN VERSION OF THE PL1 PROGRAM WALKUERE WRITTEN BY H. OWER (IKF).
C	TRANSLATED AND MODIFIED BY J. GERL.
C
c       modified by Th. Kroell:
c       12.10.94     maximum number of states increased 
c                    this number is now given by the parameter NSTATMAX
c                    32767 -> 10000
c                      
C	INPUT TO THE PROGRAM:
C
C		THE FILE  C L X . D A T
C
C	THIS FILE CONTAINS THE INPUT DATASET. THE MEANING OF THE INPUT
C	CARDS IS AS FOLLOWS.
C	( FREE FORMAT IS USED UNLESS OTHERWISE NOTED )
C
C	CARD #	  CONTENTS
C
C	  1	  TITLE ( MAX. 40 CHARACTERS )
C
C	  2	  I/O CONTROL ( FORMAT: 8I1 )
C		  MEANING OF THE CONTROL BITS:
C		   #	MEANING
C		   1	OUCNTL
C		   2	OUPSI
C		   3	OUAMP
C		   4	OUPROW
C		   5	OUMEM
C		   6	OUTENS
C		   7	OUPN
C		   8	OUDCY
C
C	  3	  NMAX, NMPOL, [EMMAX (DEF:10.0)]
C	
C	  4	  [ACCUR (DEF:0.0001)] [,XIMAX (DEF:6.0)] [,DIPOLE (DEF:0.005)]
C		  E1-POLARIZATION ACCORDING TO "ALDER WINTHER, ELECTROMAGNETIC
C		  EXCITATION (APPENDIX J)"
C	  
C	  5	  ZP, AP
C
C	  6	  ZT, AT
c
c                 01.07.01 TK   new card introduced
c         6a      IEXC  1: projectile excitation 2: target excitation 
C
C	  7	  EP ( MeV )
C
C	  8	  THETA1 [,THETA2 [,DELTATH]] ( CM - SYSTEM )
C
C	  9	  N, SPIN(N), E(N) ( MeV), PARITY(N) (1=+,-1=-), K(N)
C		  A CARD #9 IS REQUIRED FOR ALL STATES N=1,NMAX
C
C	 10	  N, M, ME(LAM,N,M), LAMDA (NEG. VALUE = MAGNETIC TRANSITION)
C		  A CARD #10 IS REQUIRED FOR EACH MATRIXELEMENT. IT IS
C		  SUFFICIENT TO USE ONLY THE MATRIXELEMENT ME(N->M).
C		  THE CORRESPONDING MATRIXELEMENT ME(M->N) IS INCLUDED
C		  BY THE PROGRAM. THE RPT PHASE CONVENTION IS USED ( SEE:
C		  "ALDER WINTHER, ELECTROMAGNETIC EXCITATION" ).
C
C	 11	  BLANK CARD
C
        parameter (nstatmax=100)
	CHARACTER*1 TITLE(40)
	REAL LLMAX,MEM(NSTATMAX,NSTATMAX,5),MEMC(NSTATMAX,20,5)
        real ME,R_QLM(5,5),I_QLM(5,5)
	REAL*8 FCT(99),FACT(99)
	REAL*8 R_AMP(2000,8),I_AMP(2000,8),R_AMPDOT(2000,8),I_AMPDOT(2000,8)
	REAL*8 R_F(2000,4,8),I_F(2000,4,8),R_AMPP(2000,8),I_AMPP(2000,8)
	REAL*8 R_Q1(2000,8),I_Q1(2000,8),R_RC,I_RC,R_RC1,I_RC1,I_TE,I_CSUM
	REAL*8 R_TE,R_CSUM
	REAL I_EX,N_FAC
	REAL R_RHOB(NSTATMAX,3,5),I_RHOB(NSTATMAX,3,5),R_RHOC(NSTATMAX,3,5)
        real I_RHOC(NSTATMAX,3,5)
	INTEGER iFOUND,ZP,ZT,SSTOP(NSTATMAX),SSTART(NSTATMAX)
        integer S_HALT(NSTATMAX),DWFLAG
	INTEGER OK,OUCNTL,OUPSI,OUAMP,OUPROW,OUMEM,OUTENS,OUPN,OUDCY,R_CPHZ
	DIMENSION SPIN(NSTATMAX),EN(NSTATMAX),IPAR(NSTATMAX),RK(NSTATMAX)
        dimension LAMDA(5),ETAN(NSTATMAX),CPSI(5)
	DIMENSION IFAC(NSTATMAX),LDNUM(5,NSTATMAX),LEAD(NSTATMAX,20,5)
        dimension CAT(2000),PSI(NSTATMAX,20,5)
	DIMENSION MRANGE(2,100000),ZETA(100000),PROB(100000)
        dimension P(NSTATMAX),MECNT(5)
C
	COMMON /QLM/ W,EPS,EROOT,LAMMAX,LAMDA,R_QLM,I_QLM,RALFA,ZPOL
	COMMON /AMP/ NMAX,SSTART,S_HALT,LEAD,LDNUM,ETAN
	COMMON /FCTRLS/ FCT,FACT
	COMMON /PAGE/ TITLE,ID
C
C	ASSIGN DEFAULT VALUES TO CERTAIN VARIABLES
C
	PI=3.14159265
	ID=1
	CALL fFCTRLS
C
C	READ IN DATA FROM FILE  C L X
C
	OPEN ( UNIT=1,FILE='clx.dat',STATUS='OLD')
	READ ( 1,1 ) TITLE
1	FORMAT ( 40A1 )
	READ ( 1,2 ) OUCNTL,OUPSI,OUAMP,OUPROW,OUMEM,OUTENS,OUPN,OUDCY
2	FORMAT ( 8I1 )
	OUCNTL=1
	OUPN=1
	OUDCY=1
	READ ( 1,* ) NMAX,NMPOL,EMMAX1
	IF ( EMMAX1 .EQ. 0. ) EMMAX1=10.
	READ ( 1,* ) ACCUR,XIMAX,DIPOL
	IF ( ACCUR .EQ. 0. ) ACCUR=0.0001
	IF ( XIMAX .EQ. 0. ) XIMAX=6.
	IF ( DIPOL .EQ. 0. ) DIPOL=0.005
	READ ( 1,* ) ZP,AP
	READ ( 1,* ) ZT,AT

        read(1,*) iexc
        if ((iexc.ne.1).and.(iexc.ne.2)) then
          write(6,90)
 90     format( ' ***MISTAKE*** CARD #6A NO NUCLEUS TO EXCITE CHOSEN' )
          stop
        endif   

	READ ( 1,* ) EP
	READ ( 1,* ) THETA1,THETA2,DTHETA
	IF ( THETA2 .EQ. 0 ) THETA2=THETA1
	IF ( DTHETA .EQ. 0. ) DTHETA=THETA2-THETA1
	DO I=1,NMAX
	  READ ( 1,* ) N,SPIN(N),EN(N),IPAR(N),RK(N)
	  OK=OK+I-N
	END DO
	IF ( OK .NE. 0 ) THEN
	 WRITE ( 6,3 )
3	 FORMAT ( ' ***MISTAKE*** CARDS #9 ARE NOT CONSISTENT WITH NMAX' )
	 STOP
	END IF
	IF ( SPIN(1) .GT. 7. ) THEN
	 WRITE ( 6,* ) ' ***FATAL*** GROUNDSTATE SPIN EXCEEDS 14/2'
	 STOP
	END IF
	IF ( EN(1) .NE. 0. ) WRITE ( 6,4 )
4	FORMAT ( ' ***WARNING*** GROUNDSTATE ENERGY IS NOT ZERO' )
	DO WHILE ( OK .EQ. 0 )
	  READ ( 1,*,END=100 ) N,M,ME,LAM
	  IF ( LAM .LT. 0 ) THEN
	   LAM=5
	   MAGN=1
	   END IF
	  IF ( N .GT. M ) THEN
	   NT=N
	   N=M
	   M=NT
	   IF ( MOD(ABS(INT(SPIN(N)-SPIN(M))),2) .NE. 0 ) ME=-ME
	  END IF
	  IFOUND=0
	    DO I1=1,LAMMAX
	      IF ( LAM .EQ. LAMDA(I1) ) IFOUND=I1
	    END DO
	  IF ( IFOUND .EQ. 0. ) THEN
	   LAMMAX=LAMMAX+1
	   LAMDA(LAMMAX)=LAM
	   IFOUND=LAMMAX
	  END IF
	  MEM(N,M,IFOUND)=ME
	  GO TO 110
100	  OK=1
110	  CONTINUE
	END DO
	CLOSE ( UNIT=1 )
C
C	CALCULATION STARTS HERE
C
C
C	CHECK AND SYMMETRIZE MEM-MATRIX
C
	DO I1=1,LAMMAX
	  LA=LAMDA(I1)
	  IF ( LA .EQ. 5 ) LA=1
	  DO N=1,NMAX
	    DO M=N,NMAX
	      IF ( MEM(N,M,I1) .NE. 0. ) THEN
	       MECNT(I1)=MECNT(I1)+1
	       IF (ABS(SPIN(N)-LA).GT.SPIN(M).OR.SPIN(M).GT.SPIN(N)+LA) THEN
		WRITE ( 6,5 )
5	FORMAT ( ' ***MISTAKE*** ANGULAR MOMENTUM CONSERVATION VIOLATED' )
		STOP
	       END IF
	       IDP=-1
	       IF ( IPAR(N) .EQ. IPAR(M) ) IDP=1
	       IF ( LAMDA(I1) .EQ. 5 ) IDP=-IDP
	       IF ( IDP .NE. (-1)**LA ) THEN
	        WRITE ( 6,6 )
6	FORMAT ( ' ***MISTAKE*** PARITY CONSERVATION VIOLATED' )
		STOP
	       END IF
	       IF ( M .NE. N ) THEN
		MECNT(I1)=MECNT(I1)+1
		FAC=(-1)**INT(ABS(SPIN(N)-SPIN(M))+.01)
		IF ( MEM(N,M,I1) .NE. 0. ) MEM(M,N,I1)=MEM(N,M,I1)*FAC
	       END IF
	      END IF
	    END DO
	  END DO
	END DO
C
C	COMPUTATION OF IFAC-ARRAY
C
	DO N=1,NMAX
	  IDPAR=0
	  IF ( IPAR(N) .NE. IPAR(1) ) IDPAR=1
	  IFAC(N)=(-1)**(IDPAR+INT(SPIN(1)-SPIN(N)))
	END DO
C
C	ACCURACY CONTROL
C
	ACC050=ACCUR/50.
C
C	NORMALIZATION FACTORS
C
	N_FAC=2.*SPIN(1)+1.
	I_SPIN=MOD(INT(2.*SPIN(1)),2)
	LLMAX=2.*(SPIN(1)+1.)
	CE3=1./(2.*SPIN(1)+1.)
C
C	PRINTOUT OF INPUT DATA
C
	CALL NEWPAGE
	WRITE ( 6,7 )
7	FORMAT ( ' PARAMETERS OF THIS CALCULATION :'// )
	WRITE ( 6,8 ) ZP,AP,EP,ZT,AT
8	FORMAT ( ' PROJECTILE : CHARGE = ',I2,' , MASS = ',F7.3,', '
	1,' ENERGY (LAB) = ',F7.2,'MeV'//' TARGET : CHARGE = '
	2,I2,' , MASS = ',F7.3/ )

        if (iexc.eq.2) write(6,91)
 91     format(' TARGET EXCITATION'/)
        if (iexc.eq.1) write(6,92)
 92     format(' PROJECTILE EXCITATION'/)

	WRITE ( 6,9 ) THETA1,THETA2,DTHETA,NMAX
9	FORMAT ( ' SCATTERING ANGLES : ',F6.2,' DEGREES TO ',F6.2,' DEGREES BY '
	1,F6.2,' DEGREES '//' NUMBER OF STATES : ',I2,//' MULTIPOLARITIES : ')
	DO LA=1,LAMMAX
	  IF ( LAMDA(LA) .EQ. 5 ) THEN
	   WRITE ( 6,10 )
10	   FORMAT ( 19X,'M1, ' )
	  ELSE
	   WRITE ( 6,11 ) LAMDA(LA)
11	   FORMAT ( 19X,'E',I1,', ' )
	  END IF
	END DO
	WRITE ( 6,111 ) DIPOL
111	FORMAT ( /' E1-POLARIZATION = ',F7.5 )
	WRITE ( 6,12 ) OUCNTL,OUPSI,OUAMP,OUPROW,OUMEM,OUTENS,OUPN,OUDCY
12	FORMAT(/' OUTPUT CONTROLS : OUCNTL = ',I1,', OUPSI = ',I1,', OUAMP = '
	1,I1,', OUPROW = ',I1/17X,'  OUMEM = ',I1,', OUTENS = ',I1,', OUPN = '
	1,I1,', OUDCY = ',I1 / )
	WRITE ( 6,13 ) EMMAX1
13	FORMAT ( ' MAX. ABSOLUTE VALUE OF M-QUANTUM NUMBER : ',F4.1/ )
	WRITE ( 6,14 ) ACCUR
14	FORMAT ( ' ACCURACY OF INTEGRATION : ',F8.6 / )
	DISTA=2*.71995*(1.+AP/AT)*ZP*ZT/EP
	VINF=4.6337*SQRT(EP/AP)
	WRITE ( 6,15 ) DISTA,VINF
15	FORMAT ( ' SMALLEST DISTANCE BETWEEN NUCLEAR CENTERS IN HEAD-ON'
	1/' COLLISIONS = ',F6.2,' fm, INITIAL PROJECTILE VELOCITY =',F4.1,'% c'/)
C
C	OUTPUT FOR DECAY CALCULATION
C
	WRITE(6,*)'OUDCY.NE.0 , 1 '
	IF ( OUDCY .NE. 0 ) THEN
	WRITE(6,*)'OUDCY.NE.0 , 2 '
	 OPEN ( UNIT=2,FILE='clxdcy.dat',STATUS='NEW' )
C	 WRITE ( 2,1 ) TITLE
	 WRITE ( 2,* ) ' '
	 WRITE ( 2,* ) ZP,AP,ZT,AT,iexc,EP,DIPOL,THETA1,THETA2,DTHETA
	 WRITE ( 2,* ) NMAX,LAMMAX
	 DO I=1,NMAX
	   WRITE ( 2,* ) SPIN(I),EN(I),IPAR(I),RK(I)
	 END DO
	 DO I1=1,LAMMAX
	   WRITE ( 2,* ) LAMDA(I1),MECNT(I1)
	   DO N=1,NMAX
	     DO M=1,NMAX
	       IF ( MEM(N,M,I1) .NE. 0. ) WRITE ( 2,16 ) N,M,MEM(N,M,I1)
16	       FORMAT ( 3(I2,1X,I2,2X,G12.5,2X) )
	     END DO
	   END DO
	 END DO
	END IF
C
C	E1-POLARIZATION ( SEE ALDER WINTHER,ELECTROMAGNETIC EXCITATION APP. J )
C
	ZPOL=DIPOL*EP*AT/(ZT*ZT*(1.+AP/AT))
C
C	COMPUTATION OF LEAD-MATRIX AND LDNUM-MATRIX
C
	DO I5=1,LAMMAX
	  DO N=1,NMAX
	    LDNUM(I5,N)=0
	    DO M=1,NMAX
	      IF ( MEM(N,M,I5) .NE. 0 ) LDNUM(I5,N)=LDNUM(I5,N)+1
	    END DO
	  END DO
	END DO
C
C	DETERMINE MAXIMUM LDNUM=LDMX
C
	LDMX=0
	DO I5=1,LAMMAX
	  DO N=1,NMAX
	    IF ( LDNUM(I5,N) .GT. LDMX ) LDMX=LDNUM(I5,N)
	  END DO
	END DO
	IF ( LDMX .GT. 20 ) THEN
	 WRITE ( 6,17 )
17	 FORMAT ( ' ***FATAL*** NUMBER OF ME''S CONNECTING A PAIR OF STATES'
	1' EXCEEDS 20' )
	 STOP
	END IF
	DO I5=1,LAMMAX
	  DO N=1,NMAX
	    LD_CNTR=0
	    DO M=1,NMAX
	      IF ( MEM(N,M,I5) .NE. 0. ) THEN
	       LD_CNTR=LD_CNTR+1
	       LEAD(N,LD_CNTR,I5)=M
	       MEMC(N,LD_CNTR,I5)=MEM(N,M,I5)
	      END IF
	    END DO
	  END DO
	END DO
C
C	ETAN - ARRAY
C
	ETA=ZP*ZT*SQRT(AP/EP)/6.34977
	DO M=1,NMAX
	  ETAN(M)=ETA/SQRT(1-(1.+AP/AT)*EN(M)/EP)
	END DO
	WRITE ( 6,18 ) ETA
18	FORMAT ( ' SOMMERFELD PARAMETER ETA = ',F6.2 )
C
C	DETERMINE LARGEST XI-VALUE THAT IS .LT. XIMAX
C
	XIM=0.
	DO I1=1,LAMMAX
	  DO N=1,NMAX
	    DO L=1,LDNUM(I1,N)
	      M=LEAD(N,L,I1)
	      XI=ETAN(N)-ETAN(M)
	      IF ( XI .GT. XIM ) XIM=XI
	    END DO
	  END DO
	END DO
	WRITE ( 6,19 ) XIMAX,XIM
19	FORMAT (/' MAXIMUM ALLOWABLE XI = ',G11.4,', LARGEST XI FOUND = ',G11.4)
	XIM=MIN(XIM,XIMAX)
	CALL NEWPAGE
	WRITE ( 6,20 )
20	FORMAT ( ' STATES KNOWN IN THIS CALCULATION'//'  N   ENERGY [MeV]'
	1,'   SPIN    PARITY     K'/ )
	DO N=1,NMAX
	  WRITE ( 6,21 ) N,EN(N),SPIN(N),IPAR(N),RK(N)
21	  FORMAT ( 1X,I2,6X,F6.4,6X,F4.1,5X,I2,7X,F3.1 )
	END DO
C
C	THE MATRIX ELEMENTS ARE PRINTED
C
	IF ( OUMEM .NE. 0 ) THEN
	 DO I5=1,LAMMAX
	   LINE_CNT=NSTATMAX
	   CALL NEWPAGE
	   LAM=LAMDA(I5)
	   IF ( LAM .EQ. 5 ) THEN
	    WRITE ( 6,22 )
22	    FORMAT ( ' MATRIX MM1(N,M) IN mu'/ )
	   ELSE
	    WRITE ( 6,23 ) LAM,LAM
23	    FORMAT ( ' MATRIX ME',I1,'(N,M) IN e*b**',I1/ )
	   END IF
	   DO M=1,NMAX,8
	     L=MIN(8,NMAX-M)
	     M1=M+1
	     M2=M+2
	     M3=M+3
	     M4=M+4
	     M5=M+5
	     M6=M+6
	     M7=M+7
	     IF ( LINE_CNT-NMAX .LT. 1 ) THEN
	      LINE_CNT=NSTATMAX
	      CALL NEWPAGE
	     END IF
	     WRITE ( 6,24 ) M,M1,M2,M3,M4,M5,M6,M7
24	     FORMAT ( '   M',3X,8(I2,6X)/' N' )
	     DO I1=1,NMAX
	       WRITE ( 6,25 ) I1,(MEM(I1,MM,I5),MM=M,M+L)
25	       FORMAT ( 1X,I2,8(1X,F7.4) )
	       LINE_CNT=LINE_CNT-1
	     END DO
	   END DO
	 END DO
	END IF
C
C	PSI - MATRIX
C
	AAZZ=1./(1.+AP/AT)/(ZP*ZT)
	CPSI(1)=5.169286*AAZZ
	CPSI(2)=14.359366*AAZZ*AAZZ
	CPSI(3)=56.982577*AAZZ**3.
	CPSI(4)=263.812653*AAZZ**4.
	CPSI(5)=0.002519*AAZZ

cccc    zz=zp target excitation;  zz=zt projectile excitation
        if (iexc.eq.1) zz=zt
        if (iexc.eq.2) zz=zp

	DO I1=1,LAMMAX
	  LAM=LAMDA(I1)
	  DO N=1,NMAX
	    PP1=(EP-(1.+AP/AT)*EN(N))**.25
	    DO L=1,LDNUM(I1,N)
	      M=LEAD(N,L,I1)
	      PP2=(EP-(1.+AP/AT)*EN(M))**.25
	      IF ( LAM .EQ. 5 ) THEN
	       PSI(N,L,I1)=CPSI(LAM)*zz*(PP1*PP2)**2.*MEMC(N,L,I1)
	      ELSE
	       PSI(N,L,I1)=CPSI(LAM)*zz*SQRT(AP)*(PP1*PP2)**(2*LAM-1)*
	1MEMC(N,L,I1)
	      END IF
	    END DO
	  END DO
	END DO
C
C	CATALOG OF MAGNETIC SUBSTATES
C	SSTART-,SSTOP- AND S_HALT-ARRAY
C
	IS=0
	DO N=1,NMAX
	  IS=IS+2*MIN(SPIN(N),EMMAX1)+1
	END DO
	ISMAX=IS
	IF ( ISMAX .GT. 2000 ) THEN
	 WRITE ( 6,* ) ' ***FATAL*** MAGNETIC SUBSTATES EXCEED 2000!'
	 STOP
	END IF
	IS=1
	DO N=1,NMAX
	  SSTART(N)=IS
	  QUAN=MIN(EMMAX1,SPIN(N))
	  DO RM=-QUAN,QUAN
	    CAT(IS)=RM
	    IS=IS+1
	  END DO
	  SSTOP(N)=IS-1
	END DO
	IF ( SPIN(1) .EQ. 0. .AND. MAGN .NE. 1 ) THEN
	 DO N=1,NMAX
	   S_HALT(N)=SSTART(N)+INT(MIN(SPIN(N),EMMAX1)+.51)
	 END DO
	ELSE
	 DO N=1,NMAX
	   S_HALT(N)=SSTOP(N)
	 END DO
	END IF
	WRITE ( 6,26 ) ISMAX
26	FORMAT ( /' TOTAL NUMBER OF MAGNETIC SUBSTATES = ',I5 )
C
C	LMAX
C
	LMAX=INT(SPIN(1)+1.01)
C
C	COUNT THE ELEMENTS IN THE MRANGE-ARRAY
C
	MR=0
	DO I1=1,LAMMAX
	  DO N=1,NMAX
	    MR=MR+LDNUM(I1,N)*(S_HALT(N)-SSTART(N)+1)
	  END DO
	END DO
C
C	COMPUTE THE MRANGE-ARRAY AND COUNT THE MEMBERS OF THE ZETA-ARRAY
C
	MRMAX=MR
	IF ( MRMAX .GT. 100000 ) THEN
	 WRITE ( 6,*) ' ***FATAL*** MRANGE EXCEEDS 100000!'
	 STOP
	END IF
	MR=0
	NZMAX=0
	DO I1=1,LAMMAX
	  LAM=LAMDA(I1)
	  IF ( LAM .EQ. 5 ) LAM=1
	  DO N=1,NMAX
	    DO I2=1,LDNUM(I1,N)
	      M=LEAD(N,I2,I1)
	      IS1=SSTART(M)
	      DO IR=SSTART(N),S_HALT(N)
	        MR=MR+1
	        ISMIN=INT(MIN(0.,CAT(IR)-CAT(IS1)-REAL(LAM)))
	        ISPLUS=INT(MAX(0.,CAT(IR)-CAT(IS1)-REAL(LAM)))
	        IS2=IS1+ISPLUS-1
	        MRANGE(1,MR)=MIN(2*LAM+1+ISMIN,SSTOP(M)-IS2)
	        MRANGE(2,MR)=IS2
	        NZMAX=NZMAX+MAX(0,MRANGE(1,MR))
	      END DO
	    END DO
	  END DO
	END DO
	IF ( NZMAX .GT. 100000 ) THEN
	 WRITE ( 6,*) ' ***FATAL*** MORE THAN 100000 ELEMENTS IN ZETA ARRAY!'
	 STOP
	END IF
	WRITE ( 6,27 ) NZMAX
27	FORMAT ( /' TOTAL NUMBER OF ELEMENTS IN ZETA ARRAY = ',I5 )
C
C	COMPUTE THE ZETA - MATRIX
C
	NZ=0
	MR=0
	DO I1=1,LAMMAX
	  RLAM=LAMDA(I1)
	  IF ( RLAM .GT. 4. ) RLAM=1.
	  SSQRT=SQRT(2.*RLAM+1.)
	  DO N=1,NMAX
	    DO I2=1,LDNUM(I1,N)
	      M=LEAD(N,I2,I1)
	      PSI_1=PSI(N,I2,I1)*SSQRT
	      DO IR=SSTART(N),S_HALT(N)
	        MR=MR+1
		DO I3=1,MRANGE(1,MR)
		  NZ=NZ+1
		  IS=MRANGE(2,MR)+I3
		  G2=CAT(IS)-CAT(IR)
		  IIEX=SPIN(M)-CAT(IS)
		  IF ( MOD(IIEX,2) .EQ. 0 ) THEN
		   PHZ=1.
		  ELSE
		   PHZ=-1.
		  END IF
		  ZETA(NZ)=PHZ*PSI_1*
	1THREEJ(SPIN(M),-CAT(IS),RLAM,G2,SPIN(N),CAT(IR))
	        END DO
	      END DO
	    END DO
	  END DO
	END DO
C
C	THE PSI - MATRIX IS PRINTED
C
	IF ( OUPSI .NE. 0 ) THEN
	 DO I1=1,LAMMAX
	   LINE_CNT=NSTATMAX
	   LAM=LAMDA(I1)
	   DO N=1,NMAX
	     DO L=1,LDNUM(I1,N)
	       IF ( LINE_CNT .GT. 49 ) THEN
	        CALL NEWPAGE
	        LINE_CNT=0
	        IF ( LAM .EQ. 5 ) THEN 
	         WRITE ( 6,28 )
28	FORMAT ( '  PSI-MATRIX FOR M1'/ )
	        ELSE
	         WRITE ( 6,29 ) LAM
29	FORMAT ( '  PSI-MATRIX FOR E',I1/ )
	        END IF
	        WRITE ( 6,30 )
30	FORMAT ( '  N      M         PSI(N,M)'/ )
	       END IF
	       M=LEAD(N,L,I1)
	       WRITE ( 6,31 ) N,M,PSI(N,L,I1)
31	FORMAT ( 1X,I2,5X,I2,9X,G11.4 )
	       LINE_CNT=LINE_CNT+1
	     END DO
	   END DO
	 END DO
	END IF
C
C	LOOP OVER THE SCATTERING ANGLES
C
	TH=THETA1-DTHETA
	DO WHILE ( TH .LT. THETA2 )
	TH=TH+DTHETA
C	IF ( OUCNTL .NE. 0 ) THEN
C	 CALL NEWPAGE
C	 WRITE ( 6,32 ) TH
C32	FORMAT ( ' THETA = ',F6.2,' DEGREES' )
C	END IF
	THETA=TH*PI/180.
C
C	RANGE AND STEP WIDTH OF INTEGRATION
C
	EPS=1./SIN(THETA/2.)
	EROOT=SQRT(EPS*EPS-1.)
	CLOSE=DISTA*(1.+EPS)/2.
	UP=LOG(1./(EPS*SQRT(ACCUR)))
	ABW=0.
	DW=40.*(ACCUR**.2)/(10.+48.*XIM+16.*XIM*EPS)
	ISTEP=UP/(DW*8.)+1.
	ISTEP=ISTEP*8
	DW=UP/(REAL(ISTEP)-.25)
	UP=DW*ISTEP
	W=-UP
	D2W=DW+DW
C	IF ( OUCNTL .NE. 0 ) THEN
C	 WRITE ( 6,33 ) EPS,CLOSE
C33	FORMAT ( /' ECCENTRICITY OF HYPERBOLA =',F7.3//
C	1' DISTANCE OF CLOSEST APPROACH IS ',F8.2,' fm' )
C	 WRITE ( 6,34 ) UP,UP,D2W
C34	FORMAT ( /' RANGE OF INTEGRATION : W = -',F6.2,' TO +',F6.2//
C	1' INITIAL STEP WIDTH, D2W = ',F8.5 )
C	END IF
C
C	------- INTEGRATION STARTS HERE -------
C
C	INITIAL VALUES OF AMPLITUDES AMP(W=-UP) : ALL AMPLITUDES
C	ARE 0, EXCEPT THE REAL PART OF THE GROUND STATE AMPLITUDE
C	WHICH IS = 1
C
	DO IR=1,ISMAX
	  DO L=1,LMAX
	    R_AMP(IR,L)=0.
	    I_AMP(IR,L)=0.
	  END DO
	END DO
	DO L=1,LMAX
	  R_AMP(L,L)=1.
	END DO
C
C	THR RUNGE-KUTTA-GILL INTEGRATION PROCEDURE 
C
C	COMPUTE STARTING VALUES OF DERIVATIVES
C
	 DO WHILE ( W .LE. UP )
	   DWFLAG=0
	   CALL Q
	   CALL A (R_AMPDOT,I_AMPDOT,R_AMP,I_AMP,MRANGE,CAT,ZETA,
	1LMAX,ISMAX,MRMAX,NZMAX)
	   DO L=1,LMAX
	     DO N=1,NMAX
	       DO IR=SSTART(N),S_HALT(N)
	         R_F(IR,1,L)=R_AMPDOT(IR,L)
	         I_F(IR,1,L)=I_AMPDOT(IR,L)
	       END DO
	     END DO
	   END DO
C
C	COMPUTE THREE MORE STARTING VALUES FOR ADAMS/MOULTON
C
	 DO NAM=2,4
	   DO L=1,LMAX
	     DO N=1,NMAX
	       DO IR=SSTART(N),S_HALT(N)
	         R_Q1(IR,L)=DW*R_AMPDOT(IR,L)
	         I_Q1(IR,L)=DW*I_AMPDOT(IR,L)
	         R_AMP(IR,L)=R_AMP(IR,L)+R_Q1(IR,L)
	         I_AMP(IR,L)=I_AMP(IR,L)+I_Q1(IR,L)
	       END DO
	     END DO
	   END DO
	   IF ( SPIN(1) .EQ. 0. .AND. MAGN .NE. 1 ) THEN
C
C	USE SYMMETRY RELATION
C
	    DO L=1,LMAX
	      DO N=1,NMAX
	        DO IR=SSTART(N),S_HALT(N)
	          IR1=IR-2*CAT(IR)
	          R_AMP(IR1,L)=IFAC(N)*R_AMP(IR,L)
	          I_AMP(IR1,L)=IFAC(N)*I_AMP(IR,L)
	        END DO
	      END DO
	    END DO
	   END IF
	   W=W+DW
	   CALL Q
	   CALL A (R_AMPDOT,I_AMPDOT,R_AMP,I_AMP,MRANGE,CAT,ZETA,
	1LMAX,ISMAX,MRMAX,NZMAX)
	   DO L=1,LMAX
	   DO N=1,NMAX
	     DO IR=SSTART(N),S_HALT(N)
	       R_AMP(IR,L)=R_AMP(IR,L)+.5857864*(DW*R_AMPDOT(IR,L)-R_Q1(IR,L))
	       I_AMP(IR,L)=I_AMP(IR,L)+.5857864*(DW*I_AMPDOT(IR,L)-I_Q1(IR,L))
	       R_Q1(IR,L)=.5857864*DW*R_AMPDOT(IR,L)+.1213204*R_Q1(IR,L)
	       I_Q1(IR,L)=.5857864*DW*I_AMPDOT(IR,L)+.1213204*I_Q1(IR,L)
	     END DO
	   END DO
	   END DO
	   IF ( SPIN(1) .EQ. 0. .AND. MAGN .NE. 1 ) THEN
C
C	USE SYMMETRY RELATION
C
	    DO L=1,LMAX
	      DO N=1,NMAX
	        DO IR=SSTART(N),S_HALT(N)
	          IR1=IR-2*CAT(IR)
	          R_AMP(IR1,L)=IFAC(N)*R_AMP(IR,L)
	          I_AMP(IR1,L)=IFAC(N)*I_AMP(IR,L)
	        END DO
	      END DO
	    END DO
	   END IF
	   CALL A (R_AMPDOT,I_AMPDOT,R_AMP,I_AMP,MRANGE,CAT,ZETA,
	1LMAX,ISMAX,MRMAX,NZMAX)
	   DO L=1,LMAX
	   DO N=1,NMAX
	     DO IR=SSTART(N),S_HALT(N)
	       R_AMP(IR,L)=R_AMP(IR,L)+3.414214*(DW*R_AMPDOT(IR,L)-R_Q1(IR,L))
	       I_AMP(IR,L)=I_AMP(IR,L)+3.414214*(DW*I_AMPDOT(IR,L)-I_Q1(IR,L))
	       R_Q1(IR,L)=3.414214*DW*R_AMPDOT(IR,L)-4.1213204*R_Q1(IR,L)
	       I_Q1(IR,L)=3.414214*DW*I_AMPDOT(IR,L)-4.1213204*I_Q1(IR,L)
	     END DO
	   END DO
	   END DO
	   IF ( SPIN(1) .EQ. 0. .AND. MAGN .NE. 1 ) THEN
C
C	USE SYMMETRY RELATION
C
	    DO L=1,LMAX
	      DO N=1,NMAX
	        DO IR=SSTART(N),S_HALT(N)
	          IR1=IR-2*CAT(IR)
	          R_AMP(IR1,L)=IFAC(N)*R_AMP(IR,L)
	          I_AMP(IR1,L)=IFAC(N)*I_AMP(IR,L)
	        END DO
	      END DO
	    END DO
	   END IF
	   W=W+DW
	   CALL Q
	   CALL A (R_AMPDOT,I_AMPDOT,R_AMP,I_AMP,MRANGE,CAT,ZETA,
	1LMAX,ISMAX,MRMAX,NZMAX)
	   DO L=1,LMAX
	   DO N=1,NMAX
	     DO IR=SSTART(N),S_HALT(N)
	       R_AMP(IR,L)=R_AMP(IR,L)+DW*R_AMPDOT(IR,L)/3.-2.*R_Q1(IR,L)/3.
	       I_AMP(IR,L)=I_AMP(IR,L)+DW*I_AMPDOT(IR,L)/3.-2.*I_Q1(IR,L)/3.
	     END DO
	   END DO
	   END DO
	   IF ( SPIN(1) .EQ. 0. .AND. MAGN .NE. 1 ) THEN
C
C	USE SYMMETRY RELATION
C
	    DO L=1,LMAX
	      DO N=1,NMAX
	        DO IR=SSTART(N),S_HALT(N)
	          IR1=IR-2*CAT(IR)
	          R_AMP(IR1,L)=IFAC(N)*R_AMP(IR,L)
	          I_AMP(IR1,L)=IFAC(N)*I_AMP(IR,L)
	        END DO
	      END DO
	    END DO
	   END IF
	   CALL A (R_AMPDOT,I_AMPDOT,R_AMP,I_AMP,MRANGE,CAT,ZETA,
	1LMAX,ISMAX,MRMAX,NZMAX)
	   DO L=1,LMAX
	   DO N=1,NMAX
	     DO IR=SSTART(N),S_HALT(N)
	       R_F(IR,NAM,L)=R_AMPDOT(IR,L)
	       I_F(IR,NAM,L)=I_AMPDOT(IR,L)
	     END DO
	   END DO
	   END DO
	 END DO
C
C	WE NOW HAVE THE STARTING VALUES FOR THE DERIVATIVES
C	AND CAN PROCEED BY THE FASTER ADAMS-MOULTON ROUTINE
C
C	THE ADAMS-MOULTON ROUTINE FOR SPIN(1)=0
C
	   DO WHILE ( W .LE. UP .AND. DWFLAG .NE. 1 )
	     DO L=1,LMAX
	     DO N=1,NMAX
	       DO IR=SSTART(N),S_HALT(N)
	       R_AMPP(IR,L)=R_AMP(IR,L)+DW/12.*(55.*R_F(IR,4,L)-59.*R_F(IR,3,L)
	1+37.*R_F(IR,2,L)-9.*R_F(IR,1,L))
	       I_AMPP(IR,L)=I_AMP(IR,L)+DW/12.*(55.*I_F(IR,4,L)-59.*I_F(IR,3,L)
	1+37.*I_F(IR,2,L)-9.*I_F(IR,1,L))
	       END DO
	     END DO
	     END DO
	     IF ( SPIN(1) .EQ. 0. .AND. MAGN .NE. 1 ) THEN
C
C	USE SYMMETRY RELATION
C
	      DO L=1,LMAX
	        DO N=1,NMAX
		  DO IR=SSTART(N),S_HALT(N)
	            IR1=IR-2*CAT(IR)
		    R_AMPP(IR1,1)=IFAC(N)*R_AMPP(IR,1)
		    I_AMPP(IR1,1)=IFAC(N)*I_AMPP(IR,1)
		  END DO
	        END DO
	      END DO
	     END IF
	     W=W+DW+DW
	     CALL Q
	     NZ=0
	     MR=0
	     DO L=1,LMAX
	     DO IR=1,ISMAX
	       R_AMPDOT(IR,L)=0.
	       I_AMPDOT(IR,L)=0.
	     END DO
	     END DO
	     DO I1=1,LAMMAX
	       IF ( LAMDA(I1) .NE. 5 ) THEN
	       DO N=1,NMAX
		 DO I2=1,LDNUM(I1,N)
		   M=LEAD(N,I2,I1)
		   XI=ETAN(N)-ETAN(M)
		   R_EX=COS(XI*RALFA)
		   I_EX=SIN(XI*RALFA)
		   DO IR=SSTART(N),S_HALT(N)
		     MR=MR+1
		     DO I3=1,MRANGE(1,MR)
		       IS=MRANGE(2,MR)+I3
		       NZ=NZ+1
		       RMU=CAT(IS)-CAT(IR)
		       MUA=INT(ABS(RMU)+1.01)
		       R_RC1=(R_QLM(I1,MUA)*R_EX-I_QLM(I1,MUA)*I_EX)*ZETA(NZ)
		       I_RC1=(R_QLM(I1,MUA)*I_EX+I_QLM(I1,MUA)*R_EX)*ZETA(NZ)
		       DO L=1,LMAX
		         R_RC=R_RC1*R_AMPP(IS,L)-I_RC1*I_AMPP(IS,L)
		         I_RC=R_RC1*I_AMPP(IS,L)+I_RC1*R_AMPP(IS,L)
		         R_AMPDOT(IR,L)=R_AMPDOT(IR,L)+I_RC
		         I_AMPDOT(IR,L)=I_AMPDOT(IR,L)-R_RC
		       END DO
		     END DO
		   END DO
		 END DO
	       END DO
	       ELSE
	       DO N=1,NMAX
	         DO I2=1,LDNUM(I1,N)
		   M=LEAD(N,I2,I1)
		   XI=ETAN(N)-ETAN(M)
		   R_EX=COS(XI*RALFA)
		   I_EX=SIN(XI*RALFA)
		   DO IR=SSTART(N),S_HALT(N)
		     MR=MR+1
		     DO I3=1,MRANGE(1,MR)
		       IS=MRANGE(2,MR)+I3
		       NZ=NZ+1
		       MU=INT(CAT(IS)-CAT(IR)+.01)
		       MUA=ABS(MU)+1
		       R_RC1=(R_QLM(I1,MUA)*R_EX-I_QLM(I1,MUA)*I_EX)*ZETA(NZ)
		       I_RC1=(R_QLM(I1,MUA)*I_EX+I_QLM(I1,MUA)*R_EX)*ZETA(NZ)
		       IF ( MU .LT. 0 ) THEN
		        R_RC1=-R_RC1
		        I_RC1=-I_RC1
		       END IF
		       DO L=1,LMAX
			 R_RC=R_RC1*R_AMPP(IS,L)-I_RC1*I_AMPP(IS,L)
			 I_RC=R_RC1*I_AMPP(IS,L)+I_RC1*R_AMPP(IS,L)
			 R_AMPDOT(IR,L)=R_AMPDOT(IR,L)+I_RC
			 I_AMPDOT(IR,L)=I_AMPDOT(IR,L)-R_RC
		       END DO
	 	     END DO
		   END DO
		 END DO
	       END DO
	       END IF
	     END DO
	     DO L=1,LMAX
	     DO N=1,NMAX
	       DO IR=SSTART(N),S_HALT(N)
		 R_AMP(IR,L)=R_AMP(IR,L)+DW/12.*
	1(9.*R_AMPDOT(IR,L)+19.*R_F(IR,4,L)-5.*R_F(IR,3,L)+R_F(IR,2,L))
		 I_AMP(IR,L)=I_AMP(IR,L)+DW/12.*
	1(9.*I_AMPDOT(IR,L)+19.*I_F(IR,4,L)-5.*I_F(IR,3,L)+I_F(IR,2,L))
	       END DO
	     END DO
	     END DO
	     IF ( SPIN(1) .EQ. 0. .AND. MAGN .NE. 1 ) THEN
C
C	USE SYMMETRY RELATION
C
	      DO L=1,LMAX
	        DO N=1,NMAX
	          DO IR=SSTART(N),S_HALT(N)
		    IR1=IR-2*CAT(IR)
		    R_AMP(IR1,1)=IFAC(N)*R_AMP(IR,1)
		    I_AMP(IR1,1)=IFAC(N)*I_AMP(IR,1)
		  END DO
	        END DO
	      END DO
	     END IF
	     CALL A (R_AMPDOT,I_AMPDOT,R_AMP,I_AMP,MRANGE,CAT,ZETA,
	1LMAX,ISMAX,MRMAX,NZMAX)
	     DO L=1,LMAX
	     DO N=1,NMAX
	       DO IR=SSTART(N),S_HALT(N)
		 R_F(IR,1,L)=R_F(IR,2,L)
		 I_F(IR,1,L)=I_F(IR,2,L)
		 R_F(IR,2,L)=R_F(IR,3,L)
		 I_F(IR,2,L)=I_F(IR,3,L)
		 R_F(IR,3,L)=R_F(IR,4,L)
		 I_F(IR,3,L)=I_F(IR,4,L)
		 R_F(IR,4,L)=R_AMPDOT(IR,L)
		 I_F(IR,4,L)=I_AMPDOT(IR,L)
	       END DO
	     END DO
	     END DO
	     IF ( W+DW .LE. UP ) THEN
C
C	FIND LARGEST AMPP - AMP
C
	      FF=0.
	      DO L=1,LMAX
	      DO IR=1,ISMAX
	        FZR=R_AMPP(IR,L)-R_AMP(IR,L)
	        FZI=I_AMPP(IR,L)-I_AMP(IR,L)
	        FZ=SQRT(FZR*FZR+FZI*FZI)/14.
	        IF ( FZ .GT. FF ) FF=FZ
	      END DO
	      END DO
	      IF ( FF .LE. ACC050 ) THEN
	       DWFLAG=1
	       DW=2.*DW
	       D2W=2.*DW
C	       IF ( OUCNTL .NE. 0 ) THEN
C	        WRITE ( 6,35 ) W,D2W
C35	FORMAT ( ' AT W = ',F7.3,' STEPWIDTH WAS DOUBLED TO BE ',F8.5 )
C	       END IF
	      END IF
	      IF ( FF .GT. ACCUR ) THEN
	       DWFLAG=1
	       DW=DW/2.
	       D2W=2.*DW
C	       IF ( OUCNTL .NE. 0 ) THEN
C	        WRITE ( 6,36 ) W,D2W
C36	FORMAT ( ' AT W = ',F7.3,' STEPWIDTH WAS HALVED TO BE ',F8.5 )
C	       END IF
	      END IF
	     END IF
C
C	THE EXCITATION PROBABILITIES DURING INTEGRATION
C
	     DO IR=1,ISMAX
	       PROB(IR)=0.
	     END DO
	     DO IR=1,ISMAX
	       DO L=1,LMAX
	       PROB(IR)=PROB(IR)+2.*(R_AMP(IR,L)**2.+I_AMP(IR,L)**2.)/N_FAC
	       END DO
	       IF ( I_SPIN .EQ. 0 ) THEN
	        PROB(IR)=PROB(IR)-(R_AMP(IR,LMAX)**2.+I_AMP(IR,LMAX)**2.)/N_FAC
	       END IF
	     END DO
	     DO N=1,NMAX
	       P(N)=0.
	     END DO
	     DO N=1,NMAX
	       DO IR=SSTART(N),SSTOP(N)
		 P(N)=P(N)+PROB(IR)
	       END DO
	     END DO
C
C	TOTAL EXCITATION PROBABILITIES
C
	     PTOT=0.
	     DO N=1,NMAX
	       PTOT=PTOT+P(N)
	     END DO
C
C	CHECK SUM OF TOTAL EXCITATION PROBABILITIES
C
	     IF ( ABS(PTOT-1.) .GT. ABS(ABW) ) ABW=PTOT-1.
	     IF ( OUPROW .NE. 0 ) THEN
C
C	PRINTOUT OF EXCITATION PROBABILITIES
C	FOR CURRENT VALUE OF W
C
	      WRITE ( 6,37 ) W
37	FORMAT ( /' EXCITATION PROBABILITIES AT W = ',F8.4/ )
	      DO N=1,NMAX,8
	 	M=N+7
	        WRITE ( 6,38 ) N,M,(P(NN),NN=N,N+7)
38	FORMAT ( ' P(N),N=',I2,',',I2,8(2X,F8.6) )
	      END DO
	      WRITE ( 6,39 ) PTOT
39	FORMAT ( ' PTOT = ',F8.6 )
	     END IF
	   END DO
	 END DO
C
C	INTEGRATION COMPLETE
C
	IF ( OUCNTL .NE. 0 ) THEN
	 PTOT=PTOT-1.
C	 WRITE ( 6,40 ) ABW,PTOT
C40	 FORMAT ( ' INTEGRATION COMPLETE'//' LARGEST DEVIATION FROM 1 = ',
C	1G11.4/' FINAL DEVIATION FROM 1 = ',G11.4 )
	END IF
C
C	PRINT-OUT OF THE FINAL AMPLITUDES
C
	IF ( OUAMP .NE. 0 ) THEN
	 IF ( SPIN(1) .NE. 0. ) THEN
	  DO L=1,LMAX
	    IR1=SSTART(1)
	    IRL=IR1+L-1
	    CALL NEWPAGE
	    WRITE ( 6,41 ) CAT(IRL)
41	FORMAT ( ' GROUNDSTATE M = ',F4.1,' AMPLITUDES'//' SPIN     M    REAL 
	1AMPLITUDE  IMAG AMPLITUDE        POPULATION'/)
	    DO N=1,NMAX
	      DO IZR=SSTART(N),SSTOP(N)
		POP=R_AMP(IZR,L)**2.+I_AMP(IZR,L)**2.
		WRITE ( 6,42 ) SPIN(N),CAT(IZR),R_AMP(IZR,L),I_AMP(IZR,L),POP
42	FORMAT ( F5.1,3X,F5.1,3X,G11.4,5X,G11.4,5X,G11.4 )
	      END DO
	    END DO
	  END DO
	 ELSE 
	  IR1=SSTART(1)
	  IRL=IR1
	  CALL NEWPAGE
	  WRITE ( 6,43 )
43	FORMAT ( ' FINAL AMPLITUDES'//' SPIN     M    REAL AMPLITUDE  IMAG 
	1AMPLITUDE  POPULATION'/ )
	  DO N=1,NMAX
	    DO IZR=SSTART(N),SSTOP(N)
	      POP=R_AMP(IZR,1)**2.+I_AMP(IZR,1)**2.
	      WRITE ( 6,42 ) SPIN(N),CAT(IZR),R_AMP(IZR,1),I_AMP(IZR,1),POP
	    END DO
	  END DO
	 END IF
	END IF
C
C	COMPUTATION OF THE DIFFERENTIAL CROSS-SECTIONS
C
	IF ( OUPN .NE. 0 ) THEN
C	 CALL NEWPAGE
	 WRITE ( 6,44 ) TH
44	FORMAT ( /,F6.2,$ )
	 DO N=1,NMAX
	  DSIG=P(N)*.000625*SQRT(EP/(EP-(1.+AP/AT)*EN(N)))*DISTA**2.*EPS**4.
	  WRITE ( 6,45 ) DSIG
C	  WRITE ( 6,45 ) P(N)
C         WRITE ( 6,45 ) TH,N,P(N),DSIG
45	FORMAT ( 2X,G11.4,$ )
	 END DO
	END IF
C
C	------- THE ANGULAR DISTRIBUTION TENSORS -------
C
	DO N=2,NMAX
	  CE2=SQRT(2.*SPIN(N)+1.)
	  KAMAX=MIN(INT(2.02*SPIN(N)),4)
	  DO KA=0,KAMAX,2
	    DO KAPPA=KA,0,-1
	    R_TE=0.
	    I_TE=0.
	    DO IR=SSTART(N),SSTOP(N)
	      IRP=IR-KAPPA
	      IF ( IRP .GE. SSTART(N) ) THEN
	       IEX=INT(SPIN(N)+CAT(IR)+.01)
	       FAC=THREEJ(SPIN(N),-CAT(IR),SPIN(N),CAT(IRP),REAL(KA),
	1REAL(KAPPA))
	       IF ( MOD(IEX,2) .NE. 0 ) FAC=-FAC
C
C	SUMMATION OVER GROUND STATE POLARIZATIONS
C
	       JR=INT(2.02*CAT(IR))
	       JRP=INT(2.02*CAT(IRP))
               IRPOS=IR-JR
	       IRPPOS=IRP-JRP
	       DO L=1,LMAX
		 IF ( LLMAX .GT. 2*L ) THEN
		  R_TE=R_TE+FAC*(R_AMP(IR,L)*R_AMP(IRP,L)+
	1I_AMP(IR,L)*I_AMP(IRP,L))
		  I_TE=I_TE+FAC*(R_AMP(IR,L)*I_AMP(IRP,L)-
	1I_AMP(IR,L)*R_AMP(IRP,L))
		 END IF
		 R_TE=R_TE+FAC*(R_AMP(IRPOS,L)*R_AMP(IRPPOS,L)+
	1I_AMP(IRPOS,L)*I_AMP(IRPPOS,L))
		 I_TE=I_TE+FAC*(R_AMP(IRPOS,L)*I_AMP(IRPPOS,L)-
	1I_AMP(IRPOS,L)*R_AMP(IRPPOS,L))
	       END DO
	      END IF
	    END DO
	    R_RHOB(N,INT(KA/2+1),KAPPA+1)=CE2*CE3*R_TE
	    I_RHOB(N,INT(KA/2+1),KAPPA+1)=CE2*CE3*I_TE
	    END DO
	  END DO
	END DO
C
C	PRINT-OUT OF THE ANGULAR DISTRIBUTION TENSORS RHOB
C
	IF ( OUTENS .NE. 0 ) THEN
	CALL NEWPAGE
	WRITE ( 6,46 )
46	FORMAT ( ' THE ANGULAR DISTRIBUTION TENSORS RHOB(N,KA,KAPPA)'//
	1' N    KA  KAPPA    REAL RHOB    IMAG RHOB' )
	DO N=2,NMAX
	  WRITE ( 6,47 )
47	FORMAT ( 2X )
	  KAMAX=2.02*SPIN(N)
	  IF ( KAMAX .GT. 4 ) KAMAX=4
	   DO KA=0,KAMAX,2
	     DO KAPPA=KA,0,-1
	       WRITE ( 6,48 ) N,KA,KAPPA,R_RHOB(N,INT(KA/2+1),KAPPA+1),
	1I_RHOB(N,INT(KA/2+1),KAPPA+1)
48	FORMAT ( I2,4X,I1,5X,I1,5X,G11.4,2X,G11.4 )
	     END DO
	   END DO
	 END DO
	END IF
C
C	ANGULAR DISTRIBUTION TENSORS RHOC(N,KA,KAPPA)
C
	DFARG=(PI+THETA)/2.
C
C	COMPUTE RHOC(N,KA,KAPPA)
C
	DO N=2,NMAX
	  KAMAX=2.02*SPIN(N)
	  IF ( KAMAX .GT. 4 ) KAMAX=4
	  DO KA=0,KAMAX,2
	    DJ=KA
	    DO KAPPA=KA,0,-1
	      DMP=KAPPA
	      R_CSUM=0.
	      I_CSUM=0.
C
C	SUMMATION OVER POSITIVE INDICES
C
	      DO KPR=0,KA
		R_CPHZ=0.
		I_CPHZ=0.
		IF ( MOD(KPR,4) .EQ. 0 ) THEN
		 R_CPHZ=1
		ELSE IF ( MOD(KPR,4) .EQ. 1 ) THEN
		 I_CPHZ=1.
		ELSE IF ( MOD(KPR,4 ) .EQ. 2 ) THEN
		 R_CPHZ=-1
		ELSE IF ( MOD(KPR,4 ) .EQ. 3 ) THEN
		 I_CPHZ=-1
		END IF
	        DM=REAL(KPR)
		DJMMVAL=DJMM(DFARG,DJ,DM,DMP)
		R_CSUM=R_CSUM+(R_CPHZ*R_RHOB(N,REAL(KA/2+1),KPR+1)+
	1I_CPHZ*I_RHOB(N,REAL(KA/2+1),KPR+1))*DJMMVAL
		I_CSUM=I_CSUM+(R_CPHZ*I_RHOB(N,REAL(KA/2+1),KPR+1)-
	1I_CPHZ*R_RHOB(N,REAL(KA/2+1),KPR+1))*DJMMVAL
	      END DO
	      IF ( KA .NE. 0 ) THEN
C
C	SUMMATION OVER NEGATIVE INDICES
C
	       DO KPR=1,KA
		 IF ( MOD(KPR,2) .EQ. 0 ) THEN
		  PHZ=1.
		 ELSE
		  PHZ=-1.
		 END IF
		 R_CPHZ=0
		 I_CPHZ=0
		 IF ( MOD (KPR,4) .EQ. 0 ) THEN
		  R_CPHZ=1
		 ELSE IF ( MOD(KPR,4) .EQ. 1 ) THEN
		  I_CPHZ=1
		 ELSE IF ( MOD(KPR,4) .EQ. 2 ) THEN
		  R_CPHZ=-1
		 ELSE IF ( MOD(KPR,4) .EQ. 3 ) THEN
		  I_CPHZ=-1
		 END IF
		 DM=-REAL(KPR)
		 DJMMVAL=DJMM(DFARG,DJ,DM,DMP)
		 R_CSUM=R_CSUM+PHZ*(R_CPHZ*R_RHOB(N,INT(KA/2+1),KPR+1)+
	1I_CPHZ*I_RHOB(N,INT(KA/2+1),KPR+1))*DJMMVAL
		 I_CSUM=I_CSUM+PHZ*(I_CPHZ*R_RHOB(N,INT(KA/2+1),KPR+1)-
	1R_CPHZ*I_RHOB(N,INT(KA/2+1),KPR+1))*DJMMVAL
	       END DO
	      END IF
	      IF ( MOD(KAPPA,2) .EQ. 0 ) THEN
	       PHZ=1.
	      ELSE
	       PHZ=-1.
	      END IF
	      R_RHOC(N,INT(KA/2+1),KAPPA+1)=PHZ*R_CSUM
	      I_RHOC(N,INT(KA/2+1),KAPPA+1)=PHZ*I_CSUM
	    END DO
	  END DO
	END DO
C
C	PRINT-OUT OF THE ANGULAR DISTRIBUTION TENSORS RHOC
C
	IF ( OUTENS .NE. 0 ) THEN
C	 CALL NEWPAGE
	 WRITE ( 6,49 )
49	FORMAT ( ' THE ANGULAR DISTRIBUTION TENSORS RHOC(N,KA,KAPPA)'//
	1' N   KA   KAPPA    REAL RHOC    IMAG RHOC ' )
	 DO N=2,NMAX
	   WRITE ( 6,47 )
	   KAMMAX=2.02*SPIN(N)
	   IF ( KAMAX .GT. 4 ) KAMAX=4
	   DO KA=0,KAMAX,2
	     DO KAPPA=KA,0,-1
	       WRITE ( 6,48 ) N,KA,KAPPA,R_RHOC(N,INT(KA/2+1),KAPPA+1),
	1I_RHOC(N,INT(KA/2+1),KAPPA+1)
	     END DO
	   END DO
	 END DO
	END IF
C
C	OUTPUT FOR DECAY CALCULATIONS
C
	IF ( OUDCY .NE. 0 ) THEN
	 DO N=2,NMAX
	   DO J1=1,3
	     WRITE ( 2,* ) (R_RHOC(N,J1,J2),J2=1,2*J1-1)
	   END DO
	 END DO
	END IF
	END DO
	CLOSE ( UNIT=2 )
	STOP
	END 
