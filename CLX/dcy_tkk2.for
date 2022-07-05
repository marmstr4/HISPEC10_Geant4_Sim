c
c       modified DCY program to                     -> K.V. 13-apr-1993 <-
c	-  incLude more than one gamma detector
c       -  furthermore it'S noW aLLoWed to define onLy one WinDOW in CLX-
c          caLcuLationS (dtheta=0 => nth=1)
c
c       further modified by Th. Kroell:
c       13.10.94  - Maximum number of states is increased.
c                   This number is now given by the parameter NSTATMAX.
c                   If NSTATMAX exceeds 50 the version CLX_TK.EXE has
c                   to be taken for the excitation calculations.
c                 - The maximum number of Ge detectors is now given by
c                   the parameter NDETMAX.
c       19.10.94  - The Phi angle of the Ge detectors has to be rotated
c                   by 7.05 + 180 degrees due to the rotation of PPAC1.
c                   The windows on the PPAC must be centered around 
c                   0 degree !!
c                   IO(18) opens the possibility only to print transition
c                   between the first 10 states
c       21.10.94  - Desorientation changed -> IO(12)
c                   Abragam & Pound reimplemented (similar to Mauthofer's
c                   implementation to TRISTAN)   
c       28.10.94    io(17) changed
c       30.10.94/31.10.94   "
c       11.11.94    Error in RHOC-Input
c                   No Phi-angle rotation anymore. DCY takes Phi-ranges
c                   NOT centered with respect to 0 degree!
c       29.11.94    Special version for the calculation of gamma-gamma-
c                   correlations DCY_TKK.FOR
c       last update:  14.12.94/16.12.94/20.12.94/21.12.94/22.12.94  
c                     23.12.94
c       04.01.95    IO(19) introduced
c       25.01.95    POPULATIONSCHEME now ok
c       27.03.95    now calculation of gamma-gamma-correlations
c                   following Alder/Winther (1975) p.69-70
c                   nkorri -> nkorrx -> nkorrf
c       28.03.95
c       
c       29.03.95    Revised version DCY_TKK2.FOR with new integration
c                   over PHI of particle
c       30.03./31.03./07.04.
c
c       07.04.      IO(13) switches to "DCY_TK-Mode"  
c       10.04.      stepwidth for integration controlled via parameter
c       18.04.      error in call filly ...
c       19.05.      IO(13) and N_D are now argument of FILLY
c
c*1*    20.12.96    IO(11) used
c*1*                Only the highest state is populated by Coulomb 
c*1*                excitation, all other states ONLY by decay
c
c***    13.01.97    If conversion coefficient equals zero, also the
c***                gamma-transition is ignored!
c***                Correction: EPSQ=0 -> DELTA=0 in some conditions  
c***    14.01./15.01.
ccc
ccc     01.07.01   IEXC introduced projectile/target excitation
ccc

cccc    03.08.05   - version for gcc 3 
cccc               - new file names
cccc               ge_pos_all.dat  : all detectors in array
cccc               ge_pos_dcy.dat  : sub-set of detectors for which cross 
cccc                                 section is calculated
cccc                                 -> angular distributions/correlations 
cccc               ge_pos_korr.dat : detectors on which is gated   
cccc                                 -> angular correlations 

c 
C	THIS PROGRAM CALCULATES THE GAMMA YIELDS OF A DECAYING
C	NUCLEUS AFTER COULOMB EXCITATION.
C	FORMULAS ARE TAKEN FROM THE ARTICLE OF DE BOER AND
C	WINTHER IN ALDER,WINTHER : COULOMB EXCITATION, PP 303 FF
C	AND - AS FAR AS THE KINETICAL PROBLEMS FOR THE TRANSFER
C	FROM THE REST SYSTEM OF THE TARGET NUCLEUS TO THE LAB
C	SYSTEM ARE CONCERNED - FROM THE ARTICLE OF D. SCHWALM
C	IN NP 192 (72), PP 449 FF.
C
C	INPUT TO THE PROGRAM:
C
C		A.) THE FILE  C L X D C Y . D A T
C	
C	THIS PROGRAM COLLABORATES CLOSELY WITH THE COULOMB
C	EXCITATION PROGRAM CLX .
C	DATA THAT HAVE ALREADY BEEN ENTERED FOR  CLX AND ARE
C	ALSO NEEDED BY THIS PROGRAM ( LIKE ENERGIES OR MATRIX-
C	ELEMENTS ) ARE HANDED OVER. THOSE DATA, AS WELL AS THE
C	ANGULAR DISTRIBUTION TENSORS, WHICH ARE CALCULATED BY
C	CLX ARE READ FROM THE FILE CLX.DAT. THE USER NEED NOT
C	WORRY ABOUT THIS FILE, SINCE CLX TAKES CARE OF EVERY-
C	THING.
C
C		B.) THE FILE  D C Y . D A T
C
C	THIS FILE CONTAINS THE CONTROL DATASET. IT MUST BE SUPPLIED
C	BY THE USER. THE MEANING OF THE INPUT CARDS IS AS FOLLOWS
C	( FREE FORMAT IS USED UNLESS OTHERWISE NOTED )
C
C	CARD #       CONTENTS
C
C	  1	     TITLE. MAY CONTAIN ANY CHARACTER STRING WHICH 
C		     IS READ IN FORMAT A(40) AND PRINTET ON TOP OF
C		     EVERY OUTPUT PAGE.
C
C	  2	     I/O CONTROL. DETERMINES WHICH INTERIM RESULTS
C		     ARE PRINTET AND WHICH SPECIAL FEATURE OF THE
C		     PROGRAM IS TO BE USED. THIS CARD IS READ IN 
C		     FORMAT 14I1.
C		     MEANING OF THE CONTROL BITS :
C			
C			#	CAUSES
C			1	CONVERSION COEFFICIENTS TO BE
C				PRINTET
C			2	NOTHING
C			3	LIFETIME OF THE STATE TO BE
C				PRINTET
C			4	THE GK-MATRIX TO BE PRINTET
C			5	=1 THE FK-MATRIX TO BE PRINTET
c                               =2 THE FKK2K1-MATRIX TO BE PRINTET
c                               =3 both TO BE PRINTET
C			6	THE AK-MATRIX TO BE PRINTET
C			7	THE GAMMA CROSSECTION TO BE
C				PRINTET FOR EVERY SCATTERING
C				ANGLE USED IN THE CLX
C				CALCULATION
C			8	THE GAMMA CROSSECTION NORMALIZED
C				TO THE RUTHERFORD CROSSECTION
C				TO BE PRINTET FOR EVERY SCATTERING 
C				ANGLE USED IN THE CLX CALCULATION
C			9	THE SUM OF THE CROSSECTION OF ALL
C				WINDOWS TO BE PRINTET
C			10	CREATION OF OUTPUT FILE ISUPP.DAT
C				CONTAINING THE CALCULATED GAMMA-
C				INTENSITIES
C			11	NOTHING
c*1*                            =1 population of states only by decay of
c*1*                            the highest excited state
C			12    -	READ IN OF SOLID ANGLE CORRECTIONS
C				FOR GAMMA DETECTOR IF EQUAL 1.
C			      -	READ IN OF SOLID ANGLE CORRECTIONS
C				FOR GAMMA DETECTOR AND READ IN FROM
C				FILE DESORIN.DAT THE DESORIENTATION
C				CORRECTION COEFFICIENTS IF EQUAL 2
c                               now coefficients for ALL states are
c                               required (if none is given, 1. is assumed)
c                             - READ IN OF SOLID ANGLE CORRECTIONS
C                               FOR GAMMA DETECTOR AND THE DESORIENTATION
C                               CORRECTION COEFFICIENTS using
c                               Abragam & Pound parametrisation IF EQUAL 3
c                               global LAMDA2 and N_V are used
c                             - READ IN OF SOLID ANGLE CORRECTIONS
C                               FOR GAMMA DETECTOR AND READ IN FROM
C                               FILE DESORIN.DAT THE DESORIENTATION
C                               CORRECTION COEFFICIENTS IF EQUAL 4
c                               now coefficients for ALL states are
c                               required (if none is given, the coefficients
c                               are calculated using Abragam & Pound 
c                               parametrisation (global LAMDA2 and N_V 
c                               are used) 
C			13	=0 DCY_TKK
c                               =1 DCY_TK 
C			14	ANGLES OF GAMMA DETECTOR ARE GIVEN
C				IN THE REST SYSTEM OF THE DECAYING
C				NUCLEUS RATHER THAN IN THE LAB SYSTEM
c
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
C			15      read ge_pos_all.dat and ge_pos_dcy.dat
c				( N_DET, THE_GAM(i),PHI_GAM(i) )
c			16	the SUM of CROSS SectionS of aLL 
c				WinDOWS of aLL gamma detectorS to be 
C				printed
c			17	gamma-intenSitieS to each correSponding
c				ppac-WinDOW( = CLX Scattering angLe)        
C                             - print sum of Ge detectors for each PPAC
c                               window if equal 1
c                             - print sum of PPAC windows for each Ge
c                               if equal 2
c                               04.08.05 not implented -> io(9).eq.1
c                             - print each Ge for each PPAC window
c                               (can produce tremendous amounts of output!)
c                               if equal 3 
c
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c                       18      only transitions between the first 10
c                               states are printed (decreases the amount
c                               of output)
c                       19    - print populationmatrix if 1
c                             - print decayscheme if 2
c                             - both if 3 
c
c *********************************************************************
c
C	  3	     RECTASCENSION ANGLE OF THE GAMMA DETECTOR IN THE
C		     LAB SYSTEM WITH RESPECT TO THE BEAM, AZIMUTHAL
C		     ANGLE OF THE GAMMA DETECTOR IN THE LAB SYSTEM WITH
C		     RESPECT TO THE SCATTERED PROJECTILE, STATE WHOSE
C		     EXCITATION ENERGY IS TO BE USED FOR THE KINEMATICAL
C		     CALCULATIONS IN THE TRANSFER FROM THE REST SYSTEM
C		     OF THE NUCLEUS TO THE LAB SYSTEM
c ***                Comment by J. Gerl: For every window the range
c NOT VALID          of the azimuthal angle MUST be symmetric with
c FOR THIS VERSION   respect to 0 degree!   
C
C	  4	     ONLY IF I/O CONTROL # 12 .NE. 0 :
C		        THE FACTORS Q1, Q2, Q3 TO BE USED FOR THE
C			CORRECTION OF THE FINITE SOLID ANGLE OF THE
C			GAMMA DETECTOR
c                    if I/O control # 12 .ge. 3
c                       the parameters LAMDA2 and N_V for the 
c                       desorientation
C
C	  5	     NUMBER OF WINDOWS ON THE PARTICLE DETECTOR, INDEX
C		     OF THE INITIAL STATE FOR THE NORMALIZATION, INDEX
C		     OF THE FINAL STATE
C
C	  6	     LIMITING RECTASCENSION ANGLES AND AZIMUTHAL ANGLES
C		     OF THE PARTICLE DETECTORS WINDOWS. ( CM-SYSTEM )
C		     THE MAXIMUM NUMBER OF WINDOWS IS 20.
C
C	  7	     NUMBER OF CONVERSION COEFFICIENTS OF THE K-SHELL,
C		     NUMBER OF CONVERSION COEFFICIENTS OF THE L-SHELL,
C		     NUMBER OF CONVERSION COEFFICIENTS OF THE M-SHELL
C		     IN THE HAGER & SELZER TABLES FOR THE CONVERSION
C		     ELECTRONS
C
C	  8	     LOWEST TABULATED ENERGY OF THE K-SHELL,
C		     LOWEST TABULATED ENERGY OF THE K-SHELL,
C		     MINIMUM GAMMA ENERGY GIVEN IN THE H&S TABLE
C		     FOR THE DECAYING NUCLEUS FOR THE K-SHELL,
C		     MAXIMUM GAMMA ENERGY GIVEN IN THE H&S TABLE
C		     FOR THE DECAYING NUCLEUS FOR THE K-SHELL
C
C	  9	     -1,LOWEST TABULATED ENERGY OF THE L2-SHELL,
C		     MINIMUM GAMMA ENERGY GIVEN IN THE H&S TABLE
C		     FOR THE DECAYING NUCLEUS FOR THE TOTAL L-SHELL,
C		     MAXIMUM GAMMA ENERGY GIVEN IN THE H&S TABLE
C		     FOR THE DECAYING NUCLEUS FOR THE TOTAL L-SHELL
C
C	  10	     -1,LOWEST TABULATED ENERGY OF THE M5-SHELL,
C		     MINIMUM GAMMA ENERGY GIVEN IN THE H&S TABLE
C		     FOR THE DECAYING NUCLEUS FOR THE TOTAL M-SHELL,
C		     MAXIMUM GAMMA ENERGY GIVEN IN THE H&S TABLE
C		     FOR THE DECAYING NUCLEUS FOR THE TOTAL M-SHELL
C
C	  11	     TYPE OF TRANSITION FOR THE CONVERSION COEFFICIENT
C		     ( E.G. E2 OR M1 ETC. ) IN COLUMNS 1 AND 2
C		     ( FORMAT A1,I1 )
C
C	  12	     THE H&S CONVERSION COEFFICIENTS FOR THE K-SHELL,
C		     THE H&S CONVERSION COEFFICIENTS FOR THE TOTAL L-SHELL,
C		     THE H&S CONVERSION COEFFICIENTS FOR THE TOTAL M-SHELL
C
C	  13	     BLANK CARD TO INDICATE END OF INPUT STREAM.
C		     NOTE: THE PROGRAM TRIES TO READ CONVERSION COEFFICIENTS
C		     UNTIL A BLANK CARD IS ENCOUNTERED
C
C		
C	AUTHORS OF THE ORIGINAL PL1 VERSION OF THIS PROGRAM ARE H. OWER 
C	AND H.-J. WOLLERSHEIM. TRANSLATED INTO FORTRAN77 AND MODIFIED
C	BY J.GERL.
C      
C
   
c        parameter (nstatmax=100)
        parameter (nstatmax=15)

c       e is array in common block YLM with fixed size e(100) in
c       filly, inty, intyy, and ykkap 
c       ... dirty solution .... parameters should be moved to an .inc file 
        parameter(enstatmax=100)

        parameter (ndetmax=40)

cc      maximal number of angles theta in clxdcy.dat
c        parameter (nthclxmax=50)
        parameter (nthclxmax=300)
 
	DIMENSION ALPHAK(23,5),ALPHAL(17,5),ALPHAM(8,5)
        dimension ALPHA1(NSTATMAX,NSTATMAX,5)
c	DIMENSION ALPHA(23),DELTA(NSTATMAX,NSTATMAX,5),E(NSTATMAX)
	DIMENSION ALPHA(23),DELTA(NSTATMAX,NSTATMAX,5),E(eNSTATMAX)
        dimension EPSQ(NSTATMAX,NSTATMAX,5),LAMDA(5)
	DIMENSION H(NSTATMAX),IO(19),IPAR(NSTATMAX),RK(NSTATMAX)
        dimension GK(NSTATMAX,NSTATMAX,10),Q(3)
	DIMENSION dcyprob(nstatmax,nstatmax)
        dimension THLB(nthclxmax)  
	DIMENSION FF(NSTATMAX,NSTATMAX,10),TPFAC(5),WIN(20,4)
        dimension W(NSTATMAX,NSTATMAX,ndetmax,2),NN(NSTATMAX)
        dimension ESORT(NSTATMAX)
	DIMENSION RHOC(NSTATMAX,10,10)
        dimension AK(nthclxmax,NSTATMAX,10,10)
        dimension ffkkk(3,3,3)
        dimension DS(NSTATMAX,NSTATMAX,20,nthclxmax)
        dimension WINLB(20,2)
	dimension ds_out(NSTATMAX,NSTATMAX,20,nthclxmax)
	DIMENSION CROSS(NSTATMAX,NSTATMAX,ndetmax)
        dimension SCROSS(NSTATMAX,NSTATMAX,ndetmax)
        dimension SCROSS_SUM(NSTATMAX,NSTATMAX)
        dimension scross_ge(nstatmax,nstatmax,20) 
	DIMENSION WSUM(NSTATMAX,NSTATMAX)
	INTEGER Z1,Z2,N_DET,N_D,ideo(nstatmax),nnpop(nstatmax,nstatmax)
        integer nnpos(nstatmax),ind(ndetmax)
c       integer sc1,sc2   never used !!!
	REAL ME(NSTATMAX,NSTATMAX,5),SPIN(NSTATMAX),KA,M1,M2
	REAL*8 FCT(99),FACT(99)
	REAL C,PI,VNORM,C1,C2,C3,C4,asum
	REAL VGAMMA,yyint(2),sig(5),dw(2)
c        real the_gamk(ndetmax),phi_gamk(ndetmax)
        integer in_d(ndetmax),in_dk(ndetmax) 
	reaL THE_GAM(ndetmax),PHI_GAM(ndetmax)
        real qg(nstatmax,3),lamda2,lamda4,n_v
	LOGICAL FOUND
	CHARACTER*1 TITLE(40),TYP,WTITLE(40)

        data sig /1.,-1.,1.,-1.,1./
        
	COMMON /FCTRLS/ FCT, FACT
	COMMON /PAGE/ TITLE
        common /gamdet/ the_gam,phi_gam,in_d,in_dk,ind,ndetmaxx
	COMMON /YLM/ a1,a2,iexc,ep,e,vcm,theta,v1,v2,nint
C
	CALL fFCTRLS
	C=2.997925E+22        ! SQRT(barn)/Sec
	PI=3.141593
	VNORM=1./SQRT(4.*PI)
C
C	------- READ DATA FROM CLX DATASET -------
C
	OPEN ( UNIT=2,FILE='clxdcy.dat',STATUS='UNKNOWN' )
	READ ( 2,1 ) (WTITLE(I),I=1,40)
	READ (2,*) Z1,A1,Z2,A2,iexc,EP,QPAR,THETA1,THETA2,DTHETA
C	Z1 = CHARGE OF PROJECTILE
C	A1 = MASS OF PROJECTILE
C	Z2 = CHARGE OF TARGET
C	A2 = MASS OF TARGET
c       IEXC = PROJECTILE/TARGET EXCITATION
C	EP = LAB. ENERGY OF THE PROJECTILE IN MEV
C	QPAR = E1-POLARISATION
C	THETA1 = LARGEST CM SCATTERING ANGLE
C	THETA2 = SMALLEST CM SCATTERING ANGLE
C	DTHETA = SCATTERING ANGLE INTERVAL
	READ (2,*) NMAX,NMPOL
C	NMAX = NUMBER OF NUCLEAR STATES
C	NMPOL = NUMBER OF MULTIPOLARITIES
	DO I=1,NMAX
	  READ (2,*) SPIN(I),E(I),IPAR(I),RK(I)
	END DO
C	SPIN(N) = SPIN OF STATE N
C	E(N) = ENERGY OF STATE N IN MEV
C	IPAR(N) = PARITY OF STATE N. (1=poS., -1=neg.)
C	RK(N) = K-QUANTUM NUMBER OF STATE N
	DO I1=1,NMPOL
	  READ (2,*) LAMDA(I1),MECNT
	  DO I=1,MECNT
	    READ ( 2,* ) N,M,CMEM
	    ME(N,M,I1)=CMEM
	  END DO
	END DO
C	ME(N,M,L) = MATRIXELEMENTS BETWEEN STATES N AND M OF MULTIPOLARITY L
C
C	------- END OF PARAMETER DATA FROM CLX -------
C	------- THE TENSORS RHOC WILL BE READ LATER -------
C
	OPEN ( UNIT=1,FILE='dcyk.dat',STATUS='UNKNOWN' )
	READ ( 1,1 ) (TITLE(I),I=1,40)
1	FORMAT ( 40A1 )
	CALL NEWPAGE
	WRITE ( 6,2 ) 
2	FORMAT ( ' INPUT DATA PASSED BY  C L X'/ )
	WRITE ( 6,3 ) Z1,A1,Z2,A2,IEXC,EP,THETA1,THETA2,DTHETA
3	FORMAT ( ' CHARGE OF PROJECTILE = ',I2//' MASS OF PROJECTILE = ',F5.1//
	1' CHARGE OF TARGET = ',I2//' MASS OF TARGET = ',F5.1//
	2' EXCITATION OF ',I2//' BOMBARDING ENERGY = ',F7.2//
	3' THETA1 = ',F6.2//' THETA2 = ',F6.2//' DTHETA = ',F6.2/ )       
	WRITE ( 6,4 ) (WTITLE(I),I=1,40)
4	FORMAT ( ' TITLE OF CLX CALCULATION WAS : ',40A1/ )
	WRITE ( 6,5 )
5	FORMAT ( '   N   ENERGY [MeV]   SPIN   PARITY    K'/ )
	DO I=1,NMAX
	  WRITE ( 6,6 ) I,E(I),SPIN(I),IPAR(I),RK(I)
6	  FORMAT ( 2X,I2,5X,F6.4,6X,F4.1,5X,I2,6X,F3.1 )
	END DO
C
C	------- READ DATA FROM THE CONTROL DATASET -------
C
C
C	INPUT AND OUTPUT CONTROLS
C
	READ ( 1,7 ) ( IO(I),I=1,19 )
7	FORMAT ( 19I1 )
C
C	GEOMETRY OF THE SYSTEM
C
	if(io(15).ne.0)then
	    READ ( 1,* ) TT,VG,K1LAB

c           All detectors
c            open(uNIt=11,fiLe='ge_posgasp_dcy.dat',StatuS='oLd')
            open(uNIt=11,fiLe='ge_pos_all.dat',StatuS='oLd')
              read(11,*)N_DET_all
              if (n_det_all.gt.ndetmax) then
                write(6,*) ' Too much Ge detectors, increase NDETMAX!!'
                write(6,*) ' Only NDETMAX detectors will be used '
                n_det_all=ndetmax
              endif
              DO i=1,N_DET_all
                read(11,*)THE_GAM(i),PHI_GAM(i)
              END DO
            cLoSe(uNIt=11)
            ndetmaxx=n_det_all

c           Detectors incremented into cross section
            open(uNIt=11,fiLe='ge_pos_dcy.dat',StatuS='oLd')
	      read(11,*)N_DET
              if (n_det.gt.ndetmaxx) then
                write(6,*) ' Too much Ge detectors !!'
                write(6,*) ' Only NDETMAXX detectors will be used '
                n_det=ndetmaxx
              endif                
	      DO i=1,N_DET
		read(11,*)THE,PHI
                do j=1,ndetmaxx
                  if ((abs(the-the_gam(j)).le.0.1).and.
     +                (abs(phi-phi_gam(j)).le.0.1)) then
                    in_d(i)=j
                    ind(j)=1
                  endif
                end do
                if (in_d(i) .eq. 0) then
                  write(6,*) ' DETECTOR #',i,' IN POS NOT FOUND'
                end if  
	      END DO
	    cLoSe(uNIt=11)

c           Ge detectors on which is gated
            if (IO(13).eq.0) then
            open(uNIt=11,fiLe='ge_pos_korr.dat',StatuS='oLd')
              read(11,*)N_DETk
              if (n_detk.gt.ndetmaxx) then
                write(6,*) ' Too much Ge detectors !!'
                write(6,*) ' Only NDETMAXX detectors will be used '
                n_detk=ndetmaxx
              endif
              DO i=1,N_DETk
                read(11,*)THE,PHI
                do j=1,ndetmaxx
                  if ((abs(the-the_gam(j)) .le. 0.1).and.
     +                (abs(phi-phi_gam(j)) .le. 0.1)) then
                    in_dk(i)=j
                    ind(j)=1
                  endif
                end do
                if (in_dk(i) .eq. 0) then
                  write(6,*) ' DETECTOR #',i,' IN KORR NOT FOUND'
                end if
              END DO
            cLoSe(uNIt=11)
            else
              n_detk=0
            endif   
            
	    
	 eLSe

	    N_DEt=1
            ndetmaxx=1
	    READ ( 1,* ) TT,VG,K1LAB
            the_gam(1)=tt
            phi_gam(1)=vg  
            in_d(1)=1
            ind(1)=1                      
cc            write(6,*) ' *** NOT SUPPORTED !!!!!! ***'
cc            stop 
            if (IO(13).eq.0) then
              write(6,*) ' *** CORRELATION WITH ONE DETECTOR??? ***'
              stop
            endif 

	ENDif

        DO i=1,NDETmaxx
          the_gam(i)=the_gam(i)/180.*pi
          phi_gam(i)=phi_gam(i)/180.*pi
        END DO


c       Gamma-Gamma-Cascade
c       in the populationmatrix nkorri -> nkorrx is regarded 
c       to be the gate
        read(1,*) nkorri,nkorrx,nkorrf


        q(1)=1.0
        q(2)=1.0
        q(3)=1.0 
	IF ( IO(12) .NE. 0 ) READ ( 1,* ) Q(1),Q(2),Q(3)
        IF ( IO(12) .ge. 3 ) READ ( 1,* ) lamda2,n_v


	READ ( 1,* ) MZAHL,NORMI,NORMF,nint
        if (nint .gt. 100) nint=100        
        found=.false.
	DO LAM=1,NMPOL
	  IF ( ME(NORMI,NORMF,LAM) .NE. 0. ) FOUND=.true.
	END DO
	IF ( .NOT. FOUND ) THEN
	 WRITE ( 6,8 )
8	FORMAT ( ' ***MISTAKE*** NORMALIZATION TRANSITION DOES NOT EXIST' )
	 STOP
	END IF
	FOUND=.FALSE.
	IF ( E(NORMI) .LT. E(NORMF) ) THEN
	 NNTEMP=NORMI
	 NORMI=NORMF
	 NORMF=NNTEMP
	END IF
	DO NI=1,MZAHL
	  READ ( 1,* ) WIN(NI,1),WIN(NI,2),WIN(NI,3),WIN(NI,4)
          IF ( WIN(NI,4) .LT. WIN(NI,3) ) THEN
           WTEMP=WIN(NI,4)
           WIN(NI,4)=WIN(NI,3)
           WIN(NI,3)=WTEMP
          END IF
	  IF ( WIN(NI,1) .LT. WIN(NI,2) ) THEN
	   WTEMP=WIN(NI,1)
	   WIN(NI,1)=WIN(NI,2)
	   WIN(NI,2)=WTEMP
	  END IF
	END DO
	CALL NEWPAGE
	WRITE ( 6,9 )
9	FORMAT ( ' INFORMATION OBTAINED FROM CONTROL DATASET'/ )
	WRITE ( 6,10 ) (IO(I),I=1,19)
10	FORMAT ( ' I/O-CONTROL :',19(2X,I1)/ )
        if (io(15).eq.0) then
	WRITE ( 6,11 ) TT,VG
11	FORMAT ( ' GAMMA DETECTOR ANGLES : THETA=',F6.2,' ,PHI=',F7.2/ )
        else
        write ( 6,*) ' GAMMA DETECTOR ANGLES READ FROM GE_POS_ALL.DAT'
        write ( 6,*) ' NUMBER OF DETECTORS : ',n_det_all
        write ( 6,*) ' '
        write ( 6,*) ' GAMMA DETECTORS READ FROM GE_POS_DCY.DAT'
        write ( 6,*) ' NUMBER OF DETECTORS : ',n_det
        write ( 6,'(1x,20i3)') (in_d(i), i=1,20)
        write ( 6,'(1x,20i3)') (in_d(i), i=21,40)
        write ( 6,*) ' '
        write ( 6,*) ' GAMMA DETECTORS READ FROM GE_POS_KORR.DAT'
        write ( 6,*) ' NUMBER OF DETECTORS : ',n_detk
        write ( 6,'(1x,20i3)') (in_dk(i), i=1,20)
        write ( 6,'(1x,20i3)') (in_dk(i), i=21,40)
        write ( 6,*) ' '
        write ( 6,*) ' CASCADE # ',nkorri,' TO # ',nkorrx,
     +               ' TO # ',nkorrf 
        write ( 6,*) ' '  
        endif  
	WRITE ( 6,13 ) MZAHL
13	FORMAT ( ' NUMBER OF PARTICLE DETECTOR WINDOWS : ',I2/ )
	WRITE ( 6,14 ) NORMI,NORMF
        write(6,'(34h STEPS FOR INTEGRATION OVER PHI : ,i5)') nint 
14	FORMAT ( ' NORMALIZATION TRANSITION : STATE # ',I2,' TO STATE # ',I2/ )
	IF ( IO(12) .NE. 0 ) WRITE ( 6,15 ) Q(1),Q(2),Q(3)
15	FORMAT ( ' FINITE SOLID ANGLE CORRECTION : Q0= ',F6.4,' Q2= ',F6.4,
     +	' Q3= ',F6.4/ )
	WRITE ( 6,16 )
16	FORMAT ( ' PARTICLE DETECTOR WINDOWS :   #    THETA FROM TO',
     +	  '      PHI FROM TO ' )
	DO NI=1,MZAHL
	  WRITE ( 6,17 ) NI,WIN(NI,2),WIN(NI,1),WIN(NI,3),WIN(NI,4)
17	FORMAT ( 31X,I2,2X,F6.2,2X,F6.2,4X,F7.2,2X,F7.2 )
	  WIN(NI,3)=WIN(NI,3)*PI/180.
	  WIN(NI,4)=WIN(NI,4)*PI/180.
	END DO
C
C	INPUT FOR HAGER+SELTZER INTERPOLATION PROGRAM
C

c       NCCK,NCCL,NCCM have to be lower or equal 23,17,8 ... otherwise the
c       ALPHAx arrays are defined to small

	READ ( 1,* ) NCCK,NCCL,NCCM
	READ ( 1,* ) CC1K,CC2K,CCKMIN,CCKMAX
	READ ( 1,* ) CC1L,CC2L,CCLMIN,CCLMAX
	READ ( 1,* ) CC1M,CC2M,CCMMIN,CCMMAX


c       03.08.05 
c       if no conversion coefficients to be fitted are given in DCK.DAT
c       the subroutine FIT produces nonsense without setting the APLHAx 
c       arrays explicitely to 0.0 (obviously this happened automatically 
c       before)     

        do lll=1,5 
          do iii=1,ncck
            alphak(iii,lll)=0.0
          end do
          do iii=1,nccl
            alphal(iii,lll)=0.0
          end do 
          do iii=1,nccm
            alpham(iii,lll)=0.0
          end do
        end do 


100	READ ( 1,18 ) TYP,LAM
18	FORMAT ( A1,I1 )
c        write(6,*) ' typ  "',typ,'"' 
	IF ( TYP .EQ. ' ' ) GO TO 110
c        write(6,*) ' nicht nach 110 gesprungen' 
	IF ( TYP .EQ. 'M' ) THEN
	 IF ( LAM .NE. 1 ) THEN
	  WRITE ( 6,19 ) LAM
19	  FORMAT ( ' ****MISTAKE**** NO M',I1,' RADIATION ALLOWED')
	  STOP
	 END IF
	 LAM=5
	ELSE IF ( TYP .EQ. 'E' ) THEN
	 IF ( LAM .LT. 1 .OR. LAM .GT. 4 ) THEN
	  WRITE ( 6,20 ) LAM
20	  FORMAT ( ' ****MISTAKE**** NO E',I1,' RADIATION  ALLOWED')
	  STOP
	 END IF
	ELSE
	 WRITE ( 6,21 )
21	 FORMAT ( ' ****MISTAKE**** ILLEGAL RADIATION TYPE ENCOUNTERED')
	 STOP
	END IF
	FOUND=.FALSE.
	DO I1=1,5
	  IF ( LAMDA(I1) .EQ. LAM ) FOUND=.TRUE.
	END DO
c       for M1 lam was set to 5, but on output I want to see M1 not M5         
        lamtemp=lam
        if (lamtemp.eq.5) lamtemp=1 
	IF ( .NOT. FOUND ) WRITE ( 6,22 ) TYP,LAMtemp
22	FORMAT ( ' ****WARNING**** CONVERSION COEFFICIENTS FOUND FOR ',A1,
	1I1/', BUT NO CORRESPONDING MATRIX ELEMENTS FOUND IN CLX DATASET' )
	 IF ( ALPHAK(1,LAM) .NE. 0. ) THEN
	 WRITE ( 6,233 )
233	 FORMAT ( ' ****ERROR**** MULTIPLE DEFINITION OF CONVERSION', 
     +	  ' COEFFICIENTS' )
	 STOP
	END IF
	READ ( 1,* ) ( ALPHAK(I,LAM),I=1,NCCK )
	READ ( 1,* ) ( ALPHAL(I,LAM),I=1,NCCL )
	READ ( 1,* ) ( ALPHAM(I,LAM),I=1,NCCM )
	GO TO 100
110	CLOSE ( UNIT=1 )
C
C	------- COMPUTATION OF dW/DO STARTS HERE -------
C
C
C	THE TOTAL CONVERSION COEFFICIENTS ALPHA K,L,M FOR E1-E4 AND M1
C
	  DO 200 I1=1,NMPOL
c	  IF ( LAMDA(I1) .EQ. 0 ) GO TO 200
	    DO 200 IN=1,NMAX
	      DO 200 IM=1,NMAX
              if ( lamda(i1) .eq. 0 ) go to 200
	      IF ( ME(IN,IM,I1) .EQ. 0. ) GO TO 200
	      LAM=LAMDA(I1)
	      IF ( LAM .EQ. 5 ) LAM=1
	      IF ( ALPHAK(1,LAMDA(I1)) .NE. 0. ) THEN
C
C	UNITS IN H+S ROUTINES = keV
C
	       CCEE=(E(IN)-E(IM))*1000.
	       IF ( CCEE .LE. 0. ) GO TO 200
	       AKTOT=0.
	       IF ( CCEE .GE. CCKMIN .AND. CCEE .LE. CCKMAX ) THEN
	        DO IW=1,NCCK
	          ALPHA(IW)=ALPHAK(IW,LAMDA(I1))
c	          write(6,*) ' alphak ',iw,i1,lamda(i1),ALPHA(IW)
	        END DO
	        CALL FIT ( CC1K,CC2K,1,LAM,NCCK,ALPHA,CCEE )
c                write(6,*) ' aktot ',CC1K,CC2K,LAM,NCCK,CCEE  
c 	        DO IW=1,NCCK
c	          write(6,*) ' alpha ',ALPHA(IW)
c	        END DO                  
	        AKTOT=ALPHA(1)
	       END IF
	       ALTOT=0.
	       IF ( CCEE .GE. CCLMIN .AND. CCEE .LE. CCLMAX ) THEN
	        DO IW=1,NCCL
	          ALPHA(IW)=ALPHAL(IW,LAMDA(I1))
	        END DO
	        CALL FIT ( CC1L,CC2L,2,LAM,NCCL,ALPHA,CCEE )
	        ALTOT=ALPHA(1)
	       END IF
	       AMTOT=0.
	       IF ( CCEE .GE. CCMMIN .AND. CCEE .LE. CCMMAX ) THEN
	        DO IW=1,NCCM
	          ALPHA(IW)=ALPHAM(IW,LAMDA(I1))
	        END DO
	        CALL FIT ( CC1M,CC2M,3,LAM,NCCM,ALPHA,CCEE )
	        AMTOT=ALPHA(1)
	       END IF
c              write(6,*) ' alpha1 ',in,im,i1,AKTOT,ALTOT,AMTOT 
	       ALPHA1(IN,IM,I1)=AKTOT+ALTOT+AMTOT

c***
c                if (((in.eq.22).and.(im.eq.4)).or.
c     +             ((in.eq.24).and.(im.eq.6))) then
c                write(*,*) in,im
c                write(*,*) i1,lam,ME(IN,IM,I1),ALPHA1(IN,IM,I1)
c                write(*,*) aktot,altot,amtot
c                endif
c***
	      END IF
 200	CONTINUE
	IF ( IO(1) .NE. 0 ) THEN
	DO 300 I1=1,NMPOL
c	  IF ( LAMDA(I1) .EQ. 0 ) GO TO 300
          if ( lamda(i1) .ne. 0 ) then 
	    IF ( LAMDA(I1) .EQ. 5 ) THEN
	      TYP='M'
	      LAM=1
 	    ELSE
	      TYP='E'
	      LAM=LAMDA(I1)
	    END IF
           endif
	   LINE=50
	   DO 300 IN=1,NMAX
	     DO 300 IM=1,NMAX
               if ( lamda(i1) .eq. 0 ) go to 300 
	       IF ( ALPHA1(IN,IM,I1) .EQ. 0. ) GO TO 300
	       IF ( LINE .GT. 45 ) THEN
	        LINE=3
		CALL NEWPAGE
	        WRITE ( 6,23 ) TYP,LAM
23	FORMAT ('     CONVERSION COEFFICIENTS OF K,L, AND M SHELL FOR '
     +	,A1,I1)
		WRITE ( 6,24 )
24	FORMAT(/5X,'STATE # INITIAL SPIN  STATE # FINAL SPIN   ENERGY ',  
     +	  ' ALPHA')
	       END IF
	       LINE=LINE+1
	       CCEE=E(IN)-E(IM)
	       WRITE (6,25) IN,SPIN(IN),IM,SPIN(IM),CCEE,
     +	              ALPHA1(IN,IM,I1)
25	       FORMAT ( 7X,I2,7X,F4.1,9X,I2,6X,F4.1,7X,F6.4,2X,G11.4 )
300	 CONTINUE
	END IF
C
C	COMPUTE THE H(N)
C
C	IN ORDER TO OBTAIN THE DELTAS IN UNITS OF SQRT(1/Sec)
C	APPR0PRIATE FACTORS TPFAC(LAMDA) ARE CHOSEN.
C	SEE: ALDER,WINTHER ELECTROMAGNETIC EXCITATION p. 273
C
	TPFAC(1)=1.5902E17
	TPFAC(2)=1.2251E13
	TPFAC(3)=5.7083E8
	TPFAC(4)=1.6965E4
	TPFAC(5)=1.7584E13
	DO 500 IN=1,NMAX
	  H(IN)=0.
	  ASUM=0.0
	  DO 400 I1=1,NMPOL
c	    IF ( LAMDA(I1) .EQ. 0. ) GO TO 400
            if ( lamda(i1) .ne. 0. ) then
	      LAM=LAMDA(I1)
	      L1=LAM
	      IF ( L1 .EQ. 5 ) L1=1
            end if
	    DO 400 IM=1,NMAX
c              write(6,*) ' 400 ',i1,lamda(i1),in,im,ME(IN,IM,I1)  
              if ( lamda(i1) .eq. 0. ) go to 400
	      IF ( ME(IN,IM,I1) .EQ. 0. ) GO TO 400
	      CCEE=(E(IN)-E(IM))
	      IF ( CCEE .LE. 0. ) GO TO 400
	      DELTA(IN,IM,I1)=SQRT(TPFAC(LAM)*CCEE**(2*L1+1)/(2*SPIN(IN)+1))
	1*ME(IN,IM,I1)
c              write(6,*) ' epsq ',in,im,i1,ALPHA1(IN,IM,I1)          
	      EPSQ(IN,IM,I1)=DELTA(IN,IM,I1)**2
	1*ALPHA1(IN,IM,I1)
c              write(6,*) ' asum ',asum,in,im,i1,
c     +                   DELTA(IN,IM,I1),EPSQ(IN,IM,I1)
	      ASUM=ASUM+DELTA(IN,IM,I1)**2+EPSQ(IN,IM,I1)
c***
c                if (((in.eq.22).and.(im.eq.4)).or.
c     +             ((in.eq.24).and.(im.eq.6))) then
c                write(*,*) in,im
c                write(*,*) i1,l1,ME(IN,IM,I1),ALPHA1(IN,IM,I1)
c                write(*,*) EPSQ(iN,im,I1),delta(iN,im,I1)
c                endif
c***
 400	  CONTINUE
 500	IF ( ASUM .NE. 0. ) H(IN)=1./ASUM
	IF ( IO(3) .NE. 0 ) THEN
	 LINE=50
	 DO NI=1,NMAX
	   IF ( LINE .GT. 45 ) THEN
	    CALL NEWPAGE
	    WRITE(6,26)
26	FORMAT ( '     LIFETIMES 0F THE STATES'/
	1/'     STATE #     LIFETIME [Sec]' )
	    LINE=3
	   END IF
	   LINE=LINE+1
	   IF ( H(NI) .EQ. 0. ) WRITE(6,27) NI
27	FORMAT ( 8X,I2,8X,'DOES NOT DECAY' )
	   IF ( H(NI) .NE. 0. ) WRITE(6,28) NI,H(NI)
28	FORMAT ( 8X,I2,8X,G11.4 )
	 END DO
	END IF
C
C	COMPUTE THE GK'S
C
	DO I1=1,NMPOL
	  LAM=LAMDA(I1)
	  IF ( LAM .EQ. 5 ) LAM=1
	  DO NI=1,NMAX
	    DO NF=1,NMAX
c***            Why EPSQ (even if ALPHA1=0 the gamma-transition exists!)?
C***            change: DELTA -> EPSQ
c***	      IF ( EPSQ(NI,NF,I1) .NE. 0. ) THEN
              IF ( delta(NI,NF,I1) .NE. 0. ) THEN
	       FAC=H(NI)*SQRT((SPIN(NI)*2.+1.)*(SPIN(NF)*2.+1.))
	1*(EPSQ(NI,NF,I1)+DELTA(NI,NF,I1)**2.)
	       DO K=1,3
		 IX=SPIN(NI)+SPIN(NF)+LAM+(2*K-2)
		 IF ( MOD(IX,2) .EQ. 0 ) THEN
		  PHZ=1
		 ELSE
		  PHZ=-1
		 END IF
		 GK(NI,NF,K)=GK(NI,NF,K)+PHZ*FAC
	1*RACAH(SPIN(NF),SPIN(NF),SPIN(NI),SPIN(NI),REAL(2*K-2),REAL(LAM))
	       END DO
	      END IF
	    END DO
	  END DO
	END DO
C
C	COMPUTE THE FF'S and the ffkkk's
C
	DO I1=1,NMPOL
	  AL1=LAMDA(I1)
	  LL1=LAMDA(I1)
	  IF ( LL1 .EQ. 5 ) THEN
	   AL1=1.
	   LL1=2
	  END IF
	  DO I2=1,NMPOL
	    AL2=LAMDA(I2)
	    LL2=LAMDA(I2)
	    IF ( LL2 .EQ. 5 ) THEN
	     AL2=1.
	     LL2=2
	    END IF
	    DO NI=1,NMAX
	      DO NF=1,NMAX
c***
c                if (((ni.eq.22).and.(nf.eq.4)).or.
c     +             ((ni.eq.24).and.(nf.eq.6))) then
c                write(*,*) ni,nf,h(ni)
c                write(*,*) i1,EPSQ(NI,NF,I1),delta(NI,NF,I1)
c                write(*,*) i2,EPSQ(NI,NF,I2),delta(NI,NF,I2)
c                write(*,*) al1,al2
c                do k2=1,3
c                write(*,*) FK(SPIN(NF),AL1,AL2,SPIN(NI),REAL(2*K2-2))
c                end do
c                endif
c***
c***            Why EPSQ which is never used in this part?
c***		IF (EPSQ(NI,NF,I1).NE.0..AND.EPSQ(NI,NF,I2).NE.0.)THEN
                IF (delta(NI,NF,I1).NE.0..AND.delta(NI,NF,I2).NE.0.)THEN
		 IF ( IPAR(NI) .EQ. IPAR(NF) ) THEN
		  IF ( MOD(LL1+LL2,4) .EQ. 0 ) THEN
		   PHZ=+1.
		  ELSE IF ( MOD(LL1+LL2,2) .EQ. 0 ) THEN
		   PHZ=-1.
		  ELSE
		   WRITE ( 6,29 )
29	FORMAT ( ' ***MISTAKE*** DELTA(L`)*DELTA(L``) IMAGINARY' )
		   STOP
		  END IF
		 ELSE
		  IF ( MOD(LL1+LL2,4) .EQ. 0 ) THEN
		   PHZ=-1.
		  ELSE IF ( MOD(LL1+LL2,2) .EQ. 0 ) THEN
		   PHZ=+1.
		  ELSE
		   WRITE ( 6,29 )
		   STOP
		  END IF
		 END IF
		 DO K2=1,3
		   FF(NI,NF,K2)=FF(NI,NF,K2)
     +               +H(NI)*PHZ*DELTA(NI,NF,I1)*DELTA(NI,NF,I2)
     +                *FK(SPIN(NF),AL1,AL2,SPIN(NI),REAL(2*K2-2))
                   if ((ni.eq.nkorri).and.(nf.eq.nkorrx)) then
                     rif=SPIN(NF)
                     rii=SPIN(NI)
                     rk2=REAL(2*K2-2)
	             do k=1,3
                       rkk=REAL(2*K-2)
                       do k1=1,3
                         rk1=REAL(2*K1-2)
                         ffkkk(k1,k2,k)=ffkkk(k1,k2,k)
     +                   +H(NI)*PHZ*DELTA(NI,NF,I1)*DELTA(NI,NF,I2)
     +                    *FKk2k1(rif,AL1,AL2,rii,rk1,rk2,rkk)
                       end do
                     end do
                   end if                 
		 END DO
		END IF
	       END DO
	    END DO
	  END DO
	END DO
C
C	DESORIENTATION CORRECTION
C
        do ijk=1,nmax
          qg(ijk,1)=1.0
          qg(ijk,2)=1.0
          qg(ijk,3)=1.0
          ideo(ijk)=0
        end do   
	IF ((IO(12).EQ.2).or.(IO(12).EQ.4)) THEN
	 OPEN ( UNIT=3,FILE='DESORIN.dat',STATUS='OLD' )
	 CALL NEWPAGE
	 WRITE ( 6,30 )
30	FORMAT ( '     THE DESORIENTATION COEFFICIENTS '//
	1,4X,' N',8X,' G0',8X,'G2',8X,'G4'/ )
        do ijk=1,nmax
 600	  READ ( 3,*,ERR=700 ) NI,QG(ni,1),QG(ni,2),QG(ni,3)
	  WRITE ( 6,31 ) NI,QG(ni,1),QG(ni,2),QG(ni,3)
31     	FORMAT ( 4X,I2,7X,F6.4,4X,F6.4,4X,F6.4 )
         ideo(ni)=1
        end do

c	 DO NF=1,NMAX
c	   IF ( FF(NI,NF,1) .NE. 0. ) THEN
c	    DO K=1,3
c	      FF(NI,NF,K)=FF(NI,NF,K)*QG(K)
c              FF(NI,NF,K)=FF(NI,NF,K)*QG(K)
c	    END DO
c	   END IF
c	 END DO

c  hier zu 600 zurueckzuspringen macht wenig Sinn  .... 02.08.05
c 	 GO TO 600

 700	 CLOSE ( UNIT=3 )
         if (io(12).eq.2) then
           write(6,1234) 
           write(6,1235)
          end if 
	END IF

        IF ( IO(12) .eq. 3 ) THEN
          write(6,*) ' '
          write(6,1236)
          write(6,*) '   LAMDA2,N_V  : ',lamda2,n_v
          write(6,*) ' '
        end if 
        IF ( IO(12) .eq. 4 ) THEN
          write(6,1234) 
          WRITE(6,1236)
          write(6,*) '   LAMDA2,N_V  : ',lamda2,n_v
          write(6,*) ' '
        end if

1234    format(/' ALL OTHER STATES : ')
1235    format(' THE DESORIENTATION COEFFICIENTS ARE 1.0')
1236    format(' THE DESORIENTATION COEFFICIENTS ',
     +   'ARE CALCULATED WITH THE PARAMETRISATION OF A&P')
             
c	DO NI=1,NMAX
c	  DO NF=1,NMAX
c	    IF ( FF(NI,NF,1) .NE. 0. ) THEN
c	     DO K=2,3
c	       FF(NI,NF,K)=FF(NI,NF,K)*Q(K)
c	     END DO
c	    END IF
c	  END DO
c	END DO
C
C	PRINTOUT OF GK'S, FF'S AND FKK2K1'S
C
	IF ( IO(4) .NE. 0 ) THEN
	 LINE=50
	 DO NI=1,NMAX
	   DO NF=1,NMAX
	     IF ( GK(NI,NF,1) .NE. 0 ) THEN
	      IF ( LINE .GT. 45 ) THEN
	       CALL NEWPAGE
	       LINE=3
	       WRITE ( 6,32 )
32	FORMAT ( '     THE G MATRIX'//5X,'STATE # INITIAL SPIN  ',
     +    'STATE # ',
     +    'FINAL SPIN   ENERGY',10X,'G0',11X,'G2',11X,'G4'/ )
	      END IF
	      LINE=LINE+1
	      CCEE=E(NI)-E(NF)
	      WRITE ( 6,33 ) NI,SPIN(NI),NF,SPIN(NF),CCEE,GK(NI,NF,1),
     +              GK(NI,NF,2),GK(NI,NF,3)
33	FORMAT ( 7X,I2,7X,F4.1,8X,I2,7X,F4.1,7X,F6.4,4X,3(2X,G11.4) )
	     END IF
	   END DO
	 END DO
	END IF

	IF ( (IO(5) .Eq. 1).or.(IO(5) .Eq. 3) ) THEN
	 LINE=50
	 DO NI=1,NMAX
	   DO NF=1,NMAX
	     IF ( FF(NI,NF,1) .NE. 0 ) THEN
	      IF ( LINE .GT. 45 ) THEN
	       CALL NEWPAGE
	       LINE=3
               WRITE ( 6,34 )
34	FORMAT ( '     THE F MATRIX'//5X,'STATE # INITIAL SPIN',
     +    '  STATE #',
     +    'FINAL SPIN   ENERGY',10X,'F0',12X,'F2',11X,'F4'/ )
	      END IF
	      LINE=LINE+1
	      CCEE=E(NI)-E(NF)
	      WRITE ( 6,35 ) NI,SPIN(NI),NF,SPIN(NF),CCEE,FF(NI,NF,1),
     +              FF(NI,NF,2),FF(NI,NF,3)
35	FORMAT ( 7X,I2,7X,F4.1,8X,I2,7X,F4.1,7X,F6.4,4X,3(2X,G11.4) )
	     END IF
	   END DO
	 END DO
	END IF

        IF ( (IO(5) .Eq. 2).or.(IO(5) .Eq. 3) ) THEN
         call newpage
         WRITE ( 6,134 )
 134     FORMAT ( '     THE FKK2K1 MATRIX'//5X,'STATE # INITIAL SPIN',
     +    '  STATE # FINAL SPIN   ENERGY'/ )
         CCEE=E(Nkorri)-E(Nkorrx)
         ni=nkorri
         nf=nkorrx
         WRITE ( 6,135 ) Ni,SPIN(Ni),Nf,SPIN(NF),CCEE
 135     FORMAT ( 7X,I2,7X,F4.1,8X,I2,7X,F4.1,7X,F6.4//)

         write(6,*) ' K2  K1          FFKKK(K=0,2,4 ,K2,K1)'
         DO k2=1,3
           DO k1=1,3
             write(6,136) 2*k2-2,2*k1-2,ffkkk(1,k2,k1),
     +             ffkkk(2,k2,k1),ffkkk(3,k2,k1)
 136         format(2i4,3e15.6)
           END DO
         END DO

        END IF


C
C	COMPUTE THE NEW AK'S
C
	DO N=1,NMAX
	  NN(N)=N
	  ESORT(N)=E(N)
	END DO
	DO N=1,NMAX
	  DO M=N+1,NMAX
	    IF ( ESORT(N) .GT. ESORT(M) ) THEN
	     NNTEMP=NN(N)
	     NN(N)=NN(M)
	     NN(M)=NNTEMP
	     ETEMP=ESORT(N)
	     ESORT(N)=ESORT(M)
	     ESORT(M)=ETEMP
	    END IF
	  END DO
	END DO
	if (dtheta.eq.0.0)then
	    nth=1
	  eLSe
	    NTH=ABS((THETA1-THETA2)/DTHETA)+1
	ENDif

cc      provisorische Notloesung!!!! 10.07.01         
        if (nth.gt.nthclxmax) then
          write(6,*) ' *** TOO MANY ANGLES IN CLXDCY ***'
          write(6,*) ' nth : ',nth,' -> ',nthclxmax
          write(6,*) '      '
          nth=nthclxmax
        endif  

	VF1=((1.+A1/A2)*Z1*Z2/EP)**2
	VCM=0.04634*SQRT(EP/A1)/(1.+A2/A1)


c        Calculate the probability that a state if populated will decay
c        to a lower lying state

        do i=1,nmax
          do j=1,nmax
            dcyprob(i,j)=0.0
          enddo
        enddo
 
        do n=0,nmax-2
          do j=1,nmax-n-1
            ni=nn(nmax-n)
            nf=nn(nmax-n-j)
            prob=0.0
            DO I1=1,NMPOL
              LAM=LAMDA(I1)
              IF ( LAM .EQ. 5 ) LAM=1
c***          IF ( EPSQ(ni,nf,I1) .NE. 0. )
              IF ( delta(ni,nf,I1) .NE. 0. )
     +          prob=prob+H(NI)*(EPSQ(NI,NF,I1)+DELTA(NI,NF,I1)**2.)
            END DO
            dcyprob(ni,nf)=prob
            if (n.gt.0) then
              do jj=0,n-1
                 nii=nn(nmax-jj)
                 dcyprob(nii,nf)=dcyprob(nii,nf)
     +             +dcyprob(nii,ni)*dcyprob(ni,nf)               
              enddo
            endif
          enddo
        enddo



C
c
c****   ========  Loop over aLL Gamma-detectorS  ==========
c
        DO N_D = 1,N_DET
            TT = THE_GAM(in_d(n_d))
            VG = PHI_GAM(in_d(N_D))
            T=tt
            COST=COS(T)
            SINT=SIN(T)
            VGAMMA=PHI_GAM(in_d(N_D))

                                      
           
C
C****	------- BIG LOOP OVER ALL SCATTERING ANGLES -------
C
	DO ITH=1,NTH
	  TH=(THETA1+(ITH-1)*DTHETA)
	  THETA=TH/180.*PI

c       If IO(12)>=3 the desorientation coefficients are calculated with
c       the Abragam & Pound parametrisation
c       lamda2 and lamda4 are given in [1/psec]
c       qg(n,1)=1.0
c       qg(n,2)=1./(1.+h(n)*lamda2*(1.-cos(theta))**n_v)
c       qg(n,3)=1./(1.+h(n)*lamda4*(1.-cos(theta))**n_v)
c       lamda4=lamda2*20./6.
c

          if (io(12).ge.3) then
            lamda4=lamda2*20./6.
            do l=1,nmax
              if (ideo(l).eq.0) then
                fac=h(l)*(1.-cos(theta))**n_v
                qg(l,2)=1./(1.+lamda2/1e-12*fac)
                qg(l,3)=1./(1.+lamda4/1e-12*fac)
              endif
            end do
          end if
                                            
        

        if (n_d .eq. 1) then
C
C	THE AK'S FULFILL THE FOLLOWING RECURSION RELATION:
C	BE NN(NMAX) THE INDEX OF THE STATE WITH THE HIGHEST EXCITATION
C	ENERGY, NN(NMAX-I) THE INDEX OF A STATE BELOW NN(NMAX) THEN:
C	AK(NN(NMAX)=RHOC(NN(NMAX))
C	AK(NN(NMAX-N))=RHOC(NN(NMAX-N))+
C	  SUM(GK(NN(NMAX-J),NN(NMAX-N))*AK(NN(NMAX-J)))
C	WHERE J RUNS FROM 0 TO N-1
C
        READ ( 2,* ) (((RHOC(N,K,KA),KA=1,2*K-1),K=1,3),N=2,NMAX)

          do k=1,nmax
            nnpos(nn(k))=k
            do kk=1,nmax
              nnpop(k,kk)=0
            end do
          end do
          nnkorri=nnpos(nkorri)
          nnkorrx=nnpos(nkorrx)      
          nnkorrf=nnpos(nkorrf) 
          nnpop(nmax,nmax)=2
       
	  DO K=1,3
	    DO KAPPA=1,2*K-1
c          write(6,*) ' rhoc1 ',NN(NMAX),K,KAPPA,RHOC(NN(NMAX),K,KAPPA)
	      AK(ith,NN(NMAX),K,KAPPA)=RHOC(NN(NMAX),K,KAPPA)
c*1*
              if (io(11).eq.1) then
                do n=1,nmax-1
                  RHOC(NN(N),K,KAPPA)=0.0
                end do
              endif
c*1*
	    END DO
	  END DO
	  DO N=1,NMAX-1
	    DO K=1,3
	      DO KAPPA=1,2*K-1
	        ASUM=0.
	        DO J=0,N-1

		  IF ( GK(NN(NMAX-J),NN(NMAX-N),1) .NE. 0. ) THEN
c***              if (j.eq.j) then

c                   Excitation energy of populated state higher than 
c                   or equal to nkorri

                    if (((nmax-n).ge.nnkorri).or.(io(13).eq.1)) then
                      asumalt=asum  
		      ASUM=ASUM+AK(ith,NN(NMAX-J),K,KAPPA)*
     +                  GK(NN(NMAX-J),NN(NMAX-N),K)*qg(nn(nmax-j),k)
c***
c                      if ((NN(NMAX-N).eq.6).and.
c     +                    (NN(NMAX-J).eq.24)) then
c                        write(*,*) ' *',NN(NMAX-J),NN(NMAX-N),k,kappa
c                        write(*,*) asumalt
c     +                             ,AK(ith,NN(NMAX-J),K,KAPPA)
c     +                             ,GK(NN(NMAX-J),NN(NMAX-N),K)
c     +                             ,asum 
c                        write(*,*) qg(nn(nmax-j),k)
c                        write(*,*) asum 
c                      endif
c*** 
                      nnpop(nn(nmax-n),nn(nmax-j))=1  
c                     If jj populates j, then jj populates also n !
                      if (j.gt.0) then
                        do jj=0,j-1
                          if (nnpop(nn(nmax-n),nn(nmax-jj)).ne.1)
     +                      nnpop(nn(nmax-n),nn(nmax-jj))=
     +                      nnpop(nn(nmax-j),nn(nmax-jj))
                        end do
                      endif 
                    endif

c                   nkorri -> nkorrx
c
c                    if ((nn(nmax-j).eq.nkorri).and.
c     +                  (nn(nmax-n).eq.nkorrx)) then
c                      ASUM=ASUM+AK(ith,nkorri,K,KAPPA)*
c     +                  GK(nkorri,nkorrf,K)*qg(nkorri,k)
c                      nnpop(nkorrf,nkorri)=1
c                     If jj populates j, then jj populates also n !
c                      if (j.gt.0) then
c                        do jj=1,n-1
c                          if (nnpop(nn(nmax-n),nn(nmax-jj)).ne.1)
c     +                      nnpop(nkorrf,nn(nmax-jj))=
c     +                      nnpop(nkorri,nn(nmax-jj))
c                        end do
c                      endif
c                    endif

c                   Excitation energy of populated state lower than nkorrx
c                   Only populated by states with excitation energies lower
c                   than or equal nkorrx
     
c                    if (((nmax-n).lt.nnkorrx).and.
c     +                  ((nmax-j).le.nnkorrx)) then
c                      ASUM=ASUM+AK(ith,NN(NMAX-J),K,KAPPA)*
c     +                  GK(NN(NMAX-J),NN(NMAX-N),K)*qg(nn(nmax-j),k)
c                      nnpop(nn(nmax-n),nn(nmax-j))=1
c                     If jj populates j, then jj populates also n !
c                      if (j.gt.0) then
c                        do jj=0,j-1
c                          if (nnpop(nn(nmax-n),nn(nmax-jj)).ne.1)
c     +                      nnpop(nn(nmax-n),nn(nmax-jj))=
c     +                      nnpop(nn(nmax-j),nn(nmax-jj))
c                        end do
c                      endif
c                    endif


	 	  END IF

c***
c                  dsum=delta(NN(NMAX-J),NN(NMAX-N),1)
c                  dsum=dsum+delta(NN(NMAX-J),NN(NMAX-N),2)
c                  dsum=dsum+delta(NN(NMAX-J),NN(NMAX-N),3)
c                  dsum=dsum+delta(NN(NMAX-J),NN(NMAX-N),4)
c                  dsum=dsum+delta(NN(NMAX-J),NN(NMAX-N),5)
c
c                  IF ((GK(NN(NMAX-J),NN(NMAX-N),k).Eq.0.).and.
c     +                (dsum.ne.0.).and.(NN(NMAX-J).eq.24)) THEN
c                    write(*,*) NN(NMAX-J),NN(NMAX-N),k
c                    write(*,*) dsum
c                    write(*,*) RHOC(NN(NMAX-N),K,KAPPA)
c                    write(*,*) ith,AK(ith,NN(NMAX-N),K,KAPPA) 
c                  endif

c***
c                if ((NN(NMAX-N).eq.6).and.(asum.ne.0.)) then
c                  write(*,*) ' **',NN(NMAX-j),NN(NMAX-N),k,kappa
c                  write(*,*) RHOC(NN(NMAX-N),K,KAPPA),asum
c     +                      ,ith,AK(ith,NN(NMAX-N),K,KAPPA)
c                endif
c***

		END DO

c               Population= decay + direct population
 
c                write(6,*) ' asum ak ',asum 

                if (((nmax-n).ge.nnkorri).or.(io(13).eq.1)) then
c                  write(6,*) ' rhoc2 ',N,NN(NMAX-N),K,KAPPA,
c     +                       RHOC(NN(NMAX-N),K,KAPPA)
		  AK(ith,NN(NMAX-N),K,KAPPA)=RHOC(NN(NMAX-N),K,KAPPA)
     +              +ASUM
                  nnpop(nn(nmax-n),nn(nmax-n))=2
                else
                  AK(ith,NN(NMAX-N),K,KAPPA)=asum
                endif

c***
c                if (NN(NMAX-N).eq.6) then
c                  write(*,*) ' ***',NN(NMAX-N),k,kappa 
c                  write(*,*) RHOC(NN(NMAX-N),K,KAPPA),asum
c     +                       ,ith,AK(ith,NN(NMAX-N),K,KAPPA)
c                endif
c***

c               Further decay of NN(NMAX-N) down to NKORRI
c
c                if ((nmax-n).gt.nnkorri) then
c                  do n1=n+1,nmax-1
c                    if (((nmax-n1).ge.nnkorri).and.
c     +                  (GK(NN(NMAX-n),NN(NMAX-N1),1).ne.0.0))
c     +                  then
c                      AK(ith,NN(NMAX-n1),K,KAPPA)= 
c     +                  AK(ith,NN(NMAX-n),K,KAPPA)*
c     +                  GK(NN(NMAX-n),NN(NMAX-N1),K)*
c     +                  qg(nn(nmax-n),k)
c                      if ((nmax-n1).eq.nnkorri)
c     +                    akkorri(nn(nmax-n),nn(nmax-n1),
c     +                      k,kappa)=AK(ith,NN(NMAX-n1),K,KAPPA)
c                      do n2=n1+1,nmax-1
c                        asum2=0.0  
c                        do j2=n1,n2-1
c                          if (((nmax-n2).ge.nnkorri).and.
c     +                       (GK(NN(NMAX-j2),NN(NMAX-N2),1).ne.
c     +                       0.0)) then
c                            ASUM2=ASUM2+AK(ith,NN(NMAX-j2),K,KAPPA)*
c     +                        GK(NN(NMAX-j2),NN(NMAX-N2),K)*
c     +                        qg(nn(nmax-j2),k)
c                          endif
c                        end do
c                        AK(ith,NN(NMAX-n2),K,KAPPA)=asum2
c                        if ((nmax-n2).eq.nnkorri)
c     +                    akkorri(nn(nmax-n),nn(nmax-n1),
c     +                      k,kappa)=asum2
c                      end do
c                    endif
c                  end do
c                endif
c
c               akkorri(nn(nmax-n),nn(nmax-n1),k,kappa) is AK of
c               nkorri if nkorri is populated ONLY by the decay of 
c               n to n1
  
                                 
	      END DO
	    END DO

	  END DO

c         Eliminate all transitions which do NOT appear in coincidence
c         with gate
c
c          do n=0,nmax-1
c            if ((nmax-n).gt.nnkorri) then
c              if (nnpop(nkorri,nn(nmax-n)).eq.0) then
c                do j=1,nmax
c                  ff(nn(nmax-n),j,1)=0.0
c                  ff(j,nn(nmax-n),1)=0.0
c                  nnpop(nn(nmax-n),j)=0
c                  nnpop(j,nn(nmax-n))=0
c                end do
c              endif
c            endif
c            if (nn(nmax-n).eq.nkorri) then
c              do j=1,nmax
c                if (j.ne.nkorrx)
c     +            ff(nn(nmax-n),j,1)=0.0
c              end do
c            endif
c            if (((nmax-n).lt.nnkorri).and.
c     +          ((nmax-n).gt.nnkorrx)) then
c               do j=1,nmax
c                 ff(nn(nmax-n),j,1)=0.0
c                 ff(j,nn(nmax-n),1)=0.0
c               end do
c            endif
c            if (nn(nmax-n).eq.nkorrx) then
c              do j=1,nmax
c                if (j.ne.nkorri)
c     +            ff(j,nn(nmax-n),1)=0.0
c              end do
c            endif 
c            if ((nmax-n).lt.nnkorrx) then
c              if (nnpop(nn(nmax-n),nkorrx).eq.0) then
c                do j=1,nmax
c                  ff(nn(nmax-n),j,1)=0.0
c                  ff(j,nn(nmax-n),1)=0.0
c                  nnpop(nn(nmax-n),j)=0
c                  nnpop(j,nn(nmax-n))=0
c                end do
c              endif
c              do j=1,nmax
c                if ((nmax-j).ge.nnkorri)  
c     +            ff(nn(nmax-j),nn(nmax-n),1)=0.0
c              end do
c            endif
c          end do


	  IF ((IO(6) .NE. 0).and.(n_d .eq. 1 )) THEN
	   LINE=50
	     DO NI=2,NMAX
	     IF ( AK(ith,NI,1,1) .NE. 0. ) THEN
	      IF ( LINE .GT. 45 ) THEN
	       CALL NEWPAGE
	       LINE=3
	       WRITE ( 6,36 ) TH
36	FORMAT ( '     THE MATRIX A(N,K,KAPPA) FOR THETA = ',F6.2,' DEGREES'//
	14X,' N',3X,' K',3X,' A(N,K,0)',4X,' A(N,K,1)',4X,' A(N,K,2)',4X,' A(N,K,3)',
	24X,' A(N,K,4)'/ )
	      END IF
	      KAMAX=MIN(INT(2.02*SPIN(NI)),4)
	      DO K=0,KAMAX,2
	        KK=K/2+1
	        WRITE ( 6,37 ) NI,K,(AK(ith,NI,KK,KAPPA),KAPPA=1,K+1)
37	FORMAT ( 4X,I2,3X,I2,5(2X,G11.4) )
	        LINE=LINE+1
	      END DO
	     END IF
	   END DO
C***
c       stop
c***
	  END IF

       end if

c      End of calculating the AKs


C
C	CALCULATION OF THE GAMMA-RAY ANGULAR DISTRIBUTION W(IN->IM)
C
C	INCLUDES THE INTEGRATION OVER THE PARTICLE COUNTER IN AZIMUTHAL
C	DIRECTION AND THE TRANSFORMATION OF dW/DO(GAMMA) INTO THE LAB SYSTEM
C
C	SYMBOL IN D.SCHWALM ET AL.              VARIABLE IN
C	NP A192 (1972) 449                     T R I S T A N
C	SMALL K					KA
C	KAPPA					DSKAPPA
C	V(J)					VP2
C	COS(THETA(J))				T2COS
C	SIN(THETA(J))				T2SIN
C	COS(THETA(GAMMA,J))			TG2COS
C
C	dCOS(THREST(GAMMA)*dPHIREST(GAMMA)
C	----------------------------------	RESTLAB
C	 dCOS(THLAB(GAMMA)*dPHILAB(GAMMA)
C
C	COS(THETA(GAMMA))			COST
C	COS(THREST(GAMMA))			TGCOS
C	SIN(THREST(GAMMA))			TGSIN
C	PHIREST(GAMMA? PARTICLE1)	        X(IX)
C	PHILAB(GAMMA)				VGAMMA
C	COS(PHIREST(GAMMA)-PHIREST(1))		VICOS(IX)
C
C
C	VALUES NEEDED FOR THE INTEGRATION OVER THE AZIMUTHAL
C	APERTURE OF THE PARTICLE DETECTOR ( FUNCTION Y )
C

	  DO I=1,MZAHL
	     V1=WIN(I,3)
	     V2=WIN(I,4)

             if (io(13).eq.0) call filly(nkorri,ith,i,0,1)

	     DO NI=1,NMAX

               if (io(13).eq.1) call filly(ni,ith,i,1,n_d)

               if (io(14).ne.0) then 
c       not supported at the moment
                 write(6,*) '  **** (IO(14).NE.0) NOT SUPPORTED ****'
                 stop                        
	       ELSE
c       only this is supported
	       END IF

               line=50
              
	       DO NF=1,NMAX
      
                W(NI,NF,N_D,1)=0.
                W(NI,NF,N_D,2)=0.

c         Gamma-Gamma-Particle angular correlations                
                if (io(13) .eq. 0) then 
                 
                 if (((ni.eq.nkorri).and.(nf.eq.nkorrx)).or.
     +               ((ni.eq.nkorrx).and.(nf.eq.nkorrf))) then

                  do n_dk=1,n_detk
                    if ( in_d(n_d) .ne. in_dk(n_dk) ) then

                      dw(1)=0.0
                      dw(2)=0.0                      
                      DO K=1,3
                      kk=2*k-2
                      rkk=real(kk)
                      DO KAPPA=-kk,kk
                      rkap=real(kappa)
                      kap=abs(kappa)+1 
                      DO K1=1,3
                      kk1=2*k1-2 
                      rk1=real(kk1)
                      DO KAPPA1=-kk1,kk1
                      rkap1=real(kappa1)
                      kap1=abs(kappa1)+1
                      do k2=1,3
                      kk2=2*k2-2
                      rk2=real(kk2)
                      do kappa2=-kk2,kk2
                      rkap2=real(kappa2)
                      kap2=abs(kappa2)+1

                        ddw=ffkkk(k1,k2,k)*qg(nkorri,k1)*q(k1)
     +                     *FF(Nkorrx,NkorrF,K2)*qg(nkorrx,k2)*q(k2)
     +                     *sqrt(2.*rkk+1.)/(2.*rk2+1.)
     +                     *((-1.)**(int(rkk-rk1-rkap2)))
     +                     *cleb(rkk,rk1,rk2,rkap,rkap1,-rkap2)
     +                     *ak(ith,nkorri,k,kap)

c           Integration over PHI of particle counter
                        yyint(1)=0.
                        yyint(2)=0. 
                        if (ddw.ne.0.0) then 
                          if (ni.eq.nkorri) then
                call intyy(n_d,k1,kappa1,n_dk,k2,kappa2,ith,i,yyint)    
                          else
                call intyy(n_d,k2,kappa2,n_dk,k1,kappa1,ith,i,yyint)
                         endif
                        endif

c           Signs of AK,Y(K1,KAPPA1) and Y(K2,KAPPA2) 
                        vorz=1.
                        if (kappa .lt. 0) vorz=vorz*sig(kap)
                        if (kappa1 .lt. 0) vorz=vorz*sig(kap1)
                        if (kappa2 .lt. 0) vorz=vorz*sig(kap2)
                        ddw=vorz*ddw

c           Calculation of complex WQ
c           Remember: the RHOC and AK are real !! 
                         dw(1)=dw(1)+ddw*yyint(1)
                         dw(2)=dw(2)+ddw*yyint(2)                       

c                        dum=ddw*(abs(yyint(1))+abs(yyint(2)))
c                        if (dum .ne. 0.0) then
c                         if (line .gt. 45) then
c                           call newpage
c                           line=3
c                         endif 
c                        write(*,*) k,kappa,k1,kappa1,k2,kappa2
c                        write(*,*) ffkkk(k1,k2,k),qg(nkorri,k1),q(k1)
c                        write(*,*) FF(Nkorrx,NkorrF,K2),
c     +                             qg(nkorrx,k2),q(k2)
c                        write(*,*) cleb(rkk,rk1,rk2,rkap,rkap1,-rkap2)
c                        write(*,*) ak(ith,nkorri,k,kap)
c                        write(*,*) yyint(1),yyint(2)
c                        write(*,*) ddw,dw(1),dw(2)
c                         write(*,*) '   ' 
c                         line=line+9
c                        endif  
                      enddo
                      enddo
                      enddo
                      ENDDO
                      enddo
                      enddo

c           Should be real at this point 
                      w(ni,nf,n_d,1)=w(ni,nf,n_d,1)+dw(1)*vnorm*vnorm
                      w(ni,nf,n_d,2)=w(ni,nf,n_d,2)+dw(2)*vnorm*vnorm

                    endif
	          END DO

	         END IF
                endif

c        Gamma-Particle angular correlations 
                if (io(13) .eq. 1) then

                  dw(1)=0.0
                  dw(2)=0.0
                  DO K=1,3
                  kk=2*k-2
                  rkk=real(kk)
                  DO KAPPA=-kk,kk
                  rkap=real(kappa)
                  kap=abs(kappa)+1

c                  write(6,*) ' ff ',ni,nf,k,FF(Ni,NF,K) 
c                  write(6,*) ' qg q ',qg(ni,k),q(k)
c                  write(6,*) ' ak ',ith,kap,ak(ith,ni,k,kap)  

                    ddw=FF(Ni,NF,K)*qg(ni,k)*q(k)
     +                  *ak(ith,ni,k,kap)

c           Integration over PHI of particle counter
                    yyint(1)=0.
                    yyint(2)=0.
                    if (ddw.ne.0.0) 
     +                 call inty(n_d,k,kappa,ith,i,yyint)

c           Signs of AK and Y(K,KAPPA) should be allright, because
c           both have the same KAPPA:   SIG(KAP)**2 = 1.

c           Calculation of complex WQ
c           Remember: the RHOC and AK are real !!
c                    write(6,*) ' 1dw ',dw(1),ddw,yyint(1) 
c                    write(6,*) ' 2dw ',dw(2),ddw,yyint(2)  
                    dw(1)=dw(1)+ddw*yyint(1)
                    dw(2)=dw(2)+ddw*yyint(2)

c                  dum=ddw*(abs(yyint(1))+abs(yyint(2)))
c                  if (dum .ne. 0.0) then
c                    if ((dum .ne. 0.0).and.(ni.eq.7)
c     +                  .and.(nf.eq.6)) then
c                     if (dw(1).lt.0.0) then 
c                      if (line .gt. 45) then
c                        call newpage
c                        line=3
c                      endif 
c                      write(*,*) ni,nf,n_d,k,kappa
c                      write(*,*) FF(Ni,NF,K),qg(ni,k),q(k)
c                      write(*,*) ak(ith,ni,k,kap)  
c                      write(*,*) yyint(1),yyint(2)
c                      write(*,*) ddw,dw(1),dw(2)
c                      write(*,*) '    '
c                      line=line+6
c                   endif

                  enddo
                  enddo

c           Should be real at this point
c                  write(6,*) ' 1w ',w(ni,nf,n_d,1),dw(1)
c                  write(6,*) ' 2w ',w(ni,nf,n_d,2),dw(2)
                  w(ni,nf,n_d,1)=w(ni,nf,n_d,1)+dw(1)*vnorm
                  w(ni,nf,n_d,2)=w(ni,nf,n_d,2)+dw(2)*vnorm

                 endif
               
	       END DO

	     END DO

	     IF ( IO(7) .NE. 0 ) THEN
	     LINE=50
	     DO NI=1,NMAX
	       DO L=1,NMAX
	         IF ((W(NI,L,N_D,1) .NE. 0.).or.
     +               (W(NI,L,N_D,2) .NE. 0.)) THEN
	          IF ( LINE .GT. 45 ) THEN
	           CALL NEWPAGE
	           WRITE ( 6,38 ) TH,I
38	FORMAT ( '     GAMMA DIFFERENTIAL CROSS-SECTION FOR THETA = ',            
     +     F6.2,' DEGREES, WI#:',I2//5X,'STATE # INITIAL SPIN',
     +     '  STATE #', 
     +	   'FINAL SPIN   ENERGY   CROSS-SECT REAL   IMG'/ )
	           LINE=3
	          END IF
	          CCEE=E(NI)-E(L)
	          WRITE ( 6,39 ) NI,SPIN(NI),L,SPIN(L),CCEE,W(NI,L,N_D,1),
     +                           W(NI,L,N_D,2)
39	FORMAT ( 6X,I2,7X,F4.1,9X,I2,7X,F4.1,6X,F6.4,2X,G11.4,2x,G11.4)
	          LINE=LINE+1
	         END IF
	       END DO
	     END DO
	    END IF
C
C	-- CALCULATION OF THE DOUBLE OR TRIPLE DIFFERENTIAL CROSS SECTION --
C	------- SEE: WINTHER, DE BOER COULOMB EXCITATION p. 314 -------
c       -- TRIPLE DIFFERENTIAL CROSS SECTION (ALDER/WINTHER 1975, p.71) --   
C
	    DO NI=1,NMAX
C
C	CALCULATION OF THE RUTHERFORD CROSS SECTION (BARNS/STERADIAN)
C
	      VF2=SQRT(EP/(EP-(1.+A1/A2)*E(NI)))
	      DSR=1.2958E-3*VF1*VF2*(1./SIN(THETA/2.)**4)
	      ETT=EP-E(NI)
C
C	TRANSFORMATION OF THE RUTHERFORD CROSS-SECTION FROM THE CM
C	INTO THE LAB SYSTEM FOR TARGET-NUCLEUS
C	SEE MARION&YOUNG, TABLE 5
C
c       revised 11.07.01 TK

 
	      IF (A1.GT.A2) THEN

C	        M&Y-SYMBOLS CORRESPONDING TO TRISTAN VARIABLES :
C		M1 = A1
C		M2 = A2
C		M3 = A2 (LIGHT PRODUCT)
C		M4 = A1 (HEAVY PRODUCT)
C		E1 = EP
C		THETA = CM ANGLE OF PROJECTILE

 	        AAA=A1*A1*(EP/ETT)/((A1+A2)**2)
	        BBB=A1*A2*(EP/ETT)/((A1+A2)**2)
	        CCC=A2*A2*(1.+(A1*-E(NI)/(A2*ETT)))/((A1+A2)**2)
	        DDD=A2*A1*(1.+(A1*-E(NI)/(A2*ETT)))/((A1+A2)**2)

c               scattered projectile A1 is detected in particle detector
                tphi=theta    
                ttheta=pi-theta            
C               recoil A2 is detected in particle detector 
C               ttheta=theta
c               tphi=pi-theta 

C               light product = recoil
                E3ETT=BBB+DDD+2.*SQRT(AAA*CCC)*COS(ttheta) 
c               heavy product = scattered projectile
                E4ETT=AAA+CCC+2.*SQRT(AAA*CCC)*COS(tphi)

c               Lab-angle of the light product = recoil 
c               SIN(PSI)
                SINpsi=SQRT(DDD/E3ETT)*SIN(ttheta)

c               Lab-angle of the heavy product = scattered projectile
C	        SIN(ZETA)
                SINzeta=SQRT(A2*(1.0-e4ett)/(A1*E4ETT))*sinpsi

c               dSigma/dOmega is given as function of CM-theta of projectile,
c               but the integration over dTheta will be done over 
c               ThetaRecoilLab!!!!
c               sigma(projectile,cm) -> sigma(recoil,lab)
C		SIGLAB = SIGMA(PSI)/SIGMA(THETA) 
                thsin=SINpsi
        	SIGLAB=E3ETT/SQRT(AAA*CCC*((DDD/BBB)-SINpsi**2))	

	      ELSE 

C		M1 = A1
C		M2 = A2
C		M3 = A1 (LIGHT PRODUCT)
C		M4 = A2 (HEAVY PRODUCT)
C		E1 = EP
C		THETA = CM ANGLE OF PROJECTILE

	        AAA=A1*A2*(EP/ETT)/((A1+A2)**2)
	        BBB=A1*A1*(EP/ETT)/((A1+A2)**2)
	        CCC=A2*A1*(1.+(A1*-E(NI)/(A2*ETT)))/((A1+A2)**2)
	        DDD=A2*A2*(1.+(A1*-E(NI)/(A2*ETT)))/((A1+A2)**2)

c               A1 is detected in particle detector
                ttheta=theta  
                tphi=pi-theta              
C               A2 is detected in particle detector
C               ttheta=pi-theta
c               tphi=theta  

C               light product = projectile
                E3ETT=BBB+DDD+2.*SQRT(AAA*CCC)*COS(tTHETA) 
c               heavy product = recoil
                E4ETT=AAA+CCC+2.*SQRT(AAA*CCC)*COS(tphi)

c               Lab-angle of the light product = scattered projectile 
c               SIN(PSI)
                SINpsi=SQRT(DDD/E3ETT)*SIN(tTHETA)

c               Lab-angle of the heavy product = recoil
C	        SIN(ZETA)
                SINzeta=SQRT(A1*(1.0-e4ett)/(A2*E4ETT))*sinpsi

c               dSigma/dOmega is given as function of theta projectile CM
c               but the integration over dTheta will be done over 
c               ThetaRecoilLab!!!!
c               SIGLAB = SIGMA(ZETA)/SIGMA(PHI)
c               sigma(projectile,cm) -> sigma(recoil,lab)
                thsin=sinzeta
                SIGLAB=E4ETT/SQRT(AAA*CCC*((CCC/AAA)-sinzeta**2))	

	       END IF

c**            Why the kinematical correction was done partly with the
c**            individual energies E(NI) and partly with E(K1LAB)????
c**            Now it's done consistently with E(NI) ...
c**            ... which is "correct" if direct population by coulex is 
c**            larger than population by decay from above ....
C**            K1LAB.NE.0 -> old version;  K1LAB.EQ.0 -> new version     

c**            thsin = sin(ThetaLab) 
c**            dOmegaLab=sin(ThetaLab)*dTheta*dPhi
c**            integration over dPhi is already done!
c**            integration over dTheta will be done later ....
 
	       IF ( NI .EQ. K1LAB ) SINLAB=THSIN
	       DO L=1,NMAX
	         IF (W(NI,L,N_D,1).NE.0.) then
                   DS(NI,L,I,ITH)=W(NI,L,N_D,1)*DSR*SIGLAB
	           if (K1LAB.eq.0)
     +               DS(NI,L,I,ITH)=DS(NI,L,I,ITH)*THSIN
                 end if  
	       END DO
	     END DO

             if (K1LAB.ne.0) then
	       DO NI=1,NMAX
	         DO L=1,NMAX
	           IF (W(NI,L,N_D,1).NE.0.)
     +               DS(NI,L,I,ITH)=DS(NI,L,I,ITH)*SINLAB
	         END DO
	       END DO
             end if

C
C	NOTE THAT DS(NI,NF) HAS BEEN TRANSFERRED TO THE LAB SYSTEM !
C
	     IF ( IO(8) .NE. 0 ) THEN
	     LINE=50
	      DO NI=1,NMAX
	        DO L=1,NMAX
	          IF ( W(NI,L,N_D,1) .NE. 0. ) THEN
	           IF ( LINE .GT. 45 ) THEN
c	            CALL NEWPAGE
	            if(n_d.eq.n_det)WRITE ( 6,40 ) TH,I
	            LINE=3
	           END IF
	           CCEE=E(NI)-E(L)
		   
		   ds_out(ni,l,i,ith)=ds_out(ni,l,i,ith)+DS(NI,L,I,ITH)

	           if(n_d.eq.n_det)
     +		        WRITE ( 6,41 )NI,SPIN(NI),
     +			L,SPIN(L),CCEE,DS_out(NI,L,I,ITH)/float(n_det)
	           LINE=LINE+1
	          END IF
	        END DO
	      END DO
	     END IF
	    
	  END DO
	END DO

40      FORMAT (//'     GAMMA DIFFERENTIAL CROSS-SECTION',
     +   ' NORMALIZED TO THE'/
     +   '     RUTHERFORD CROSS-SECTION FOR THETA = ',
     +   F6.2,' DEGREES, WI#',I2//5X,'STATE # INITIAL SPIN  STATE #',
     +   'FINAL SPIN   ENERGY   CROSS-SECTION'/ )
41      FORMAT ( 7X,I2,7X,F4.1,9X,I2,7X,F4.1,7X,F6.4,4X,G11.4)

c****       end of loop over PPAC-windows
C
C	TRANSFORM THE SCATTERING ANGLES USED BY CLX AND THE
C	PARTICLE DETECTOR WINDOWS FROM THE CM TO THE LAB SYSTEM
c       Theta(Projectile,CM) -> Theta(Recoil,Lab)!!!!
C
	DO ITH=1,NTH
	  THLB(ITH)=(180.-((ITH-1)*DTHETA+THETA1))/2.*PI/180.
	END DO
	DO I=1,MZAHL
	  WINLB(I,1)=(180.-WIN(I,1))/2.*PI/180.
	  WINLB(I,2)=(180.-WIN(I,2))/2.*PI/180.
	END DO
	IF ( IO(10) .NE. 0 ) THEN
	  OPEN ( UNIT=3, FILE='ISUPP', STATUS='UNKNOWN' )
	  WRITE ( 3,42 ) NMAX,MZAHL
	END IF

C****   second loop over PPAC-windows

	DO I=1,MZAHL
	  IF ( IO(10) .NE. 0 ) WRITE ( 3,43 ) WIN(I,2),WIN(I,1)
C
C	COMPUTE THE INTEGRATION LIMITS CORRESPONDING TO THE WINDOW
C
	  WUU=(WIN(I,1)-THETA1)/DTHETA+1.
	  WLL=(WIN(I,2)-THETA1)/DTHETA+1.
	  IF ( DTHETA .LT. 0. ) THEN
	    INCR=1
	    IUL=INT(WUU+.99999)
	    IUU=IUL-1
	    ILL=INT(WLL)
	    ILU=ILL+1
	  ELSE
	    INCR=-1
	    IUL=INT(WUU)
	    IUU=IUL+1
	    ILL=INT(WLL+.99999)
	    ILU=ILL-1
	  END IF
	  IF ( IUL .GE. 1 .AND. IUU .GE. 1 .AND. ILL .GE. 1 .AND. ILU .GE. 1
     1      .AND. IUL .LE. NTH .AND. IUU .LE. NTH .AND. ILL .LE. NTH
     2      .AND. ILU .LE. NTH ) THEN
	    DO NI=1,NMAX
	      DO L=1,NMAX
	        IF ( W(NI,L,N_D,1) .NE. 0. ) THEN
	          CROSS(NI,L,N_D)=0.
	          DO ITH=IUL,ILL-INCR,INCR
		    CROSS(NI,L,N_D)=CROSS(NI,L,N_D)+
     1              (DS(NI,L,I,ITH)+DS(NI,L,I,ITH+INCR))/2.
     2              *ABS(THLB(ITH)-THLB(ITH+INCR))
	          END DO
C
C	INTERPOLATE AT THE BORDERS
c
C  does not work if WINLB(1) and WINLB(2) are within the same intervall
c  THLB(iul=ilu) und THLB(iuu=ill) 24.07.01
c
c	          M1=(DS(NI,L,I,IUL)-DS(NI,L,I,IUU))/(THLB(IUL)-THLB(IUU))
c	          B1=DS(NI,L,I,IUL)-M1*THLB(IUL)
c	          CROSS(NI,L,N_D)=CROSS(NI,L,N_D)+(M1*WINLB(I,1)+B1+
c     1            DS(NI,L,I,IUL))/2.*ABS(WINLB(I,1)-THLB(IUL))
c	          M2=(DS(NI,L,I,ILL)-DS(NI,L,I,ILU))/(THLB(ILL)-THLB(ILU))
c	          B2=DS(NI,L,I,ILL)-M2*THLB(ILL)
c	          CROSS(NI,L,N_D)=CROSS(NI,L,N_D)+(M2*WINLB(I,2)+B2+
c     1            DS(NI,L,I,ILL))/2.*ABS(WINLB(I,2)-THLB(ILL))
c
c  new version 24.07.01

                  M1=(DS(NI,L,I,IUL)-DS(NI,L,I,IUU))
                  m1=m1/(THLB(IUL)-THLB(IUU))
	          B1=DS(NI,L,I,IUL)-M1*THLB(IUL)
                  b1=b1+M1*WINLB(I,1)

	          M2=(DS(NI,L,I,ILL)-DS(NI,L,I,ILU))
                  m2=m2/(THLB(ILL)-THLB(ILU))
	          B2=DS(NI,L,I,ILL)-M2*THLB(ILL)
                  b2=b2+M2*WINLB(I,2)

                  if (iuu.ne.ill) then
                    b1=(b1+DS(NI,L,I,IUL))/2.
                    b1=b1*ABS(WINLB(I,1)-THLB(IUL))
                    b2=(b2+DS(NI,L,I,ILL))/2.
                    b2=b2*ABS(WINLB(I,2)-THLB(ILL))
                    CROSS(NI,L,N_D)=CROSS(NI,L,N_D)+b1+b2
                  else
                    b=(b1+b2)/2.*ABS(WINLB(I,2)-WINLB(I,1))
                    CROSS(NI,L,N_D)=CROSS(NI,L,N_D)+b
                  endif   

	        END IF
	      END DO
	    END DO
C
C	COMPUTE THE NORMALIZED CROSS-SECTION
C
	    CN=CROSS(NORMI,NORMF,N_D)
c            write(6,*) '  **  ',NORMI,NORMF,N_D,CROSS(NORMI,NORMF,N_D) 
            if (cn.eq.0.) cn=1.0 
	    LINE=50
	    DO NI=1,NMAX
	      DO L=1,NMAX
	        IF ( W(NI,L,N_D,1) .NE. 0. ) THEN
 	          IF ( LINE .GT. 45 ) THEN
	            IF ( IO(17) .eq. 3 ) then
	              call newpage
                      WRITE ( 6,44 ) n_d,I,WIN(I,2),WIN(I,1)
	              WRITE ( 6,45 )
                    end if 
	            LINE=3
	          END IF
	          CCEE=E(NI)-E(L)
	          CROSSNM=CROSS(NI,L,N_D)/CN
	          SCROSS(NI,L,N_D)=SCROSS(NI,L,N_D)+CROSS(NI,L,N_D)
                  if (n_d.eq.1) scross_ge(ni,l,i)=0.0
                  scross_ge(ni,l,i)=scross_ge(ni,l,i)+CROSS(NI,L,N_D)
	          IF ( IO(17) .eq. 3 ) then
                    if (((io(18).eq.1).and.(ni.lt.11).and.(l.lt.11)).or.
     +                (io(18).eq.0)) then
c                      write(6,*) '  *  ',cn,ni,l,n_d,W(NI,L,N_D,1)
         	      WRITE(6,46)NI,SPIN(NI),L,SPIN(L),CCEE,
     +                  CROSS(NI,L,N_D),CROSSNM
                      line=line+1
                    endif
	            IF ( IO(10) .NE. 0 ) then 
                      WRITE ( 3,47 ) NI,L,CROSS(NI,L,N_D)
	              LINE=LINE+1
                    endif
	          END IF
                end if
	      END DO
	    END DO
	  ELSE 
	   WRITE ( 6,48 ) I
	  END IF
	END DO

c***    end of second loop over PPAC-windows

42      FORMAT ( I2,2X,I2 )
43      FORMAT ( F6.2,2X,F6.2 )
44      FORMAT ( ' THE GAMMA-INT. OF DETECTOR #',i2,' IN', 
     +   ' WINDOW #',I2,' FROM ',F6.2,' TO ',F6.2,' DEGREES'/ )
45      FORMAT ( 5X,'STATE # INITIAL SPIN  STATE # ',
     +    'FINAL SPIN   ENERGY   CROSS-SECTION  NORM. CROSS-SECTION'/ )
46      FORMAT ( 7X,I2,7X,F4.1,9X,I2,7X,F4.1,7X,F6.4,4X,G11.4,5X,G11.4 )
47      FORMAT ( I2,2X,I2,2X,G12.5 )
48      FORMAT ( ' ****MISTAKE**** WINDOW #',I2,' IS OUTSIDE THE',
     +   'RANGE OF THE'/' CLX CALCULATIONS - WINDOW IS IGNORED' )


C	PRINTOUT OF SUM OF CROSS-SECTIONS
C
	IF ( IO(9) .NE. 0 ) THEN 
          CN=SCROSS(NORMI,NORMF,N_D)
          if(cn.eq.0.) cn=1.0  
	  LINE=50
	  DO NI=1,NMAX
	    DO L=1,NMAX
	      IF ( W(NI,L,N_D,1) .NE. 0. ) THEN
	        IF ( LINE .GT. 45 ) THEN
	          CALL NEWPAGE
	          WRITE ( 6,49 )N_D
	          WRITE ( 6,45 )
	          LINE=3
	        END IF
	        CCEE=E(NI)-E(L)
	        CROSSNM=SCROSS(NI,L,N_D)/CN
                if (((io(18).eq.1).and.(ni.lt.11).and.(l.lt.11)).or.
     +             (io(18).eq.0)) then
     	             WRITE(6,46)NI,SPIN(NI),L,SPIN(L),CCEE,
     +               SCROSS(NI,L,N_D),CROSSNM
	             LINE=LINE+1
                end if
	      END IF
	    END DO
	  END DO
	END IF
c
	END DO

c****       end of loop over Ge detectors

49      FORMAT(' THE GAMMA-INTENSITIES for',
     +    ' gamma-detector #',i2,' SUMMED OVER ALL PPAC WINDOWS'/)


c
c    Sum over the Ge detectors
c
   
	IF (((IO(16).NE.0).or.(io(17).eq.1)).AND.(IO(15).NE.0)) THEN

	  CNSUM = 0.0
	  DO NI=1,NMAX
	    DO L=1,NMAX
              SCROSS_SUM(NI,L)=0.0	
              WSUM(NI,L) = 0.0
	    END DO
	  END DO

	  DO N_D=1,N_DET
	    CNSUM = CNSUM + SCROSS(NORMI,NORMF,N_D)
	    DO NI=1,NMAX
	      DO L=1,NMAX
c                write(6,*) ' w ',ni,l,n_d,W(NI,L,N_D,1)
	        IF ( W(NI,L,N_D,1) .NE. 0. ) THEN
	          CCEE=E(NI)-E(L)
	          SCROSS_SUM(NI,L)=SCROSS_SUM(NI,L)+SCROSS(NI,L,N_D)
        	  WSUM(NI,L) = WSUM(NI,L) + W(NI,L,N_D,1)
C		      Write(6,*)SCROSS_SUM(NI,L),WSUM(NI,L)
C	              Write(6,*)THE_GAM(N_D),PHI_GAM(N_D)
C	              Write(6,*)N_D,N_DET,NI,Spin(NI),L,Spin(L)
C	              Write(6,*)CCEE,SCROSS(NI,L,N_D),CROSSNM
                end if 
	      END DO
	    END DO
	  END DO

        end if	 

c
C   Print sum of all Ge detectors for each PPAC-window
c
        IF ((IO(17).eq.1).AND.(IO(15).NE.0)) THEN
          do i=1,mzahl
            LINE=50
            cnge=scross_ge(NORMI,NORMF,i)
            if (cnge.eq.0.) cnge=1.0 
            DO NI=1,NMAX
              DO L=1,NMAX
                if (((io(18).eq.1).and.(ni.lt.11).and.(l.lt.11)).or.
     +             (io(18).eq.0)) then
c                write(6,*) ' wsum ',i,ni,l,WSUM(NI,L),SCROSS_ge(NI,L,i) 
                  IF ( WSUM(NI,L) .NE. 0. ) THEN
                    IF ( LINE .GT. 45 ) THEN
                      CALL NEWPAGE
                      WRITE ( 6,441 )i,win(i,2),win(i,1)
                      WRITE ( 6,45 )
                      LINE=3
                    END IF
                    CCEE=E(NI)-E(L)
                    CROSSNM=SCROSS_ge(NI,L,i)/CNge
                    WRITE(6,46)NI,SPIN(NI),L,SPIN(L),CCEE,
     +                   SCROSS_ge(NI,L,i),CROSSNM
                    LINE=LINE+1
                  END IF
                end if
              END DO
            END DO
          end do
        end if

 441    format('    THE SUM OF THE GAMMA-INTENSITIES IN WINDOW #',I2,
     +    ' FROM ',F6.2,' TO ',F6.2,' DEGREES'/) 


c
c   Print sum of all Ge detectors and of all PPAC-windows
c
        if (cnsum.eq.0.) cnsum=1.0
        if ((io(16).eq.1).and.(io(15).ne.0)) then 
          LINE=50
          DO NI=1,NMAX
	    DO L=1,NMAX
              if (((io(18).eq.1).and.(ni.lt.11).and.(l.lt.11)).or.
     +           (io(18).eq.0)) then 
                IF ( WSUM(NI,L) .NE. 0. ) THEN
	          IF ( LINE .GT. 45 ) THEN
	            CALL NEWPAGE
	            WRITE ( 6,149 )
	            WRITE ( 6,45 )
	            LINE=3
	          END IF
	          CCEE=E(NI)-E(L)
	          CROSSNM=SCROSS_SUM(NI,L)/CNSUM
	          WRITE(6,46)NI,SPIN(NI),L,SPIN(L),CCEE,
     +                SCROSS_SUM(NI,L),CROSSNM
	          LINE=LINE+1
	        END IF
              end if 
	    END DO
	  END DO
	END IF

149     FORMAT(' THE SUM OF THE GAMMA-INTENSITIES FOR ALL',
     +        ' GAMMA DETECTORS AND PPAC WINDOWS'/)


        if ((io(19).eq.1).or.(io(19).eq.3)) then

        call newpage
        write(6,*) '  POPULATIONMATRIX '
        WRITE(6,*) ' '
        write(6,*) '          1         2         3         4',
     +    '         5         6         7'
        write(6,*) ' 1        0         0         0         0',
     +    '         0         0         0'           
        do i=1,70
            write(6,'(2x,70i1)') (nnpop(i,j), j=1,70)
        end do

        endif


        if ((io(19).eq.2).or.(io(19).eq.3)) then

        line=50
        do i=1,nmax
          do j=1,nmax
            if (line .gt. 45) then
              call newpage
              write(6,*) '  DECAYSCHEME '
              write(6,*) ' '
              line=3
            endif
            if (dcyprob(i,j).ne.0) then
              write(6,*) i,j,dcyprob(i,j)
              line=line+1
            endif
          enddo
        enddo

        endif
               
 
	IF ( IO(10) .NE. 0 ) CLOSE ( UNIT=3 )
	STOP
	END 	
