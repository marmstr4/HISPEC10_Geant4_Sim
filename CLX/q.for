	SUBROUTINE Q
C
C	COMPUTE QLM FUNCTIONS
C
	DIMENSION LAMDA(5)
	REAL R_QLM(5,5), I_QLM(5,5)
	COMMON /QLM/W,EPS,EROOT,LAMMAX,LAMDA,R_QLM,I_QLM,RALFA,ZPOL
C
	EW=EXP(W)
	COSHY=.5*(EW+1./EW)
	SINHY=.5*(EW-1./EW)
	DEN=EPS*COSHY+1.
	DEN1=DEN*DEN
	DEN2=DEN1*DEN1
	SH1=EROOT*SINHY
	SH2=SH1*SH1
	SH3=SH2*SH1
	CH1=COSHY+EPS
	CH2=CH1*CH1
	DO I1=1,LAMMAX
C
C	E1
C
	  IF ( LAMDA(I1) .EQ. 1 ) THEN
	   R_QLM(I1,1)=.5*CH1/DEN1
	   I_QLM(I1,1)=0.
	   R_QLM(I1,2)=0.
	   I_QLM(I1,2)=-.3535534*SH1/DEN1
C
C	E2
C
	  ELSE IF ( LAMDA(I1) .EQ. 2 ) THEN
C
C	POL ACCOUNTS FOR E1-POLARIZATION
C
	   POL=1.0-ZPOL/DEN
	   R_QLM(I1,1)=.75*(2.*CH2-SH2)/DEN2*POL
	   I_QLM(I1,1)=0.
	   R_QLM(I1,2)=0.
	   I_QLM(I1,2)=-1.8371173*CH1*SH1/DEN2*POL
	   R_QLM(I1,3)=-.9185587*SH2/DEN2*POL
	   I_QLM(I1,3)=0.
C
C	E3
C
	  ELSE IF ( LAMDA(I1) .EQ. 3 ) THEN
	   DEN3=DEN2*DEN1
	   R_QLM(I1,1)=1.875*CH1*(2.*CH2-3.*SH2)/DEN3
	   I_QLM(I1,1)=0.
	   R_QLM(I1,2)=0.
	   I_QLM(I1,2)=-1.6237976*(4.*CH2-SH2)*SH1/DEN3
	   R_QLM(I1,3)=-5.1348989*CH1*SH2/DEN3
	   I_QLM(I1,3)=0.
	   R_QLM(I1,4)=0.
	   I_QLM(I1,4)=2.0963137*SH3/DEN3
C
C	E4
C
	  ELSE IF ( LAMDA(I1) .EQ. 4 ) THEN
	   DEN4=DEN2*DEN2
	   CH4=CH2*CH2
	   SH4=SH2*SH2
	   R_QLM(I1,1)=1.09375*(8.*CH4-24.*CH2*SH2+3.*SH4)/DEN4
	   I_QLM(I1,1)=0.
	   R_QLM(I1,2)=0.
	   I_QLM(I1,2)=-4.8913987*CH1*(4.*CH2-3.*SH2)*SH1/DEN4
	   R_QLM(I1,3)=-3.4587411*(6.*CH2-SH2)*SH2/DEN4
	   I_QLM(I1,3)=0.
	   R_QLM(I1,4)=0.
	   I_QLM(I1,4)=12.9414244*CH1*SH3/DEN4
	   R_QLM(I1,5)=4.5754844*SH4/DEN4
	   I_QLM(I1,5)=0.
C
C	M1
C
	  ELSE IF ( LAMDA(I1) .EQ. 5 ) THEN
	   R_QLM(I1,1)=0.
	   I_QLM(I1,1)=0.
	   R_QLM(I1,2)=-.3535534*EROOT/DEN1
	   I_QLM(I1,2)=0.
	  END IF
	END DO
	RALFA=EPS*SINHY+W
	RETURN
	END
