	FUNCTION DJMM(BETA,RJ,RM,RMP)
C
	REAL*8 FCT(99),FACT(99)	
	COMMON/FCTRLS/FCT,FACT
C
	J=2.01*RJ
	M=2.01*RM
	MP=2.01*RMP
	IF ( J .LT. ABS(M) .OR. J .LT. ABS(MP) .OR. J .LT. 0 ) THEN
	 WRITE ( 6,1 ) 
1	FORMAT ( ' ***MISTAKE*** ILLEGAL ARGUMENTS IN CALL FOR DJMM' )
	 RETURN
	END IF
	JA=(J+MP)/2+1
	JB=(J-MP)/2+1
	JC=(J+M)/2+1
	JD=(J-M)/2+1
	B1=FACT(JA)+FACT(JB)+FACT(JC)+FACT(JD)
	MINSIG=-MIN(0,M+MP)
	MAXSIG=J-MAX(MP,M)
	ISIG=MINSIG
	SUM=0.
	DO WHILE ( ISIG .LE. MAXSIG )
	  JA=ISIG/2+1
	  JB=(J-M-ISIG)/2+1
	  JC=(J-MP-ISIG)/2+1
	  JD=(M+MP+ISIG)/2+1
	  IF ( MOD((J-MP-ISIG)/2,2) .EQ. 0 ) THEN
	   FASE=1.
	  ELSE
	   FASE=-1.
	  END IF
	  IC=ISIG+(M+MP)/2
	  IS=J-ISIG-(M+MP)/2
	  B2=FACT(JA)+FACT(JB)+FACT(JC)+FACT(JD)
	  SUM=SUM+FASE*COS(BETA/2.)**IC*SIN(BETA/2)**IS*
	1EXP(B1)/(EXP(B2)**2.)
	  ISIG=ISIG+2
	END DO
	DJMM=SUM
	RETURN
	END
