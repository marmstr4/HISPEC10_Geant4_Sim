	FUNCTION FK(AJF,AL1,AL2,AJI,AK)
C
C	THE GAMMA-GAMMA CORRELATION COEFFICIENTS FK(JF,L1,L2,JI)
C
	FK=((-1.)**(INT(AJF-AJI-1.00000)))
	FK=FK*SQRT((2.*AL1+1.)*(2.*AL2+1.)*(2.*AJI+1.))
	FK=FK*CLEB(AL1,AL2,AK,1.,-1.,0.)*RACAH(AJI,AJI,AL1,AL2,AK,AJF)
	RETURN
	END
