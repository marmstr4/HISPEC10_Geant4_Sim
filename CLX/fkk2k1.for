	FUNCTION FKk2k1(AJF,AL1,AL2,AJI,AK,ak2,ak1)
C
C	THE GAMMA-GAMMA CORRELATION COEFFICIENTS FKK2K1(JF,L1,L2,JI)
C       Alder/Winther (1975) p.334
c       TK 27.3.95
C
	FKkk=((-1.)**(INT(AK2+AK1+1.0-AL1)))
	FKkk=FKkk*SQRT((2.*AJI+1.)*(2.*AJF+1.))
        fkkk=FKkk*SQRT((2.*AL1+1.)*(2.*AL2+1.))
        fkkk=FKkk*SQRT((2.*Ak1+1.)*(2.*Ak2+1.))
	FKkk=FKkk*CLEB(AL1,AL2,AK,1.,-1.,0.)
        fkkk=fkkk*wninej(AJf,Al1,Aji,Ajf,Al2,AJi,ak2,ak,ak1)
        fkk2k1=fkkk 
	RETURN
	END
