	SUBROUTINE SPLINE(N)
C
	REAL*8 A(23),B(4),C(23),D(4,23),E(23),F(23),H(23),S(23),T(23)
	REAL*8 X,Y,ST(23),AA,Z,ZP
	COMMON / C1 / A,B,C,D,E,F,H,S,T,X,Y,ST
C
	M=N-1
	  DO 100 I=1,M
	  E(I)=C(I+1)-C(I)
100	  F(I)=A(I+1)-A(I)
	MS=N-2
	  DO 200 I=1,MS
	  D(1,I)=E(I+1)
	  D(2,I)=2.*(E(I)+E(I+1))
	  D(3,I)=E(I)
200	  D(4,I)=3.*(E(I+1)*F(I)/E(I)+E(I)*F(I+1)/E(I+1))
	S(1)=-1.
	T(1)=2.*F(1)/E(1)
	  DO 300 I=2,M
	  L=I-1
	  S(I)=-(D(1,L)+D(2,L)*S(L))/(D(3,L)*S(L))
300	  T(I)=(-T(L)*(D(2,L)+D(3,L)*S(I))+D(4,L))/D(3,L)
	AA=2.*F(M)/E(M)
	H(N)=(T(M)+AA*S(M))/(S(M)+1.)
	  DO 400 I=1,M
	  K=N-I
400	  H(K)=(H(K+1)-T(K))/S(K)
	  DO 500 K=1,M
500	  IF ( X .LT. C(K+1) ) GO TO 600
600	K1=K
	Z=(C(K1+1)-C(K1))/2.
	B(3)=(H(K1+1)-H(K1))/(4.*Z)
	B(4)=.25*((H(K1+1)+H(K1))*Z-A(K1+1)+A(K1))/Z**3
	B(1)=(A(K1+1)+A(K1)-2.*B(3)*Z**2)/2.
	B(2)=(H(K1+1)+H(K1)-6.*B(4)*Z**2)/2.
	ZP=X-(C(K1+1)+C(K1))/2.
	Y=B(1)+ZP*(B(2)+ZP*(B(3)+ZP*B(4)))
	RETURN
	END
