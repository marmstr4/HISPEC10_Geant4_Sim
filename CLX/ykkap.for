	subroutine YKKAP(ni,n_d,K,KAPPA,theta,phi,y)
C	TRANSFORMATION OF THE ANGLES FROM THE LAB INTO THE REST
C	SYSTEM (NP A192 (1972) 449)
C
C	TRANSFORMATION OF THE CROSS SECTION FROM THE REST INTO
C	THE LAB SYSTEM
C
C	RESTLAB = FUNCTION OF PARTICLE ANGLE ( AZIMUTHAL-ANGLE = X )
C
C	CALCULATION OF THE SPHERICAL HARMONIC FUNCTION Y(K,KAPPA)
C	INCLUDING THE TRANSFORMATION OF THE CROSS-SECTIONS FROM
C	THE REST INTO THE LAB SYSTEM
C
c       TK  29.3.95/31.03.
c           06.04.   now "right" KAPPAs and complex Ys are used
c                    K still means: K=1,2,3 <=> k=0,2,4 
c                    signs corresponding to Alder/Winther 1975
c                    SUBROUTINE
c                    Y(1) real part of Y
c                    Y(2) imaginary part of Y
c           07.04./10.04./18.04.
c
c       TK 01.07.01
c       projectile excitation: iexc=1; target excitation iexc=2


        real e(100),ka,y(2)
        integer in_d(40),in_dk(40),ndetmaxx,ind(40)
        reaL TH_GAM(40),PH_GAM(40)

        COMMON /YLM/ a1,a2,iexc,ep,e,vcm,thet,v1,v2,nint   
        common /gamdet/ th_gam,ph_gam,in_d,in_dk,ind,ndetmaxx

        sint=sin(th_gam(n_d))
        cost=cos(th_gam(n_d))
        vgamma=ph_gam(n_d)

        KA=(1.+A1/A2)*E(NI)/EP
        if (iexc.eq.2) then 
          DSKAPPA=SQRT(1.-KA)
          VP2=VCM*SQRT(1.+(DSKAPPA**2)-2.*DSKAPPA*COS(THETA))
          T2COS=(VCM/VP2)*(1.-DSKAPPA*COS(THETA))
c        T2SIN=SQRT(1.-T2COS*T2COS)
        endif
        if (iexc.eq.1) then 
          DSKAPPA=A2/A1*SQRT(1.-KA)
          VP2=VCM*SQRT(1.+(DSKAPPA**2)+2.*DSKAPPA*COS(THETA))
          T2COS=(VCM/VP2)*(1.+DSKAPPA*COS(THETA))
c        T2SIN=SQRT(1.-T2COS*T2COS)
        endif
        vorz=1.
        if (t2cos .lt. 0.0) vorz=-1.
        if (abs(t2cos) .gt. 1.0) t2cos=vorz
        T2SIN=sin(acos(T2COS))

        COSVG=COS(VGAMMA-phi)
        TG2COS=COST*T2COS-T2SIN*SINT*COSVG
        RESTLAB=ABS((1.-VP2*VP2)/((1.-VP2*TG2COS)**2))
        TGCOS=(COST*SQRT(1.-VP2*VP2)-
     +     T2COS*VP2+T2COS*TG2COS*(1.-SQRT(1.-VP2*VP2)))
     +     /(1.-VP2*TG2COS)
        vorz=1.
        if (tgcos .lt. 0.0) vorz=-1. 
        if (abs(tgcos) .gt. 1.0) tgcos=vorz
        TGSIN=sin(acos(TGCOS))

        VICOS=(-1.)*(TG2COS-T2COS*TGCOS-
     +     VP2*(1.-T2COS*TGCOS*TG2COS))
     +     /(TGSIN*T2SIN*(1.-VP2*TG2COS))
c        vialt=vicos
        if (iexc.eq.1) vicos=-vicos
        
        vorz=1.
        if (vicos .lt. 0.0) vorz=-1.
        if (abs(vicos) .gt. 1.0) vicos=vorz
c        write(*,*) th_gam(n_d),ph_gam(n_d)
c        write(*,*) thet,ith,v1,v2,ni
c        write(*,*) ni,n_d,K,KAPPA,theta,phi
c        write(*,*) vicos,vicos+1.,-vicos/1.,vialt,vialt+1.,-vialt/1.
c        write(*,*) vialt-vicos
c        write(*,*) acos(vicos)
c        write(*,*) sin(acos(vicos))
c        write(*,*) ' ***********************************'
        visin=sin(acos(viCOS)) 

        kap=abs(kappa)
        sig=1.
        if (kappa .lt. 0) then 
          sig=-1.
          restlab=restlab*((-1.)**(kap))
        end if

	IF ( KAP .EQ. 0 .AND. K .EQ. 1 ) THEN    ! Y(0,0)
	 y(1)=RESTLAB*.282095
         y(2)=0.
	 RETURN
	END IF

	IF ( K .EQ. 2 ) THEN
	 IF ( KAP .EQ. 0 ) THEN                  ! Y(2,0)
	  y(1)=RESTLAB*.630783*(TGCOS*TGCOS-.5*tgsin*tgsin)
          y(2)=0.
	  RETURN
	 ELSE IF ( KAP .EQ. 1 ) THEN             ! Y(2,1), Y(2,-1)
	  y(1)=RESTLAB*(-0.772548)*TGCOS*TGSIN
          y(2)=sig*y(1)
          y(1)=y(1)*vicos
          y(2)=y(2)*visin
	  RETURN
	 ELSE IF ( KAP .EQ. 2 ) THEN             ! Y(2,2), Y(2,-2)
	  y(1)=RESTLAB*0.386274*TGSIN*TGSIN
          y(2)=sig*y(1)
          y(1)=y(1)*(vicos*vicos-visin*visin)
          y(2)=y(2)*visin*vicos
	  RETURN
	 END IF
	END IF

	IF ( K .EQ. 3 ) THEN
	 IF ( KAP .EQ. 0 ) THEN                  ! Y(4,0)
	  y(1)=RESTLAB*.105786*(35.*TGCOS*TGCOS*TGCOS*
     +                          TGCOS-30.*TGCOS*TGCOS+3.)
          y(2)=0.
	  RETURN
	 ELSE IF ( KAP .EQ. 1 ) THEN             ! Y(4,1), Y(4,-1)
	  y(1)=RESTLAB*(-.473087)*(7.*TGCOS*TGCOS-3.)*TGSIN*
     +         TGCOS
          y(2)=sig*y(1)
          y(1)=y(1)*vicos
          y(2)=y(2)*visin
	  RETURN
	 ELSE IF ( KAP .EQ. 2 ) THEN             ! Y(4,2), Y(4,-2)
	  y(1)=RESTLAB*.3345248*(7.*TGCOS*TGCOS-1.)*TGSIN*
     +         TGSIN
          y(2)=sig*y(1)
          y(1)=y(1)*(vicos*vicos-visin*visin)
          y(2)=y(2)*visin*vicos
	  RETURN
	 ELSE IF ( KAP .EQ. 3 ) THEN             ! Y(4,3), Y(4,-3)
	  y(1)=RESTLAB*(-1.251671)*TGSIN*TGSIN*TGSIN*TGCOS
          y(2)=sig*y(1)
          y(1)=y(1)*(4.*vicos*vicos*vicos-3.*vicos)
          y(2)=y(2)*(3.*visin-4.*visin*visin*visin)
	  RETURN
	 ELSE IF ( KAP .EQ. 4 ) THEN             ! Y(4,4), Y(4,-4)
	  y(1)=RESTLAB*.442533*TGSIN*TGSIN*TGSIN*TGSIN
          y(2)=sig*y(1)
          y(1)=y(1)*(8.*vicos*vicos*vicos*vicos
     +               -8.*vicos*vicos+1.)
          y(2)=y(2)*(8.*vicos*vicos*vicos*visin
     +                   -4.*vicos*visin)
	  RETURN
	 END IF
	END IF

        y(1)=0.
        y(2)=0.

	RETURN
	END
