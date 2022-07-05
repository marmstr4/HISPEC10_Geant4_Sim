        subroutine intyy(n1,k1,kappa1,n2,k2,kappa2,ith,iph,yyint)
c       Integration of Y(K1,KAPPA1)*Y(K2,KAPPA2) over PHI of particle
c       TK 29.03.95/31.03./18.04.
 
        real aykkap(40,3,5,101,2),e(100),yyint(2)
        integer in_d(40),in_dk(40),ndetmaxx,ind(40)
        reaL THE_GAM(40),PHI_GAM(40)
 
        common / y / aykkap
        common /gamdet/ the_gam,phi_gam,in_d,in_dk,ind,ndetmaxx
        common /ylm/ a1,a2,iexc,ep,e,vcm,theta,v1,v2,nint
        
        dv=(v2-v1)/(real(nint))

        kap1=abs(kappa1)+1
        kap2=abs(kappa2)+1
        vorz1=1.
        vorz2=1.
        if (kappa1 .lt. 0) vorz1=-1.
        if (kappa2 .lt. 0) vorz2=-1.

        yyint(1)=0.
        yyint(2)=0.
        y1r=aykkap(in_d(n1),k1,kap1,1,1)
        y1i=vorz1*aykkap(in_d(n1),k1,kap1,1,2)
        y2r=aykkap(in_dk(n2),k2,kap2,1,1)
        y2i=vorz2*aykkap(in_dk(n2),k2,kap2,1,2)
        yaltr=y1r*y2r-y1i*y2i
        yalti=y1r*y2i+y1i*y2r
        do i=2,nint+1
          y1r=aykkap(in_d(n1),k1,kap1,i,1)
          y1i=vorz1*aykkap(in_d(n1),k1,kap1,i,2)
          y2r=aykkap(in_dk(n2),k2,kap2,i,1)
          y2i=vorz2*aykkap(in_dk(n2),k2,kap2,i,2)
          yneur=y1r*y2r-y1i*y2i
          yneui=y1r*y2i+y1i*y2r 
          yyint(1)=yyint(1)+(yaltr+yneur)/2.
          yyint(2)=yyint(2)+(yalti+yneui)/2.
          yaltr=yneur
          yalti=yneui
        end do
        yyint(1)=yyint(1)*dv
        yyint(2)=yyint(2)*dv 

	RETURN
	END
