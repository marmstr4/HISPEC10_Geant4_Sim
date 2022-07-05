        subroutine inty(n,k,kappa,ith,iph,yyint)
c       Integration of Y(K,KAPPA) over PHI of particle
c       TK 07.04.95/10.04./18.04.
 
        real aykkap(40,3,5,101,2),e(100),yyint(2)
        integer in_d(40),in_dk(40),ndetmaxx,ind(40)
        reaL THE_GAM(40),PHI_GAM(40)
 
        common / y / aykkap
        common /gamdet/ the_gam,phi_gam,in_d,in_dk,ind,ndetmaxx
        common /ylm/ a1,a2,iexc,ep,e,vcm,theta,v1,v2,nint
        
        dv=(v2-v1)/(real(nint))

        kap=abs(kappa)+1
        vorz=1.
        if (kappa .lt. 0) vorz=-1.

        yyint(1)=0.
        yyint(2)=0.
        y1r=aykkap(in_d(n),k,kap,1,1)
        y1i=vorz*aykkap(in_d(n),k,kap,1,2)
        yaltr=y1r
        yalti=y1i
        do i=2,nint+1
          y1r=aykkap(in_d(n),k,kap,i,1)
          y1i=vorz*aykkap(in_d(n),k,kap,i,2)
          yneur=y1r
          yneui=y1i 
          yyint(1)=yyint(1)+(yaltr+yneur)/2.
          yyint(2)=yyint(2)+(yalti+yneui)/2.
          yaltr=yneur
          yalti=yneui
        end do
        yyint(1)=yyint(1)*dv
        yyint(2)=yyint(2)*dv 

	RETURN
	END
