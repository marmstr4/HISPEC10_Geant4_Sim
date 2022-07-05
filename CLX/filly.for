	 subroutine filly(ni,ith,iph,io13,n_d)
c        Fills array with values of Y(K,KAPPA)
c        In order to save memory only Y(K,KAPPA) with positive KAPPA
c        are stored. If the array is used, the values have to be
c        multiplied by (-1.)**(KAPPA) if KAPPA < 0 to obtain the right
c        signs used in Alder/Winther 1975
c        Furthermore if KAPPA < 0 the imaginary part must be multiplied
c        with -1.
c        The index KAP used in array equals KAPPA+1 
c        TK 31.03.95/07.04./10.04./18.04./19.05.


        real aykkap(40,3,5,101,2),e(100),y(2) 
        integer in_d(40),in_dk(40),ndetmaxx,ind(40)
        reaL THE_GAM(40),PHI_GAM(40)

        common / y / aykkap
        common / gamdet / the_gam,phi_gam,in_d,in_dk,ind,ndetmaxx
        COMMON /YLM/ a1,a2,iexc,ep,e,vcm,theta,v1,v2,nint

        dv=(v2-v1)/(real(nint))

        ndetx=ndetmaxx
        if (io13 .eq. 1) ndetx=1
   
        do i=1,ndetx
          ii=i
          if (io13 .eq. 1) ii=in_d(n_d)
          if (ind(ii) .eq. 1) then
            do k=1,3
              do kappa=0,2*k-2
                do j=1,nint+1
                  phi=v1+(j-1)*dv
                  call ykkap(ni,ii,k,kappa,theta,phi,y)
                  aykkap(ii,k,kappa+1,j,1)=y(1)
                  aykkap(ii,k,kappa+1,j,2)=y(2)
                enddo
              enddo
            enddo
          endif
        enddo 
     
        return
        END
