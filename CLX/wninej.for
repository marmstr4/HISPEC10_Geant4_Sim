	FUNCTION wninej(A1,A2,A3,A4,A5,A6,a7,a8,a9)
C
C	THE WIGNER 9j SYMBOL
C       deShalit/Feshbach Theoret. Nucl. Phys. Vol.1 p.931 (2.107)
c       p.928 (2.86,2.87,2.88) 
c       TK 24.3.95
C

        j1=int(2.*abs(a4-a8))
        j2=int(2.*abs(a2-a6))
        j3=int(2.*abs(a1-a9))
        jmin=max0(j1,j2,j3)
        j1=int(2.*(a4+a8))
        j2=int(2.*(a2+a6))
        j3=int(2.*(a1+a9))
        jmax=min0(j1,j2,j3)        
        sa=a1+a2+a4+a6+a8+a9
        wine=0.0
        if (jmin .gt. jmax) goto 10 
        do j=jmin,jmax,2
          rj=real(j)/2.
          w=((-1.)**(int(2.*(rj+sa))))         
          w=w*(2.*rj+1.)
          w=w*racah(a1,a4,a9,a8,a7,rj)
          w=w*racah(a8,a2,a4,a6,a5,rj)
          w=w*racah(a6,a9,a2,a1,a3,rj)
          wine=wine+w
        end do
  10    wninej=wine 
	RETURN
	END
