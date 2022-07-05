C====================================================================
C
			PROGRAM DWEIKO
C
C Coupled-channels program for nuclear scattering at high energies.
C Authors: C.A. Bertulani, C. Campbell and T. Glasmacher -  06/28/2002 
C   
C Elastic scattering, Coulomb excitation of E1, E2, E3, M1 and M2 modes
C + Nuclear excitation for E0, E1, E2. Total reaction cross sections.
C
C Units are in MeV and fm. Output cross sections in milibarns. 
C   Angular momentum notation of: 
C   A.R. Edmonds, "Angular Momentum in Quantum Mechanics",
C   Princeton University Press (1957).
C--------------------------------------------------------------------
C        SEE END OF PROGRAM FOR INSTRUCTIONS HOW TO MAKE AN INPUT
C--------------------------------------------------------------------
C   Inputs are in files 'dweiko.dim', 'dweiko.in' and 'optw.in' (input of 
C   optical potential: optional).
C---------------------------------------------------------------------------
C                     START PROGRAM
C---------------------------------------------------------------------------
C
      implicit real*8(a-h,o-z)
      include 'dweiko.dim'
C
      integer*4 mj(nstmax),ii(0:nstmax),istop(0:nmax),istat(ngrid)
      real*8 mm(nstmax),jj(nstmax),intb(nbmax),b(nbmax),mnkap,
     & mat(5,0:nmax,0:nmax),spin(nmax),fd(5),delte(0:3,nmax),
     & exx(nstmax),pcc(nmax,nbmax),ex(nmax),dens1(0:ngrid),
     & dens2(0:ngrid),dru(0:ngrid),r(0:ngrid),z(0:ngrid),tb(nbmax),
     & ru(0:ngrid),psitot(5,nstmax,nstmax),psi(5,nbmax,nstmax,nstmax),
     & xi(nstmax,nstmax),spinave(nmax,nbmax),pint(nbmax),pint_2(nbmax)
     & ,bcorr(nbmax),simp(0:ngrid),sigj(nmax),frd0(nmax),frd1(nmax),
     & frd2(nmax),frd3(nmax),theta(ngrid)
      complex*16 ei,ca(nstmax),dca(nstmax),wn(0:ngrid),wn2(0:ngrid),
     & phasen(nbmax),psinuc(0:3,0:ngrid,nstmax,nstmax),wn1(0:ngrid),
     & fact,ubm(0:3,0:ngrid),stat(nmax,0:10,-10:10),
     & phases(nbmax)
     & ,vsp(0:ngrid),cmu(-20:20,nbmax)
      real*8 iv(0:130),fak(0:130),fad(0:130)
C
      common /gf1 / iv
      common /gf2 / fak
      common /gf3 / fad
      common/redfak/fakm(100)
      common/a0/ap,zp,at,zt
      common/a1/ beta,gamma,pi,hc,gcm,betacm
      common/a2/signn,alnn
      common/a3/xi
      common/a4/mm,n
      common/a5/psitot
      common/a6/psinuc
      common/a7/bb
      common/a8/iopnuc
      common/a9/delr
      common/a10/jj,istop
      common/a11/mat,ex,sigj
      common/a12/ii
      common/a13/frd0,frd1,frd2,frd3
      common/a14/bmin1,bmax1,eta,pcmhc,b13			    
      common/a15/bmax,delb1,delb2,spin			    
      common/a16/nb1
      common/lix/phasen			    
      common/simps/ simp
C
      character*14 fln_out,fln_out2
C
      external dcadt,rkqs
C
C initialize factorial routines
      call gfv(30)
      call factor
C
      open(10,file='dweiko.in', status='old')
      open(12,file='dweiko.out', status='unknown')
C
C  'rsharps' allows for comment lines in input file
      call rsharps(10)
C
      read(10,*)ap,zp,at,zt,eca
      call rsharps(10)
      read(10,*)iw,iopm,ioelas,ioinel,iogam 
      call rsharps(10)
      read(10,*)nb,accur,bmin,iob
      if (nb.gt.nbmax) then
       write(12,*) '******    STOP: NB LARGER THAN NBMAX'
       write(12,*)' Decrease value of NB, or increase NBMAX'
       stop
      endif
C
C constants 
      zero = 0.d0
      one = 1.d0
      ei = cmplx(0.,1.)
      pi = acos(-1.)
      e2 = 1.44
      amu = 931.49432D0
      hc  = 197.3270533D0
      e2ohbc = 1./137.
      cpion = 1.41
C cm and laboratory energies and momenta for the nuclei
      elab=eca*ap
	redm=ap*at/(ap+at)
      apmass=ap*amu
      atmass=at*amu
      plab=sqrt(elab*(elab+2.*apmass))
      gamma=(eca+amu)/amu
      beta=plab/(elab+apmass)
      ecm=sqrt(atmass**2+apmass**2+2.*(elab+apmass)*atmass)        
      pcm=plab*atmass/ecm
      pcmhc=pcm/hc
      rho=ap/at
      e1=eca/amu
      gcm=(1.+e1+rho)/sqrt((1+rho)**2+2.*rho*e1)
      betacm=sqrt((gcm+1.)*(gcm-1.))/gcm
	rhog=rho*(1+rho*(1.+e1))/((1.+e1)+rho)
      a0=zp*zt*e2/(redm*amu*beta*beta)
      aclose=a0/gamma
      eta=a0*pcmhc
c for the nucleons
      plab_nn=plab/ap
      ecm_nn=sqrt(2.*amu**2+2.*(elab/ap+amu)*amu)        
      pcm_nn=plab_nn*amu/ecm_nn
      pcm_nn=pcm_nn/hc
C
C grazing impact parameter
      rp13=1.2*ap**(1./3.)
      rt13=1.2*at**(1./3.)  
      b13=rp13+rt13
      binf=1.d-10
      if(bmin.ne.0.) then
         thet_bm =  2.*atan(aclose/(bmin))*180./pi
        else
         thet_bm = 2.*atan(aclose/(b13+binf))*180./pi
      endif
C 
C mesh in impact parameter and integration factors
C bcorr is impact parameter corrected for Coulomb deflection
      bmin1=b13/3.
      bmax1=2.*b13
      bmax=200.
      nb1=int(nb/2)
      if(iv(nb1).ge.0.)nb1=nb1+1
      delb1=(bmax1-bmin1)/(nb1-1)
C separate impact parameter into two regions: close collisions
      do ib=1,nb1
       b(ib)=bmin1+(ib-1)*delb1
       bcorr(ib)=aclose+dsqrt(aclose**2+b(ib)**2)
      enddo
      nb2=nb-nb1
      delb2=(bmax-bmax1)/nb2
C distant collisions
      do ib=nb1+1,nb
       b(ib)=bmax1+(ib-nb1)*delb2
       bcorr(ib)=aclose+dsqrt(aclose**2+b(ib)**2)
      enddo
C
C mesh in radial coordinates (rmax=maximum range of the potential)
      rmin = 0.001
      rmax = 50.
      delr=(rmax-rmin)/dfloat(ngrid)
      do i=0,ngrid
	 r(i)=rmin + i*delr
	 z(i)=rmin + i*delr
      enddo
C
C factors for some of the Simpson's integrations
      simp(0)=1./3.
      simp(ngrid)=1./3.
      ig=4
      do i=1,ngrid-1
       simp(i)=ig/3.
       ig=6-ig
      enddo
C
C nucleon-nucleon cross section (in fm^2)
      call signne(eca,sigpp,sigpn)
      signn = 0.1 * (zp*zt+(ap-zp)*(at-zt))/ap/at * sigpp +
     &	    0.1 * (zp*(at-zt)+zt*(ap-zp))/ap/at * sigpn
C real part of the nucleon-nucleon scattering amplitude
      call phnne(eca,phpp,phpn,ahpp,ahpn)
      alnn = (zp*zt+(ap-zp)*(at-zt))/ap/at * phpp +
     &	    (zp*(at-zt)+zt*(ap-zp))/ap/at * phpn
C
C Index T is reserved for the excited nucleus (depends on option IW)
C
      IAp = idnint(ap)
      IZp = idnint(zp)
      IAt = idnint(at)
      IZt = idnint(zt)
C
      if(iw.eq.0) then
       dum1=ap
       dum2=zp
       ap=at
       zp=zt
       at=dum1
       zt=dum2
      endif
C
      call rsharps(10)
      read(10,*)iopw,iopnuc
      if(iopw.eq.2) open(11,file='optw.in', status='old')
C
      do i=0,ngrid
        wn(i)=0.
        vsp(i)=0.
      end do
C Woods-Saxon optical potential
      if(iopw.eq.1) then
        call rsharps(10)
        read(10,*)V0_ws,r0,d_ws,VI_ws,r0_I,dI_ws
        R_ws=r0*(ap**0.3333+at**0.3333)
        RI_ws=r0_I*(ap**0.3333+at**0.3333)
        call omp_ws(r,V0_ws,VI_ws,R_ws,RI_ws,d_ws,dI_ws,wn,iopm)
C  for protons include spin-orbit and surface term
        if(IAp.eq.1.or.IAt.eq.1) then
          call rsharps(10)
          read(10,*)VS,r0_S,dS,V_surf,d_surf
          RS=r0_S*(ap**0.3333+at**0.3333)
          do i=1,ngrid
            temp=(r(i)-RS)/dS
            if(temp.lt.30.)vsp(i)=-2.*cpion**2/r(i)*VS/dS*
     & 		  exp(temp)/(1.+exp(temp))**2  
            temp=(r(i)-RI_ws)/dI_ws
            if(temp.lt.30.)wn(i)=wn(i)+4.*ei*V_surf*d_surf/dI_ws*
     & 		  exp(temp)/(1.+exp(temp))**2  
          enddo 
          call phnuc(vsp,bcorr,nb,z,phases)
        endif
      endif
C optical potential read from file 'optw.in'
      if(iopw.eq.2)then
        call omp_read(r,wn,iopm)
      endif
C t-rho-rho optical potential
      if(iopw.eq.3)then
        call omp_den(r,rmin,rmax,wn,iopm)
      endif
C M3Y optical potential
      if(iopw.eq.4)then
        call rsharps(10)
        read(10,*)Wrat
        call omp_M3Y(r,rmin,rmax,wn,Wrat,iopm)
      endif
C
C nuclear eikonal phase and absorption factor 
      if(iopw.eq.0) then
	  ioelas=0.
        call omp_den(r,rmin,rmax,wn,0)
      endif
      call phnuc(wn,bcorr,nb,z,phasen)
      if(iopw.eq.3)call phnucF(pcm_nn,eca,b,r,delr,phasen)
      sum=0.
      do ib=1,nb
	 tb(ib)=exp(-2.*dimag(phasen(ib)))
      enddo 
C
C elastic or inelastic scattering cross section
      if(ioelas.ge.1.or.ioinel.ge.1) then
        call rsharps(10)
        read(10,*)thmax,ntheta
        if (ntheta.gt.ngrid) then
          write(12,*) '******    STOP: NTHETA LARGER THAN NGRID '
          write(12,*)' Decrease value of NTHETA, or increase NGRID'
         stop
        endif
        dtheta=thmax/ntheta
        do i=1,ntheta
          theta(i)=i*dtheta
        enddo 		             
        if(ioelas.ge.1) then
          if(IAp.gt.1.and.IAt.gt.1)
     &     call elast(b,delb1,phasen,ioelas,rhog,theta,ntheta) 
          if(IAp.eq.1.or.IAt.eq.1)
     &       call elastp(b,delb1,phasen,phases,ioleas,rhog,theta,ntheta)		  					  	         
        endif
        if(ioinel.ge.1) then
          call rsharps(10)
          read(10,*)jinel
        endif 	   
      endif
C
      write(12,*)' '
      write(12,'('' ****       Output of progam DWEIKO       ****'')')
      write(12,*)' '
      write(12,109)
109   format(' Input parameters:',/,' -----------------')
      write(12,1)IAp,IZp,IAt,IZt,eca,iw,ngrid,nb,bmin,accur,iopm,
     &	ioelas,ioinel,iogam
C     
1     format(2x,' Nuclear mass and charge numbers             ',
     + '           :'
     1 ,2(1x,2(1x,i3)),/,2x, 
     2 ' Laboratory energy/per nucleon [MeV]                    :',
     3  1x,F7.1,/,2x,
     4 ' Option for the excited nucleus [0(1) for proj (targ)]  :',
	5 1x,I2,/,2x,
     6 ' Integration parameters (NGRID,NB,BMIN,ACCUR)           :'
	7  ,1x,I4,2x,I4,2x,F5.1,2x,F9.4,/,2x,
     8 ' Output of optical potential 0(no) 1(yes)               :',
     9 1x,I2,/,2x
     & ' Output of elastic cross section 0(no) 1(cm)  2(lab.)   :',
     9 1x,I2,/,2x
     & ' Output of inelastic cross section 0(no) 1(cm)  2(lab)  :',
     & 1x,I2,/,2x
     & ' Statistical tensors(1)  gamma-rays distrib(2)  none(0) :',
     & 1x,I2)
      write(12,'(a,1x,i1)')
     & "   Nuclear excitation [0: none   1: vibrational]          : "
     & ,iopnuc  
      write(12,'(a,1x,i1)')
     & "   Optical Pot. [0:no  1:WS   2:read  3:t-rho-rho  4:M3Y] : "
     & ,iopw       
      if(iopw.eq.1)then
        write(12,'(a,6(1x,F7.2))')
     &  "   V0_real[MeV],R[fm],d[fm],V0_imag[MeV],R_I[fm],d_I[fm]  :"
     &  ,V0_ws,R_ws,d_ws,VI_ws,RI_ws,dI_ws 
        if(IAp.eq.1.or.IAt.eq.1) then
           write(12,'(a,5(1x,F7.2))')
     &  "   V_sp-orb[MeV],RS[fm],dS[fm],V_surf[MeV],d_surf[fm]     :"
     &     ,VS,RS,dS,V_surf,d_surf
        endif 
      endif
C
      call rsharps(10)
      read(10,*)nst
      if (nst.gt.nmax) then
       write(12,*) '******    STOP: NST LARGER THAN NMAX '
       write(12,*)' Decrease value of NST, or increase NMAX'
       stop
      endif
      write(12,2)nst
2     format(2x,' Number of nuclear states = ',I2,//,
     1   '  State     E [MeV]     Spin ')
C
      if (ioinel.eq.1) then
	  if(jinel.gt.nst.or.jinel.eq.0) then
          write(12,*) '******    STOP: JINEL LARGER THAN NST '
          write(12,*)'If IOINEL = 1, enter 2 <= JINEL <= NST'
          stop
        endif
      endif
C
      call rsharps(10)
      do i=1,nst
	 read(10,*)j,ex(j),spin(j)
       call rsharps(10)
	 write(12,3)j,ex(j),spin(j)
3      format(2x,I3,3x,F10.3,4x,F5.1)
      enddo
C
      write(12,4)
4     format(/,' Reduced matrix elements - units: [M/E;L = e * fm^L]',/,
     1         ' ---------------------------------------------------', 
     2       /,2x,14x,'E1',8x,'E2',9x,'E3',9x,'M1',9x,'M2')
C read matrices (E/M;L, I_i-->I_j) values
      call rsharps(10)
      do k=1,nstmax
       read(10,*)i,j,MAT(1,i,j),MAT(2,i,j),MAT(3,i,j),
     &	 MAT(4,i,j),MAT(5,i,j)
       if(i.eq.0.and.k.eq.1) then
         write(12,*)'No electromagnetic matrix elements were given.'
       endif
       if(i.eq.0) goto 51
       write(12,5)i,j,mat(1,i,j),mat(2,i,j),mat(3,i,j),
     &     mat(4,i,j),mat(5,i,j)
      enddo
5     format(i2,' --> ',i2,2X,5(f9.4,2x))    
C
C output of B-values
51    if(k.ne.1) write(12,52)
52    format(/,' B(E/M;L)-values - units: [e^2 * b^L]',/,
     1         ' ------------------------------------', 
     2        /,2x,14x,'E1',8x,'E2',9x,'E3',9x,'M1',9x,'M2')
      do i=1,nst
        do j=1,nst
          be1=mat(1,i,j)**2/(2.*spin(1)+1.)/100
          be2=mat(2,i,j)**2/(2.*spin(1)+1.)/10**4
          be3=mat(3,i,j)**2/(2.*spin(1)+1.)/10**6
          bm1=mat(4,i,j)**2/(2.*spin(1)+1.)/100.
          bm2=mat(5,i,j)**2/(2.*spin(1)+1.)/10**4
      if(be1.ne.0.or.be2.ne.0.or.be3.ne.0.or.bm1.ne.0.or.bm2.ne.0)         
     &		  write(12,5)i,j,be1,be2,be3,bm1,bm2
        enddo
      enddo
C
C construct magnetic quantum numbers and catalogue states
C mj=mag. mom. degeneracy;  jj=ang. mom.;  mm= z-component  
      j=1
      ii(0)=0
      istop(0)=0
      do i=1,nst       
        mj(i)=nint(2.*spin(i)+1.)
        do m=0,mj(i)-1
          ii(j)=i
          jj(j)=spin(i)
          mm(j)=-spin(i)+dfloat(m)
          exx(j)=ex(i)
          j=j+1
        enddo
        istop(i)=istop(i-1)+mj(i)
      enddo  
      n=j-1
      if(n.gt.nstmax) then
       write(12,*)' WARNING - number of states is too big'
       write(12,*)' Increase parameter NSTMAX in file DWEIKO.DIM'
       write(12,*)' Try NSTMAX ~ (2*J_max+1)*NST
     & where J_max = largest spin state'
	 go to 2000
      endif  
C
C Calculates the deformations lengths for nuclear excitation 
C See Ref.  C.A. Bertulani et al., Phys.Rev. C53 (1996) 224 
C Eqs. (42) and (45)
      if(iopnuc.eq.1) then 
C
        write(12,6)
6       format(/,'  Deform. Par. DELTE0,',
     1         ' DELTE1 [fm],',' DELTE2 [fm] and ','DELTE2 [fm]',
     2     /,1x,
	3         ' ----------------------------------------------',
	4         '---------------------'
     5         ,/,1x,' Transit.',11x,'ALPHA0',6x,'DELTE1',6x,
     6         'DELTE2',6x,'DELTE3')
C
        do j=2,nst 
          call rsharps(10)
          read(10,*)idum,frd0(j),frd1(j),frd2(j),frd3(j)
        end do
        call deform(nst,r,delr,dnp,xr2x,delte)
        do j=2,nst 
          write(12,7)j,delte(0,j),delte(1,j),delte(2,j),delte(3,j)
7         format(1x,'1  --> ',i2,6X,4(f10.3,2x))    
        end do
C
        write(12,*)"   "
        write(12,*)"Parameters of the excited nucleus"
        write(12,*)"---------------------------------"
        write(12,'(A19,2x,I3,2x,I3,2x,I3)')
     & " at,zt,nt       : ",idnint(at),idnint(zt),idnint(at-zt)
        write(12,'(A,f12.3)')" <r^2>     [fm^2]:",xr2x 
        write(12,'(A,f12.3)')" sqrt[<r^2>] [fm]:",dsqrt(xr2x) 
C
      endif
C
      write(12,*)"  "
      write(12,'(a)')
     &"  Useful quantities"
      write(12,'(a)')
     &"  -----------------"
      write(12,'(a,f10.3)')
     &"    beta   (v/c)                        :",beta 
      write(12,'(a,f10.3)')
     &"    gamma  (Lorentz factor)             :",gamma
      write(12,'(a,f10.3)')
     &"    eta (Sommerfeld parameter)          :",eta
      write(12,'(a,f10.3)')
     &"    k (wavenumber)                [1/fm]:",pcmhc
      write(12,'(a,f10.3)')
     &"    Minimum impact parameter        [fm]:",bmin
      write(12,'(a,f10.3)')
     &"    Maximum Coulomb scatter. angle [deg]:",thet_bm
      write(12,'(a,i3)')
     &"    Number of coupled channels          :",n
C
      if(iopnuc.eq.1)then
C Nuclear potentials in the Bohr-Mottelson model
        call derivative(wn,wn1,wn2,delr,ngrid)
        do i=0,ngrid
	    ubm(0,i)=3.*wn(i)+r(i)*wn1(i)
	    ubm(1,i)=-1.5*dnp*(wn1(i)+r(i)*wn2(i)/3.)
	    ubm(2,i)=wn1(i)
	    ubm(3,i)=wn1(i)
        enddo
      endif
C     
C computation of psi tensor and nuclear excitation matrices
C notation: 1=E1, 2=E2, 3=E2, 4=M1, 5=M2
      fac=zp*e2/beta/hc 
      fd(1)=fac*sqrt(2.*pi/3.)
      fd(2)=-fac*sqrt(2.*pi/15.)/2. 
      fd(3)=fac*sqrt(2.*pi/105.)/3.
	fd(4)=fd(1)
	fd(5)=fd(2)
      do i=1,n-1
       do j=i+1,n 
        mu=nint(mm(j)-mm(i)) 
        do l=1,5
         psi(l,ib,i,j)=0.
         ll=l
         if(l.eq.4)ll=1 
         if(l.eq.5)ll=2
         bmat=mat(l,ii(i),ii(j))
         iv1=iv(nint(jj(j)+dfloat(ll)-mm(j)+1))
         three= threej(jj(j),dfloat(ll),jj(i),-mm(j),dfloat(mu),mm(i))
         aux=bmat*iv1*three
         do ib=1,nb
          psi(l,ib,i,j)=fd(ll)*aux/b(ib)**ll 
         enddo
        enddo
C nuclear excitation part - vibrational model
C notation for 'ln': 0=Monopole, 1=(Isoscalar)Dipole, 2=Quadrupole, 3=octupole
        do ln=0,3
         bmatn=delte(ln,ii(j))
         iv2=iv(nint(jj(j)-mm(j)))
         three= threej(jj(j),dfloat(ln),jj(i),-mm(j),dfloat(mu),mm(i))
         three0=threej(jj(j),dfloat(ln),jj(i),zero,zero,zero)
         if(iopnuc.eq.1)then
           auxn=-bmatn*sqrt((2.*jj(j)+1.)/
     &          (2.*jj(i)+1.)/4./pi)*iv2*three*three0*gamma
           do ir=0,ngrid
             psinuc(ln,ir,i,j)=psinuc(ln,ir,i,j)+auxn*ubm(ln,ir)
           enddo
         endif	 
        enddo
       enddo
      enddo 
C
      if(iogam.eq.1)then
        fln_out2="dweiko_stat.out"
        open(unit=15,file=fln_out2,status='unknown')
        write(15,110)
110     format(1H1,9x,
     &     49HThe angular distribution tensors STAT(J,KA,KAPPA))
        write(15,111)
111     format(7x,1HJ,5x,2HKA,4x,5HKAPPA,9x,9HReal STAT,11x,9HImag STAT)
      endif
C
C for angular distributions of gamma-rays
      if(iogam.eq.2)then 
        call rsharps(10)
        read(10,*)iff,igg,thmin,thmax,ntheta
      endif
      if (ntheta.gt.ngrid) then
        write(12,*) '******    STOP: NTHETA LARGER THAN NGRID '
        write(12,*)' Decrease value of NTHETA, or increase NGRID'
        stop
      endif
C
C ==========================
C start integration in time for each impact parameter 
C
C impact parameter loop      
      icc=0
      do 101 ib=1,nb
       e0=hc*beta*gamma/b(ib)
       bb=b(ib)
       do i=1,n
        do j=1,n
         xi(i,j)=(exx(i)-exx(j))/e0
         if(j.lt.nst)spinave(j,ib)=0.
         do l=1,5
          psitot(l,i,j)=psi(l,ib,i,j)
         enddo
        enddo 	   
       enddo
C initialize statistical tensors
       do ka=0,4,2
        do kappa=-ka,ka
         do j=2,nst
          stat(j,ka,kappa)=0.
         enddo
        enddo
       enddo          
C
C initialize amplitudes and average over initial polarization
       mj1=1 
1001   do k=1,n
        ca(k)=0.
        dca(k)=0.
       enddo
       ca(mj1)=cmplx(1.,0.)
C
C ************** integrate coupled-channels equations
C initial guess for the time mesh
       tmax=15. 
       nt=50
       dtau=2.*tmax/nt    
       t1=-tmax
       t2=tmax
       nvar=n 
       h1=dtau
       hmin=0.
C accuracy
       eps=accur
C Numerical Recipes routine for adaptive Runge-Kutta integration 
C
       call odeint(ca,nvar,t1,t2,eps,h1,hmin,nok,nbad,dcadt,rkqs)
C
C ************ end time integration
C
C excitation probabilities summed over final spin
       do j=1,nst
         sum=0.
         do i=istop(j-1)+1,istop(j)
           sum=sum+abs(ca(i))**2 
           if(j.eq.jinel) then
              mu=nint(mm(i)-mm(mj1)) 
              cmu(mu,ib)=ca(i)
           endif
         enddo
         pcc(j,ib)=sum*tb(ib)
       enddo
C
C check if norm is conserved 
       xn=0.
       do j=1,nst
        xn=xn+pcc(j,ib)
       enddo
       xn=xn+1.-tb(ib)
       eps10=eps*10.
       if(abs(xn-1.).gt.eps10) then
        write(12,*)' Error in probability exceeds 10 x ACCUR'
        go to 2000
       endif          
C
C sum excitation probabilities over initial spin
       do j=1,nst
        spinave(j,ib)=spinave(j,ib)+pcc(j,ib)/mj(1)
       enddo
C
C statistical tensors    
       if(iogam.ne.1)goto 113
       do ka=0,4,2
        do kappa=-ka,ka
         do i=1,nst
          jini=istop(ii(i-1))+1
          jfin=istop(ii(i))
          do j=jini,jfin
           amnkap=mm(j)+kappa
           iv3=iv(nint(abs(jj(j)+amnkap)))
           three= threej(jj(j),jj(j),dfloat(ka),-amnkap,
     &                             mm(j),dfloat(kappa))
           jkap=j+kappa
           if(jkap.ge.jini.and.jkap.le.jfin) then
            stat(i,ka,kappa)=stat(i,ka,kappa)+iv3*three*
     &            sqrt(2.*jj(j)+1.)*conjg(ca(j))*ca(jkap)/mj(1)
	     endif 
          enddo
         enddo
        enddo
       enddo
C
113    mj1=mj1+1
       if(mj1.le.mj(1))goto 1001 
C end sum over intial orientation
C
C output of statistical tensors for each impact parameter
       if(iogam.eq.1) then 
        fln_out2="dweiko_stat.out"
        open(unit=18,file=fln_out2,status='unknown')
        write(18,*)'b=',b(ib)
        do j=2,nst
         do ka=0,4,2
          do kappa=0,ka
           if(stat(j,ka,kappa).ne.0.) then
       	  write(18,114)j,ka,kappa,stat(j,ka,kappa)
           endif
114        format(1x,3I7,2E22.4)
          enddo
         enddo
        enddo
       endif
C
101   continue
C end impact parameter loop 
C =========================
C	        
      if(iob.eq.1) then
       write(12,116)
116    format(//,
     &   ' Impact parameter, state, and occupation probabilities',/)
       do ib=1,nb
        do j=2,nst
         if(j.eq.2) then
          write(12,118)b(ib),j,spinave(j,ib)
          else
          write(12,119)j,spinave(j,ib)
         endif
118      format(2x,g12.4,i2,2x,g12.4)
119      format(14x,i2,2x,g12.4)
        enddo
       enddo
      endif
C
C  inelastic scattering cross sections
      if(ioinel.ge.1) then
        mumax=nint(mm(istop(jinel))-mm(1))
        call inelast(b,cmu,mumax,phasen,ioinel,rhog,theta,ntheta) 
      endif
C
C calculates nuclear reaction cross section in mb
      do ib=1,nb
        pint(ib)=20.*pi*b(ib)*(1.-tb(ib))
      enddo      
      big=1.d30
      call spline(b,pint,nb,big,big,pint_2)
      call qsimp(b,pint,pint_2,nb,bmin1,bmax,sigr)
151   write(12,42)sigr
42    format(/,3x,' Total nuclear reaction cross section (mb) = ',
     &       F8.2)
C
C **** Compute total cross sections 
      write(12,160)
160   format(/,' Excitation cross sections in mb',/,
     1         ' -------------------------------')
C
      sum1=0.
      do k=2,nst
C interpolate and integrate - recoil correction included
        do ib=1,nb
         pint(ib)=20.*pi*b(ib)*spinave(k,ib)
        enddo      
        big=1.d30
        call spline(b,pint,nb,big,big,pint_2)
        if(bmin.eq.0.)bmin=bmin1
        call qsimp(b,pint,pint_2,nb,bmin,bmax,cross)     
        sigj(k)=cross
C output
        write(12,'('' State '',i3,''  = '','//
     &    'e12.3)')k,cross
        sum1 = sum1 + cross
      enddo
      write(12,'('' Total '',3x,''  = '','//
     &    'e12.3)')sum1
C **** end cross section 
C
C angular distribution of gamma-rays
      if(iogam.eq.2)
     &      call gamdis(iw,iff,igg,thmin,thmax,ntheta)
C
      goto 2000
1500  write(12,1501)
1501  format(26h0ERROR: NST EXCEEDS NSTMAX)
2000  stop
      end   
C--------------------------------------------------------------------
C                  END MAIN PROGRAM
C--------------------------------------------------------------------
C
C====================================================================
C
      subroutine dcadt(tau,ca,dca)  
C
C time derivative of the occupation amplitudes 
C          
      include 'dweiko.dim'
      implicit real*8(a-h,o-z)
      complex*16 ei,ca(nstmax),dca(nstmax),expt,vmat(nstmax,nstmax),
     &           vnuc
      real*8 xi(nstmax,nstmax),mm(nstmax)
      common /a1/ beta,gamma,pi,hc,gcm,betacm
      common/a3/xi
      common/a4/mm,n
      ei=cmplx(0.,1.)
      call vint(tau,vmat)      
      do i=1,n
       dca(i)=0.
       do j=1,n
        expt=exp(ei*tau*xi(j,i))
        dca(i)=dca(i)-ei*vmat(i,j)*expt*ca(j)
       enddo
      enddo
      return
      end
C
C=====================================================================
C
      subroutine vint(t,vmat)
C
C time-dependent fields for E1, E2, E3, M1 and M2 transitions
C vei = Q-functions in Ref. 
C C.A. Bertulani, Comp. Phys. Comm. 116 (1999) 345 
C
      implicit real*8(a-h,o-z)
      include 'dweiko.dim'
      complex*16 ei,q,ve1,ve2,ve3,vm1,vm2,vmat(nstmax,nstmax),Ylm,
     &                  psinuc(0:3,0:ngrid,nstmax,nstmax)
      real*8 mm(nstmax)
      real*8 psitot(5,nstmax,nstmax),xi(nstmax,nstmax)
      common /a1/ beta,gamma,pi,hc,gcm,betacm
      common/a3/xi
      common/a4/mm,n
      common/a5/psitot
      common/a6/psinuc
      common/a7/bb
      common/a8/iopnuc
      common/a9/delr
C useful constants
      ei=cmplx(0.,1.)
      t2=t*t
      phi=1./sqrt(1.+t2)
      phi3=phi**3
      phi5=phi**5
      phi7=phi**7            
      beta2=beta*beta
      gamma2=gamma*gamma
C loop in transitions state i --> j
      do i=1,n-1
       do j=i+1,n
        mu=nint(mm(j)-mm(i)) 
        vmat(i,j)=0.
        x=xi(j,i)
	  musign=1.
        if(mu.ne.0.)musign=mu/iabs(mu)
C electric dipole (E1) t.d. field
        if(mu.eq.0) then
         q=-sqrt(2.)/gamma*t*phi3
	  else
         q=-musign*phi3
        endif
        ve1=psitot(1,i,j)*q
C electric quadrupole (E2) t.d. field
        if(iabs(mu).eq.2) then
         q=3.*phi5
         else
         if(iabs(mu).eq.1) then
           q=musign*3.*gamma*(2.-beta2)*t*phi5
          else
           q=sqrt(6.)*(2.*t2-1.)*phi5
         endif
        endif
        ve2=psitot(2,i,j)*q      
C electric octupole (E3) t.d. field
        if(iabs(mu).eq.3) then
         q=-musign*15./2.*sqrt(3./2.)*phi7
         else
         if(iabs(mu).eq.2) then
	     q=-15./2.*gamma*(3.-beta2)*t*phi7
          else
	    if(iabs(mu).eq.1) then
	     q=-musign*gamma2*3./sqrt(40.)*(15.-11.*beta2)*(4.*t2-1.)*phi7
	     else
           q=-3.*sqrt(3./10.)*gamma*(5.-beta2)*t*(2.*t2-3.)*phi7
          endif
         endif
        endif
        ve3=psitot(3,i,j)*q      
C magnetic dipole (M1) t.d. field
        q=ei*iabs(mu)*beta*phi3
        vm1=psitot(4,i,j)*q
C magnetic quadrupole (M2) t.d. field
        if(iabs(mu).eq.2)q=-musign*3.*ei*beta*phi5
        if(iabs(mu).eq.1)q=-3.*ei*beta*gamma*t*phi5
        vm2=psitot(5,i,j)*q
C
C add E1, E2, E3, M1 and M2
        vmat(i,j)=ve1+ve2+ve3+vm1+vm2
        if(iopnuc.ne.1)goto 20
C
C add nuclear interaction
       rr=bb/phi
       ir=int(rr/delr)
       h=1.-(rr-delr*ir)/delr
       theta=t*phi
       zero=0.d0
       vnuc=0.
       if(ir+1.lt.ngrid.and.ir.gt.0.) then
        do ln=0,3
         if(abs(mu).le.ln)vnuc=vnuc+(h*psinuc(ln,ir,i,j)+
     &	   (1.-h)*psinuc(ln,ir+1,i,j))*Ylm(ln,mu,theta,zero)
        enddo
        vmat(i,j)=vmat(i,j)+vnuc
       endif
C
C interaction is hermitian
20      vmat(j,i)=conjg(vmat(i,j))
       enddo
       vmat(i,i)=0.
      enddo
      return
      end	          
C
C=====================================================================
C
      subroutine omp_ws(r,V0_ws,VI_ws,R_ws,RI_ws,d_ws,dI_ws,wn,iout)
C
C Woods-Saxon optical model potential. Output in 'dweiko_omp.out'.
C
      include 'dweiko.dim'
      implicit real*8(a-h,o-z)
      real*8 r(0:ngrid)
      complex*16 ei,wn(0:ngrid)
C
      ei=cmplx(0.,1.)
      do i=0,ngrid
        vdum1 = -V0_ws/(1.+dexp((r(i)-R_ws)/d_ws))
        vdum2 = -VI_ws/(1.+dexp((r(i)-RI_ws-1.)/dI_ws))
        wn(i) = vdum1+ei*vdum2
      end do     
C Out the OMP
      if(iout.gt.0)then
        open(unit=71,file="dweiko_omp.out",status="unknown")
        write(71,'(A4,1x,A12,2x,(A12,1x,A12))')
     &    "#  I","       r[fm]"," Re{U} [MeV]"," Im{U} [MeV]"
        do i=0,ngrid
          write(71,'(I4,1x,F12.5,2x,(e12.5,1x,e12.5))')i,r(i),wn(i)
        end do
        close(unit=71)
      endif
      return 
      end
C
C=====================================================================
C
      subroutine omp_read(r,wn,iout)
C
C reads optical potential from file 'optw.in'. Writes in 'dweiko_omp.out'. 
C
      include 'dweiko.dim'
      implicit real*8(a-h,o-z)
      real*8 r(0:ngrid),ru(0:ngrid),dru(0:ngrid)
      complex*16 ei,wn(0:ngrid),u(0:ngrid)
C
      open(unit=11,file='optw.in',status='unknown')
      ei=cmplx(0.,1.)
      read(11,*)nr
      if (nr.ge.ngrid) then
       write(12,*) '******    STOP: NR LARGER THAN NGRID '
       stop
      endif
      ru(0)=0.
      do i=1,nr
       read(11,*)ru(i),vdum1,vdum2
       u(i)=vdum1+ei*vdum2
       dru(i)=ru(i)-ru(i-1)
      enddo
      u(0)=u(1) 
      close(unit=11)
C interpolate optical potential
      ir=0
      ncount=0
  1   ir=ir+1
      if(ir.gt.nr-1)goto 11
      do i=ncount,ngrid
       wn(i)=0. 
       aux=ru(ir)-r(i)
       if(aux.ge.0.) then
        al=1.-aux/dru(ir) 
        wn(i)=al*u(ir)+(1.-al)*u(ir-1)
        ncount=ncount+1
        else
 	  goto 1
       endif
      enddo
C  
   11 continue
C
C Out the OMP
      if(iout.gt.0)then
        open(unit=71,file="dweiko_omp.out",status="unknown")
        write(71,'(A4,1x,A12,2x,(A12,1x,A12))')
     &  "#  I","       r[fm]"," Re{U} [MeV]"," Im{U} [MeV]"
        do i=0,ngrid
          write(71,'(I4,1x,F12.5,2x,(e12.5,1x,e12.5))')i,r(i),wn(i)
        end do
        close(unit=71)
      endif
C
      return
      end
C
C=====================================================================
C
      subroutine omp_den(r,rmin,rmax,wn,iout)
C
C optical model potential from nuclear densities and nucleon-nucleon
C scattering cross section and forward scatering phase (t-rho-rho model)
C M.S.Hussein, R.A.Rego and C.A. Bertulani, Phys. Rep. 201 (1991) 279 
C output in 'dweiko_omp.out'
C
      implicit real*8(a-h,o-z)
      include 'dweiko.dim'
      common/a0/ap,zp,at,zt
      common/a1/ beta,gamma,pi,hc,gcm,betacm
      common/a2/signn,alnn
      real*8 r(0:ngrid),ru(0:ngrid),dru(0:ngrid)
      real*8 dens1(0:ngrid),dens2(0:ngrid)
      complex*16 ei,fact,wn(0:ngrid),u(0:ngrid)
      ei=cmplx(0.,1.)
C compute liquid-drop densities
      delru = (rmax-rmin)/ngrid
      do ir=0,ngrid
        ru(ir) = rmin + ir*delru
      end do
      do ir=1,ngrid
        dru(ir) = ru(ir)-ru(ir-1)
      end do
      do ir=0,ngrid
       dens1(ir)=rhonp(ru(ir),ap,zp)+rhopp(ru(ir),ap,zp)
       dens2(ir)=rhonp(ru(ir),at,zt)+rhopp(ru(ir),at,zt)
      enddo
C
      fact=-(alnn+ei)*signn*hc*beta/2.
      call twofold(fact,ru,rmax,delru,dens1,dens2,u)
C interpolate u 
      ir=0
      ncount=0
30    ir=ir+1
      if(ir.gt.ngrid-1)goto 31
      do i=ncount,ngrid
       wn(i)=0. 
       aux=ru(ir)-r(i)
       if(aux.ge.0.) then
        al=1.-aux/dru(ir) 
        wn(i)=al*u(ir)+(1.-al)*u(ir-1)
        ncount=ncount+1
        else
        goto 30
       endif
      enddo
31    continue
C Out the OMP
      if(iout.gt.0)then
        open(unit=71,file="dweiko_omp.out",status="unknown")
        write(71,'(A4,1x,A12,2x,(A12,1x,A12))')
     &  "#  I","       r[fm]"," Re{U} [MeV]"," Im{U} [MeV]"
        do i=0,ngrid
          write(71,'(I4,1x,F12.5,2x,(e12.5,1x,e12.5))')i,r(i),wn(i)
        end do
        close(unit=71)
      endif      
      return
      end
C
C=====================================================================
C
      subroutine omp_M3Y(r,rmin,rmax,wn,Wrat,iout)
C
C M3Y optical model potential 
C G.F. Bertsch et al., Nucl. Phys. A284, 399 (1977)
C A.M. Kobos et al., Nucl. Phys. A425, 205 (1984)
C output in 'dweiko_omp.out'
C
      implicit real*8(a-h,o-z)
      include 'dweiko.dim'
      common/a0/ap,zp,at,zt
      common/a1/ beta,gamma,pi,hc,gcm,betacm
      real*8 r(0:ngrid),ru(0:ngrid),dru(0:ngrid),t(0:ngrid)
      real*8 dens1(0:ngrid),dens2(0:ngrid),td3(0:ngrid)
      complex*16 ei,fact,wn(0:ngrid),u(0:ngrid),td2(0:ngrid)         	
C compute liquid-drop densities
      ei=cmplx(0.,1.)
      delru = (rmax-rmin)/ngrid
      do ir=0,ngrid
        ru(ir) = rmin + ir*delru
      end do
      do ir=1,ngrid
        dru(ir) = ru(ir)-ru(ir-1)
      end do
      do ir=0,ngrid
       dens1(ir)=rhonp(ru(ir),ap,zp)+rhopp(ru(ir),ap,zp)
       dens2(ir)=rhonp(ru(ir),at,zt)+rhopp(ru(ir),at,zt)
      enddo
C M3Y on a grid
      do ir=0,ngrid
       s=ru(ir)
       t(ir)=tM3Y(s)
      enddo
C  get folding of tM3Y and dens2
      fact=cmplx(1.d0,0.d0)
      call twofold(fact,ru,rmax,delru,t,dens2,td2)
      C3Y=-276.
      do ir=0,ngrid
       td2(ir)=td2(ir)+C3Y*dens2(ir)
      enddo
C  fold with the other density
      do ir=0,ngrid
       td3(ir)=dble(td2(ir))
      enddo
      fact=cmplx(1.d0,Wrat) 
      call twofold(fact,ru,rmax,delru,td3,dens1,u)
C interpolate u 
      ir=0
      ncount=0
30    ir=ir+1
      if(ir.gt.ngrid-1)goto 31
      do i=ncount,ngrid
       wn(i)=0. 
       aux=ru(ir)-r(i)
       if(aux.ge.0.) then
        al=1.-aux/dru(ir) 
        wn(i)=al*u(ir)+(1.-al)*u(ir-1)
        ncount=ncount+1
        else
        goto 30
       endif
      enddo
31    continue
C Out the OMP
      if(iout.gt.0)then
        open(unit=71,file="dweiko_omp.out",status="unknown")
        write(71,'(A4,1x,A12,2x,(A12,1x,A12))')
     &  "#  I","       r[fm]"," Re{U} [MeV]"," Im{U} [MeV]"
        do i=0,ngrid
          write(71,'(I4,1x,F12.5,2x,(e12.5,1x,e12.5))')i,r(i),wn(i)
        end do
        close(unit=71)
      endif      
      return
      end
C
C=====================================================================
C
      real*8 FUNCTION tM3Y(s)
      implicit real*8 (a-h,o-z)
C
C    M3Y effective interation - delta-function treated separately
C
      A3Y=7999.
      B3Y=-2134.
      beta1=4.
      beta2=2.5
C
      s1=s*beta1
      s2=s*beta2
      tM3Y= A3Y*exp(-s1)/s1+B3Y*exp(-s2)/s2
      return
      end
C
C=====================================================================
C
      subroutine deform(nst,r,delr,dnp,xr2x,delte)
C
C deformation parameters for nuclear excitation in the collective
C (vibrational, or Bohr-Mottelson) model.
C G.R. Satchler, Nucl. Phys. A472, 215 (1987) 
C
      implicit real*8 (a-h,o-z)
      include 'dweiko.dim'
      real*8 delte(0:3,nmax),ex(nmax),dens1(0:ngrid),
     & mat(5,0:nmax,0:nmax),dens2(0:ngrid),r(0:ngrid),aux1(0:ngrid),
     & sigj(nmax),simp(0:ngrid),frd0(nmax),frd1(nmax),
     & frd2(nmax),frd3(nmax)
      complex*16 ei
      common/a0/ap,zp,at,zt
      common/a1/ beta,gamma,pi,hc,gcm,betacm
      common/a11/mat,ex,sigj
      common/a13/frd0,frd1,frd2,frd3
      common/simps/ simp
      amu = 931.49432D0
      ei=cmplx(0.,1.)
C average value of r^2 and differences in neutron and proton radii
      do j=2,nst 
        do i=0,ngrid
          dens1(i) = rhonp(r(i),at,zt) + rhopp(r(i),at,zt)
          aux1(i) = 4.*pi*dens1(i)*(r(i)**2) * (r(i)**2)
        end do
        xr2x= 0.
        do i=0,ngrid
          xr2x = xr2x + simp(i)*delr*aux1(i)
        enddo
        xr2x = xr2x/at
C difference between neutron and proton radii
        do i=0,ngrid
          aux1(i) = 4.*pi*rhonp(r(i),at,zt)*(r(i)**2) * (r(i)**2)
        end do
        xrn2x= 0.
        do i=0,ngrid
          xrn2x = xrn2x + simp(i)*delr*aux1(i)
        enddo
        xrn2x = xrn2x/(at-zt)
        do i=0,ngrid
          aux1(i) = 4.*pi*rhopp(r(i),at,zt)*(r(i)**2) * (r(i)**2)
        end do
        xrp2x= 0.
        do i=0,ngrid
          xrp2x = xrp2x + simp(i)*delr*aux1(i)
        enddo
        xrp2x = xrp2x/zt
C neutron skin: to be used later in nuclear dipole potential
        dnp=(sqrt(xrn2x)-sqrt(xrp2x))/sqrt(xr2x)
C
C monopole deformation parameter
        xdel02 = 2.*pi*(hc**2)/(amu*xr2x)*1./(at*ex(j))
        delte(0,j) = dsqrt(frd0(j)*xdel02)
C dipole deformation parameter
        xdel12 = pi/2.*(hc**2)/amu*at/((at-zt)*zt)*(1./ex(j))
        delte(1,j) = dsqrt(frd1(j)*xdel12)
C quadrupole deformation parameter
        xdel22 = (2.*pi/3.)*((hc**2)/amu)*10.*(1./(at*ex(j)))
        delte(2,j) = dsqrt(frd2(j)*xdel22)
C octupole deformation parameter
        xdel23 = (2.*pi/3.)*((hc**2)/amu)*21.*(1./(at*ex(j)))
        delte(3,j) = dsqrt(frd3(j)*xdel23)
      enddo
      return
      end
C
C=====================================================================
C
      real*8 FUNCTION rhopp(rg,ap,zp)
C
C    Droplet densities for protons -  See  Myers and Swiatecki, Ann. of
C                                     Phys. A55 (1969) 395 and
C                                     Ann. of Phys. A84 (1974) 186
      implicit real*8 (a-h,o-z)
	pi=3.141597265
	tpp = 2.4
	bpp = 0.413 * tpp
	deltap = ((ap-2.*zp)/ap + 8.076d-3*zp*ap**(-.66666))
     #         / ( 1. + 4.871 * ap**(-.33333) )
	epsp = -0.1724*ap**(-.33333) + 3.051d-3 * zp**2 * ap**(-1.33333)
     #       + .4166666 * deltap**2
	rp = 1.18 * ap**0.33333 * ( 1. + epsp )
	dp = 0.66666 * ((ap-2.*zp)/ap - deltap ) * rp
	rpp = rp - ((ap-zp)/ap) * dp
	cpp = rpp * ( 1. - (bpp/rpp)**2 )
	rho0pp = 3. * zp / ( 12.5664 * cpp**3 *
     #         ( 1. + pi**2 * tpp**2 / (19.36 * cpp**2)))
	rhopp = rho0pp / ( 1. + exp((rg-cpp)/(tpp/4.4)) )
      return
      end
C
C=====================================================================
C
      real*8 FUNCTION rhonp(rg,ap,zp)
C
C    Droplet densities for neutrons - From Myers and Swiatecki, Ann. of
C                                     Phys. A55 (1969) 395 and
C                                     Ann. of Phys. A84 (1974) 186
      implicit real*8 (a-h,o-z)
	pi=3.141597265
	tnp = 2.4
	bnp = 0.413 * tnp
	deltap = ((ap-2.*zp)/ap + 8.076d-3*zp*ap**(-.66666))
     #         / ( 1. + 4.871 * ap**(-.33333) )
	epsp = -0.1724*ap**(-.33333) + 3.051d-3 * zp**2 * ap**(-1.33333)
     #       + .4166666 * deltap**2
	rp = 1.18 * ap**0.33333 * ( 1. + epsp )
	dp = 0.66666 * ((ap-2.*zp)/ap - deltap ) * rp
	rnp = rp + (zp/ap) * dp
	cnp = rnp * ( 1. - (bnp/rnp)**2 )
	rho0np = 3. * ( ap - zp )/ ( 12.5664 * cnp**3 *
     #         ( 1. + pi**2 * tnp**2 / (19.36 * cnp**2)))
	rhonp = rho0np / ( 1. + exp((rg-cnp)/(tnp/4.4)) )
      return
      end
C
C=====================================================================
C
      SUBROUTINE phnuc(wn,b,nb,z,phasen)
C
C  nuclear eikonal phase
C
      implicit real*8 (a-h,o-z)
      include 'dweiko.dim' 
      real*8 b(nbmax),z(0:ngrid), simp(0:ngrid)
      complex*16 ei,phasen(nbmax),vnucl,wn(0:ngrid)
      common /a1/ beta,gamma,pi,hc,gcm,betacm
      common/a9/ delr
      common/simps/ simp
      ei=cmplx(0.,1.)
      do k=1,nb
	phasen(k)=0.
C  integrate in z coordinate - Simpson's rule
	do i=0,ngrid
       vnuc=0.
	 rv=sqrt(b(k)**2+z(i)**2)
	 l=int(rv/delr)
	 al=1.-(rv-delr*l)/delr
	 if(l+1.lt.ngrid.and.l.gt.0)vnucl=al*wn(l)+(1.-al)*wn(l+1)
	 phasen(k)=phasen(k)+simp(i)*vnucl
	enddo
	phasen(k)=-2./hc/beta*delr*phasen(k)
      enddo
      return
      end
C
C=====================================================================
C
      subroutine gamdis(iw,iff,igg,thmin,thmax,nth)
C
C calculates angular distribution of gamma rays
C
      include 'dweiko.dim'
      implicit real*8(a-h,o-z)
      integer*4 mj(nstmax),ii(0:nstmax),istop(0:nmax)
      real*8 mm(nstmax),jj(nstmax),mat(5,0:nmax,0:nmax),
     & ex(nmax),simp(0:ngrid),gamd(0:ngrid),sigj(nmax)
      real*8 iv(0:130),fak(0:130),fad(0:130)
C
      character*14 fln_out
C
      common /gf1 / iv
      common /gf2 / fak
      common /gf3 / fad
      common/a1/ beta,gamma,pi,hc,gcm,betacm
      common/a4/mm,n
      common/a10/jj,istop
      common/a11/mat,ex,sigj
      common/a12/ii
      common/simps/ simp
      call gfv(30)
C
      pi=acos(-1.)
      zero=1.d0
      one=1.d0
      tr=pi/180.
      ajf=jj(iff)
      ajg=jj(igg)
      ntheta=ngrid-1
      dt0=180./ntheta
      itheta=0
      ifmin=istop(ii(iff-1))+1
      ifmax=istop(ii(iff))
      i1max=istop(ii(1))
C
      fln_out="dweiko_gam.out"
      open(unit=65,file=fln_out,status='unknown')
      write(65,*)'  Gamma-ray distribution (laboratory system)'
      write(65,*)'    theta       W(theta) [1/sr]'
C
      do theta=0.,180.,dt0
        sumka=0.
        do ka=0,4,2	   
         sumi=0.
         do i1=1,i1max
          sumf=0.
          do iif=ifmin,ifmax
           amf=mm(iif)
           iv1=iv(nint(abs(amf)))
           three1= threej(ajf,ajf,dfloat(ka),amf,-amf,zero)
           suml=0.
	     do l=1,5
            ll=l
            if(l.eq.4)ll=1
            if(l.eq.5)ll=2
            bmatl=mat(l,iff,igg)
            isl=ll+1 
            if(l.lt.4)isl=ll
            ivsl=iv(isl)             
            redl=ivsl*sqrt(8.*pi*(ll+1)/ll/(fad(2*ll+1))**2*
     &           ((ex(iff)-ex(igg))/hc)**(2*ll+1)/(2*ajg+1))
     &           *bmatl
            sumlp=0.
            itest=min(iabs(l-ka),5)
            if(itest.lt.1)itest=1		   		  	 
            do lp=itest,5
             llp=lp
             if(lp.eq.4)llp=1
             if(lp.eq.5)llp=2
             islp=llp+1 
             if(lp.lt.4)islp=llp
             ivslp=iv(islp)             
             three2=threej(dfloat(ll),dfloat(llp),dfloat(ka),
     &            			one,-one,zero)
             rac=racah(dfloat(ll),dfloat(llp),dfloat(ka),
     &             ajf,ajf,ajg)
             part1=sqrt((2.*ll+1.)*(2.*llp+1.)*(2.*ajf+1.)*(2.*ka+1.))
             part2=iv(nint(abs(ajf-ajg-1.)))
		   fka=part1*part2*three2*rac
             bmatlp=mat(lp,iff,igg)
             redlp=ivslp*sqrt(8.*pi*(llp+1)/llp/(fad(2*llp+1))**2*
     &           ((ex(iff)-ex(igg))/hc)**(2*llp+1)/(2*ajg+1))
     &           *bmatlp
             sumlp=sumlp+fka*redlp
            enddo
            suml=suml+sumlp*redl
           enddo
           sumf=sumf+suml	 
          enddo          
          sumi=sumi+sumf 
         enddo
         call legendc(ka,theta,pleg)
         sumka=sumka+sqrt(2.*ka+1.)*pleg*sigj(iff)
        enddo
        gamd(itheta)=sumka        
        sumt=sumka+sin(theta*tr)*simp(itheta)*sumka 	  	  	    	    
        itheta=itheta+1
      enddo       
      sumt=sumt*2.*pi*(dt0*tr)
C normalize angular distribution
      do itheta=0,ntheta
        gamd(itheta)=gamd(itheta)/sumt
      enddo 
C interpolate and write output
      ntheta=nth
      if(thmin.eq.0.)thmin=0.01
      dt=(thmax-thmin)/ntheta
      do theta=thmin,thmax,dt
        itheta=theta/dt0
        theta0=itheta*dt0
        al=1.-(theta-theta0)/dt0 
        gam=al*gamd(itheta)+(1.-al)*gamd(itheta+1)
C convert to laboratory system (if iw=0, projectile excitation)
        fac=1.
        thetal=theta
        if(iw.eq.0)fac=gcm**2*(1.+betacm*cos(theta*tr))**2
        gam=fac*gam
        fac=sin(theta*tr)/gcm/(cos(theta*tr)+betacm)
        if(iw.eq.0)thetal=atan(fac)/tr 
        write(65,162)thetal,gam
      enddo
162   format(2x,g12.4,2x,g12.4)
      return
      end
C
C=====================================================================
C
      subroutine signne(e0,sigpp,sigpn)
C
C free nucleon-nucleon cross sections (in mb)
C 
      implicit real*8(a-h,o-z)
      e=e0
      if(e.gt.650.) e=650.
      ei=1./e
      if(e.gt.111.511) then
        sigpn=28.1451+ei*(1431.74+ei*297627.)
       else
        if(e.gt.38.22) then
          sigpn=4.923+ei*(5673.25+ei*113412.)
         else
          sigpn=-1.86+e*(0.09415+e*1.306e-4)
          sigpn=3./(1.206*e+sigpn*sigpn)
          sigpn=sigpn+1./(1.206*e+(.4233+e*0.13)**2)
          sigpn=3141.59*sigpn
         endif
       endif
      if(e.gt.158.555) then
        sigpp=11.7386+0.0189*e+1362.11*ei
       else
        if(e.gt.42.738) then
          sigpp=17.8465+ei*(454.414+ei*65760.5)
         else
          if(e.gt.10.54) then
            sigpp=-2.0331+ei*(2690.66+ei*6498.86)
           else
            ei=1./(e+0.62)
            sigpp=-229.51+ei*(6596.46-ei*3920.61)
            sigpp=e*ei*sigpp
           endif
         endif
       endif
      return
      end
C
C
C=====================================================================
C
      subroutine phnne(e0,phpp,phpn,ahpp,ahpn)
C
C phase of nucleon-nucleon scattering amplitude
C from table 1 of Hussein, Rego and Bertulani, Phys. Rep. 201 (1991) 279
C 
      implicit real*8(a-h,o-z)
	real*8 e(10),php(10),phn(10),ahp(10),ahn(10),y2a(10)
	data e/100.,150.,200.,325.,425.,550.,650.,800.,1000.,2200./
	data php/1.87,1.53,1.15,0.45,0.47,0.32,0.16,0.06,-0.09,-0.17/
	data phn/1.,0.96,0.71,0.16,0.25,-0.24,-0.35,-0.2,-0.46,-0.5/
	data ahp/0.66,0.57,0.56,0.26,0.21,0.04,0.07,0.09,0.09,0.12/
	data ahn/0.36,0.58,0.68,0.36,0.27,0.085,0.09,0.12,0.12,0.14/
C  interpolate
      big=1.d30
      call spline(e,php,10,big,big,y2a)
      call splint(e,php,y2a,10,e0,phpp)      
C
      call spline(e,phn,10,big,big,y2a)
      call splint(e,phn,y2a,10,e0,phpn)      
C
      call spline(e,ahn,10,big,big,y2a)
      call splint(e,ahn,y2a,10,e0,ahpn)      
C
      call spline(e,ahp,10,big,big,y2a)
      call splint(e,ahp,y2a,10,e0,ahpp)      
C
      return
      end
C
C=====================================================================
C
      SUBROUTINE twofold(fact,r,rmax,delr,dens1,dens2,v)
C
C  Folding of densities
C
      implicit real*8 (a-h,o-z)
      include 'dweiko.dim'
      real*8 dens1(0:ngrid),dens2(0:ngrid),r(0:ngrid),x(0:ngrid),
     &       rint(0:ngrid), simp(0:ngrid)
      complex*16 v(0:ngrid),fact
      common/simps/ simp
C
      pi=3.141597265
      delx=2./ngrid
      do i=0,ngrid
       x(i)=-1.+i*delx
       rint(i)=r(i)
      enddo
      do i=0,ngrid
        sum2=0.
        do j=0,ngrid
	    sum1=0.
	    do k=0,ngrid
	      aux=sqrt(abs(r(i)**2+rint(j)**2-2.*r(i)*rint(j)*x(k)))
	      ifrac=int(aux/delr)
	      al=1.-(aux-delr*ifrac)/delr
	      if(ifrac+1.lt.ngrid.and.ifrac.gt.0)
     &            sum1=sum1+simp(k)*(al*dens1(ifrac)+
     &                  (1.-al)*dens1(ifrac+1))
	    enddo
	    ifrac=int(rint(j)/delr+1.)
	    al=1.-(rint(j)-delr*ifrac)/delr
	    if(ifrac+1.lt.ngrid.and.ifrac.gt.0)
     &           sum2=sum2+simp(j)*rint(j)**2*(al*dens2(ifrac)
     &             +(1.-al)*dens2(ifrac+1))*delx*sum1
	  enddo
	  v(i)=2.*pi*fact*sum2*delr
      enddo
      return
      end
C
C=====================================================================
C      
      subroutine phnucF(pcm_nn,eca,b,r,delr,chi)
C
C  returns eikonal phase in 't-rho-rho' model
C
      implicit real*8 (a-h,o-z)
      include 'dweiko.dim'
      parameter(nbig=5000)	  
      real*8 r(0:ngrid),q(0:nbig),b(nbmax),rho1(0:ngrid),
     &       simp(0:ngrid),fq1(0:nbig),fq2(0:nbig),rho2(0:ngrid)		    
      complex*16 chi(nbmax),sum,fn,ei 
      common/a0/ap,zp,at,zt
      common/a14/bmin1,bmax1,eta,pcmhc			    
      pi=acos(-1.)
      ei=cmplx(0.,1.)
      delq=50./bmax1/nbig
      do iq=0,nbig
        q(iq)=iq*delq
      enddo
      do ir=0,ngrid
       rho1(ir)=rhonp(r(ir),ap,zp)+rhopp(r(ir),ap,zp)
       rho2(ir)=rhonp(r(ir),at,zt)+rhopp(r(ir),at,zt)
      enddo
      call  fourier(r,delr,rho1,nbig,q,fq1)
      call  fourier(r,delr,rho2,nbig,q,fq2) 
      call signne(eca,sigpp,sigpn)
      signn = 0.1 * (zp*zt+(ap-zp)*(at-zt))/ap/at * sigpp +
     &	    0.1 * (zp*(at-zt)+zt*(ap-zp))/ap/at * sigpn
      call phnne(eca,phpp,phpn,ahpp,ahpn)
      alnn = (zp*zt+(ap-zp)*(at-zt))/ap/at * phpp +
     &	    (zp*(at-zt)+zt*(ap-zp))/ap/at * phpn
      xinn = (zp*zt+(ap-zp)*(at-zt))/ap/at * ahpp +
     &	    (zp*(at-zt)+zt*(ap-zp))/ap/at * ahpn
      do ib=1,nbmax	      
        sum=0.
        ig=4
        do iq=0,nbig 
          xinn=0.
          fn=pcm_nn/4/pi*signn*(ei+alnn)*exp(-xinn*q(iq)*q(iq))
          sum=sum+ig*q(iq)*
     &    		fq1(iq)*fq2(iq)*fn*bessj0(q(iq)*b(ib))
         ig=6-ig
         enddo
         chi(ib)=sum*delq/3.
      enddo
      return
      end
C
C=====================================================================
C      
      subroutine fourier(r,delr,fr,nbig,q,fq)
C
C  returns Fourier transform of densities
C
      implicit real*8 (a-h,o-z)
      include 'dweiko.dim'  
      real*8 r(0:ngrid),fr(0:ngrid),q(0:nbig),fq(0:nbig),
     &       simp(0:ngrid)		    
      common/simps/ simp
      pi=acos(-1.)
      do iq=0,nbig  
        sum=0.
        do ir=1,ngrid	      
          if(abs(q(iq)).le.1.d-20)go to 25
          fac=r(ir)*sin(r(ir)*q(iq))/q(iq)*fr(ir)
	    go to 26
25        continue
          fac=r(ir)**2*fr(ir)
26        continue
          sum=sum+simp(ir)*fac
         enddo
         fq(iq)=4.*pi*sum*delr
      enddo
      return
      end
C
C=====================================================================
C      
      subroutine elast(b,delb1,phasen,ioelas,rhog,theta,ntheta)
C
C  elastic scatering cross section / Rutherford cross section
C  theta in degrees
C
      implicit real*8 (a-h,o-z)
      include 'dweiko.dim'
      parameter (nbig=10000)
C NOTE: keep 'nbig' large (~ 10000) so that scattering of large
C charged particles is accurate.
      integer ig(0:nbig)		  
      real*8 b(nbmax),theta(ntheta)
      complex*16 ei,phasen(0:ngrid),pint,sumc,amp(ntheta),fc(ntheta),
     &           phasni(nbig)
      common/a1/ beta,gamma,pi,hc,gcm,betacm
      common/a14/bmin1,bmax1,eta,pcmhc,b13
      open(unit=81,file="dweiko_elast.out",status="unknown")
      pi=acos(-1.)
      ei=cmplx(0.,1.)
      tr=pi/180.
C initialize Coulomb amplitude suboutine
      call fcoul(ei,theta,ntheta,fc) 
C interpolate nuclear phase
      bstep=bmax1/nbig
      ig(0)=4 
      do k=1,nbig        
        phasni(k)=0.
        bi=k*bstep
        if(bi.le.(bmin1+delb1)) then
          phasni(k)=phasen(1)
         else   
          ir=(bi-bmin1)/delb1
          al=1.-(bi-(ir*delb1+bmin1))/delb1
          phasni(k)=al*phasen(ir)+(1.-al)*phasen(ir+1)
        endif
        ig(k)=6-ig(k-1)
      enddo
21    ig(nbig)=1
C        	  	   
C elastic amplitudes 
      do i=1,ntheta
        th=tr*theta(i)/2.
        qt=2.*pcmhc*sin(th)
C impact parameter integration
        amp(i)=0.	      
        sumc=0.
        do k=1,nbig
         bi=k*bstep
         phasec=2.*eta*dlog(pcmhc*bi)
         dum1=qt*bi
         pint=bi*bessj0(dum1)*exp(ei*phasec)
     &             *(1.-exp(ei*phasni(k)))
         sumc=sumc+ig(k)*pint
        enddo      
        amp(i)=ei*pcmhc*bstep*sumc/3.
      enddo  
C  end integration
150   format(2x,g8.3,5x,g8.3,5x,g8.3)
      if(ioelas.eq.2)goto 200
C  elastic scattering cross section in the center of mass
      write(81,'(''   Angle       Cross sect.   Ratio'')')
      write(81,'('' (degrees)      (mb/sr) '')')
      do i=1,ntheta
        sigrut=abs(fc(i))**2  
        sigma=abs(amp(i)+fc(i))**2/sigrut
        write(81,150)theta(i),sigrut,sigma
      enddo
      goto 300
C  elastic scattering cross section in the laboratory
200   write(81,'(''   Angle       Cross sect.   Ratio'')')
      write(81,'('' (degrees)      (mb/sr) '')')
      do i=1,ntheta
        fac=(gcm**2*(cos(theta(i)*tr)+rhog)**2+
     &		  sin(theta(i)*tr)**2)**1.5/
     &    	  gcm/abs(1.+rhog*cos(theta(i)*tr))
        sigrut=abs(fc(i))**2  
        sigrut=fac*sigrut
        sigma=fac*abs(amp(i)+fc(i))**2/sigrut
        fac=sin(theta(i)*tr)/gcm/(cos(theta(i)*tr)+rhog)
        theta(i)=atan(fac)/tr 
        if(i.gt.2.and.theta(i).lt.theta(i-1))goto 300
        write(81,150)theta(i),sigrut,sigma
      enddo   
300   return 
      end 
C
C
C=====================================================================
C      
      subroutine elastp(b,delb1,phasen,phases,ioelas,rhog,theta,ntheta)
C
C  elastic scatering cross section for protons / Rutherford cross section
C  theta in degrees
C
      implicit real*8 (a-h,o-z)
      include 'dweiko.dim'
      parameter (nbig=3000)
      integer ig(0:nbig)		  
      real*8 b(nbmax),theta(ntheta)
      complex*16 ei,phasen(0:ngrid),pint1,sumc1,amp1(ntheta),fc(ntheta),
     &           phasni(nbig),phases(0:ngrid),phassi(nbig),sum2,pint2,
     &           amp2(ntheta)
      common/a1/ beta,gamma,pi,hc,gcm,betacm
      common/a14/bmin1,bmax1,eta,pcmhc,b13
      open(unit=81,file="dweiko_elast.out",status="unknown")
      pi=acos(-1.)
      ei=cmplx(0.,1.)
      tr=pi/180.
C initialize Coulomb amplitude suboutine
      call fcoul(ei,theta,ntheta,fc) 
C interpolate nuclear phase
      bstep=bmax1/nbig
      ig(0)=4 
      do k=1,nbig        
        phasni(k)=0.
        phassi(k)=0.
        bi=k*bstep
        if(bi.le.(bmin1+delb1)) then
          phasni(k)=phasen(1)
          phassi(k)=phases(1)
         else   
          ir=(bi-bmin1)/delb1
          al=1.-(bi-(ir*delb1+bmin1))/delb1
          phasni(k)=al*phasen(ir)+(1.-al)*phasen(ir+1)
          phassi(k)=al*phases(ir)+(1.-al)*phases(ir+1)
        endif
        ig(k)=6-ig(k-1)
      enddo
21    ig(nbig)=1
C        	  	   
C elastic amplitudes 
      do i=1,ntheta
        th=tr*theta(i)/2.
        qt=2.*pcmhc*sin(th)
C impact parameter integration
        amp1(i)=0.
        amp2(i)=0.	  	      
        sumc1=0.
        sumc2=0.
        do k=1,nbig
         bi=k*bstep
         phasec=2.*eta*dlog(pcmhc*bi)
         dum1=qt*bi
         pint1=bi*bessj0(dum1)*exp(ei*phasec)
     &          *(1.-exp(ei*phasni(k)))*cos(dble(pcmhc*bi*phassi(k)))
         pint2=bi*bessj1(dum1)*exp(ei*(phasec+phasni(k)))*
     & 	   sin(dble(pcmhc*bi*phassi(k)))
         sumc1=sumc1+ig(k)*pint1
         sumc2=sumc2+ig(k)*pint2 
        enddo      
        amp1(i)=ei*pcmhc*bstep*sumc1/3.
        amp2(i)=ei*pcmhc*bstep*sumc2/3.
      enddo  
C  end integration
150   format(2x,g8.3,5x,g8.3,5x,g8.3)
      if(ioelas.eq.2)goto 200
C  elastic scattering cross section in the center of mass
      write(81,'(''   Angle       Cross Sect.   Ratio'')')
      write(81,'('' (degrees)      (mb/sr) '')')
      do i=1,ntheta
        sigrut=abs(fc(i))**2  
        sigma=(abs(amp1(i)+fc(i))**2+abs(amp2(i))**2)
        sigmar=sigma/sigrut        
        write(81,150)theta(i),sigma,sigmar
      enddo
      goto 300
C  elastic scattering cross section in the laboratory
200   write(81,'(''   Angle       Cross sect.   Ratio'')')
      write(81,'('' (degrees)      (mb/sr) '')')
      do i=1,ntheta
        fac=(gcm**2*(cos(theta(i)*tr)+rhog)**2+
     &		  sin(theta(i)*tr)**2)**1.5/
     &    	  gcm/abs(1.+rhog*cos(theta(i)*tr))
        sigrut=abs(fc(i))**2  
        sigrut=fac*sigrut
        sigma=(abs(amp1(i)+fc(i))**2+abs(amp2(i))**2)
        sigmar=sigma/sigrut
        fac=sin(theta(i)*tr)/gcm/(cos(theta(i)*tr)+rhog)
        theta(i)=atan(fac)/tr 
        if(i.gt.2.and.theta(i).lt.theta(i-1))goto 300
        write(81,150)theta(i),sigma,sigmar
      enddo   
300   return 
      end 
C
C=====================================================================
C      
      subroutine 
     & inelast(b,cmu,mumax,phasen,ioinel,rhog,theta,ntheta)
C
C  inelastic scatering cross section
C  theta in degrees
C
      implicit real*8 (a-h,o-z)
      include 'dweiko.dim'
      parameter (nbig=5000)
      integer ig(0:nbig),ih(0:nbig)		  
      real*8 b(nbmax),theta(ntheta),spin(nmax)
      complex*16 ei,phasen(0:ngrid),pint,sumc,amp(ntheta),fc(ntheta),
     &      phasni(nbig),cmu(-20:20,nbmax),sumn,finel(-20:20,ntheta),
     &      cmuni(-20:20,nbig)
      common/a1/ beta,gamma,pi,hc,gcm,betacm
      common/a14/bmin1,bmax1,eta,pcmhc,b13
      common/a15/bmax,delb1,delb2,spin			    
      common/a16/nb1			    			    
      open(82,file='dweiko_inel.out', status='unknown')
      pi=acos(-1.)
      ei=cmplx(0.,1.)
      tr=pi/180.
C interpolate nuclear phase
      bstep=bmax/nbig
      ig(0)=4 
      do k=1,nbig        
        phasni(k)=0.
        bi=k*bstep
        if(bi.le.(bmin1+delb1))phasni(k)=phasen(1)
        if(bi.lt.bmax1.and.bi.gt.bmin1+delb1) then   
          ir=(bi-bmin1)/delb1
          al=1.-(bi-(ir*delb1+bmin1))/delb1
          phasni(k)=al*phasen(ir)+(1.-al)*phasen(ir+1)
        endif
        if(bi.gt.bmax1)phasni(k)=0.
        ig(k)=6-ig(k-1)
      enddo
      ig(1)=1
      ig(nbig)=1	 
C interpolate nuclear excitation amplitudes
      do mu=-mumax,mumax
        do k=1,nbig        
          cmuni(mu,k)=0.
          bi=k*bstep
          if(bi.le.(bmin1+delb1))cmuni(mu,k)=cmu(mu,1)
          if(bi.lt.bmax1.and.bi.gt.bmin1+delb1) then   
            ir=(bi-bmin1)/delb1
            al=1.-(bi-(ir*delb1+bmin1))/delb1
     		  cmuni(mu,k)=al*cmu(mu,ir)+(1.-al)*cmu(mu,ir+1)
          endif
	    if(bi.ge.bmax1) then
            ir=(bi-bmax1)/delb2
            al=1.-(bi-(ir*delb2+bmax1))/delb2
            if(ir+nb1.gt.0.and.ir+nb1.lt.nbmax)
     & 		  cmuni(mu,k)=al*cmu(mu,ir+nb1)+(1.-al)*cmu(mu,ir+nb1+1)
          endif
        enddo
      enddo 
C        	  	   
C inelastic scattering amplitudes 
      do i=1,ntheta
        th=tr*theta(i)/2.
        qt=2.*pcmhc*sin(th)
C impact parameter integration
        do mu=-mumax,mumax
          sumn=0.
          do k=1,nbig
            bi=k*bstep
            phasec=2.*eta*dlog(pcmhc*bi)
            dum1=qt*bi
            muabs=abs(mu)
            if(muabs.eq.0)bj=bessj0(dum1)
            if(muabs.eq.1)bj=bessj1(dum1)
            if(muabs.gt.1)bj=bessj(muabs,dum1)
            sumn=sumn+ig(k)*bi*bj*bstep*
     &        		exp(ei*(phasec+phasni(k)))*cmuni(mu,k)
          enddo
          finel(mu,i)=sumn*ei*pcmhc/3.
        enddo
      enddo  
C  end integration
150   format(2x,g8.3,10x,g8.3)
      if(ioinel.eq.2)goto 200
C  inelastic scattering cross section in the center of mass
      write(82,'('' Angle (degrees)   Inelastic cross section 
     & (mb/sr)'')')
      do i=1,ntheta
C average over initial and sum over final spins
        sigma=0. 
        do mu=-mumax,mumax
          sigma=sigma+abs(finel(mu,i))**2
        enddo
        sigma=10.*sigma/(2.*spin(1)+1.)
        write(82,150)theta(i),sigma
      enddo
      goto 300
C  inelastic scattering cross section in the laboratory
200   write(82,'('' Angle (degrees)   Inelastic cross section 
     & (mb/sr)'')')
      do i=1,ntheta
        fac=(gcm**2*(cos(theta(i)*tr)+rhog)**2+
     &		  sin(theta(i)*tr)**2)**1.5/
     &    	  gcm/abs(1.+rhog*cos(theta(i)*tr))
C average over initial and sum over final spins
        sigma=0. 
        do mu=-20,20
          sigma=sigma+abs(finel(mu,i))**2
        enddo
        sigma=10.*fac*sigma/(2.*spin(1)+1.)
        fac=sin(theta(i)*tr)/gcm/(cos(theta(i)*tr)+rhog)
        theta(i)=atan(fac)/tr 
        if(i.gt.2.and.theta(i).lt.theta(i-1))goto 300
        write(82,150)theta(i),sigma
      enddo   
300   return 
      end 
C
C=====================================================================
C      
      subroutine fcoul(y,theta,ntheta,fc)
C 
C  elastic Coulomb amplitude
C
      implicit real*8 (a-h,o-z)
      parameter(neuler=5000)
C NOTE: keep 'neuler' big (~ 5000) so that Coulomb phase for large
C charges is accurate.
      complex*16  y, fc(ntheta)
      real*8 theta(ntheta),delc(neuler),dif(neuler)
      common/a14/bmin1,bmax1,eta,pcmhc
      euler=0.5772156649
      pi=acos(-1.)
      tr=pi/180.
C Coulomb phase
	delc(1)=atan(eta)
	do 49 l=2,neuler
	  delc(l)=delc(l-1)+atan(eta/float(l))
	  fase=eta*log(l+0.5)
	  dif(l)=delc(l)-fase
	  dife=dif(l)-dif(l-1)
	  if(abs(dife).lt.1.0e-04)go to 490
   49 continue
  490 argam=-dif(l)
C       
      do i=1,ntheta
        th=tr*theta(i)/2.
        fact=-2.*eta*dlog(sin(th))+2.*argam
        sqrut=eta/(2.*pcmhc*sin(th)**2)
        fc(i)=-exp(y*fact)*sqrut
      enddo
      return
      end
C
C=====================================================================
C      
      subroutine rsharps(iunit)
C
C allows comments in input files
C
      character*1 isharp
    1 read(iunit,'(A1)')isharp
      if(isharp.eq."#")goto 1
      backspace iunit
      return
      end
C
C=====================================================================
C
      real*8 function threej(a1,a2,a3,b1,b2,b3)
C
C     Returns WIGNER 3J-COEFFICIENTS (for Clebsh coeficients change
C     name to clebsh and delete lines where b3=-b3)
C     a1,a2,a3 are j1,j2,j3 and b1,b2,b3 are m1,m2,m3.
C     array f(n+1)=n!/a**n, where a is an arbitrary number
C     that cancels in the calculation of clebsh; its purpose
C     is to prevent under/overflows in array f
C      
      implicit real*8 (a-h,o-z) 
      dimension f(85)
      save f
      data ifirst /1/
C
      if (ifirst .eq. 1) then
	f(1)=1.0
	do 10 i=1,84
10      f(i+1)=f(i)*dfloat(i)/15.
	ifirst=0
      endif
C
C  delete this line to run this program for clebsh coefficients      
      b3=-b3 
C
      clebsh=0. 
      threej=0.
      if(abs(b1)-a1 .gt. 0.01)  goto 100
      if(abs(b2)-a2 .gt. 0.01)  goto 100
      if(abs(b3)-a3 .gt. 0.01)  goto 100
      if(abs(b1+b2-b3).gt.0.01) goto 100
      if(a3-a1-a2 .gt. 0.01)    goto 100 
      if(abs(a1-a2)-a3 .gt. 0.01)  goto 100
      l1=a1+a2-a3+1.01
      l2=a1+b1+1.01
      l3=a2+b2+1.01
      l4=a2-b2+1.01
      l5=a3+a2-a1+1.01
      l6=a1-b1+1.01
      l7=a3+b3+1.01
      l8=a3-b3+1.01
      l9=a3+a1-a2+1.01
      l10=a3+a1+a2+1.01
      m1=l1
      m2=l6
      m3=l3
      m4=(l7+l8-l3-l4+l2-l6)/2+1
      m5=(l7+l8-l2-l6-l3+l4)/2+1
      m6=1
      zw=sqrt(f(l1)*f(l2)*f(l6))
      tw=sqrt(f(l9)/f(l10)*(2.0*a3+1.)/float(l10))
      sw=tw*sqrt(f(l3)*f(l7)*f(l8)*f(l4)*f(l5))
      min=max0(-m4,-m5,-m6)+2
      max=min0(m1,m2,m3)
      fr=-(-1.)**min
      if(max.lt.min)  return
      do m=min,max
	 mw=m-1
	 clebsh=clebsh+sw*fr/f(m4+mw)/(f(m2-mw)/zw)/f(m3-mw)/
     1                   f(m1-mw)/f(m5+mw)/f(m6+mw)
	 fr=-fr
      enddo
      threej=(-1)**nint(a1-a2+b3)/sqrt(2.*a3+1.)*clebsh
100   continue
C
C  delete this line to run this program for clebsh coefficients      
      b3=-b3 
C
      return
      end
C
C=====================================================================
C      
      Real*8 FUNCTION RACAH(A1,A2,A3,B1,B2,B3)
C
C     6-J SYMBOLS, RACAH FORMULA,  A.MESSIAH II (1962), APENDIX C
C     FAKM(N) = GAMMA(N)*10**(-N+1)
C
      IMPLICIT REAL*8 (A-H,O-Z)
      implicit integer(i-n)
      COMMON/REDFAK/FAKM(100)
      external delta
C
      RACAH=0.0D0
      DEL1=DELTA(A1,A2,A3)
      IF(DEL1.EQ.0.0D0) RETURN
      DEL2=DELTA(A1,B2,B3)
      IF(DEL2.EQ.0.0D0) RETURN
      DEL3=DELTA(B1,A2,B3)
      IF(DEL3.EQ.0.0D0) RETURN
      DEL4=DELTA(B1,B2,A3)
      IF(DEL4.EQ.0.0D0) RETURN
      K1=A1+A2+A3+0.001D0
      K2=A1+B2+B3+0.001D0
      K3=B1+A2+B3+0.001D0
      K4=B1+B2+A3+0.001D0
      K5=A1+A2+B1+B2+2.001D0
      K6=A2+A3+B2+B3+2.001D0
      K7=A3+A1+B3+B1+2.001D0
      MU=MAX0(K1,K2,K3,K4)+1
      MO=MIN0(K5,K6,K7)-1
      IF(MO.LT.MU) RETURN
      SU=0.0D0
      DO 1 K=MU,MO
      J1=K-K1
      J2=K-K2
      J3=K-K3
      J4=K-K4
      J5=K5-K
      J6=K6-K
      J7=K7-K
1     SU=-SU+FAKM(K+1)/(FAKM(J1)*FAKM(J2)*FAKM(J3)*FAKM(J4)*FAKM(J5)
     1*FAKM(J6)*FAKM(J7))
      SS=-1.0D0
      IF(MOD(MO,2).EQ.1) SS=1.0D0
      RACAH=SS*SU*DSQRT(DEL1*DEL2*DEL3*DEL4)*10.0D0
      IF(DABS(RACAH).LT.1.0D-6)RACAH=0.0D0
      RETURN
      END
C
C=====================================================================
C      
      Real*8 FUNCTION DELTA(A,B,C)
C
C     SUBPROGRAM OF RACAH-ROUTINE
C     FAKM(N) = GAMMA(N)*10**(-N+1)
C
      IMPLICIT REAL*8 (A-H,O-Z)
      implicit integer(i-n)
      COMMON/REDFAK/FAKM(100)
      DELTA=0.0D0
      TEST=A+B+C-DINT(A+B+C+0.001d0)+0.001D0 ! Integer relation
      IF(TEST.GT.0.1D0) RETURN
      J1=A+B-C+100.01D0          ! triangular relation
      J2=B+C-A+100.01D0
      J3=C+A-B+100.01D0
      J4=A+B+C+2.01D0
      J1=J1-99
      J2=J2-99
      J3=J3-99
      IF(MIN0(J1,J2,J3 ).LE.0) RETURN
      DELTA=0.1d0*FAKM(J1)*FAKM(J2)/FAKM(J4)*FAKM(J3)
      RETURN
      END
C
C=====================================================================
C      
      SUBROUTINE FACTOR
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/REDFAK/FAKM(100)
      FAKM(1)=1.0D0
      DO 1 I=1,40
        AI=dfloat(I)
        FAKM(I+1)=FAKM(I)*0.1D0*AI
1     CONTINUE
      RETURN
      END
C
C=====================================================================
C      
      SUBROUTINE Derivative(Wf,Wf1,Wf2,Dx,Npnts)
C
C   Calculates first (WF1) and second derivative (WF2) of function WF
C   Initialize Ncount=0 in calling program
C
      implicit real*8(a-h,o-z)
      real*8 A(6,6),B(6,6)
      complex*16 Wf(0:Npnts),Wf1(0:Npnts),Wf2(0:Npnts)
      common/coeff/ Ncount
      if(Ncount.gt.0) go to 1
      Ncount=1
      ds=120.0d0*dx
      A(1,1)=-274.0d0/ds
      A(1,2)=600.0d0/ds
      A(1,3)=-600.0d0/ds
      A(1,4)=400.0d0/ds
      A(1,5)=-150.0d0/ds
      A(1,6)=24.0d0/ds
      A(2,1)=-24.0d0/ds
      A(2,2)=-130.0d0/ds
      A(2,3)=240.0d0/ds
      A(2,4)=-120.0d0/ds
      A(2,5)=40.0d0/ds
      A(2,6)=-6.0d0/ds
      A(3,1)=6.0d0/ds
      A(3,2)=-60.0d0/ds
      A(3,3)=-40.0d0/ds
      A(3,4)=120.0d0/ds
      A(3,5)=-30.0d0/ds
      A(3,6)=4.0d0/ds
      A(4,1)=-4.0d0/ds
      A(4,2)=30.0d0/ds
      A(4,3)=-120.0d0/ds
      A(4,4)=40.0d0/ds
      A(4,5)=60.0d0/ds
      A(4,6)=-6.0d0/ds
      A(5,1)=6.0d0/ds
      A(5,2)=-40.0d0/ds
      A(5,3)=120.0d0/ds
      A(5,4)=-240.0d0/ds
      A(5,5)=130.0d0/ds
      A(5,6)=24.0d0/ds
      A(6,1)=-24.0d0/ds
      A(6,2)=150.0d0/ds
      A(6,3)=-400.0d0/ds
      A(6,4)=600.0d0/ds
      A(6,5)=-600.0d0/ds
      A(6,6)=274.0d0/ds
C
      ds=60.0d0*Dx**2
C
      B(1,1)=225.0d0/ds
      B(1,2)=-770.0d0/ds
      B(1,3)=1070.0d0/ds
      B(1,4)=-780.0d0/ds
      B(1,5)=305.0d0/ds
      B(1,6)=-50.0d0/ds
      B(2,1)=50.0d0/ds
      B(2,2)=-75.0d0/ds
      B(2,3)=-20.0d0/ds
      B(2,4)=70.0d0/ds
      B(2,5)=-30.0d0/ds
      B(2,6)=5.0d0/ds
      B(3,1)=-5.0d0/ds
      B(3,2)=80.0d0/ds
      B(3,3)=-150.0d0/ds
      B(3,4)=80.0d0/ds
      B(3,5)=-5.0d0/ds
      B(3,6)=0.0d0
      B(4,1)=0.0d0
      B(4,2)=-5.0d0/ds
      B(4,3)=80.0d0/ds
      B(4,4)=-150.0d0/ds
      B(4,5)=80.0d0/ds
      B(4,6)=-5.0d0/ds
      B(5,1)=5.0d0/ds
      B(5,2)=-30.0d0/ds
      B(5,3)=70.0d0/ds
      B(5,4)=-20.0d0/ds
      B(5,5)=-75.0d0/ds
      B(5,6)=50.0d0/ds
      B(6,1)=-50.0d0/ds
      B(6,2)=305.0d0/ds
      B(6,3)=-780.0d0/ds
      B(6,4)=1070.0d0/ds
      B(6,5)=-770.0d0/ds
      B(6,6)=225.0d0/ds
    1 continue
C
      do 2 Ir=3,NPNTS-3
	WF1(Ir)=A(3,1)*Wf(Ir-2)+A(3,2)*Wf(Ir-1)+A(3,3)*Wf(Ir)
     &    +A(3,4)*Wf(Ir+1)+A(3,5)*Wf(Ir+2)+A(3,6)*Wf(Ir+3)
	WF2(Ir)=B(3,1)*Wf(Ir-2)+B(3,2)*Wf(Ir-1)+B(3,3)*Wf(Ir)
     &    +B(3,4)*Wf(Ir+1)+B(3,5)*Wf(Ir+2)+B(3,6)*Wf(Ir+3)
    2 continue
      Ir=1
      WF1(IR)=A(1,1)*Wf(Ir)+A(1,2)*Wf(Ir+1)+A(1,3)*Wf(Ir+2)
     &    +A(1,4)*Wf(Ir+3)+A(1,5)*Wf(Ir+4)+A(1,6)*Wf(Ir+5)
C
      WF2(IR)=B(1,1)*Wf(Ir)+B(1,2)*Wf(Ir+1)+B(1,3)*Wf(Ir+2)
     &    +B(1,4)*Wf(Ir+3)+B(1,5)*Wf(Ir+4)+B(1,6)*Wf(Ir+5)
C
      Ir=2
      WF1(IR)=A(2,1)*Wf(Ir-1)+A(2,2)*Wf(Ir)+A(2,3)*Wf(Ir+1)
     &    +A(2,4)*Wf(Ir+2)+A(2,5)*Wf(Ir+3)+A(2,6)*Wf(Ir+4)
C
      WF2(IR)=B(2,1)*Wf(Ir-1)+B(2,2)*Wf(Ir)+B(2,3)*Wf(Ir+1)
     &    +B(2,4)*Wf(Ir+2)+B(2,5)*Wf(Ir+3)+B(2,6)*Wf(Ir+4)
C
      Ir=NPNTS-2
      WF1(IR)=A(4,1)*Wf(Ir-3)+A(4,2)*Wf(Ir-2)+A(4,3)*Wf(Ir-1)
     &    +A(4,4)*Wf(Ir)+A(4,5)*Wf(Ir+1)+A(4,6)*Wf(Ir+2)
      WF2(IR)=B(4,1)*Wf(Ir-3)+B(4,2)*Wf(Ir-2)+B(4,3)*Wf(Ir-1)
     &    +B(4,4)*Wf(Ir)+B(4,5)*Wf(Ir+1)+B(4,6)*Wf(Ir+2)
C
      Ir=NPNTS-1
      WF1(IR)=A(5,1)*Wf(Ir-4)+A(5,2)*Wf(Ir-3)+A(5,3)*Wf(Ir-2)
     &    +A(5,4)*Wf(Ir-1)+A(5,5)*Wf(Ir)+A(5,6)*Wf(Ir+1)
      WF2(IR)=B(5,1)*Wf(Ir-4)+B(5,2)*Wf(Ir-3)+B(5,3)*Wf(Ir-2)
     &    +B(5,4)*Wf(Ir-1)+B(5,5)*Wf(Ir)+B(5,6)*Wf(Ir+1)
C
      Ir=NPNTS
      WF1(IR)=A(6,1)*Wf(Ir-5)+A(6,2)*Wf(Ir-4)+A(6,3)*Wf(Ir-3)
     &    +A(6,4)*Wf(Ir-2)+A(6,5)*Wf(Ir-1)+A(6,6)*Wf(Ir)
      WF2(IR)=B(6,1)*Wf(Ir-5)+B(6,2)*Wf(Ir-4)+B(6,3)*Wf(Ir-3)
     &    +B(6,4)*Wf(Ir-2)+B(6,5)*Wf(Ir-1)+B(6,6)*Wf(Ir)
C
      return
      end
C
C=====================================================================
C      
      complex*16 function  Ylm(l,mm,theta,phi)
C
C  This routine computes spherical harmonics (theta, phi in radians)
C
      implicit real*8 (a-h,o-z)
      complex*16 cephi
      real*8 fak(0:130)
      common /gf2 / fak
      m=abs(mm)
	pi=3.141597265
      x=cos(theta)  
      sinp = sin(phi)
      cosp = cos(phi)
      ci=cmplx(0.,1.)
      cephi = cosp + ci*sinp    
      if(m.gt.l.or.abs(x).gt.1.) then
       write(12,*)'bad arguments in Ylm'
       stop
      endif	  
      pmm=1.
      if(m.gt.0) then
        somx2=sqrt((1.-x)*(1.+x))
        fact=1.
        do 11 i=1,m
          pmm=-pmm*fact*somx2
          fact=fact+2.
11      continue
      endif
      if(l.eq.m) then
        plgndr=pmm
      else
        pmmp1=x*(2*m+1)*pmm
        if(l.eq.m+1) then
          plgndr=pmmp1
        else
          do 12 ll=m+2,l
            pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
            pmm=pmmp1
            pmmp1=pll
12        continue
          plgndr=pll
        endif
      endif
      Ylm=sqrt((2.*l+1.)/4./pi*fak(l-m)/fak(l+m))*plgndr*cephi**m
      if(mm.lt.0)Ylm=(-1.)**m*conjg(Ylm) 
      return
      end
C
C=====================================================================
C      
      subroutine legendc(l,theta,plc)
c
c  legendre polynomials of degree l of cos(theta)
c
c  NOTE: theta in degrees
c
      real*8 theta,plc,x,tr,pl6,pl7,pl8,pl9
      if (l.lt.0) l=0
      if (l.gt.14) l=14 
      tr= 0.017453293*theta
      x=cos(tr)
      if (l.eq.0) plc=1
      if (l.eq.1) plc=x
      if (l.eq.2) plc=1.5*(x**2)-0.5
      if (l.eq.3) plc=2.5*(x**3)-1.5*x
      if (l.eq.4) plc=4.375*(x**4)-3.75*(x**2)+0.375 
      if (l.eq.5) plc=(1./8.)*(63.*(x**5)-70.*(x**3)+15.*x)
      if (l.eq.6) plc=(1./16.)*(231.*(x**6)-315.*(x**4)+105.*(x**2)-5.)
      if (l.eq.7) plc=(1./16.)*(429.*(x**7)-693.*(x**5)+315.*(x**3)
     &                 -35.*x)
      if (l.gt.7) then
          pl6 =(1./16.)*(231.*(x**6)-315.*(x**4)+105.*(x**2)-5.)
          pl7 =(1./16.)*(429.*(x**7)-693.*(x**5)+315.*(x**3)
     &          -35.*x)
          pl8 =(15./8.)*x*pl7-(7./8.)*pl6
          pl9 =(17./9.)*x*pl8-(8./9.)*pl7
          pl10=(19./10.)*x*pl9-(9./10.)*pl8
          pl11=(21./11.)*x*pl10-(10./11.)*pl9
          pl12=(23./12.)*x*pl11-(11./12.)*pl10
          pl13=(25./13.)*x*pl12-(12./13.)*pl11
          pl14=(27./14.)*x*pl13-(13./14.)*pl12
          If (l.eq.8)  plc=pl8
          If (l.eq.9)  plc=pl9
          If (l.eq.10) plc=pl10
          If (l.eq.11) plc=pl11
          If (l.eq.12) plc=pl12
          If (l.eq.13) plc=pl13
          If (l.eq.14) plc=pl14
      endif
      return
      end	
C
C=====================================================================
C      
      real*8 function bessj(n,x)
C
C  Bessel function of order n
C
      implicit real*8(a-h,o-z)
C
      INTEGER*4 n,IACC
      REAL*8 x,BIGNO,BIGNI
      PARAMETER (IACC=40,BIGNO=1.e10,BIGNI=1.e-10)
CU    USES bessj0,bessj1
      INTEGER*4 j,jsum,m
      REAL*8 ax,bj,bjm,bjp,sum,tox,bessj0,bessj1
      if(n.lt.2)pause 'bad argument n in bessj'
      ax=abs(x)
      if(ax.eq.0.)then
        bessj=0.d0
      else if(ax.gt.float(n))then
        tox=2./ax
        bjm=bessj0(ax)
        bj=bessj1(ax)
        do 11 j=1,n-1
          bjp=j*tox*bj-bjm
          bjm=bj
          bj=bjp
11      continue
        bessj=bj
      else
        tox=2./ax
        m=2*((n+int(dsqrt(dfloat(IACC*n))))/2)
        bessj=0.d0
        jsum=0
        sum=0.d0
        bjp=0.d0
        bj=1.
        do 12 j=m,1,-1
          bjm=j*tox*bj-bjp
          bjp=bj
          bj=bjm
          if(abs(bj).gt.BIGNO)then
            bj=bj*BIGNI
            bjp=bjp*BIGNI
            bessj=bessj*BIGNI
            sum=sum*BIGNI
          endif
          if(jsum.ne.0)sum=sum+bj
          jsum=1-jsum
          if(j.eq.n)bessj=bjp
12      continue
        sum=2.*sum-bj
        bessj=bessj/sum
      endif
      if(x.lt.0..and.mod(n,2).eq.1)bessj=-bessj
      return
      END
C
C=====================================================================
C      
      real*8 function bessj0(x)
C
C  Bessel function of order 0
C
      implicit real*8(a-h,o-z)
      data p1,p2,p3,p4,p5/1.d0,-.1098628627d-2,.2734510407d-4,
     * -.2073370639d-5,.2093887211d-6/,q1,q2,q3,q4,q5/-.1562499995d-1,
     * .1430488765d-3,-.6911147651d-5,.7621095161d-6,-.934945152d-7/
      data r1,r2,r3,r4,r5,r6/57568490574.d0,-13362590354.d0,
     * 651619640.7d0,-11214424.18d0,77392.33017d0,-184.9052456d0/,
     * s1,s2,s3,s4,s5,s6/57568490411.d0,1029532985.d0,
     * 9494680.718d0,59272.64853d0,267.8532712d0,1.d0/
      if(abs(x).lt.8.)then
      y=x**2
      bessj0=(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))
     * /(s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6)))))
      else
      ax=abs(x)
      z=8./ax
      y=z**2
      xx=ax-.785398164
      bessj0=sqrt(.636619772/ax)*(cos(xx)*(p1+y*(p2+y*(p3+y*(p4+
     * y*p5))))-z*sin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))
      endif
      return
      end
C
C=====================================================================
C      
      real*8 function bessj1(x)
C
C  Bessel function of order 1
C
      implicit real*8(a-h,o-z)
      real*8 x
      real*8 ax,xx,z
      real*8  p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,
     *s1,s2,s3,s4,s5,s6,y
      save p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,s1,s2,s3,s4,
     *s5,s6
      data r1,r2,r3,r4,r5,r6/72362614232.d0,-7895059235.d0,
     *242396853.1d0,-2972611.439d0,15704.48260d0,-30.16036606d0/,s1,s2,
     *s3,s4,s5,s6/144725228442.d0,2300535178.d0,18583304.74d0,
     *99447.43394d0,376.9991397d0,1.d0/
      data p1,p2,p3,p4,p5/1.d0,.183105d-2,-.3516396496d-4,
     *.2457520174d-5,-.240337019d-6/, q1,q2,q3,q4,q5/.04687499995d0,
     *-.2002690873d-3,.8449199096d-5,-.88228987d-6,.105787412d-6/
      if(abs(x).lt.8.)then
        y=x**2
        bessj1=x*(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))/(s1+y*(s2+y*(s3+
     *y*(s4+y*(s5+y*s6)))))
      else
        ax=abs(x)
        z=8./ax
        y=z**2
        xx=ax-2.356194491
        bessj1=sqrt(.636619772/ax)*(cos(xx)*(p1+y*(p2+y*(p3+y*(p4+y*
     *p5))))-z*sin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))*sign(1.d0,x)
      endif
      return
      end
C
C=====================================================================
C      
      subroutine gfv(ie)
C
C  Calculates sign and factorials  of integers and half int.
C
C     iv(n)  =  (-1)**n
C     fak(n) =  n!
C     fad(n) =  n!!
C
      implicit real*8 (a-h,o-z)
      real*8 iv(0:130),fak(0:130),fad(0:130)
      common /gf1 / iv
      common /gf2 / fak
      common /gf3 / fad
C
      if (ie.gt.130) then
       write(12,*) '******    STOP IN GFV: IE LARGER 130 '
       stop
      endif
      iv(0)  = +1
      fak(0) =  1.d0
      fad(0)=1.
      fad(1)=1.
      do 10 i=1,ie
         iv(i)  = -iv(i-1)
10       fak(i) = i*fak(i-1)
      do 11 i=3,ie,2
         fad(i) = i*fad(i-2)
11    continue 
      do 12 i=2,ie,2
         fad(i) = i*fad(i-2)
12    continue 
      return
	end
C
C
C=====================================================================
C **** Routines for Runge-Kutta integration begin here
C
      SUBROUTINE odeint(ystart,nvar,x1,x2,eps,h1,hmin,nok,nbad,derivs,
     *rkqs)
C
C Complex Runge-Kutta integration routine with adaptive stepsize control     
C
      implicit real*8(a-h,o-z)
      COMPLEX*16 ystart(nvar)
      PARAMETER (MAXSTP=10000,KMAXX=200,TINY=1.d-30)
      include 'dweiko.dim'
      COMPLEX*16 dydx(NSTMAX),y(NSTMAX),yp(NSTMAX,KMAXX),
     &     	yscal(NSTMAX)
      EXTERNAL derivs,rkqs
      REAL*8 xp(KMAXX)
      COMMON /path/ kmax,kount,dxsav,xp,yp
      x=x1
      h=sign(h1,x2-x1)
      nok=0
      nbad=0
      kount=0
      do 11 i=1,nvar
        y(i)=ystart(i)
11    continue
      if (kmax.gt.0) xsav=x-2.*dxsav
      do 16 nstp=1,MAXSTP
        call derivs(x,y,dydx)
        do 12 i=1,nvar
          yscal(i)=abs(y(i))+abs(h*dydx(i))+TINY
12      continue
        if(kmax.gt.0)then
          if(abs(x-xsav).gt.abs(dxsav)) then
            if(kount.lt.kmax-1)then
              kount=kount+1
              xp(kount)=x
              do 13 i=1,nvar
                yp(i,kount)=y(i)
13            continue
              xsav=x
            endif
          endif
        endif
        if((x+h-x2)*(x+h-x1).gt.0.) h=x2-x
        call rkqs(y,dydx,nvar,x,h,eps,yscal,hdid,hnext,derivs)
        if(hdid.eq.h)then
          nok=nok+1
        else
          nbad=nbad+1
        endif
        if((x-x2)*(x2-x1).ge.0.)then
          do 14 i=1,nvar
            ystart(i)=y(i)
14        continue
          if(kmax.ne.0)then
            kount=kount+1
            xp(kount)=x
            do 15 i=1,nvar
              yp(i,kount)=y(i)
15          continue
          endif
          return
        endif
        if(abs(hnext).lt.hmin) then
         write(12,*)'******    stepsize smaller than minimum in odeint'
         stop
        endif 
        h=hnext
16    continue
      write(12,*)'******    too many steps in odeint'
      stop
      return
      END
C
C=====================================================================
C
      SUBROUTINE rkqs(y,dydx,n,x,htry,eps,yscal,hdid,hnext,derivs)
C
CU    USES derivs,rkck
      implicit real*8(a-h,o-z)
      COMPLEX*16 dydx(n),y(n),yscal(n)
      EXTERNAL derivs
      include 'dweiko.dim'
      COMPLEX*16 yerr(NSTMAX),ytemp(NSTMAX)
      PARAMETER (SAFETY=0.9,PGROW=-.2,PSHRNK=-.25,ERRCON=1.89d-4)
      h=htry
1     call rkck(y,dydx,n,x,h,ytemp,yerr,derivs)
      errmax=0.
      do 11 i=1,n
        errmax=max(errmax,abs(yerr(i)/yscal(i)))
11    continue
      errmax=errmax/eps
      if(errmax.gt.1.)then
        h=SAFETY*h*(errmax**PSHRNK)
        if(h.lt.0.1*h)then
          h=.1*h
        endif
        xnew=x+h
        if(xnew.eq.x) then 
         write(12,*)'******    stepsize underflow in rkqs'
         stop
        endif
        goto 1
      else
        if(errmax.gt.ERRCON)then
          hnext=SAFETY*h*(errmax**PGROW)
        else
          hnext=5.*h
        endif
        hdid=h
        x=x+h
        do 12 i=1,n
          y(i)=ytemp(i)
12      continue
        return
      endif
      END
C
C=====================================================================
C
      SUBROUTINE rkck(y,dydr,n,x,h,yout,yerr,derivs)
C
      implicit real*8(a-h,o-z)
      COMPLEX*16 dydr(n),y(n),yerr(n),yout(n)
      EXTERNAL derivs
      include 'dweiko.dim'
CU    USES derivs
      COMPLEX*16 ak2(NSTMAX),ak3(NSTMAX),ak4(NSTMAX),ak5(NSTMAX),
     & 	ak6(NSTMAX),ytemp(NSTMAX)
      PARAMETER (A2=.2,A3=.3,A4=.6,A5=1.,A6=.875,B21=.2,B31=3./40.,
     *B32=9./40.,B41=.3,B42=-.9,B43=1.2,B51=-11./54.,B52=2.5,
     *B53=-70./27.,B54=35./27.,B61=1631./55296.,B62=175./512.,
     *B63=575./13824.,B64=44275./110592.,B65=253./4096.,C1=37./378.,
     *C3=250./621.,C4=125./594.,C6=512./1771.,DC1=C1-2825./27648.,
     *DC3=C3-18575./48384.,DC4=C4-13525./55296.,DC5=-277./14336.,
     *DC6=C6-.25)
      do 11 i=1,n
        ytemp(i)=y(i)+B21*h*dydr(i)
11    continue
      call derivs(x+A2*h,ytemp,ak2)
      do 12 i=1,n
        ytemp(i)=y(i)+h*(B31*dydr(i)+B32*ak2(i))
12    continue
      call derivs(x+A3*h,ytemp,ak3)
      do 13 i=1,n
        ytemp(i)=y(i)+h*(B41*dydr(i)+B42*ak2(i)+B43*ak3(i))
13    continue
      call derivs(x+A4*h,ytemp,ak4)
      do 14 i=1,n
        ytemp(i)=y(i)+h*(B51*dydr(i)+B52*ak2(i)+B53*ak3(i)+B54*ak4(i))
14    continue
      call derivs(x+A5*h,ytemp,ak5)
      do 15 i=1,n
        ytemp(i)=y(i)+h*(B61*dydr(i)+B62*ak2(i)+B63*ak3(i)+B64*ak4(i)+
     *B65*ak5(i))
15    continue
      call derivs(x+A6*h,ytemp,ak6)
      do 16 i=1,n
        yout(i)=y(i)+h*(C1*dydr(i)+C3*ak3(i)+C4*ak4(i)+C6*ak6(i))
16    continue
      do 17 i=1,n
        yerr(i)=h*(DC1*dydr(i)+DC3*ak3(i)+DC4*ak4(i)+DC5*ak5(i)+DC6*
     *ak6(i))
17    continue
      return
      END
C
C=====================================================================
C
      SUBROUTINE rk4(y,dydx,n,x,h,yout,derivs)
C
      implicit real*8(a-h,o-z)
      COMPLEX*16 dydx(n),y(n),yout(n)
      EXTERNAL derivs
      include 'dweiko.dim'
      COMPLEX*16 dym(NSTMAX),dyt(NSTMAX),yt(NSTMAX)
      hh=h*0.5
      h6=h/6.
      xh=x+hh
	call derivs(x,y,dydx)
      do 11 i=1,n
        yt(i)=y(i)+hh*dydx(i)
11    continue
      call derivs(xh,yt,dyt)
      do 12 i=1,n
        yt(i)=y(i)+hh*dyt(i)
12    continue
      call derivs(xh,yt,dym)
      do 13 i=1,n
        yt(i)=y(i)+h*dym(i)
        dym(i)=dyt(i)+dym(i)
13    continue
      call derivs(x+h,yt,dyt)
      do 14 i=1,n
        yout(i)=y(i)+h6*(dydx(i)+dyt(i)+2.*dym(i))
14    continue
      return
      END
C
C End routines for Runge-Kutta integration
C
C=====================================================================
C *** Routines for spline interpolation begin here
C      
      SUBROUTINE spline(x,y,n,yp1,ypn,y2)
C
      implicit real*8(a-h,o-z)
      REAL*8 x(n),y(n),y2(n)
      PARAMETER (NMAX=500)
      REAL*8 u(NMAX)
      if (yp1.gt..99d30) then
        y2(1)=0.
        u(1)=0.
      else
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+
     *1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*
     *u(i-1))/p
11    continue
      if (ypn.gt..99d30) then
        qn=0.
        un=0.
      else
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      return
      END
C
C=====================================================================
C      
      SUBROUTINE splint(xa,ya,y2a,n,x,y)
C  spline interpolation
C
      implicit real*8(a-h,o-z)
      REAL*8 xa(n),y2a(n),ya(n)
      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.) then
       write(12,*)'******    bad xa input in splint'
       stop
      endif
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**
     *2)/6.
      return
      END
C
C End routines for spline interpolation 
C=======================================================================
C *** Routines for Simpson integration begin here
C
      SUBROUTINE trapzd(X,Y,Y2,Nu,a,b,s,n)
c
c Integration of a real function by trapezoidal rule.  
C N = integration points, A (B) is lower (upper) limit, Y (Y2) is
c a function (its second derivative) calculated at Nu array points X.
C Integral result = S  
c
      implicit real*8(a-h,o-z)       
      REAL*8 x(nu),y(nu),Y2(nu)
      IF (n.EQ.1) THEN
        CALL splint(X,Y,Y2,nu,a,funca)
        CALL splint(X,Y,Y2,nu,b,funcb)
        s = .5*(b-a) * (funca + funcb)
      ELSE
        it = 2**(n-2)
        tnm = it
        del = (b-a)/tnm
        xx = a+.5*del
        sum=.0
        DO j= 1,it
          CALL splint(x,y,y2,nu,xx,funcx)
          sum = sum + funcx
          xx = xx + del
        ENDDO
        s = .5 * (s + (b-a)*sum/tnm)
      ENDIF
      RETURN
      END
c
C=======================================================================
c
      SUBROUTINE qsimp(x,y,y2,n,a,b,s)
c
c Integration of a real function by Simpson's rule. Use routine TRAPZD.  
C N = integration points, A (B) is lower (upper) limit, Y (Y2) is
c a function (its second derivative) calculated at N array points X.
C Integral result = S  
c
      implicit real*8(a-h,o-z)   
      REAL*8 x(n),y(n),y2(n)
      PARAMETER( eps = 1.d-3,jmax = 20)
      ost = -1.d30
      os  = -1.d30
      DO j = 1,jmax
        CALL trapzd(x,y,y2,n,a,b,st,j)
        s = (4.*st-ost)/3.
        IF (abs(s-os) .le. eps*abs(os)) RETURN
        os = s
        ost = st
      ENDDO
      WRITE(12,*)'******     too many steps in qsimp'
      STOP
      END	           
C End routines for Simpson's integration
C
C--------------------------------------------------------------------
C
C      HOW TO MAKE AN INPUT FOR THIS PROGRAM?
C--------------------------------------------
C 
C             'dweiko.dim' (dimension of vectors and matrices)
C             -----------
C   NMAX = maximum number of levels.  
C   NBMAX = maximum # of impact parameters. 
C   NSTMAX = maximum # of total magnetic substates 
C    A good estimate is NSTMAX = (2*J_max+1)*NST, where J_max is the
C    maximum angular momentum of the states. 
C   NGRID= maximum # of points in optical potential grid (EVEN).
C
C              'dweiko.in' (unit 10)
C              ---------- 
C A comment line can be introduced by using a '#' at the front of 
C each line.
C
C ****   ap,zp,at,zt,eca    ****
C input of system charges and masses (AP,ZP,AT,ZT), bombarding energy
C per nucleon (ECA) in MeV.
C
C ****   iw,iopm,ioelas,ioinel, iogam   **** 
C IOPM=1(0) for output (none) of optical model potentials. 
C IOELAS=(0)[1]{2} for (no) [center of mass] {laboratory} elastic scattering 
C        cross section.
C IOINEL=(0)[1]{2} for (no) [center of mass] {laboratory} inelastic scattering 
C        cross section.
C IOGAM=(0)[1]{2} for (no) output [output of statistical tensors] {output of 
C gamma-ray angular distributions}
C 
C ****   nb,accur,bmin,iob   ***
C NB=number of points in impact parameter mesh (<=NBMAX)
C ACCUR=accuracy required for time integration at each impact parameter.
C BMIN=minimum impact parameter (enter 0 for default)
C IOB=1(0) prints (does not print) out impact parameter probabilities.  
C 
C ****   iopw,iopnuc   ***
C  IOPW=0 (no optical potential,IOELAS =0)
C  IOPW=1 (Woods-Saxon potential) 
C  IOPW=2 (read from file 'optw.in'.)
C  IOPW=3 ('t-rho-rho' folding potential) 
C  IOPW=4 ('M3Y' folding potential)
C If optical potential is provided (IOPW=2), it should  be stored in 
C 'optw.in' in rows of R x Real[U(R)] x Imag[U(R)]. 
C The first line in this file is the number of rows (maximum=NGRID). 
C The program makes an interpolation for intermediate values. 
C  IOPNUC=1(0) is an option to compute (or not) nuclear excitation.
C
C****   V0_ws,r0,d_ws,VI_ws,r0_I,dI_ws   ***   
C If IOPW=1, enter V0_ws [VI_ws] = real part [imaginary] (> 0, both) of Woods-Saxon.
C                  r0 [r0_I] = radius parameter (R_ws = r0 * (ap^1/3 + at^1/3)).
C                  d_ws [dI_ws] = diffuseness.
C If IOPW is not equal to 1, place a '#' sign at the front of this line, or delete it.
C
C****   VS, r0_S, dS, V_surf, d_surf   **** 
C If IOPW=1 and projectile, or target is a proton, include spin-orbit and surface
C parameters. If not, place a '#' sign at the front of this line, or delete it.
C          VS = spin orbit potential depth parameter (> 0).
C          r0_S = radius parameter.
C          dS = diffuseness.
C          V_surf = surface potential depth parameter (> 0).
C          d_surf = difuseness.
C
C ****   Wrat  ***
C If IOPW=4, enter Wrat = ratio of imaginary to real part of M3Y interaction.
C If IOPW is not equal to 4, place a '#' sign at the beginning of this line, or delete it.
C
C ****  thmax, ntheta  *** 
C If IOELAS=1,2 or IOINEL=1,2 enter here THMAX, maximum angle (in degrees and in the
C center of mass), and NTHETA, the number of points in scatering angle (<=NGRID).
C If IOELAST or IOINEL are not 1, or 2, place a '#' sign at the beginning of this line,
C or delete it.
C
C ****   jinel    *** 
C If IOINEL=1 enter the state (JINEL) for the inelastic angular distribution (1 < JINEL <= NST). 
C If IOINEL is not 1, or 2, place a '#' sign at the beginning of this line, or delete it.
C
C ****   nst   *** 
C NST is the number of nuclear levels (<=NMAX).
C
C ****   i, ex(i), spin(i)   ***
C Input of state label (I), energy (EX), and  spin (SPIN). 
C I ranges from 1 to NST.
C
C ****   i,j,matrixE1,matrixE2,matrixE3,matrixM1,matrixM2   ***
C Reduced matrix elements for E1, E2, E3, M1 and M2 excitations:
C    <I_j||O(E/M;L)||I_i>,      j > i ,
C for the electromagnetic transitions.
C Add a row of zeros at the end of this list. If no electromagnetic excitation
C is wanted, just enter a row of zeros.
C 
C ****   j,f0(j),f1(j),f2(j),f3(j)   ***
C If IOPNUC=1 enter sum rule fraction of the nuclear deformation parameters for 
C monopole, dipole, quadrupole nuclear excitations (DELTE0,DELTE1,DELTE2,DELTE3) 
C for each excited state J: DELTE_i = f_i * (sum rules).
C If IOPNUC=0 insert a comment card ('#') in front of each entry row, or 
C delete them.
C
C **** iff,igg,thmin,thmax,ntheta   ***
C If IOGAM=2, enter here the initial and final states (iff > igg) for the 
C gamma transition, the minimum and maximum values of theta, and the number 
C of theta points (< NGRID).
C_______________________________________________
C  SPECIAL CASE: EXCITATION OF GIANT RESONANCES:
C  A good guess for the electric reduced matrix elements is obtained from 
C
C  |<I_j||O(E1)||I_i>|^2 = (2I_i+1) .
C                       9/(4 pi) . hbar^2/(2 m_N) . e^2 . NZ/A / E_x  
C
C  where  E_x = E_j - E_i .  
C
C    |<I_j||O(E2)||I_i>|^2  = (2I_i+1) .
C        15/(4 pi) . hbar^2/m_N . e^2 . R^2 . Z^2/A / E_x (for isoscalar)
C        15/(4 pi) . hbar^2/m_N . e^2 . R^2 . NZ/A / E_x (for isovector)
C
C  The excitation of GR's exhaust a fraction of these sum rules, which is 
C  close to unity.
C--------------------------------------------------------------------
C                         END PROGRAM
C--------------------------------------------------------------------
