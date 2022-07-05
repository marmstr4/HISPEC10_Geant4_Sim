# make for clx,clx_tk,dcy_tkk2
#FFLAGS=-B/usr/lib/gcc/x86_64-linux-gnu/3.4.6/
#FFLAGS=-B/usr/lib/x86_64-linux-gnu/
FFLAGS=-B/usr/lib/x86_64-linux-gnu/

ALLOBJ=newpage.o ffctrls.o cleb.o
CLXOBJ=djmm.o q.o threej.o
DCYOBJ=fit.o fk.o racah.o intyy.o ykkap.o spline.o wninej.o \
       fkk2k1.o filly.o inty.o

#F77=/usr/bin/gfortran
FFLAGS+=-fPIE
F77=/usr/bin/g77
#F77=/usr/bin/fort77

#g77
####
#F77=gfortran
#FFLAGS+=-std=legacy

# old version clx (max. 50 levels)

clx: clx.o ampder_50.o $(ALLOBJ) $(CLXOBJ)
	$(F77) $(FFLAGS) -o clx clx.o ampder_50.o $(ALLOBJ) $(CLXOBJ)

# my version clx_tk (max. 100 levels, number controlled by PARAMETER)
# my means in this case Thorsten
#clx_tk: clx_tk.o ampder.o $(ALLOBJ) $(CLXOBJ)
#	$(F77) $(FFLAGS) -o clx_tk clx_tk.o ampder.o $(ALLOBJ) $(CLXOBJ)


# mymy version (K. Wimmer), modified output, more useful for plotting, nothing else changed
#clx_kw: clx_kw.o ampder.o $(ALLOBJ) $(CLXOBJ)
#	$(F77) $(FFLAGS) -o clx_kw clx_kw.o ampder.o $(ALLOBJ) $(CLXOBJ)

# my dcy_tkk2 (Thorsten)
#dcy_tkk2: dcy_tkk2.o $(ALLOBJ) $(DCYOBJ)
#	$(F77) $(FFLAGS) -o dcy_tkk2 dcy_tkk2.o $(ALLOBJ) $(DCYOBJ)

####

clx.o: clx.for
	$(F77) -c $(FFLAGS) clx.for

clx_tk.o: clx_tk.for
	$(F77) -c $(FFLAGS) clx_tk.for

clx_kw.o: clx_kw.for
	$(F77) -c $(FFLAGS) clx_kw.for

dcy_tkk2.o: dcy_tkk2.for
	$(F77)  -c $(FFLAGS) dcy_tkk2.for

ampder_50.o: ampder_50.for
	$(F77) -c $(FFLAGS) ampder_50.for

ampder.o: ampder.for
	$(F77) -c $(FFLAGS) ampder.for

newpage.o: newpage.for
	$(F77) -c $(FFLAGS) newpage.for

djmm.o: djmm.for
	$(F77) -c $(FFLAGS) djmm.for

ffctrls.o: ffctrls.for
	$(F77) -c $(FFLAGS) ffctrls.for

q.o: q.for
	$(F77) -c $(FFLAGS) q.for

threej.o: threej.for
	$(F77) -c $(FFLAGS) threej.for

cleb.o: cleb.for
	$(F77) -c $(FFLAGS) cleb.for

fit.o: fit.for
	$(F77) -c $(FFLAGS) fit.for

fk.o: fk.for
	$(F77) -c $(FFLAGS) fk.for

racah.o: racah.for
	$(F77) -c $(FFLAGS) racah.for

intyy.o: intyy.for
	$(F77) -c $(FFLAGS) intyy.for

ykkap.o: ykkap.for
	$(F77) -c $(FFLAGS) ykkap.for

spline.o: spline.for
	$(F77) -c $(FFLAGS) spline.for

wninej.o: wninej.for
	$(F77) -c $(FFLAGS) wninej.for

fkk2k1.o: fkk2k1.for
	$(F77) -c $(FFLAGS) fkk2k1.for

filly.o: filly.for
	$(F77) -c $(FFLAGS) filly.for

inty.o: inty.for
	$(F77) -c $(FFLAGS) inty.for
