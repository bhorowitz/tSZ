FC	= g95
FFLAGS	= -O3 
OBJS	= multipole.o cosmo.o nfw.o da.o linearpk.o sigma.o spline.o hunt.o chebyshev.o growth.o dverk.o integrand.o qromb.o qgaus1.o rombint.o rombint1.o mvir2mdel.o zbrent.o mf_tinker_redshift.o
.SUFFIXES: .f90
.f90.o: 
	$(FC) $(FFLAGS) -c $<
default: szfast 
szfast: $(OBJS) szfast.o
	$(FC) $(FFLAGS) -o $@ $(OBJS) $@.o $(LDFLAGS)
clean:
	-rm -f *.o *.mod
tidy: clean
	-rm -f szfast
