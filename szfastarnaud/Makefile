FC	= ftn
FFLAGS	= -O3 
OBJS_ps	= multipole.o cosmo.o nfw.o da.o linearpk.o sigma.o spline.o hunt.o chebyshev.o growth.o dverk.o integrand_ps.o qromb.o qgaus1.o rombint.o rombint1.o rombint2.o mvir2mdel.o zbrent.o mf_tinker_redshift.o
OBJS_ts	= multipole.o cosmo.o nfw.o da.o linearpk.o sigma.o spline.o hunt.o chebyshev.o growth.o dverk.o integrand_ts.o qromb.o qgaus1.o rombint.o rombint1.o rombint2.o mvir2mdel.o zbrent.o mf_tinker_redshift.o

.SUFFIXES: .f90
.f90.o: 
	$(FC) $(FFLAGS) -c $<
default: szfast_ps
szfast_ps: $(OBJS_ps) szfast_ps.o
	$(FC) $(FFLAGS) -o $@ $(OBJS_ps) $@.o $(LDFLAGS)
szfast_ts:$(OBJS_ts) szfast_ts.o
	$(FC) $(FFLAGS) -o $@ $(OBJS_ts) $@.o $(LDFLAGS)
clean:
	-rm -f *.o *.mod
tidy: clean
	-rm -f szfast_ps
