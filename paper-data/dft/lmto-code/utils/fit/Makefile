# ---- Machine-specific compiler flags ---
include Make.inc


PROGS =  sfit lfit nlfit prony

# --- Make rules ---
.f.o: 
	$(FC) $(FFLAGS) -c $<

progs:	nlfitlm pfit $(PROGS)

nlfitlm:	%:	%.o
	$(LK) $(LKFLAGS1) $@.o $(LIBES) $(LKFLAGS2) -o $@

pfit:	%:	%.o p2.o powell.o
	$(LK) $(LKFLAGS1) $@.o p2.o powell.o $(LIBES) $(LKFLAGS2) -o $@

$(PROGS):	%:	%.o p2.o
	$(LK) $(LKFLAGS1) $@.o p2.o $(LIBES) $(LKFLAGS2) -o $@

mvbin:	
	mv -f $(PROGS) $(BIN)

clean:	
	rm -f $(PROGS)

tar:	
	rm -r -f fit; mkdir fit
	tar cf - Make.inc Makefile pfit.f lfit.f nlfit.f nlfitlm.f sfit.f prony.f p2.f powell.f | ( cd fit; tar xvf - )
	mkdir fit/examples
	(cd examples; tar cf - . | ( cd ../fit/examples ; tar xvf - ))
	tar cvf - fit | gzip >FIT.tar.gz
	rm -r -f fit
