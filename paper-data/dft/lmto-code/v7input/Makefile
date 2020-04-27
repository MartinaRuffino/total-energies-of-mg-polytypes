RANLIB = /usr/bin/ranlib

# --- compiler switches ---
include Make.inc
 
# --- archive library ---
SUBS = $(TOPDIR)/v7input/subs.a
 
ALL:	m_gtv.mod m_toksw.mod m_rdctrl.mod $(SUBS)
# ALL:	modules 

.PRECIOUS:	$(SUBS)

subs:	$(SUBS)
 
subs.a:	$(SUBS)
 
# --- Routines that depend on structures.h ---
$(SUBS)(rdctrl.o):	structures.h

init:
	if test -r "$(INAME)"; then \
	env PATH=$(PATH):$(CCOMPDIR) /bin/csh -f $(INAME) ; \
	fi
 
# objects = m_toksw.mod m_gtv.mod m_rdctrl.mod 

m_gtv.mod :	m_gtv.f
	$(FC) $(FFLAGS) -c $< 

m_gtv.o m_gtv.mod :	m_gtv.f
	$(FC) $(FFLAGS) -c $< ; ar rv subs.a  m_gtv.o

m_toksw.o m_toksw.mod:	m_toksw.f m_gtv.mod
	$(FC) $(FFLAGS) -c $< ; ar rv subs.a  m_toksw.o

m_rdctrl.o m_rdctrl.mod:	m_rdctrl.f m_gtv.mod m_toksw.mod
	$(FC) $(FFLAGS_NONE) -c $< ; ar rv subs.a  m_rdctrl.o

m_rdctrlchk.o m_rdctrlchk.mod :	m_rdctrlchk.f
	$(FC) $(FFLAGS_NONE) -c $< ; ar rv subs.a  m_rdctrlchk.o

# $(SUBS)(m_toksw.o): 
# 	$(FC) $(FFLAGS_NONE) -c m_toksw.f
# 	ar rv $(SUBS) m_toksw.o
# 	rm -f m_toksw.o

$(SUBS)(m_rdctrl.o): 
	$(FC) $(FFLAGS_NONE) -c m_rdctrl.f
	ar rv $(SUBS) m_rdctrl.o
	rm -f m_rdctrl.o

$(SUBS)(m_rdctrlchk.o): 
	$(FC) $(FFLAGS_NONE) -c m_rdctrlchk.f
	ar rv $(SUBS) m_rdctrlchk.o
	rm -f m_rdctrlchk.o

$(SUBS):	$(SUBS)(m_gtv.o) \
	$(SUBS)(m_toksw.o) \
	$(SUBS)(m_rdctrl.o) \
	$(SUBS)(rdctrl.o)\
	$(SUBS)(m_rdctrlchk.o) \
	$(SUBS)(rdctrlchk.o)

clean: 
	-rm -f $(SUBS) Names m_gtv.mod  m_rdctrl.mod  m_rdctrlchk.mod  m_toksw.mod
