this= .
PROG = MainHART_hres
O = .

SRCS =  READ_HART_cosmo.f cooling.f trans3.f trans2.f trans.f getrow.f \
	findplane.f findcent.f components.f velcent.f Vcir_comp.f \
	vradrd_cosmo.f barprop_cosmo.f dispQ_cosmo.f fouriermod_cosmo.f \
        Main.f

OBJS =  $(O)/READ_HART_cosmo.o $(O)/trans3.o $(O)/trans2.o \
	$(O)/trans.o $(O)/getrow.o $(O)/cooling.o $(O)/findplane.o \
	$(O)/findcent.o $(O)/components.o $(O)/velcent.o $(O)/Vcir_comp.o \
	$(O)/vradrd_cosmo.o $(O)/barprop_cosmo.o $(O)/dispQ_cosmo.o  \
        $(O)/fouriermod_cosmo.o  $(O)/Main.o

LIBS =

FC = ifort

#FFLAGS  = -O3 -i-dynamic -convert big_endian -mcmodel=medium
#LDFLAGS = -O3 -i-dynamic -convert big_endian -mcmodel=medium

#FFLAGS  = -O3 -i-dynamic -convert big_endian -mcmodel=large -auto-scalar
#LDFLAGS = -O3 -i-dynamic -convert big_endian -mcmodel=large -auto-scalar

FFLAGS  = -O3  -convert big_endian
LDFLAGS = -O3 -convert big_endian

all: $(PROG)

$(PROG): $(OBJS)
	$(FC) $(LDFLAGS) -o $(this)/$(PROG) $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS)


LIB =

$(O)/Main.o: Main.f parameters.h $(depall)
	$(FC) $(FFLAGS) -c Main.f $(LIB)

$(O)/READ_HART_cosmo.o: READ_HART_cosmo.f parameters.h $(depall)
	$(FC) $(FFLAGS) -c READ_HART_cosmo.f $(LIB)

$(O)/cooling.o: cooling.f parameters.h $(depall)
	$(FC) $(FFLAGS) -c cooling.f $(LIB)

$(O)/trans3.o: trans3.f parameters.h $(depall)
	$(FC) $(FFLAGS) -c trans3.f $(LIB)

$(O)/trans2.o: trans2.f parameters.h $(depall)
	$(FC) $(FFLAGS) -c trans2.f $(LIB)

$(O)/trans.o: trans.f parameters.h $(depall)
	$(FC) $(FFLAGS) -c trans.f $(LIB)

$(O)/getrow.o: getrow.f parameters.h $(depall)
	$(FC) $(FFLAGS) -c getrow.f $(LIB)

$(O)/findplane.o: findplane.f parameters.h $(depall)
	$(FC) $(FFLAGS) -c findplane.f $(LIB)

$(O)/findcent.o: findcent.f parameters.h $(depall)
	$(FC) $(FFLAGS) -c findcent.f $(LIB)

$(O)/velcent.o: velcent.f parameters.h $(depall)
	$(FC) $(FFLAGS) -c velcent.f $(LIB)

$(O)/components.o: components.f parameters.h $(depall)
	$(FC) $(FFLAGS) -c components.f $(LIB)

$(O)/Vcir_comp.o: Vcir_comp.f parameters.h $(depall)
	$(FC) $(FFLAGS) -c Vcir_comp.f $(LIB)

$(O)/vradrd_cosmo.o: vradrd_cosmo.f parameters.h $(depall)
	$(FC) $(FFLAGS) -c vradrd_cosmo.f $(LIB)

$(O)/barprop_cosmo.o: barprop_cosmo.f parameters.h $(depall)
	$(FC) $(FFLAGS) -c barprop_cosmo.f $(LIB)

$(O)/dispQ_cosmo.o: dispQ_cosmo.f parameters.h $(depall)
	$(FC) $(FFLAGS) -c dispQ_cosmo.f $(LIB)

$(O)/fouriermod_cosmo.o: fouriermod_cosmo.f parameters.h $(depall)
	$(FC) $(FFLAGS) -c fouriermod_cosmo.f $(LIB)
