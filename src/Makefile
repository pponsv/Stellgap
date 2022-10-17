LIBSTELL  = $${HOME}/bin/libstell_dir/
METRIC    = metric_element_create_ver8.46
NOSOUND   = stellgap_ver5
SOUND     = stellgap_soundwave_lagrng_ver7
SERIAL    = SERIAL # PARALLEL or SERIAL

MKLFLAGS  = -lmkl_gf_lp64 -lmkl_core -lmkl_sequential -lpthread -lm
PREPFLAGS = -cpp -P -C -D$(SERIAL)
GFLAGS    = -O2 -ffixed-form
OFILES    = fitpack.o Fourier_lib_convolve.o fourier_lib.o

BUILDPATH = ../bin/

all: makebin $(BUILDPATH)xmetric $(BUILDPATH)xstgap $(BUILDPATH)xstgap_snd

makebin:
	mkdir -p $(BUILDPATH)

clean:
	rm -f *.mod *.o
	rm -f ../bin/*

kind_spec.o: kind_spec.f
	gfortran $(GFLAGS) -c kind_spec.f

fourier_lib.o: fourier_lib.f kind_spec.o
	gfortran $(GFLAGS) -c fourier_lib.f

fitpack.o: fitpack.f kind_spec.o
	gfortran $(GFLAGS) -c fitpack.f

Fourier_lib_convolve.o: Fourier_lib_convolve.f fourier_lib.o fitpack.o
	gfortran $(GFLAGS) -c Fourier_lib_convolve.f

metric.o: $(METRIC).f
	gfortran $(GFLAGS) -c $(METRIC).f -I$(LIBSTELL) -o metric.o

$(BUILDPATH)xmetric: metric.o
	gfortran -O2  metric.o $(LIBSTELL)libstell.so -lnetcdf -llapack -lblas \
		-o $(BUILDPATH)xmetric

$(BUILDPATH)xstgap_snd: $(OFILES) $(SOUND).f
	gfortran $(PREPFLAGS) $(GFLAGS) $(SOUND).f -I. $(MKLFLAGS) $(OFILES) \
		-o $(BUILDPATH)xstgap_snd
		
$(BUILDPATH)xstgap: $(OFILES) $(NOSOUND).f
	gfortran $(PREPFLAGS) $(GFLAGS) $(NOSOUND).f -I. $(MKLFLAGS) $(OFILES) \
		-o $(BUILDPATH)xstgap