LIBSTELL = $${HOME}/bin/libstell_dir/
METRIC   = metric_element_create_ver8.46
NOSOUND  = stellgap_ver5
SOUND    = stellgap_soundwave_lagrng_ver7

MKLFLAGS = -lmkl_gf_lp64 -lmkl_core -lmkl_sequential -lpthread -lm

all: clean nosound sound metric 

auxfiles:
	gfortran -O2 -ffixed-form -c fitpack.f
	gfortran -O2 -ffixed-form -c Fourier_lib_convolve.f
	
metric:
	@printf "\nBuilding metric... \n\n"
	gfortran -O2 -c $(METRIC).f -I$(LIBSTELL)
	gfortran -O2 -o xmetric $(METRIC).o -I$(LIBSTELL) $(LIBSTELL)libstell.so -lnetcdf -llapack -lblas
		
nosound:
	@printf "\nBuilding nosound... \n\n"
	gfortran -cpp -E -P -C -DSERIAL $(NOSOUND).f > stellgap.f
	gfortran -O2 -ffixed-form -c stellgap.f
	$(MAKE) auxfiles
	gfortran -O2 -ffixed-form -o xstgap -I. stellgap.o $(MKLFLAGS) fitpack.o Fourier_lib_convolve.o
	rm stellgap.f
	
sound: 
	@printf "\nBuilding sound... \n\n"
	gfortran -cpp -E -P -C -DSERIAL $(SOUND).f > stellgap_snd.f
	gfortran -O2 -ffixed-form -c stellgap_snd.f
	$(MAKE) auxfiles
	gfortran -O2 -ffixed-form -o xstgap_snd -I. stellgap_snd.o $(MKLFLAGS) fitpack.o Fourier_lib_convolve.o
	rm stellgap_snd.f

clean:
	@printf "\nCleaning... \n\n"
	-rm -r *.o *.mod x*

silo:
	# It does not work and I don't need it so I won't even try
	@printf "Building silo... \n\n"
