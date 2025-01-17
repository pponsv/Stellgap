#	Paths
BIN      = ./bin
SRC      = ./src
BLD      = ./bld
LIBSTELL = ${HOME}/bin/libstell_dir/

#	Filenames
METRIC   = metric_element_create_ver8.46
NOSOUND  = stellgap_ver5
SOUND    = stellgap_soundwave_lagrng_ver7
SERIAL   = SERIAL # PARALLEL or SERIAL
OFILES_i = fitpack.o Fourier_lib_convolve.o fourier_lib.o post_process.o
OFILES   = $(patsubst %, $(BLD)/%, $(OFILES_i))

#	Flags
# MKLFLAGS    = -lmkl_gf_lp64 -lmkl_core -lmkl_sequential -lpthread -lm
MKLFLAGS    = -llapack -lm -lblas -lpthread
METRICFLAGS = $(LIBSTELL)libstell.so -lnetcdf -llapack -lblas # To compile xmetric
PREPFLAGS   = -cpp -P -C -D$(SERIAL) # For preprocessing
GFLAGS      = -O2 -ffixed-form -J$(BLD)

.PHONY: all build link clean debug rmlink

all: build link

build : $(BIN) $(BLD) $(BIN)/xmetric $(BIN)/xstgap $(BIN)/xstgap_snd

debug: clean all

clean: rmlink
	rm -f ./bld/* ./bin/*
	
link:
	mkdir -p $${HOME}/.local/bin/
	ln -r -sf ./bin/x* $${HOME}/.local/bin/
	
rmlink:
	rm -f ~/.local/bin/xstgap
	rm -f ~/.local/bin/xstgap_snd
	rm -f ~/.local/bin/xmetric

#	Compilation
$(BIN)/xmetric: $(SRC)/$(METRIC).f
	gfortran $(GFLAGS) -I $(LIBSTELL) $^ $(METRICFLAGS) -o $@

$(BIN)/xstgap_snd: $(SRC)/$(SOUND).f $(OFILES)
	gfortran $(PREPFLAGS) $(GFLAGS) $^ $(MKLFLAGS) -o $@
		
$(BIN)/xstgap: $(SRC)/$(NOSOUND).f $(OFILES)
	gfortran $(PREPFLAGS) $(GFLAGS) $^ $(MKLFLAGS) -o $@

$(BLD)/%.o : $(SRC)/%.f
	gfortran $(GFLAGS) -c $< -o $@

#	Dependencies
$(BLD)/fourier_lib.o: $(BLD)/kind_spec.o

$(BLD)/fitpack.o: $(BLD)/kind_spec.o

$(BLD)/post_process.o: $(BLD)/kind_spec.o

$(BLD)/Fourier_lib_convolve.o: $(BLD)/fourier_lib.o $(BLD)/fitpack.o

#	Folders
$(BIN):
	@mkdir -p $(BIN)

$(BLD):
	@mkdir -p $(BLD)