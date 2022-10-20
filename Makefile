.PHONY: all build link clean debug rmlink

all: build link

debug: clean all

clean: rmlink
	@$(MAKE) -C ./src/ clean

build:
	@$(MAKE) -C ./src/ all
	
link:
	mkdir -p $${HOME}/.local/bin/
	ln -r -sf ./bin/x* $${HOME}/.local/bin/
	
rmlink:
	rm ~/.local/bin/xstgap
	rm ~/.local/bin/xstgap_snd
	rm ~/.local/bin/xmetric
