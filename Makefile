.PHONY : all clean doc compile link rmlink configure

all: install

bld:
	meson bld

configure: bld
	meson setup --wipe bld

install: configure
	meson install -C bld
	$(MAKE) link

clean: 
	rm -rf bld
	rm -rf bin/*
	$(MAKE) rmlink

link:
	mkdir -p $${HOME}/.local/bin/
	ln -r -sf ./bin/x* $${HOME}/.local/bin/

rmlink:
	rm -f ~/.local/bin/xstgap
	rm -f ~/.local/bin/xstgap_snd
	rm -f ~/.local/bin/xmetric
