.PHONY : build

all: clean build rmlink link

link:
	mkdir -p $${HOME}/.local/bin/
	ln -r -sf ./bin/x* $${HOME}/.local/bin/
	
rmlink:
	rm -f ~/.local/bin/xstgap
	rm -f ~/.local/bin/xstgap_snd
	rm -f ~/.local/bin/xmetric

build:
	mkdir -p ./build
	cmake --build ./build --target xstgap_new

clean:
	cmake --build ./build --target clean

clean_bld:
	rm -rf bld/*

build_all:
	cmake --build ./build --target all