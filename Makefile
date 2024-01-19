.PHONY : build

all: clean configure build_all rmlink link

link:
	mkdir -p $${HOME}/.local/bin/
	ln -r -sf ./bin/x* $${HOME}/.local/bin/
	
rmlink:
	rm -f ~/.local/bin/xstgap
	rm -f ~/.local/bin/xstgap_snd
	rm -f ~/.local/bin/xmetric

configure:
	cmake -S . -B ./bld
build:
	cmake --build ./bld --target xstgap

clean:
	rm -rf bld/ bin/

build_all:
	cmake --build ./bld --target all