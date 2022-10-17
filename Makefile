all: build link

clean:
	cd ./src && $(MAKE) clean

build:
	cd ./src && $(MAKE) all
	
link:
	mkdir -p $${HOME}/.local/bin/
	ln -r -sf ./bin/x* $${HOME}/.local/bin/
	
rmlink:
	rm ~/.local/bin/xstgap
	rm ~/.local/bin/xstgap_snd
	rm ~/.local/bin/xmetric
