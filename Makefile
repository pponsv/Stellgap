all:
	cd ./src && $(MAKE) all
	mkdir -p ./bin/
	mv ./src/x* ./bin/
	$(MAKE) link

clean:
	cd ./src && $(MAKE) clean
	
link:
	mkdir -p $${HOME}/.local/bin/
	ln -r -sf ./bin/* $${HOME}/.local/bin/
	
rmlink:
	rm ~/.local/bin/xstgap
	rm ~/.local/bin/xstgap_snd
	rm ~/.local/bin/xmetric
