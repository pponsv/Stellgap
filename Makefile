all:
	cd ./src && $(MAKE) all
	mv ./src/x* ./bin/
	$(MAKE) link

clean:
	cd ./src && $(MAKE) clean
	
link:
	# cd ./bin && ln -f * ~/.local/bin/
	ln -r -sf ./bin/* $${HOME}/.local/bin/
	
rmlink:
	rm ~/.local/bin/xstgap
	rm ~/.local/bin/xstgap_snd
	rm ~/.local/bin/xmetric

tests:
	echo $${HOME}
