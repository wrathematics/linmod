all: mk_linmod mk_pkg 


mk_linmod:
	cd linmod && make && cd ..


mk_pkg:
	mkdir -p build;
	cd build; \
	cmake ..; \
	make install


mk_pkg: mk_linmod


clean:
	clean:
	rm -rf ./linmod/build

