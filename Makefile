all:
	mkdir -p build;
	cd build; \
	cmake ..; \
	make

install:
	mkdir -p build;
	cd build; \
	cmake ..; \
	make install

clean:
	rm -rf ./build

