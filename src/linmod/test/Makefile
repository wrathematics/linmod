FC = gfortran

all: distributions transpose

distributions:
	$(FC) -o distributions distributions.f90 -I../build/include ../build/liblinmod_static.a -lgomp

transpose:
	$(FC) -o transpose transpose.f90 -I../build/include ../build/liblinmod_static.a -lgomp


clean:
	rm -rf ./transpose
	rm -rf ./distributions

