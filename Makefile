LIB2= -Wl,--start-group /opt/intel/mkl/lib/intel64/libmkl_intel_lp64.a /opt/intel/mkl/lib/intel64/libmkl_core.a /opt/intel/mkl/lib/intel64/libmkl_gnu_thread.a -Wl,--end-group -lpthread -lm -ldl -lgfortran

all: clean dlasru.o dbdsqru.o CGS2.o DGEBRDG_4_BISIDE.o DGEBRDG_K.o DGEBRDG_LP1.o ERR.o RESGKL.o doqds1.o doqds3.o doqds3.o dlartg2.o dlartg3.o dlartg4.o dlartg5.o dlartg6.o dlartg7.o dfma0.o RESGKL_MAIN.o fileinput_gkl

%.o: src/%.f90
	gfortran -fopenmp -Wall -O3 -mtune=native -march=native -mcmodel=medium -c -o src/$@ $< 

%.o: lib/%.f
	gfortran -fopenmp -Wall -O3 -mtune=native -march=native -mcmodel=medium  -c -o lib/$@ $<

dlartg7.o:
	gfortran -fopenmp -Wall -mcmodel=medium -O2 -c -o lib/dlartg7.o lib/dlartg7.f

dfma0.o:
	gcc -O3 -fopenmp -Wall -mtune=native -march=native -mcmodel=medium -c -o lib/dfma0.o lib/dfma0.c

fileinput_gkl:
	gcc -fopenmp -Wall -O3 -mtune=native -march=native -mcmodel=medium -o fileinput_gkl fileinput_gkl.c lib/dlasru.o lib/dbdsqru.o src/CGS2.o src/DGEBRDG_4_BISIDE.o src/DGEBRDG_K.o src/DGEBRDG_LP1.o src/ERR.o src/RESGKL.o lib/doqds1.o lib/doqds3.o lib/dlartg2.o lib/dlartg3.o lib/dlartg4.o lib/dlartg5.o lib/dlartg6.o lib/dlartg7.o lib/dfma0.o src/RESGKL_MAIN.o ${LIB2}

clean:
	rm -rf lib/*.o src/*.o fileinput_gkl


#%.o: src/%.f90
#	${FC} ${FCFLAG} -c -o src/$@ $< ${LIB}
#%.o: lib/%.f
#	${FC} ${FCFLAG} -c -o lib/$@ $<
#	gfortran -O3 -mtune=native -march=native -c -o lib/$@ $<
#-L/usr/lib/gcc/x86_64-linux-gnu/4.8/libgfortran.a

#DRV=./Driver
#TriEVD=./srcTriEVD
#RGKL=srcRGKL
#MKMAT=SVD_MKMAT
#ifdef debug
#   CCFLAG=-Wall -g -mcmodel=medium -shared-intel -Wall
#   FCFLAG=-warn all -CB -traceback -g -mcmodel=medium -shared-intel
#else
#   CCFLAG=-fp-model precise -O3 -ipo -xHOST -mcmodel=medium -shared-intel -Wall
#   FCFLAG=-fp-model precise -O3 -ipo -xHOST -mcmodel=medium -shared-intel
#endif
##LIB=-mkl -L/opt/intel/lib/intel64 -lifcore -lifcore_pic -lifcoremt -lifport -lchkpwrap -qopenmp


