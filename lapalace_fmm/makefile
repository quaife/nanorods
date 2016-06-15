#MEX=/Applications/MATLAB_R2015a.app/bin/mex
MEX = /usr/local/matlab/bin/mex

FORTRAN = gfortran
#FORTRAN = ifort
#FFLAGS = -c -fpic -openmp -O3
FFLAGS = -c -fpic -O3 -fopenmp

#MEXLIBDIR = /share/apps/intel/Compiler/12.1/composer_xe_2011_sp1.7.256/compiler/lib/intel64/
#for Ronaldo
#MEXLIB = -lifcore -lirc -lsvml -limf -lgfortran -liomp5 
MEXLIB = -lgomp

OBJS=\
lfmm2drouts.o\
cfmm2dpart.o\
rfmm2dpart.o\
d2mtreeplot.o\
d2tstrcr_omp.o\
l2dterms.o\
laprouts2d.o\
prini.o\
laplaceDLP.o

all: cfmm2dpart d2mtreeplot d2tstrcr_omp l2dterms laprouts2d lfmm2drouts prini rfmm2dpart laplaceDLP mex

mex: laplaceDLPfmm.F $(OBJS) 
	$(MEX) -O laplaceDLPfmm.F $(OBJS) $(MEXLIB)

cfmm2dpart: cfmm2dpart.f
	$(FORTRAN) $(FFLAGS) cfmm2dpart.f

d2mtreeplot: d2mtreeplot.f
	$(FORTRAN) $(FFLAGS) d2mtreeplot.f

d2tstrcr_omp: d2tstrcr_omp.f
	$(FORTRAN) $(FFLAGS) d2tstrcr_omp.f

l2dterms: l2dterms.f
	$(FORTRAN) $(FFLAGS) l2dterms.f

laprouts2d: laprouts2d.f
	$(FORTRAN) $(FFLAGS) laprouts2d.f

lfmm2drouts: lfmm2drouts.f
	$(FORTRAN) $(FFLAGS) lfmm2drouts.f

prini: prini.f
	$(FORTRAN) $(FFLAGS) prini.f

rfmm2dpart: rfmm2dpart.f
	$(FORTRAN) $(FFLAGS) rfmm2dpart.f

laplaceDLP: laplaceDLP.f
	$(FORTRAN) $(FFLAGS) laplaceDLP.f

clean: 
	rm *.o
	cp *.mex* ../src
	rm *.mex*
