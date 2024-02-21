#-----------------------------------------------------------------------
#                     S P A R S K I T  makefile 
#-----------------------------------------------------------------------
#
# There are three types of source files in SPARSKIT:
# 1. library source files
# 2. non-library and unsupported source files
# 3. test programs
#
# Simply type "make" to make the SPARSKIT library.
#
# "make all" will additionally make all the unsupported objects.
#
# To make the test programs, use the makefiles in each subdirectory,
# or see the dotests script.
# There are no references to test programs in this makefile.
# Some test programs use the non-library objects from other directories,
# so do a "make all" before making test programs.
#
# Last Updated:  Nov. 19, 2009 

#This make file has been update by Jaad Tannous of Clemson University
#on Feb 21st 2024

SHELL = /bin/sh
F90 = gfortran

# F90 = F90 
AR = ar -rv 
RANLIB = ranlib	

OPT = -c -O -Wall 
LIB = libskit.a

#
# library objects
#

OBJ =	SPAGged/blassm.o	\
		SPAGged/matvec.o	\
		SPAGged/formats.o	\
		SPAGged/unary.o		\
		SPAGged/infofun.o	\
		SPAGged/inout.o		\
		SPAGged/ilut.o		\
		SPAGged/iters.o		\
		SPAGged/genmat.o	\
		SPAGged/elmtlib2.o	\
		SPAGged/femgen.o	\
		SPAGged/meshes.o	\
		SPAGged/sobel.o		\
		SPAGged/zlatev.o	\
		SPAGged/ccn.o		\
		SPAGged/color.o		\
		SPAGged/dsepart.o

#
# non-library and unsupported objects
#
OBJ2 =	ITSOL/itaux.o		\
	MATGEN/FDIF/functns.o	\
	MATGEN/FEM/functns2.o	\
	UNSUPP/BLAS1/blas1.o	\
	UNSUPP/MATEXP/exppro.o	\
	UNSUPP/MATEXP/phipro.o	\
	UNSUPP/PLOTS/psgrd.o	\
	UNSUPP/PLOTS/texgrid1.o	\
	UNSUPP/PLOTS/texplt1.o

DIRS =	.			\
	BLASSM			\
	FORMATS			\
	INFO			\
	INOUT			\
	ITSOL			\
	MATGEN/FDIF		\
	MATGEN/FEM		\
	SPAGged		\
	SPAGged		\
	UNSUPP/BLAS1		\
	UNSUPP/MATEXP		\
	UNSUPP/PLOTS

$(LIB): $(OBJ) 
	$(AR) $@ $(OBJ) 
	$(RANLIB) $@

# do not ranlib on some architectures

clean:
	@for dir in $(DIRS) ;\
          do \
          echo cleaning $$dir ;\
          (cd $$dir; rm -f *~ *.o *.ex .,* fort.* ftn?? *.mat core, *.a) ;\
          done

tarit: 
	(cd ..; tar cvf - SPARSKIT2 | gzip -c > SPARSKIT2.tar.gz) 

all: $(OBJ) $(OBJ2) libskit.a

SPAGged/blassm.o: SPAGged/blassm.f90
	(cd SPAGged ; $(F90)  $(OPT) blassm.f90)
SPAGged/matvec.o: SPAGged/matvec.f90
	(cd SPAGged ; $(F90)  $(OPT) matvec.f90)
SPAGged/formats.o: SPAGged/formats.f90
	(cd SPAGged ; $(F90)  $(OPT) formats.f90)
SPAGged/unary.o: SPAGged/unary.f90
	(cd SPAGged ; $(F90)  $(OPT) unary.f90)
SPAGged/infofun.o: SPAGged/infofun.f90
	(cd SPAGged ; $(F90)  $(OPT) infofun.f90)
SPAGged/inout.o: SPAGged/inout.f90
	(cd SPAGged; $(F90)  $(OPT) inout.f90)
SPAGged/ilut.o: SPAGged/ilut.f90
	(cd SPAGged; $(F90)  $(OPT) ilut.f90)
SPAGged/iters.o: SPAGged/iters.f90
	(cd SPAGged; $(F90)  $(OPT) iters.f90)
SPAGged/itaux.o: SPAGged/itaux.f90
	(cd SPAGged; $(F90)  $(OPT) itaux.f90)
SPAGged/genmat.o: SPAGged/genmat.f90
	(cd SPAGged ; $(F90)  $(OPT) genmat.f90)
SPAGged/functns.o: SPAGged/functns.f90
	(cd SPAGged ; $(F90)  $(OPT) functns.f90)
SPAGged/elmtlib2.o: SPAGged/elmtlib2.f90
	(cd SPAGged ; $(F90)  $(OPT) elmtlib2.f90)
SPAGged/femgen.o: SPAGged/femgen.f90
	(cd SPAGged ; $(F90)  $(OPT) femgen.f90)
SPAGged/functns2.o : SPAGged/functns2.f90 
	(cd SPAGged ; $(F90)  $(OPT) functns2.f90)
SPAGged/meshes.o: SPAGged/meshes.f90
	(cd SPAGged ; $(F90)  $(OPT) meshes.f90)
SPAGged/sobel.o: SPAGged/sobel.f90
	(cd SPAGged ; $(F90)  $(OPT) sobel.f90)
SPAGged/zlatev.o: SPAGged/zlatev.f90
	(cd SPAGged ; $(F90)  $(OPT) zlatev.f90)
SPAGged/ccn.o: SPAGged/ccn.f90
	(cd SPAGged ; $(F90)  $(OPT) ccn.f90)
SPAGged/color.o: SPAGged/color.f90
	(cd SPAGged ; $(F90)  $(OPT) color.f90)
SPAGged/dsepart.o: SPAGged/dsepart.f90
	(cd SPAGged ; $(F90)  $(OPT) dsepart.f90)
UNSUPP/BLAS1/blas1.o: UNSUPP/BLAS1/blas1.f
	(cd UNSUPP/BLAS1 ; $(F90)  $(OPT) blas1.f)
UNSUPP/MATEXP/exppro.o: UNSUPP/MATEXP/exppro.f
	(cd UNSUPP/MATEXP ; $(F90)  $(OPT) exppro.f)
UNSUPP/MATEXP/phipro.o: UNSUPP/MATEXP/phipro.f
	(cd UNSUPP/MATEXP ; $(F90)  $(OPT) phipro.f)
UNSUPP/PLOTS/psgrd.o : UNSUPP/PLOTS/psgrd.f 
	(cd UNSUPP/PLOTS ; $(F90) $(OPT) psgrd.f)
UNSUPP/PLOTS/texgrid1.o : UNSUPP/PLOTS/texgrid1.f 
	(cd UNSUPP/PLOTS ; $(F90) $(OPT) texgrid1.f)
UNSUPP/PLOTS/texplt1.o : UNSUPP/PLOTS/texplt1.f 
	(cd UNSUPP/PLOTS ; $(F90) $(OPT) texplt1.f)
