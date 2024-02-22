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

OBJ =	src/blassm.o	\
		src/matvec.o	\
		src/formats.o	\
		src/unary.o		\
		src/infofun.o	\
		src/inout.o		\
		src/ilut.o		\
		src/iters.o		\
		src/genmat.o	\
		src/elmtlib2.o	\
		src/femgen.o	\
		src/meshes.o	\
		src/sobel.o		\
		src/zlatev.o	\
		src/ccn.o		\
		src/color.o		\
		src/dsepart.o

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
	src		\
	src		\
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

src/blassm.o: src/blassm.f90
	(cd src ; $(F90)  $(OPT) blassm.f90)
src/matvec.o: src/matvec.f90
	(cd src ; $(F90)  $(OPT) matvec.f90)
src/formats.o: src/formats.f90
	(cd src ; $(F90)  $(OPT) formats.f90)
src/unary.o: src/unary.f90
	(cd src ; $(F90)  $(OPT) unary.f90)
src/infofun.o: src/infofun.f90
	(cd src ; $(F90)  $(OPT) infofun.f90)
src/inout.o: src/inout.f90
	(cd src; $(F90)  $(OPT) inout.f90)
src/ilut.o: src/ilut.f90
	(cd src; $(F90)  $(OPT) ilut.f90)
src/iters.o: src/iters.f90
	(cd src; $(F90)  $(OPT) iters.f90)
src/itaux.o: src/itaux.f90
	(cd src; $(F90)  $(OPT) itaux.f90)
src/genmat.o: src/genmat.f90
	(cd src ; $(F90)  $(OPT) genmat.f90)
src/functns.o: src/functns.f90
	(cd src ; $(F90)  $(OPT) functns.f90)
src/elmtlib2.o: src/elmtlib2.f90
	(cd src ; $(F90)  $(OPT) elmtlib2.f90)
src/femgen.o: src/femgen.f90
	(cd src ; $(F90)  $(OPT) femgen.f90)
src/functns2.o : src/functns2.f90 
	(cd src ; $(F90)  $(OPT) functns2.f90)
src/meshes.o: src/meshes.f90
	(cd src ; $(F90)  $(OPT) meshes.f90)
src/sobel.o: src/sobel.f90
	(cd src ; $(F90)  $(OPT) sobel.f90)
src/zlatev.o: src/zlatev.f90
	(cd src ; $(F90)  $(OPT) zlatev.f90)
src/ccn.o: src/ccn.f90
	(cd src ; $(F90)  $(OPT) ccn.f90)
src/color.o: src/color.f90
	(cd src ; $(F90)  $(OPT) color.f90)
src/dsepart.o: src/dsepart.f90
	(cd src ; $(F90)  $(OPT) dsepart.f)
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
