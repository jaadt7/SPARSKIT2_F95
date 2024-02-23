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
F95 = gfortran
AR = ar -rv 
RANLIB = ranlib	

OPT = -c -O -Wall -w
LIB = libskit.a

#
# library objects
#

OBJ =	src/library/blassm.o	\
		src/library/matvec.o	\
		src/library/formats.o	\
		src/library/unary.o		\
		src/library/infofun.o	\
		src/library/inout.o		\
		src/library/ilut.o		\
		src/library/iters.o		\
		src/library/genmat.o	\
		src/library/elmtlib2.o	\
		src/library/femgen.o	\
		src/library/meshes.o	\
		src/library/sobel.o		\
		src/library/zlatev.o	\
		src/library/ccn.o		\
		src/library/color.o		\
		src/library/dsepart.o

#
# non-library
#

DIRS =	.				\
		src/library 	\
		src/non-library	\
		test_programs/	\

$(LIB): $(OBJ) 
	$(AR) $@ $(OBJ) 
	$(RANLIB) $@

# do not ranlib on some architectures

clean:
	@for dir in $(DIRS) ;\
          do \
          echo cleaning $$dir ;\
          (cd $$dir; rm -f *~ *.o *.ex *.fort *.mat *.a *.ps *.pic fort.* *. *.hb) ;\
          done
clean_obj:
	@for dir in $(DIRS) ;\
          do \
          echo cleaning $$dir ;\
          (cd $$dir; rm -f *.o) ;\
          done

tarit: 
	(cd ..; tar cvf - SPARSKIT2_F95 | gzip -c > SPARSKIT2_F95.tar.gz) 


src/library/blassm.o: src/library/blassm.f95
	(cd src/library ; $(F95)  $(OPT) blassm.f95)
src/library/matvec.o: src/library/matvec.f95
	(cd src/library ; $(F95)  $(OPT) matvec.f95)
src/library/formats.o: src/library/formats.f95
	(cd src/library ; $(F95)  $(OPT) formats.f95)
src/library/unary.o: src/library/unary.f95
	(cd src/library ; $(F95)  $(OPT) unary.f95)
src/library/infofun.o: src/library/infofun.f95
	(cd src/library ; $(F95)  $(OPT) infofun.f95)
src/library/inout.o: src/library/inout.f95
	(cd src/library; $(F95)  $(OPT) inout.f95)
src/library/ilut.o: src/library/ilut.f95
	(cd src/library; $(F95)  $(OPT) ilut.f95)
src/library/iters.o: src/library/iters.f95
	(cd src/library; $(F95)  $(OPT) iters.f95)
src/library/genmat.o: src/library/genmat.f95
	(cd src/library ; $(F95)  $(OPT) genmat.f95)
src/library/elmtlib2.o: src/library/elmtlib2.f95
	(cd src/library ; $(F95)  $(OPT) elmtlib2.f95)
src/library/femgen.o: src/library/femgen.f95
	(cd src/library ; $(F95)  $(OPT) femgen.f95)
src/library/meshes.o: src/library/meshes.f95
	(cd src/library ; $(F95)  $(OPT) meshes.f95)
src/library/sobel.o: src/library/sobel.f95
	(cd src/library ; $(F95)  $(OPT) sobel.f95)
src/library/zlatev.o: src/library/zlatev.f95
	(cd src/library ; $(F95)  $(OPT) zlatev.f95)
src/library/ccn.o: src/library/ccn.f95
	(cd src/library ; $(F95)  $(OPT) ccn.f95)
src/library/color.o: src/library/color.f95
	(cd src/library ; $(F95)  $(OPT) color.f95)
src/library/dsepart.o: src/library/dsepart.f95
	(cd src/library ; $(F95)  $(OPT) dsepart.f95)
