#-----------------------------------------------------------------------
#                     S P A R S K I T  makefile 
#-----------------------------------------------------------------------
#
# There are three types of source files in SPARSKIT:
# 1. library source files
# 2. non-library and unsupported source files
# 3. test programs
#
# Simply type "make" to compile the SPARSKIT static library.
#
# To compile and run all the test files, simply execute the bash script 
# by typing ./run_tests.sh;  Execution permissions might need to be chaged
#
#This make file has been rewritten  by Jaad Tannous of Clemson University

#Compilation options
SHELL = /bin/sh
F95 = gfortran
AR = ar -rv 
RANLIB = ranlib	
OPT = -c -O -Wall -w

#Library name
LIB = libskit.a

#Source directory
SRC_DIR = src/library/

#temprorary directory
TMP_DIR = tmp


#List of source files
SRCS := $(wildcard $(SRC_DIR)/*.f95)

#Object files
OBJS := $(patsubst $(SRC_DIR)/%.f95,$(TMP_DIR)/%.o,$(SRCS))


$(LIB): $(OBJS) 
	$(AR) $@ $(OBJS) 
	$(RANLIB) $@
	rm -rf $(TMP_DIR)


.PHONY: clean clean_all tarball

$(TMP_DIR)/%.o: $(SRC_DIR)/%.f95 | $(TMP_DIR)
	$(F95) $(FFLAGS) -c $< -o $@

$(TMP_DIR):
	mkdir -p $(TMP_DIR)


clean:
	rm $(LIB)

clean_all: clean
	(cd test_programs/; rm -f *~ *.o *.ex *.fort *.mat *.a *.ps *.pic fort.* *. *.hb)

tarball:
	(cd ..; tar cvf - SPARSKIT2_F95 | gzip -c > SPARSKIT2_F95.tar.gz)
