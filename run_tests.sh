#/bin/sh
#
#SPARSKIT script for making and running all test programs.
#Updated by Jaad A. Tannous to run with the modern fortran versions

make all

cd test_programs

echo compiling Blassm Test

gfortran -o blassm_test.ex blassm_tester.f95  ../src/non-library/functns.f95 -L ../. -lskit -Wl,--unresolved-symbols=ignore-all

echo Running Blassm Test

./blassm_test.ex

echo compiling matvec Test

gfortran -o matvec_test.ex matvec_tester.f95  ../src/non-library/functns.f95 -L ../. -lskit -Wl,--unresolved-symbols=ignore-all

echo Running matvec Test

./matvec_test.ex
