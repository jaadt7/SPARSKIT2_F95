#/bin/sh
#
#SPARSKIT script for making and running all test programs.
#Updated by Jaad A. Tannous to run with the modern fortran versions

make

cd test_programs

echo Compiling Blassm Test

gfortran -o blassm_test.ex blassm_tester.f95  ../src/non-library/functns.f95 -L ../. -lskit -Wl,--unresolved-symbols=ignore-all

echo Running Blassm Test

./blassm_test.ex

echo Compiling matvec Test

gfortran -o matvec_test.ex matvec_tester.f95  ../src/non-library/functns.f95 -L ../. -lskit -Wl,--unresolved-symbols=ignore-all

echo Running matvec Test

./matvec_test.ex

echo Compiling format tester

gfortran -o formats_test.ex format_tester.f95  ../src/non-library/functns.f95 -L ../. -lskit -Wl,--unresolved-symbols=ignore-all

echo Running format tester

./formats_test.ex
grep ERROR *.mat
ls -s *.mat
echo Compliling unary tester

gfortran -o unary_test.ex unary_tester.f95  ../src/non-library/functns.f95 -L ../. -lskit -Wl,--unresolved-symbols=ignore-all

echo Running unary tester

./unary_test.ex
cat unary.mat

echo Compiling Variable Block Tester

gfortran -o var_block_test.ex var_block_tester.f95  ../src/non-library/functns.f95 -L ../. -lskit -Wl,--unresolved-symbols=ignore-all

echo Running Varaiable Block Tester
./var_block_test.ex

echo Compiling Info Tester

gfortran -o info_test.ex info_tester.f95 info_test_supp.f95  ../src/non-library/functns.f95 -L ../. -lskit -Wl,--unresolved-symbols=ignore-all

echo Running Info Tester
./info_test.ex < test_supp_files/saylr1

echo Compiling I/O test

gfortran -o io_test.ex io_tester.f95 ../src/non-library/functns.f95 -L ../. -lskit -Wl,--unresolved-symbols=ignore-all

echo Running I/O test
./io_test.ex

echo Compiling Harwell-Boeing test

gfortran -o harwell_boeing_test.ex harwell_boeing_tester.f95 -L ../. -lskit -Wl,--unresolved-symbols=ignore-all

echo Running Harwell-Boeing test

./harwell_boeing_test.ex < test_supp_files/saylr1 > saylr1.ps

echo Compiling Harwell-Boeing Read test

gfortran -o harwell_boeing_read_test.ex harwell_boeing_read_tester.f95 -L ../. -lskit -Wl,--unresolved-symbols=ignore-all

echo Running Harwell-Boeing Read test

./harwell_boeing_read_test.ex < test_supp_files/saylr1 > saylr1.pic

echo Compiling Basic Iterative Solver test

gfortran -o basic_iter_test.ex basic_iter_tester.f95 ../src/non-library/blas1.f95 ../src/non-library/itaux.f95 -L ../. -lskit -Wl,--unresolved-symbols=ignore-all

echo Running Basic Iterative Solve test
./basic_iter_test.ex

echo Compiling Preconditioned GMRES Test

gfortran -o precon_gmres_test.ex precon_gmres_tester.f95 ../src/non-library/blas1.f95 ../src/non-library/itaux.f95 -L ../. -lskit -Wl,--unresolved-symbols=ignore-all

echo Running Preconditioned GMRES test

./precon_gmres_test.ex

echo Compiling Iterative solver with Harwell-Boeing matrix Test

gfortran -o harwell_boeing_iter_test.ex harwell_boeing_iter_tester.f95 ../src/non-library/blas1.f95 ../src/non-library/itaux.f95 -L ../. -lskit -Wl,--unresolved-symbols=ignore-all

echo Iterative solver with Harwell-Boeing matrix Test

./harwell_boeing_iter_test.ex < test_supp_files/saylr1