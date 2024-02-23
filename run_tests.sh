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

echo Compiling 5 and 7 point test in Harwell-Boeing format

gfortran -o 5_pt_harwell_boeing_test.ex 5_pt_harwell_boeing_tester.f95 ../src/non-library/functns.f95 -L ../. -lskit -Wl,--unresolved-symbols=ignore-all

echo Running 5 and 7 point test

./5_pt_harwell_boeing_test.ex << \EOF
10 10 1
test_pt.mat
EOF

echo Compiling Block Matrix Harwell Boeing test
gfortran -o block_mat_test.ex block_mat_test.f95 ../src/non-library/functns.f95 -L ../. -lskit -Wl,--unresolved-symbols=ignore-all

echo Running Block Matrix Harwell Boeing test
./block_mat_test.ex << \EOF
10 10 1
4
test_bl.mat
EOF

echo Compiling convective-diffusive test
gfortran -o conv_test.ex conv_diff_test.f95 ../src/non-library/functns2.f95 ../src/non-library/psgrd.f95 -L ../. -lskit -Wl,--unresolved-symbols=ignore-all

echo Running convection test
./conv_test.ex << \EOF
2
2
EOF
cat mat.hb

echo Compiling Sobel test

gfortran -o sobel_test.ex sobel_test.f95 -L ../. -lskit -Wl,--unresolved-symbols=ignore-all

echo Running Sobel test

./sobel_test.ex << \EOF
10
EOF

echo Compiling Zlatev test

gfortran -o zlatev_test.ex zlatev_test.f95 -L ../. -lskit -Wl,--unresolved-symbols=ignore-all

echo Running Zlatev test

./zlatev.ex
cat zlatev1.mat
cat zlatev2.mat
cat zlatev3.mat

echo Compiling markov test

gfortran -o markov_test.ex markov_test.f95 -L ../. -lskit -Wl,--unresolved-symbols=ignore-all

echo Running markov tests

./markov.ex << \EOF
10
EOF
cat markov.mat

echo Compiling exponential propagator test

gfortran -o prop_test.ex exp_prop_test.f95 ../src/non-library/exppro.f95

echo Running exponential propagator test

./prop_test.ex << \EOF
0.1
0.00001
10
EOF

echo Compiling phi approximation

gfortran -o phi_approx.ex phi_approx.f95 ../src/non-library/phipro.f95

echo Running Phi Approximation

./phi_approx.ex << \EOF
0.1
0.00001
10
EOF
