#!/bin/sh

echo "Test reading in parameters"

make test_read_par

cd tests/

chmod u+x test_par_class.out

./test_par_class.out parameters_test.dat

cd ..


echo "Test constancy of parameters"

make test_const_par

cd tests/

chmod u+x test_const_par.out

./test_const_par.out parameters_test.dat

cd ..

make clean_test
