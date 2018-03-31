CXX=c++
CXXFLAGS=-Wall -pedantic -std=c++14 -DNDEBUG -O3
CXXTESTFLAGS=-Wall -pedantic -std=c++14 -DNDEBUG

.PHONY=all

all : CEvo_sim.cpp global_fct.cpp par_class.cpp 
	$(CXX) $(CXXFLAGS) $^ -o CEvo_sim.out

.PHONY=test_read_par

test_read_par : tests/test_read_pars.o par_class.o
	$(CXX) $(CXXTESTFLAGS) $^ -o tests/test_par_class.out

.PHONY=test_const_par

test_const_par : tests/test_par_const.o par_class.o
	$(CXX) $(CXXTESTFLAGS) $^ -o tests/test_const_par.out
	

tests/test_read_pars.o : tests/test_read_pars.cpp
	$(CXX) $(CXXFLAGS) -c $^ -o $@
tests/test_par_const.o : tests/test_par_const.cpp
	$(CXX) $(CXXFLAGS) -c $^ -o $@
par_class.o : par_class.cpp par_class.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@
global_fct.o : global_fct.cpp global_fct.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@
CEvo_sim.o : CEvo_sim.cpp	
	$(CXX) $(CXXFLAGS) -c $< -o $@

.PHONY=clean

clean : 
	rm -f tests/*.out
	rm -f tests/*.o
	rm -f *.out
	rm -f *.o

.PHONY=clean_test

clean_test : 
	rm -f tests/*.out
	rm -f tests/*.o
