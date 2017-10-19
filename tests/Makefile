
CXXFLAGS += -Wall -Wextra -g
OPTFLAGS += -O2	# -DNDEBUG
CXX11FLAGS += -std=c++11


test_%: test_%.o catch.o
	$(CXX) $(CXXFLAGS) $^ -o $@

test_%.o: test_%.cpp ../%.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

catch.o: catch.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

bench_query.o: bench_query.cpp ../geotree.hpp ../geodensity.hpp ../geodensity_bounded.hpp
	$(CXX) $(CXXFLAGS) $(CXX11FLAGS) $(OPTFLAGS) -c $< -o $@

bench_query: bench_query.o
	$(CXX) $(CXXFLAGS) $(CXX11FLAGS) $(OPTFLAGS) $^ -o $@


.PRECIOUS: test_%.o
.PHONY: clean

clean:
	-rm -rf *.o *.exe
