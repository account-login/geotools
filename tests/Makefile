
CXXFLAGS += -Wall -Wextra -g


test_%: test_%.o catch.o
	$(CXX) $(CXXFLAGS) $^ -o $@

test_%.o: test_%.cpp ../%.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

catch.o: catch.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@


.PRECIOUS: test_%.o
.PHONY: clean

clean:
	-rm -rf *.o
	-rm test_ttlmap
