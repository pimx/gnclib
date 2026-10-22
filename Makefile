# gnclib -- статическая библиотека GNC (g++)
CXX = g++
CXXFLAGS = -std=c++17 -O2 -Wall -Iinclude -fno-exceptions -fno-rtti
SRC = $(wildcard src/*.cpp)
OBJ = $(SRC:src/%.cpp=build/%.o)

all: build/libgnc.a

build:
	mkdir -p build

build/%.o: src/%.cpp | build
	$(CXX) $(CXXFLAGS) -c $< -o $@

build/libgnc.a: $(OBJ)
	ar rcs $@ $(OBJ)

clean:
	rm -f build/*.o build/libgnc.a

# gnclib 1.1.0: unit tests of the integration and the rotating Earth;
# 1.6.0: utilities, tables, observers, gusts, wind profile, orbit (test_common);
# 1.7.0: propellant slosh (test_slosh; CSV of the check -> build/slosh_out);
# 1.7.1: attitude order of IntegrateRK4 (test_integrate, check f);
# 1.8.0: atmosphere GOST 4401-81 to 250 km against the standard's table (test_atmosphere);
# 1.9.0: its Smooth mode (test_atmosphere, f, g);
# 1.10.0: upper atmosphere GOST R 25645.166-2004 against the standard's tables (test_atm_upper)
test: build/libgnc.a
	$(CXX) $(CXXFLAGS) tests/test_integrate.cpp -o build/test_integrate
	./build/test_integrate
	$(CXX) $(CXXFLAGS) tests/test_ini.cpp build/libgnc.a -o build/test_ini
	./build/test_ini
	$(CXX) $(CXXFLAGS) tests/test_point.cpp -o build/test_point
	./build/test_point
	$(CXX) $(CXXFLAGS) tests/test_atmosphere.cpp -o build/test_atmosphere
	./build/test_atmosphere
	$(CXX) $(CXXFLAGS) tests/test_atm_upper.cpp -o build/test_atm_upper
	./build/test_atm_upper
	$(CXX) $(CXXFLAGS) tests/test_common.cpp build/libgnc.a -o build/test_common
	./build/test_common
	mkdir -p build/slosh_out
	$(CXX) $(CXXFLAGS) tests/test_slosh.cpp build/libgnc.a -o build/test_slosh
	./build/test_slosh build/slosh_out
