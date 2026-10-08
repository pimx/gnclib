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

# gnclib 1.1.0: unit tests of the integration and the rotating Earth
test: build
	$(CXX) $(CXXFLAGS) tests/test_integrate.cpp -o build/test_integrate
	./build/test_integrate
	$(CXX) $(CXXFLAGS) tests/test_ini.cpp build/libgnc.a -o build/test_ini
	./build/test_ini
	$(CXX) $(CXXFLAGS) tests/test_point.cpp -o build/test_point
	./build/test_point

