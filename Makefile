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
