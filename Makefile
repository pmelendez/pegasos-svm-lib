CXX=g++
CXXFLAGS=-std=c++0x -I.

all:
	$(CXX) $(CXXFLAGS) main.cpp

