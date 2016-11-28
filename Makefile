CC = gcc
CXX = g++
CXXFLAGS = -std=c++11 
INCLUDES = -I ~/Software/bamtools/include/ -I /opt/gurobi651/linux64/include/
LDFLAGS = -L ~/Software/bamtools/lib/ -L /opt/gurobi651/linux64/lib/
LDLIBS = -lbamtools -lgurobi_c++ -lgurobi65 -lz

SRCS = src/main.cpp src/ReadRec.cpp src/SegmentGraph.cpp src/WriteIO.cpp src/Config.cpp

all: squid

squid:
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(LDFLAGS) -o bin/squid $(SRCS) $(LDLIBS)
