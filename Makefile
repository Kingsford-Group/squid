BAMTOOLS = /opt/local/stow/bamtools-2.4.0
GUROBI = /opt/local/stow/gurobi702
CC = gcc
CXX = g++
CXXFLAGS = -std=c++11 -D_GLIBCXX_USE_CXX11_ABI=0
INCLUDES = -I ${BAMTOOLS}/include -I ${GUROBI}/include
LDFLAGS = -L ${BAMTOOLS}/lib/bamtools -L ${GUROBI}/lib/
LDLIBS = -lbamtools -lgurobi_c++ -lgurobi70 -lz
RPATH =  ${BAMTOOLS}/lib/bamtools:${GUROBI}/lib/

SRCS = src/main.cpp src/ReadRec.cpp src/SegmentGraph.cpp src/WriteIO.cpp src/Config.cpp

all: squid

squid:
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(LDFLAGS) -o bin/squid $(SRCS) $(LDLIBS) -Wl,-rpath,${RPATH}
