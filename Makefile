BOOST = /opt/local/boost # decompressed boost path
BAMTOOLS = /opt/local/bamtools-2.4.0 # bamtools path
GUROBI = /opt/local/gurobi702	# gurobi path

CC = gcc
CXX = g++
CXXFLAGS = -std=c++11 -D_GLIBCXX_USE_CXX11_ABI=0
INCLUDES = -I ${BAMTOOLS}/include -I ${GUROBI}/include -I ${BOOST}
LDFLAGS = -L ${BAMTOOLS}/lib/ -L ${GUROBI}/lib/
LDLIBS = -lbamtools -lgurobi_c++ -lgurobi70 -lz
RPATH =  ${BAMTOOLS}/lib/bamtools:${GUROBI}/lib/

SRCS = src/main.cpp src/ReadRec.cpp src/SegmentGraph.cpp src/WriteIO.cpp src/Config.cpp

all: squid

squid:
	mkdir -p bin
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(LDFLAGS) -o bin/squid $(SRCS) $(LDLIBS) -Wl,-rpath,${RPATH}
