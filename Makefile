BOOST = /opt/local/stow/boost-1.55.0-gcc
BAMTOOLS = /opt/local/stow/bamtools-2.4.0
GUROBI = /opt/local/stow/gurobi702

CC = gcc
CXX = g++
CXXFLAGS = -std=c++11 -D_GLIBCXX_USE_CXX11_ABI=0 $(INCLUDES)
INCLUDES = -I $(BAMTOOLS)/include -I $(GUROBI)/include -I $(BOOST)
LDFLAGS = -L $(BAMTOOLS)/lib/bamtools -L $(GUROBI)/lib/
LDLIBS = -lbamtools -lgurobi_c++ -lgurobi70 -lz
RPATH =  $(BAMTOOLS)/lib/bamtools:$(GUROBI)/lib/

SRCS = src/main.cpp src/ReadRec.cpp src/SegmentGraph.cpp src/WriteIO.cpp src/Config.cpp

all: bin/squid

bin/squid: $(subst .cpp,.o,$(SRCS))
        mkdir -p bin
        $(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS) -Wl,-rpath,$(RPATH)

clean:
        rm -f bin/squid src/*.o
