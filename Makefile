BOOST = /opt/local/stow/boost-1.55.0-gcc
BAMTOOLS = /opt/local/stow/bamtools-2.4.0
GLPK = /opt/local/stow/glpk-4.62

CC = gcc
CXX = g++
CXXFLAGS = -std=c++11 $(INCLUDES)
INCLUDES = -I $(BAMTOOLS)/include/bamtools -I $(GLPK)/include -I $(BOOST)
LDADD = $(BAMTOOLS)/lib/libbamtools.a $(GLPK)/lib/libglpk.a
LDLIBS = -lz -lm
#RPATH =  $(BAMTOOLS)/lib/:$(GLPK)/lib/

SRCS = src/main.cpp src/ReadRec.cpp src/SegmentGraph.cpp src/WriteIO.cpp src/Config.cpp

all: bin/squid

bin/squid: $(subst .cpp,.o,$(SRCS))
	mkdir -p bin
	$(CXX) -o $@ $^ $(LDADD) $(LDLIBS)

clean:
	rm -f bin/squid src/*.o
