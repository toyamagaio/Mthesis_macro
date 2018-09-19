CC=gcc
CXX=g++
#CC=icc
#CXX=icpc
#CFLAGS  = -parallel -par-report3 -par-threshold0 -O3

CFLAGS  = -O2

BINDIR = ./bin
LIBDIR = ./lib

ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs)
ROOTGLIBS = $(shell root-config --glibs)
CXXFLAGS = -Wall -O2 $(ROOTFLAGS) 
CXXLIBS = $(ROOTLIBS)

TARGET1=     H3L_life_data
OBJS1=       H3L_life_data.o
TARGET2=     Recoil_mom
OBJS2=       Recoil_mom.o
TARGET3=     cherenkov
OBJS3=       cherenkov.o
TARGET4=     scintiemit_MPPCeff
OBJS4=       scintiemit_MPPCeff.o
TARGET5=     mppc_spec
OBJS5=       mppc_spec.o
TARGET6=     dedx
OBJS6=       dedx.o
TARGET7=     beta_energy
OBJS7=       beta_energy.o
TARGET8=     draw_wave_spice
OBJS8=       draw_wave_spice.o
TARGET9=     proton_range
OBJS9=       proton_range.o

all: $(TARGET1) \
     $(TARGET2) \
     $(TARGET3) \
     $(TARGET4) \
     $(TARGET5) \
     $(TARGET6) \
     $(TARGET7) \
     $(TARGET8) \
     $(TARGET9) \

$(LIBDIR)/%.o : %.cc
	$(CXX) $(CFLAGS) -c -o $@ $< $(CXXFLAGS)

$(TARGET1): $(patsubst %,$(LIBDIR)/%,$(OBJS1))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $< $(CXXLIBS) $(CXXFLAGS)

$(TARGET2): $(patsubst %,$(LIBDIR)/%,$(OBJS2))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $< $(CXXLIBS) $(CXXFLAGS)

$(TARGET3): $(patsubst %,$(LIBDIR)/%,$(OBJS3))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $< $(CXXLIBS) $(CXXFLAGS)

$(TARGET4): $(patsubst %,$(LIBDIR)/%,$(OBJS4))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $< $(CXXLIBS) $(CXXFLAGS)

$(TARGET5): $(patsubst %,$(LIBDIR)/%,$(OBJS5))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $< $(CXXLIBS) $(CXXFLAGS)

$(TARGET6): $(patsubst %,$(LIBDIR)/%,$(OBJS6))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $< $(CXXLIBS) $(CXXFLAGS)

$(TARGET7): $(patsubst %,$(LIBDIR)/%,$(OBJS7))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $< $(CXXLIBS) $(CXXFLAGS)

$(TARGET8): $(patsubst %,$(LIBDIR)/%,$(OBJS8))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $< $(CXXLIBS) $(CXXFLAGS)

$(TARGET9): $(patsubst %,$(LIBDIR)/%,$(OBJS9))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $< $(CXXLIBS) $(CXXFLAGS)

.PHONY: clean
clean:
	rm -f $(LIBDIR)/*.o core $(BINDIR)/*

