TARGET := transpose.exe
SRC := shm.cc mpp.cc utils.cc transpose.cc

COMM := _none_
ALGO := _RINGONE_

CXX := sst++
CC := sstcc
CXXFLAGS := --std=c++11 -fPIC -D_PACKING_ -D_NOPAD_ -D_CELL_ -D$(COMM) -D$(ALGO) -D_SHM_ -D_EXTERNAL_ -D_STAT_ -D_COMPUTATION_ -D_OUTPUT_
CXXFLAGS := --std=c++11 -fPIC -D_PACKING_ -D_NOPAD_ -D_CELL_ -D$(COMM) -D$(ALGO) -D_SHM_ -D_EXTERNAL_ -D_NOMEM_
CPPFLAGS := -I. -I/usr/include/libxml2 -I$(FFTW)/include -I$(BOOST)/include
LIBDIR :=
PREFIX :=
LDFLAGS := -L$(FFTW)/lib -lfftw3 -lm -lxml2 -lrt

OBJ := $(SRC:.cc=.o)
OBJ := $(OBJ:.cpp=.o)
OBJ := $(OBJ:.c=.o)

.PHONY: clean install

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) -o $@ $+ $(LDFLAGS) $(LIBS)  $(CXXFLAGS)

%.o: %.cc
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

%.o: %.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(TARGET) $(OBJ)

install: $(TARGET)
	cp $< $(PREFIX)/bin
