TARGET := shm_mgr.x
SRC := shm.cc mpp.cc utils.cc shm_mgr.cc

CXX :=    smpicxx
CC :=     smpicc
CXXFLAGS := -D_CELL_
CPPFLAGS := -I. -I/usr/include/libxml2 -I$(FFTW)/include
LIBDIR :=  
PREFIX := 
LDFLAGS := -L$(FFTW)/lib -lfftw3 -lm -lxml2 -lrt

OBJ := $(SRC:.cc=.o) 
OBJ := $(OBJ:.cpp=.o)
OBJ := $(OBJ:.c=.o)

.PHONY: clean install 

all: $(TARGET)

shm_chk.x: shm_chk.o 
	$(CXX) -o $@ shm_chk.o $(LDFLAGS) $(LIBS)  $(CXXFLAGS)
$(TARGET): $(OBJ) 
	$(CXX) -o $@ $+ $(LDFLAGS) $(LIBS)  $(CXXFLAGS)

%.o: %.cc 
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

%.o: %.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

clean: 
	rm -f $(TARGET) $(OBJ) shm_chk.o shm_chk.x 

install: $(TARGET)
	cp $< $(PREFIX)/bin

