
IDIR =./include
CC=mpic++
#CFLAGS=-I$(IDIR) -I/usr/include/gdal 
CFLAGS=-I$(IDIR) -w -g -fpermissive -I/usr/include/gdal 
#CFLAGS=-I$(IDIR) -O3 -fpermissive
ADIR=./apps
SDIR=./src

ODIR=./obj
LDIR =./lib

#LIBS=-lm -lgdal -lcurl -lboost_program_options -lboost_iostreams
LIBS=-lm -lboost_program_options -lboost_iostreams -lgdal


_APP_SRC = pp2g.cpp
APP_SRC = $(patsubst %,$(ADIR)/%,$(_APP_SRC))
APP_OBJ=  $(patsubst %.cpp,$(ODIR)/%.o,$(_APP_SRC))


_LIB_SRC = GridFile.cpp GridMap.cpp InCoreInterp.cpp Interpolation.cpp OutCoreInterp.cpp MpiInterp.cpp
LIB_SRC = $(patsubst %,$(SDIR)/%,$(_LIB_SRC))
LIB_OBJ=  $(patsubst %.cpp,$(ODIR)/%.o,$(_LIB_SRC))

SRC = $(APP_SRC) $(LIB_SRC)

#depend: .depend
#.depend: $(SRC)
#	rm -f ./.depend
#	$(CC) $(CFLAGS) -MM $^ -MF  ./.depend;
#include .depend


p_points2grid: $(APP_OBJ) $(LIB_OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)
	mkdir -p bin
	cp $@ bin
#cp $@ ~/p_points2grid_test

$(ODIR)/%.o: $(ADIR)/%.cpp
	$(CC) -c -o $@ $< $(CFLAGS)

$(ODIR)/%.o: $(SDIR)/%.cpp
	$(CC) -c -o $@ $< $(CFLAGS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o ./.depend pp2g core p_points2grid bin/p_points2grid
	rm -rf bin


