
IDIR =./include
CC=mpic++
#CFLAGS=-I$(IDIR) -I/usr/include/gdal 
CFLAGS=-I$(IDIR) -g
ADIR=./apps
SDIR=./src

ODIR=./obj
LDIR =./lib

#LIBS=-lm -lgdal -lcurl -lboost_program_options -lboost_iostreams
LIBS=-lm -lboost_program_options -lboost_iostreams


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


pp2g: $(APP_OBJ) $(LIB_OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

$(ODIR)/%.o: $(ADIR)/%.cpp
	$(CC) -c -o $@ $< $(CFLAGS)

$(ODIR)/%.o: $(SDIR)/%.cpp
	$(CC) -c -o $@ $< $(CFLAGS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o ./.depend pp2g core points2grid
	
