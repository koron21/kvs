#============================================================================
# General Option 
#=============================================================================

OUTDIR = ./bin
CC = g++
INCLUDE = -I$(KVS_DIR)/include -I../include
LIB = -L$(KVS_DIR)/lib -L../lib

#=============================================================================
# Flags, Definitions
#=============================================================================

CFLAGS = -c -Wall 
DEFINITIOINS = -DKVS_SUPPORT_GLEW -DKVS_SUPPORT_GLUT
ifdef DEBUG
CFLAGS += -O0 -ggdb3 
DEFINITIOINS += -D_DEBUG -DDEBUG -DKVS_ENABLE_DEBUG
else
CFLAGS += -O3
DEFINITIOINS += -DNDEBUG 
endif
LDFLAGS = -lCellLocator -lkvsSupportGLUT -lglut -lkvsSupportGLEW -lGLEW -lkvsCore -lGLU -lGL -lm -lpthread 

#=============================================================================
# Project Related Files
#=============================================================================

SOURCES := $(wildcard *.cpp) 
OBJECTS := $(SOURCES:.cpp=.o)
OBJECTS_WIN := $(SOURCES:.cpp=.obj)
EXECUTABLE = HyperStreamline
SLIB = lib$(EXECUTABLE).a
DLIB = lib$(EXECUTABLE).so
SLIB_WIN = $(EXECUTABLE).dll
DLIB_WIN = $(EXECUTABLE).lib

#=============================================================================
# Rules
#=============================================================================

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LIB) $(OBJECTS) -o $@  $(LDFLAGS)

.cpp.o: 
	$(CC) $(INCLUDE) $(CFLAGS) $(DEFINITIOINS) $< -o $@ 

clean:
	rm -rf $(OBJECTS) $(OBJECTS_WIN) $(EXECUTABLE) $(OUTDIR) $(SLIB) $(DLIB) $(SLIB_WIN) $(DLIB_WIN)

rebuild:
	make clean
	make

install:
	@if ! [ -d $(OUTDIR) ]; \
		then echo "making directory $(OUTDIR) for binary files" &&\
		mkdir $(OUTDIR); \
		fi
	mv $(OBJECTS) $(EXECUTABLE) $(OUTDIR)

lib: $(OBJECTS)
	$(AR) $(ARFLAGS) $(SLIB) $(OBJECTS)
