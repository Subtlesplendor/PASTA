#Makefile for PTA library

#Choose C++ compiler:
CXX = g++
#Linking flags:
LDFLAGS=-lgsl -lgslcblas -lm -larb -lflint

#Debugging:
DEBUG = -g -pg
#Optimiziation options:
OPTIMIZATION = -march=native -O3

#Compiler flags
CXXFLAGS = -Wall -pedantic $(OPTIMIZATION) $(DEBUG)
#directory where user written programs are located:
VPATH = scripts

#---------------Handles the linking of the library-------------------
OBJDIR=lib
MODOBJDIR=$(OBJDIR)
LIBDIR=$(OBJDIR)
SRCDIR=bin
MODDIR=models
SOURCES= $(addprefix $(SRCDIR)/,*.cpp)
MODELS= $(addprefix $(MODDIR)/,*.cpp)
SCRIPTS = $(addprefix $(VPATH)/,*.cpp)
OBJECTS :=  $(notdir $(patsubst %.cpp,%.o,$(wildcard $(SOURCES))))
FINALS :=  $(notdir $(patsubst %.cpp,%,$(wildcard $(SCRIPTS))))
MODOBJECTS :=  $(notdir $(patsubst %.cpp,%.o,$(wildcard $(MODELS))))
LIB=libPTA.a
LDFLAGS+=-L./$(LIBDIR) -lPTA
LIBS=
#--------------------------------------------------------------------

all: cleanlib lib 

lib: $(addprefix $(LIBDIR)/, $(LIB))

$(addprefix $(LIBDIR)/, $(LIB)): $(addprefix $(OBJDIR)/, $(OBJECTS)) $(addprefix $(MODOBJDIR)/, $(MODOBJECTS))
	@ echo "Making library $@"
	@ ar rcs $@ $(addprefix $(OBJDIR)/, $(OBJECTS)) $(addprefix $(MODOBJDIR)/, $(MODOBJECTS))

$(OBJDIR)/%.o : $(addprefix $(SRCDIR)/, %.cpp) $(addprefix $(SRCDIR)/, %.h) $(addprefix $(SRCDIR)/, common.h)
	$(CXX) $(CXXFLAGS) -c $< -o $@  

$(MODOBJDIR)/%.o : $(addprefix $(MODDIR)/, %.cpp) $(addprefix $(MODDIR)/, %.h)
	$(CXX) $(CXXFLAGS) -c $< -o $@

%: %.cpp $(addprefix $(LIBDIR)/, $(LIB))
	@ $(CXX) $< -Isrc $(CXXFLAGS) $(LDFLAGS) $(addprefix $(LIBDIR)/, $(LIBS)) -o $@
    
cleanlib:
	@ echo "Cleaning library"
	@ rm -f $(addprefix $(OBJDIR)/, *.o)
	@ rm -f $(addprefix $(LIBDIR)/, $(LIB)) 


clean: cleanlib
	@ echo "Cleaning scripts"	
	@ rm -rf *.dSYM
	@ rm -f $(addprefix ./, $(FINALS))