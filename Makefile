
HOSTNAME = $(shell hostname)
USERNAME = $(shell whoami)

EXEC = ./kmst -f data/g01.dat -m scf -k 5

DEBUG = 0

CPLEX_VERSION = 1210
CPLEX_DIR = /Applications/CPLEX_Studio${CPLEX_VERSION}

SRCDIR = src
OBJDIR = obj
ARCH = x86-64
GPP = g++

CPPFLAGS = -DIL_STD \
	-isystem $(CPLEX_DIR)/cplex/include \
	-isystem $(CPLEX_DIR)/concert/include

CXXFLAGS += -Wall -Wno-non-virtual-dtor -pipe

LDFLAGS = -L$(CPLEX_DIR)/cplex/lib/$(ARCH)_osx/static_pic \
	-L$(CPLEX_DIR)/concert/lib/$(ARCH)_osx/static_pic

LDFLAGS += -lilocplex -lcplex -lconcert -lpthread -ldl

ifeq ($(DEBUG), 1)
  CXXFLAGS += -g -p
else
  CXXFLAGS += -O3
endif

STARTUP_SOURCE = $(SRCDIR)/Main.cpp

CPP_SOURCES = \
	src/Instance.cpp \
	src/kMST_ILP.cpp \
	src/Tools.cpp \
	src/CutCallback.cpp \
	src/Maxflow.cpp \

# $< the name of the related file that caused the action.
# $* the prefix shared by target and dependent files.
# $? is the names of the changed dependents.
# $@ is the name of the file to be made.
# $^ dependencies above


# ----- object files ---------------------------------------------------------------

OBJ_FILES = $(addprefix $(OBJDIR)/, $(patsubst %.cpp,%.o, \
	$(patsubst src/%, %, $(CPP_SOURCES) ) ) )
STARTUP_OBJ = $(addprefix $(OBJDIR)/, $(patsubst %.cpp,%.o, \
	$(patsubst src/%, %,$(STARTUP_SOURCE) ) ) )


all: kmst

depend:
	@echo 
	@echo "creating dependencies ..."
	$(GPP) -MM $(CPPFLAGS) $(CPP_SOURCES) $(SINGLE_FILE_SOURCES) \
	$(STARTUP_SOURCE) $(LD_FLAGS) \
	| sed -e "s/.*:/$(OBJDIR)\/&/" > depend.in

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp $(SRCDIR)/%.h
	@echo
	@echo "compiling $<"
	$(GPP) $(CPPFLAGS) $(CXXFLAGS) -o $@ -c $< 

# for startup (no header file available)
$(OBJDIR)/Main.o: $(SRCDIR)/Main.cpp
	@echo 
	@echo "compiling $<"
	$(GPP) $(CPPFLAGS) $(CXXFLAGS) -o $@ -c $< 

# ----- linking --------------------------------------------------------------------

kmst: $(STARTUP_OBJ) $(OBJ_FILES)	
	@echo 
	@echo "linking ..."
	@echo
	$(GPP) $(CPPFLAGS) $(CXXFLAGS) -o kmst $(OBJ_FILES) $(STARTUP_OBJ) $(LDFLAGS)

# ----- debugging and profiling ----------------------------------------------------

gdb: all
	gdb --args $(EXEC)

test: all
	$(EXEC)

clean:
	rm -rf obj/*.o kmst gmon.out

