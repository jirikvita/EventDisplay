######################################################
# A generic Makefile which does dependencies         #
# autmoatically                                      #
#                                                    #
# Original by Kurt Rinnert                           #
# Modified by Carl Gwilliam                          #
#                                                    #
# 1. Place in your top dir                           #
# 2. Change HEADERPAT and SRCPAT to your extension   #
#    naming (e.g. C cpp cxx h hh)                    #
# 3. Change LIBNAME to the wanted library name       #
# 4. Do 'gmake setup' to make dir structure and      #
#    copy files there:                               #
#    a. headers in include                           #
#    b. source files in src                          #
#    c. program files (containing main) in prg       #
#                                                    #
######################################################

#directories
LIBDIR = lib
BINDIR = bin
INCDIR = include
BOOSTDIR = /usr/local/boost_1_54_0/
SRCDIR = src
PRGDIR = prg
TMPDIR = tmp
BINOBJSDIR = $(TMPDIR)/bin
LIBOBJSDIR = $(TMPDIR)/lib
DEPSDIR    = $(TMPDIR)/deps

#paterns
HEADERPAT = hpp
SOURCEPAT = cpp

#source files
INCLUDEFILES = $(wildcard $(INCDIR)/*.$(HEADERPAT))
LIBCPPFILES  = $(wildcard $(SRCDIR)/*.$(SOURCEPAT))
PRGCPPFILES  = $(wildcard $(PRGDIR)/*.$(SOURCEPAT))

#targets
LIBNAME = EventDisplay
LIB     = $(LIBDIR)/lib$(LIBNAME).a
LIBOBJS = $(subst $(SRCDIR),$(LIBOBJSDIR),$(subst .$(SOURCEPAT),.o,$(LIBCPPFILES)))
BINS    = $(subst $(PRGDIR),$(BINDIR),$(subst .$(SOURCEPAT),,$(PRGCPPFILES)))
BINOBJS = $(subst $(PRGDIR),$(BINOBJSDIR),$(subst .$(SOURCEPAT),.o,$(PRGCPPFILES)))
LIBDEPS = $(subst $(SRCDIR),$(DEPSDIR),$(LIBCPPFILES:.$(SOURCEPAT)=.d))
BINDEPS = $(subst $(PRGDIR),$(DEPSDIR),$(PRGCPPFILES:.$(SOURCEPAT)=.d))

#tools
CPP = g++ -std=c++11
LD  = g++ -std=c++11
AR  = ar
DEP = gcc -std=c++11

#flags
#CPPBASEFLAGS = -c -Wall -Werror -ansi -pedantic -I$(INCDIR)
CPPBASEFLAGS =  -c -Wall -I$(INCDIR) -I$(BOOSTDIR)
DEPFLAGS = -MM -I$(INCDIR) -I$(shell root-config --incdir)

DEBUG=yes

ifeq ($(DEBUG),yes)
CPPDBGFLAGS = -g
else
CPPDBGFLAGS =
endif

ifeq ($(OPTIMIZE),yes)
CPPOPTFLAGS = -O2
else
CPPOPTFLAGS =
endif

ifeq ($(PROFILING),yes)
CPPPROFLAGS = -pg
else
CPPPROFLAGS =
endif

# orig:
# ROOTCPPFLAGS = $(shell root-config --cflags) -Wno-long-long 
ROOTCPPFLAGS = $(shell root-config --cflags) -Wno-long-long 

# default:
ROOTLDFLAGS  = $(shell root-config --libs) 

#EXTRALIBSPATH=${TestArea}/DataQuality/GoodRunsLists/StandAlone/GoodRunsListsLib.so
#EXTRALIBSPATH=./StandAlone/GoodRunsListsLib.so
#EXTRALIBSPATH=""

CPPFLAGS = $(CPPDBGFLAGS) $(CPPPROFLAGS) $(CPPOPTFLAGS) $(CPPBASEFLAGS) $(ROOTCPPFLAGS) -std=c++11

ARFLAGS = -rc

LDFLAGS = $(ROOTLDFLAGS) -lEG -lMinuit

ifeq ($(DEBUG),yes)
override LDFLAGS += -g 
endif
ifeq ($(PROFILING),yes)
override LDFLAGS += -pg 
endif

LOADLIBS = -L$(LIBDIR) -l$(LIBNAME)

# All depends on bin -> will "make" it if not up-to-date
all: bin
	echo "* all done."

# bin depends on static library (.a) for classes + binary objects for prgs (.o) + binarys for prgs
bin: $(LIB) $(BINOBJS) $(BINS) $(EXTRALIBSPATH)
# lib dependes on static library for classes (LIB).  Will "make" if not up-to-date
lib: $(LIB)

$(BINS):$(BINDIR)%:$(BINOBJSDIR)/%.o $(LIB)
	echo "* linking: $(@F)"
	$(LD) $(LDFLAGS)  $(EXTRALIBSPATH) -o $@ $< $(LOADLIBS)


# Static library (LIB) depends on library objects (.o) for each class (LIBOBJS).  
$(LIB):$(LIBOBJS)
	echo "* building library: $(@F)"
	$(AR) $(ARFLAGS) $@ $(LIBOBJS)

# Library object (.o) for each class (LIBOBJS) depends on corresponding source file
# (in srcdir with .o replaced by .SOURCEPAT) and will compile if not up-to-date  
$(LIBOBJS):$(LIBOBJSDIR)/%.o:$(SRCDIR)/%.$(SOURCEPAT) 
	echo "* compiling: $(<F)"
	$(CPP) $(CPPFLAGS) -o $@ $<

# Binary object for each prg (BINOBJS) depends on corresponding main file 
# (in prgdir with .o replaced by .SOURCEPAT) and will compile if no up-to-date 
$(BINOBJS):$(BINOBJSDIR)/%.o:$(PRGDIR)/%.$(SOURCEPAT)
	echo "* compiling: $(<F)"
	$(CPP) $(CPPFLAGS) -o $@ $<

$(LIBDEPS):$(DEPSDIR)/%.d:$(SRCDIR)/%.$(SOURCEPAT)
	echo "* creating dependencies: $(<F)"
	set -e; $(DEP) $(DEPFLAGS) $< | sed 's!\w.*\.o[ :]*!$(LIBOBJSDIR)/$*.o $@ : !' > $@;\
	[ -s $@ ] || rm -f $@

$(BINDEPS):$(DEPSDIR)/%.d:$(PRGDIR)/%.$(SOURCEPAT)
	echo "* creating dependencies: $(<F)"
	set -e; $(DEP) $(DEPFLAGS) $< | sed 's!\w.*\.o[ :]*!$(BINOBJSDIR)/$*.o $@ : !' > $@;\
	[ -s $@ ] || rm -f $@

# Clean up by removing all binaries, binary objects, libaries, library objects 
# and dependencies for both
clean:
	echo "* removing all targets"
	rm -f $(BINS) $(BINOBJS) $(LIB) $(LIBOBJS) $(BINDEPS) $(LIBDEPS)
	echo "* clean done"

# Make correct directory structure and copy files into it.  Files ending in:
# - .HEADERPAT go to INCDIR
# - .SRCPAT go to
#   - PRGDIR if contain main()
#   - SRCDIR otherwise
setup:
	echo "* creating directory structure"
	mkdir -p $(SRCDIR) $(PRGDIR) $(BINDIR) $(LIBDIR) $(INCDIR) $(TMPDIR)/deps $(TMPDIR)/lib $(TMPDIR)/bin
	echo "* Moving program files to $(PRGDIR)"
	mv $(shell grep -H "main()" *.$(SOURCEPAT) | cut -d : -f 1) $(PRGDIR)/.
	echo "* Moving header files to $(INDIR)"
	mv *.$(HEADERPAT) $(INCDIR)/.
	echo "* Moving source files to $(SRCDIR)"
	mv *.$(SOURCEPAT) $(SRCDIR)/.

dust:
	echo "* cough, cough"

ifneq ($(VERBOSE),yes)
.SILENT:
endif

.PHONY: all clean dust lib bin setup

# include the dependency files if the target does not contain "clean", "setup"
# the "-" suppresses warnings if files do not exist
ifeq (,$(findstring clean,$(MAKECMDGOALS)))
ifeq (,$(findstring setup,$(MAKECMDGOALS)))
-include $(BINDEPS)
-include $(LIBDEPS)
endif
endif


# MACROS:
#$@ Full name of the current target. 
#$? A list of files for current dependency which are out-of-date. 
#$< The source file of the current (single) dependency.

#* creating dependencies: IPBase.C
#* creating dependencies: IPAnal.C
#* creating dependencies: Main.C
#* creating dependencies: Fit.C
#echo "* compiling: IPAnal.C"
#g++    -c -Wall -Iinclude -pthread -m32 -I/batchsoft/atlas/athena/13.0.30_PCACHE/sw/lcg/external/root/5.14.00h/slc3_ia32_gcc323/root/include -Wno-long-long -o tmp/lib/IPAnal.o src/IPAnal.C
#echo "* compiling: IPBase.C"
#g++    -c -Wall -Iinclude -pthread -m32 -I/batchsoft/atlas/athena/13.0.30_PCACHE/sw/lcg/external/root/5.14.00h/slc3_ia32_gcc323/root/include -Wno-long-long -o tmp/lib/IPBase.o src/IPBase.C
#echo "* building library: libIP.a"
#ar -rc lib/libIP.a tmp/lib/IPAnal.o tmp/lib/IPBase.o
#echo "* compiling: Fit.C"
#g++    -c -Wall -Iinclude -pthread -m32 -I/batchsoft/atlas/athena/13.0.30_PCACHE/sw/lcg/external/root/5.14.00h/slc3_ia32_gcc323/root/include -Wno-long-long -o tmp/bin/Fit.o prg/Fit.C
#echo "* compiling: Main.C"
#g++    -c -Wall -Iinclude -pthread -m32 -I/batchsoft/atlas/athena/13.0.30_PCACHE/sw/lcg/external/root/5.14.00h/slc3_ia32_gcc323/root/include -Wno-long-long -o tmp/bin/Main.o prg/Main.C
#gmake: *** No rule to make target `tmp/bin//Fit.o', needed by `bin/Fit'.  Stop.
