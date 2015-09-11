#******************************************************************************
#
# H O P R   M A K E F I L E
#
#******************************************************************************

MAINDIR = ../
include ../Makefile.defs

SRC  = $(shell ls `cat srcfiles.mk`)
OBJf= $(SRC:.f90=.o)
OBJ = $(OBJf:.c=.o)

all: builddeps
	$(MAKE) $(PREPROC_LIB)
	$(FC) $(FCFLAGS) $(INCDIRS) -c hopr.f90 -o hopr.o
	@echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
	@echo ' HOPR COMPILED, NOW LINKING...'
	@echo '----------------------------------------------------------------------------------------'
	$(FC) $(FLFLAGS) hopr.o ../lib/$(PREPROC_LIB) -o ../bin/hopr \
	$(LIBDIRS) \
	$(EXTERNAL_LIBS)
	cp ../bin/hopr ../bin/preproctool
	@if [ -f $(HDF5_BASE)/bin/h5diff ] ; then \
	  cp $(HDF5_BASE)/bin/h5diff ../bin/. ; \
	fi
	@echo '----------------------------------------------------------------------------------------'
	@echo ' HOPR LINKED!'
	@echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

builddeps:
	@echo $(SRC) | sed -e 's/ /\n/g' > srcfiles_pp.mk
	@python ../tools/builddeps/builddeps --source=srcfiles_pp.mk deplist.mk

$(PREPROC_LIB): $(OBJ)
	$(AR) $(AR_FLAGS) ../lib/$(PREPROC_LIB) $(OBJ)

%.o: %.c
	$(CC) -c $< -o $@

%.o: %.f90
	$(FC) $(FCFLAGS) $(INCDIRS) -c $< -o $@


#------------------------------------------------------------------------------
# DEPENDENCIES
#------------------------------------------------------------------------------

include deplist.mk

#------------------------------------------------------------------------------
# UTILITY TARGETS
#------------------------------------------------------------------------------
.PHONY: clean veryclean

clean:
	rm -f $(OBJ) *.mod hopr.o *.i

veryclean: clean
	rm -f *~ */*~
