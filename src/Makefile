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
	$(FC) $(FLFLAGS) hopr.o $(PREPROC_LIB) -o ../bin/hopr \
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
	@echo "Use HDF5 in " $(HDF5_BASE)
	@echo $(SRC) | sed -e 's/ /\n/g' > srcfiles_pp.mk
	@python ../tools/builddeps/builddeps --source=srcfiles_pp.mk deplist.mk

$(PREPROC_LIB): $(OBJ)
	$(AR) $(AR_FLAGS) $(PREPROC_LIB) $(OBJ)

%.o: %.c
	$(CC) -c $< -o $@

%.o: %.f90
	$(FC) $(FCFLAGS) $(INCDIRS) -c $< -o $@

#------------------------------------------------------------------------------
# TAGS
#------------------------------------------------------------------------------
TAG = $(SRC:.f90=.tags)

tags: $(TAG)
	@rm -f ctags_unsorted.txt
	@find . -name "*.tags" -exec cat {} >> ctags_unsorted.txt \;
	@python -c "lines = open('ctags_unsorted.txt', 'r').readlines(); lines = sorted(lines); f = open('ctags_sorted.txt', 'w'); [f.write(l) for l in lines]; f.close()"
	@echo "!_TAG_FILE_FORMAT	2	/extended format; --format=1 will not append ;\" to lines/" > ctags.txt
	@echo "!_TAG_FILE_SORTED	1	/0=unsorted, 1=sorted, 2=foldcase/" >> ctags.txt
	@echo "!_TAG_PROGRAM_AUTHOR	Darren Hiebert	/dhiebert@users.sourceforge.net/" >> ctags.txt
	@echo "!_TAG_PROGRAM_NAME	Exuberant Ctags	//" >> ctags.txt
	@echo "!_TAG_PROGRAM_URL	http://ctags.sourceforge.net	/official site/" >> ctags.txt
	@echo "!_TAG_PROGRAM_VERSION	5.9~svn20110310	//" >> ctags.txt
	@cat ctags_sorted.txt >> ctags.txt
	@rm -f ctags_unsorted.txt ctags_sorted.txt
	@find . -name "*.tags" -exec rm {} \;


%.tags: %.f90
	@ctags -f $@ --language-force=fortran --fortran-kinds=+i+L $< 
	@sed -i '/^!_/d' $@



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
