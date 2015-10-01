# set paths and compiler flags
include Makefile.defs

# "make" builds all
all: hopr
	@echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
	@echo ' SUCCESS: ALL EXECUTABLES GENERATED!'
	@echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

hopr: shared
	cd src && touch deplist.mk && $(MAKE)  #create deplist.mk for builddebs to prevent errors if not existing

hoprtools:
	cd tools/yplusestimator && $(MAKE) 

shared:
	cd share && $(MAKE) 

# utility targets
.PHONY: clean veryclean cleanshare

clean:
	cd src   && $(MAKE) clean
	cd tools/yplusestimator && $(MAKE) clean

veryclean:
	cd src   && $(MAKE) veryclean
	cd tools/yplusestimator && $(MAKE) veryclean
	rm -f lib/$(PREPROC_LIB)
	rm -f *~ */*~ */*/*~

cleanshare:
	cd share && $(MAKE) clean
