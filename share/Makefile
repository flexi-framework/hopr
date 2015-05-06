# set paths and compiler flags
MAINDIR = ../
include ../Makefile.defs

all: hdf5 cgns

dir:
	@# Create the directory for the share programs
	@ if [ ! -d $(SHARECOMP)/ ] ; then \
		mkdir $(SHARECOMP)/ ; \
	fi
hdf5: dir
	@# HDF5 libraries
	@ if [ -n "$(HDF5_VERSION)" ] && [ ! -f $(SHARECOMP)/hdf5-$(HDF5_VERSION)/hdf5/lib/libhdf5.a ] ; then \
		echo "===============================================================";\
		echo "= M A K E   H D F 5                                           =";\
		echo "===============================================================";\
		if [ ! -f hdf5-$(HDF5_VERSION).tar.bz2 ] ; then \
		  wget -N http://www.hdfgroup.org/ftp/HDF5/prev-releases/hdf5-$(HDF5_VERSION)/src/hdf5-$(HDF5_VERSION).tar.bz2 ; \
		fi ; \
		if [ ! -d $(SHARECOMP)/hdf5-$(HDF5_VERSION) ] ; then \
			cd $(SHARECOMP)/ && tar -xvjf ../hdf5-$(HDF5_VERSION).tar.bz2 1>/dev/null && cd .. ; \
		fi ; \
		cd $(SHARECOMP)/hdf5-$(HDF5_VERSION) && CC='$(CC)' FC='$(FC)' ./configure --enable-fortran --disable-shared && $(MAKE) && $(MAKE) install ; \
		echo "===============================================================";\
		echo "= H D F 5   D O N E                                           =";\
		echo "===============================================================";\
	fi
cgns: dir
	@# CGNS library (if 64bit needed: set --enable-64bit in configure command and CGNS_INT=64 in Makefile.defs
	@ if [ -n "$(CGNS_VERSION)" ] && [ ! -f $(CGNS_LIBDIR)/libcgns.a ] ; then \
#	@ if [ -n "$(CGNS_VERSION)" ] && [ ! -f $(SHARECOMP)/cgnslib_$(CGNS_VERSION)/src/lib/libcgns.a ] ; then \
		echo "===============================================================";\
		echo "= M A K E   C G N S                                           =";\
		echo "===============================================================";\
		if [ ! -f cgnslib_$(CGNS_VERSION)-$(CGNS_SUBVERSION).tar.gz ] ; then \
			wget -N http://mesh.dl.sourceforge.net/project/cgns/cgnslib_$(CGNS_VERSION)/Release%20$(CGNS_SUBVERSION)/cgnslib_$(CGNS_VERSION)-$(CGNS_SUBVERSION).tar.gz ; \
		fi ; \
		if [ ! -d $(SHARECOMP)/cgnslib_$(CGNS_VERSION) ] ; then \
			cd $(SHARECOMP) && tar -xvzf ../cgnslib_$(CGNS_VERSION)-$(CGNS_SUBVERSION).tar.gz 1>/dev/null && cd .. ; \
		fi ; \
 		cd $(CGNS_DIR) && FC='$(FC)' ./configure --enable-lfs  $(CGNSF2C) && $(MAKE) ; \
		echo "===============================================================";\
		echo "= C G N S   D O N E                                           =";\
		echo "===============================================================";\
	fi
recipes: dir
	@# NUMERICAL RECIPES
	@ if [ ! -f $(SHARECOMP)/RECIPES/librecipes.a ] ; then \
    echo "===============================================================";\
    echo "= M A K E   R E C I P E S                                     =";\
    echo "===============================================================";\
		if [ ! -d $(SHARECOMP)/RECIPES ] ; then \
			cd $(SHARECOMP) && tar -xvzf ../recipes.tar.gz 1>/dev/null && cd .. ; \
		fi ; \
		cd  $(SHARECOMP)/RECIPES && $(MAKE) ; \
		echo "===============================================================";\
		echo "= R E C I P E S   D O N E                                     =";\
		echo "===============================================================";\
	fi		
	

# utility targets
.PHONY: clean
clean:  cleanhdf5 cleancgns cleannr
cleanhdf5:
	rm -f *~
	@# HDF5 libraries
	@ if [ -d $(SHARECOMP)/hdf5-$(HDF5_VERSION) ] ; then \
		cd $(SHARECOMP)/hdf5-$(HDF5_VERSION) && $(MAKE) clean && $(MAKE) uninstall ; \
	fi
cleancgns:
	@# CGNS library
	@ if [ -d $(SHARECOMP)/cgnslib_$(CGNS_VERSION) ] ; then \
#		cd $(SHARECOMP)/cgnslib_$(CGNS_VERSION)/src && $(MAKE) allclean ; \
		cd $(CGNS_DIR) && $(MAKE) allclean ; \
	fi
cleannr:
	@# Numerical Recipes
	@ if [ -d $(SHARECOMP)/RECIPES ] ; then \
		cd $(SHARECOMP)/RECIPES && $(MAKE) clean ; \
	fi
	@ rm -f $(MAINDIR)lib/librecipes.a
