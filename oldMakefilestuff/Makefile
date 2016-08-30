# Makefile to build the Software tools from R. Rietbroek

export ROOT_DIR=$(PWD)

# include general system parameters
include make.inc


all:mod libs prog #install

#make all module files	
mod:
	cd FORTtlbx;$(MAKE) mod
	cd GRDtlbx;$(MAKE) mod	
	cd SHtlbx;$(MAKE) mod
	cd GPStlbx;$(MAKE) mod

libs:
	cd FORTtlbx;$(MAKE) libs
	cd GRDtlbx;$(MAKE) libs	
	cd SHtlbx;$(MAKE) libs
	cd GPStlbx;$(MAKE) libs


prog:
	cd FORTtlbx;$(MAKE) prog
	cd GRDtlbx;$(MAKE) prog	
	cd SHtlbx;$(MAKE) prog
	cd GPStlbx;$(MAKE) prog


install:
	cd FORTtlbx;$(MAKE) install
	cd GRDtlbx;$(MAKE) install	
	cd SHtlbx;$(MAKE) install
	cd GPStlbx;$(MAKE) install
	cp $(LIBOUT) $(INSTALL)/lib/
	cp include/*mod $(INSTALL)/include/
#clean up goodies
clean:
	rm -f $(LIBOUT) include/* 
	cd FORTtlbx; $(MAKE) clean
	cd GRDtlbx; $(MAKE) clean
	cd SHtlbx; $(MAKE) clean
	cd GPStlbx; $(MAKE) clean


