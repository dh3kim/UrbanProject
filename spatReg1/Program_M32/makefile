# Makefile.

# Set some variables.
Compiler = gfortran
Option   = -u -V -dcfuns -dusty -ieee=full 
Module   = /usr/src/LAPACK95/lapack95_modules
LAPACK95 = /usr/src/LAPACK95/lapack95.a
LAPACK   = /usr/lib/liblapack.a
BLAS     = /usr/lib/libblas-3.a
TMG      = /usr/lib/libtmglib.a 

# The symbol \ dontes line continuation.
# use the symbol $ to refer to a variable.
Compile = $(Compiler) $(Option) -o $@.exe -I$(Module) $@.f95 \
          $(LAPACK95) $(LAPACK) $(BLAS) $(TMG) 
Execute = ./$@.exe > $@.check

All: GMSPanelUnbalanced00 GMSPanelUnbalanced01 \
     GMSPanelUnbalanced02 GMSPanelUnbalanced10 \
     GMSPanelUnbalanced11 GMSPanelUnbalanced12 \
     GetTexForm

GMSPanelUnbalanced00:
	$(Compile) ; $(Execute)
GMSPanelUnbalanced01:
	$(Compile) ; $(Execute)
GMSPanelUnbalanced02:
	$(Compile) ; $(Execute)
GMSPanelUnbalanced10:
	$(Compile) ; $(Execute)
GMSPanelUnbalanced11:
	$(Compile) ; $(Execute)
GMSPanelUnbalanced12:
	$(Compile) ; $(Execute)
GetTexForm:
	$(Compiler) -o $@.exe $@.f95 ; ./$@.exe

clear:
	rm -f *.exe *.check
        
	
