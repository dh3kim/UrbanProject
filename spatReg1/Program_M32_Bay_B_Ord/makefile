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

#All: Bay_SPanel_Alg1 Bay_SPanel_Alg2 LogDetCompare
All: Bay_SPanel_Alg1 

Bay_SPanel_Alg1:
	$(Compile) ; $(Execute)
Bay_SPanel_Alg2:
	$(Compile) ; $(Execute)
LogDetCompare:
	$(Compile) ; $(Execute)
clear:
	rm -f *.exe *.check
        
	
