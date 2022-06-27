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

All: bay_re_alg1 bay_re_alg2 bay_fe_alg1 \

bay_re_alg1:
	$(Compile) ; $(Execute)
bay_re_alg2:
	$(Compile) ; $(Execute)
bay_fe_alg1:
	$(Compile) ; $(Execute)
#GetTexForm:
#	$(Compiler) -o $@.exe $@.f95 ; ./$@.exe

clear:
	'rm' -f *.exe *.check
        
	
