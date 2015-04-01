############################# VARIABLES ################################

PYVERSION = 2.7

CY = cython
	# Cython compilater

CC = gcc
	# C Compiler

CFLAGS = -pthread -fno-strict-aliasing -DNDEBUG -g -fwrapv -O2 -Wall -Wstrict-prototypes -fPIC -I/usr/include/python$(PYVERSION)
	# Compilation options

LFLAGS = -pthread -shared -Wl,-O1 -Wl,-Bsymbolic-functions -Wl,-Bsymbolic-functions -Wl,-z,relro -fno-strict-aliasing -DNDEBUG -g -fwrapv -O2 -Wall -Wstrict-prototypes -D_FORTIFY_SOURCE=2 -g -fstack-protector --param=ssp-buffer-size=4 -Wformat -Werror=format-security
	# link edition options

CYFLAGS =
	# cython compiler options

#################### INSTRUCTIONS DE COMPILATION #######################
# $@ =  Target # $^ = list of dependencies # $< First dependency #

all: AdapterTrimmer.so

AdapterTrimmer.so: AdapterTrimmer.o
	#Link editing
	$(CC) $(LFLAGS) $^ -o $@
	
AdapterTrimmer.o: AdapterTrimmer.c
	#Compilation of source object
	$(CC) $(CFLAGS) -c $^ -o $@ 

AdapterTrimmer.c: AdapterTrimmer.pyx
	#Cythonizing AdapterTrimmer.pyx
	$(CY) $(CYFLAGS) $^ -o $@ 


##################### INSTRUCTIONS DE NETTOYAGE ########################

.PHONY: clean mrproper

clean:
	#Clean intermediate files
	rm -rf AdapterTrimmer.o AdapterTrimmer.c AdapterTrimmer.html ./build/ *.pyc

mrproper:
	#Clean everything
	rm -rf AdapterTrimmer.o AdapterTrimmer.c AdapterTrimmer.html AdapterTrimmer.so ./build/ *.pyc
