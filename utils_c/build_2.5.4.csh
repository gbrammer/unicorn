#!/bin/csh
setenv CC gcc-4.2
setenv MACOSX_DEPLOYMENT_TARGET 10.6
./setup_interp.py build_ext -i

#### Need to force compile with gcc-4.0 
gcc-4.2 -I/usr/stsci/pyssgx/Python-2.5.4/include -L/usr/stsci/pyssgx/Python-2.5.4/lib -O -fPIC -bundle -undefined dynamic_lookup build/temp.macosx-10.6-x86_64-2.5/interp_c.o -o interp_c.so

gcc-4.2 -I/usr/stsci/pyssgx/Python-2.5.4/include -L/usr/stsci/pyssgx/Python-2.5.4/lib -O -fPIC -bundle -undefined dynamic_lookup build/temp.macosx-10.6-x86_64-2.5/reduce_c.o -o reduce_c.so
