#!/bin/csh
#setenv CC "gcc-4.0"
#python2.5 setup_interp.py build_ext -i

cython -a interp_c.pyx

gcc-4.0 -fno-strict-aliasing -Wno-long-double -no-cpp-precomp -mno-fused-madd -DNDEBUG -I/usr/stsci/pyssg/Python-2.5.4/include -O -fPIC -I/usr/stsci/pyssg/Python-2.5.4/include/python2.5 -c interp_c.c -o build/temp.macosx-10.6-i386-2.5/interp_c.o

gcc-4.0 -I/usr/stsci/pyssg/Python-2.5.4/include -L/usr/stsci/pyssg/Python-2.5.4/lib -O -fPIC -bundle -undefined dynamic_lookup build/temp.macosx-10.6-i386-2.5/interp_c.o -o interp_c.so

cython -a reduce_c.pyx

gcc-4.0 -fno-strict-aliasing -Wno-long-double -no-cpp-precomp -mno-fused-madd -DNDEBUG -I/usr/stsci/pyssg/Python-2.5.4/include -O -fPIC -I/usr/stsci/pyssg/Python-2.5.4/include/python2.5 -c reduce_c.c -o build/temp.macosx-10.6-i386-2.5/reduce_c.o

gcc-4.0 -I/usr/stsci/pyssg/Python-2.5.4/include -L/usr/stsci/pyssg/Python-2.5.4/lib -O -fPIC -bundle -undefined dynamic_lookup build/temp.macosx-10.6-i386-2.5/reduce_c.o -o reduce_c.so

