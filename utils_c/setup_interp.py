#!/usr/stsci/pyssgx/Python-2.7/bin/python
# /usr/stsci/pyssg/Python-2.5.4/bin/python
"""
Use distutils to compile Cython routines with all of the appropriate 
CFLAGS and libraries.

(Thanks to "Rob Wolfe": http://bytes.com/topic/python/answers/839099-cython-dynamic-library-problem)

"""
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext as build_pyx

import os, sys
import glob

#### Generate HTML files for analysis
files=glob.glob('*.pyx')
for file in files:
    root='%s' %file.split('.pyx')[0]
    os.system('cython -a %s' %file)  ## to generate html diagnostic    
    setup(name = root, ext_modules=[Extension(root, [file])], cmdclass = { 'build_ext': build_pyx })
    
# os.system('cython -a reduce_c.pyx')
# 
# setup(name = 'interp_c',
# ext_modules=[Extension('interp_c', ['interp_c.pyx'])],
# cmdclass = { 'build_ext': build_pyx })
# 
# setup(name = 'reduce_c',
# ext_modules=[Extension('reduce_c', ['reduce_c.pyx'])],
# cmdclass = { 'build_ext': build_pyx })

