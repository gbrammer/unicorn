# pygrism

October 9, 2015
---------------

This is the full 3D-HST grism pipeline used to produce the products described in
Momcheva et al. (http://arxiv.org/abs/1510.02106v1).  Note that it has many relative imports that
expect the module to be called "unicorn" (don't ask), so for now the easiest installation is to clone
the repository with that name.

Many additional utilities are needed from the legacy code at https://github.com/gbrammer/threedhst.

We will be working to clean up the two modules stripping out only the most useful utilities for a 
more portable pipeline distribution.

