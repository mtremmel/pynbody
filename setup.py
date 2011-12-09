import os
from distutils.core import setup, Extension
from distutils.sysconfig import get_python_lib

import numpy
import numpy.distutils.misc_util

libraries=[ ]
extra_compile_args = ['-ftree-vectorizer-verbose=1', '-ftree-vectorize',
                      '-fno-omit-frame-pointer',
                      '-funroll-loops',
                      '-fprefetch-loop-arrays',
                      '-fstrict-aliasing',
                      '-std=c99',
                      '-Wall',
                      '-O0']

extra_link_args = []

incdir = numpy.distutils.misc_util.get_numpy_include_dirs()
incdir.append('pynbody/pkdgrav2')
incdir.append('pynbody/pkdgrav2/mdl2/null')

#os.path.join(get_python_lib(plat_specific=1), 'numpy/core/include')
kdmain = Extension('pynbody/kdmain',
                   sources = ['pynbody/kdmain.c', 'pynbody/kd.c', 
                              'pynbody/smooth.c'],
                   include_dirs=incdir,
                   undef_macros=['DEBUG'],
                   libraries=libraries,
                   extra_compile_args=extra_compile_args,
                   extra_link_args=extra_link_args)

gravity = Extension('pynbody/pkdgrav',
                    sources = ['pynbody/gravity/pkdgravlink.c',
                               'pynbody/pkdgrav2/cl.c',
                               'pynbody/pkdgrav2/cosmo.c',
                               'pynbody/pkdgrav2/ewald.c',
                               'pynbody/pkdgrav2/fio.c',
                               'pynbody/pkdgrav2/grav2.c',
                               'pynbody/pkdgrav2/ilc.c',
                               'pynbody/pkdgrav2/ilp.c',
                               'pynbody/pkdgrav2/listcomp.c',
                               'pynbody/pkdgrav2/mdl2/null/mdl.c',
                               'pynbody/pkdgrav2/moments.c',
                               'pynbody/pkdgrav2/outtype.c',
                               'pynbody/pkdgrav2/pkd.c',
                               'pynbody/pkdgrav2/psd.c',
                               'pynbody/pkdgrav2/romberg.c',
                               'pynbody/pkdgrav2/smooth.c',
                               'pynbody/pkdgrav2/smoothfcn.c',
                               'pynbody/pkdgrav2/rbtree.c',
                               'pynbody/pkdgrav2/tree.c',
                               'pynbody/pkdgrav2/walk2.c'],
                   include_dirs=incdir,
                   undef_macros=['DEBUG','INSTRUMENT'],
                   define_macros=[('HAVE_CONFIG_H',None),
                                  ('__USE_BSD',None)],
                   libraries=libraries,
                   extra_compile_args=extra_compile_args,
                   extra_link_args=extra_link_args)

dist = setup(name = 'pynbody',
             author = '',
             author_email = '',
             version = '0.13alpha',
             description = '',
             package_dir = {'pynbody/': ''},
             packages = ['pynbody', 'pynbody/analysis', 'pynbody/bc_modules', 
                         'pynbody/plot', 'pynbody/gravity', 'examples' ],
# treat weave .c files like data files since weave takes
# care of their compilation for now
# could make a separate extension for them in future
             package_data={'pynbody': ['default_config.ini', 
                                       'sph_image.c','sph_to_grid.c'],
                           'pynbody/analysis': ['cmdlum.npz',
                                                'ionfracs.npz',
                                                'interpolate.c',
                                                'interpolate3d.c',
                                                'CAMB_WMAP7'],
                           'pynbody/plot': ['tollerud2008mw'],
                           'pynbody/gravity': ['direct.c']},
             ext_modules = [kdmain]
      )

#if dist.have_run.get('install'):
#    install = dist.get_command_obj('install')

