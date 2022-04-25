#!/usr/bin/env python

# This sample setup.py can be used as a template for any project using d2to1.
# Simply copy this file and, if desired, delete all the comments.  Also remove
# the 'namespace_packages' and 'packages' arguments to setup.py if this project
# does not contain any packages belonging to a namespace package.

# This import statement attempts to import the setup() function from setuptools
# (this replaces the setup() one uses from distutils.core when using plain
# distutils).
#
# If the import fails (the user doesn't have the setuptools package) it then
# uses the ez_setup bootstrap script to install setuptools, then retries the
# import.  This is common practice for packages using setuptools.
try:
    from setuptools import setup
except ImportError:
    from ez_setup import use_setuptools
    use_setuptools()
    from setuptools import setup

# The standard setup() call.  Notice, however, that most of the arguments
# normally passed to setup() are absent.  They will instead be read from the
# setup.cfg file using d2to1.
#
# In order for this to work it is necessary to specify setup_requires=['d2to1']
# If the user does not have d2to1, this will boostrap it.  Also require
# stsci.distutils to use any of the common setup_hooks included in
# stsci.distutils (see the setup.cfg for more details).
#
# The next line, which defines namespace_packages and packages is only
# necessary if this projet contains a package belonging to the stsci namespace
# package, such as stsci.distutils or stsci.tools.  This is not necessary for
# projects with their own namespace, such as acstools or pyraf.
#
# d2to1=True is required to enable d2to1 and read the remaning project metadata
# from the setup.cfg file.
#
# use_2to3 and zip_safe are common options support by setuptools; these can
# also be placed in the setup.cfg, as will be demonstrated in a future update
# to this sample package.
setup(
    name='clusterpluck',
    version = '0.1',
    description = 'isolate and analyse Gaia DR2 data',
    long_description = 'README.rst',
    author = 'Marc Canale',
    author_email = 'marc.canale@hotmail.com',
    download_url = 'https://github.com/canalem/clusterpluck',
    classifiers=[
        'Intended Audience :: Science/Research/Education',
        'License :: None',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Astronomy',
        'Topic :: Software Development :: Libraries :: Python Modules'
    ],
    packages=['clusterpluck'],
    install_requires=[
          'numpy<1.23',
          'bokeh<2.4.0,>=2.3.0',
          'pandas',
          'astropy>=3.2.3',
          'matplotlib',
          'scipy',
          'seaborn'
          ],
    zip_safe = False,
)