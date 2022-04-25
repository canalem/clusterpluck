clusterpluck
============
Simple, user-friendly python project to isolate and analyse Gaia DR2 cluster data. I'm creating this as a hands on aid to learn Python as well as develop something that interests me and may be of use at some point.

In its current form it can be used to find, download and plot star data in the form of CSV table and CMDs, proper motion scatter plots and distance kernel density graphs. The estimated cluster distance is extracted based on a simple inversion of the parallax with a zero point error correction.

Currently working on making available as package for initial alpha release.

Packages Required
-----------------
- astropy_
- pandas
- numpy
- seaborn
- scipy
.. _astropy: https://github.com/astropy/astropy

Contents
--------

.. toctree::

   tutorials