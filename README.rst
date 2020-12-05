clusterpluck
============
Simple, user-friendly python package to isolate and analyse Gaia open cluster data.

Modules:

clusterpluck.search(ra, dec, radius)
  Return Gaia archive (https://gea.esac.esa.int/archive/) cone search (via astropy https://github.com/astropy/astropy TAP+ query) centered on inputted RA, DEC and radius
clusterpluck.search_name(name, radius)
  Return Gaia archive (https://gea.esac.esa.int/archive/) cone search (via astropy https://github.com/astropy/astropy TAP+ query) of inputted radius, centered on named cluster (cluster attributes downloaded from SIMBAD (http://simbad.u-strasbg.fr/simbad/), then filtered by distance and proper motion derived from cluster attributes
