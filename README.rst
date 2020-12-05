clusterpluck
============
Simple, user-friendly python package to isolate and analyse Gaia open cluster data.

Modules:
clusterpluck.search(ra, dec, radius)
  Return Gaia cone search centered on declared RA, DEC and of declared radius
clusterpluck.search_name(name, radius)
  Return Gaia cone search centered on named cluster (cluster attributes downloaded from SIMBAD (http://simbad.u-strasbg.fr/simbad/), search then filtered by distance and proper motion derived from cluster attributes
