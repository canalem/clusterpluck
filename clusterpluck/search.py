from astroquery.gaia import Gaia

def search(ra_coord, dec_coord, radius):
    query = """SELECT parallax, parallax_error, pmra, pmra_error, pmdec, pmdec_error, bp_rp, phot_g_mean_mag, ra, dec \
               FROM gaiadr2.gaia_source \
               WHERE CONTAINS(POINT('ICRS',gaiadr2.gaia_source.ra,gaiadr2.gaia_source.dec),CIRCLE('ICRS',{ra},{dec},{area}))=1 \
               AND bp_rp IS NOT NULL ;""".format(ra=str(ra_coord), dec=str(dec_coord), area=str(radius))
    job1 = Gaia.launch_job_async(query)
    search_results = job1.get_results()
    return search_results

