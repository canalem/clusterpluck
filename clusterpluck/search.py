from astropy import units as u
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia
from astroquery.simbad import Simbad
import os
import pandas as pd
import numpy as np
import seaborn as sns

home_dir = str(os.getcwd())
download_dir = home_dir + r'\clusterpluck\data'
if not os.path.isdir(home_dir + r'\clusterpluck\data'):
    os.mkdir(r'clusterpluck\data')


def search(ra_coord, dec_coord, radius, d_near_plx=1000, d_far_plx=0.00001, pm_ra_min=-30, pm_ra_max=30,
           pm_dec_min=-30, pm_dec_max=30):
    query = """SELECT parallax, parallax_error, pmra, pmra_error, pmdec, pmdec_error, bp_rp, phot_g_mean_mag, ra, dec \
                FROM gaiadr2.gaia_source \
                WHERE CONTAINS(POINT('ICRS',gaiadr2.gaia_source.ra,gaiadr2.gaia_source.dec),CIRCLE('ICRS',{ra},{dec},{area}))=1 \
                AND bp_rp IS NOT NULL \
                AND parallax<{ig_near} AND parallax>{ig_far} \
                AND visibility_periods_used>8 \
                AND pmra<{ra_max} AND pmra>{ra_min} \
                AND pmdec<{dec_max} AND pmdec>{dec_min};""".format(ra=str(ra_coord), dec=str(dec_coord),
                                                                   area=str(radius), ig_near=str(d_near_plx),
                                                                   ig_far=str(d_far_plx),
                                                                   ra_min=str(pm_ra_min), ra_max=str(pm_ra_max),
                                                                   dec_min=str(pm_dec_min), dec_max=str(pm_dec_max))
    job1 = Gaia.launch_job_async(query, dump_to_file=True, output_format="csv",
                                 output_file=download_dir + '\output.csv')
    r = job1.get_results()
    if len(r) == 0:
        print('No values returned')
    else:
        print('Number of stars:', len(r))
    return r


def search_name(name, radius):
    customsimbad = Simbad()
    customsimbad.add_votable_fields('pmra', 'pmdec', 'plx', 'fe_h')
    result_table = customsimbad.query_object(name)
    raval = result_table['RA']
    decval = result_table['DEC']
    pmraval = result_table['PMRA']
    pmdecval = result_table['PMDEC']
    plxval = result_table['PLX_VALUE']
    coord = SkyCoord(ra=raval, dec=decval, unit=(u.hour, u.degree), frame='icrs')
    ra_coord = coord.ra.deg[0]
    dec_coord = coord.dec.deg[0]
    ra_pm = pmraval[0]
    dec_pm = pmdecval[0]
    pm_r = 2
    plx = plxval[0]
    d_near_plx, d_far_plx = distance_range(plx, radius)
    pm_ra_min, pm_ra_max, pm_dec_min, pm_dec_max = pm_range(ra_pm, dec_pm, pm_r)
    r = search(ra_coord, dec_coord, radius, d_near_plx, d_far_plx, pm_ra_min, pm_ra_max, pm_dec_min, pm_dec_max)
    return r


def distance_range(plx, radius):
    if type(plx) != np.float64:
        est_dist = (10 / (2 * (np.tan(radius) / 2) * 0.0174533))
        d_near_plx = 1 / ((est_dist / 5) / 1000)
        d_far_plx = 1 / ((est_dist * 3) / 1000)
    else:
        est_dist = (1 / plx) * 1000
        d_near_plx = 1 / ((est_dist - (est_dist / 2)) / 1000)
        d_far_plx = 1 / ((est_dist + (est_dist / 2)) / 1000)
    return d_near_plx, d_far_plx


def pm_range(ra_pm, dec_pm, pm_r):
    if type(ra_pm) != np.float64:
        ra_in = 0
        pm_r = 30
    if type(dec_pm) != np.float64:
        dec_in = 0
        pm_r = 30
    pm_ra_min = ra_pm - pm_r
    pm_ra_max = ra_pm + pm_r
    pm_dec_min = dec_pm - pm_r
    pm_dec_max = dec_pm + pm_r
    return pm_ra_min, pm_ra_max, pm_dec_min, pm_dec_max

class Plotting:

    def load(self):
        # Load complete table abd refine by removing stars mag > 18

        data = pd.read_csv(download_dir + '\output.csv')
        data = data.sort_values(by=['phot_g_mean_mag'])
        data = data.reset_index(drop=True)

        # Adjust parallax for median offset (Lindegren, et al. 2018)

        data['parallax'] = data['parallax'] + 0.029

        # Adjust parallax error to 2 sigma confidence level

        data['parallax_error'] = data['parallax_error'] * 2

        # Find trivial distance and error.

        data['distance'] = 1 / (data['parallax'] / 1000)
        data['distance_error'] = data['distance'] - (1 / ((data['parallax'] + data['parallax_error']) / 1000))

        # Calculate M_v Tycho (Jordi, et al. 2010)

        data['m_v_tycho'] = data['phot_g_mean_mag'] - (
                (0.02157 * (data['bp_rp'] ** 3)) - (0.2346 * (data['bp_rp'] ** 2)) - (
                    0.06629 * data['bp_rp']) - 0.01842)

        # Calculate (B - V) (Jordi, et al. 2010)

        data['b_v'] = 0.0043372 * (316.97 * np.sqrt(
            10046700000000 * data['bp_rp'] ** 2 - 43399787980000 * data['bp_rp'] + 701395108514151) + 1004670000 * data[
                                       'bp_rp'] - 2169989399) ** (1 / 3) - 17506 / (316.97 * np.sqrt(
            10046700000000 * data['bp_rp'] ** 2 - 43399787980000 * data['bp_rp'] + 701395108514151) + 1004670000 * data[
                                                                                        'bp_rp'] - 2169989399) ** (
                              1 / 3) + 1.4699

        # Save to CSV

        data.to_csv(download_dir + r'\complete.csv', index=None)

        return data

    def d_plot(self):
        f = sns.kdeplot(self['distance'])
        return f

    def cmd(self):
        f = self.plot.scatter(x='bp_rp', y='phot_g_mean_mag', s=1, alpha=0.25)
        f.invert_yaxis()
        return f

    def pm_plot(self):
        f = self.plot.scatter(x='pmra', y='pmdec', s=1, alpha=0.25)
        return f

    def pm_kde(self):
        f = sns.kdeplot(x=self['pmra'], y=self['pmdec'], fill=True, thresh=0, levels=100, cmap="mako")
        return f

