from astropy import units as u
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia
from astroquery.simbad import Simbad
import os
import pandas as pd
import numpy as np
import seaborn as sns
from .cp_clustering import k_means, dp

home_dir = str(os.getcwd())
download_dir = home_dir + r'\clusterpluck\data'
if not os.path.isdir(home_dir + r'\clusterpluck\data'):
    os.mkdir(r'clusterpluck\data')


def search(ra_input, dec_input, radius=0.5, pmra=0, pmdec=0, pm_r=10, d_near=1, d_far=10000, dr=2):
    """Carry out a cone search using position, radius and optionally, proper motion and distance filters. Converts
    from HH MM SS to decimal degrees. Converts distances to parallax. Performs asynchronous query then saves table as
    a CSV. Prints search parameters as well as len(table). Arguments = RA, Dec, Radius, PM RA centroid,
    PM Dec Centroid, PM Radius, Distance Near and Distance Far."""
    if d_near < 1:
        d_near = 1
    else:
        pass
    coord = SkyCoord(ra=ra_input, dec=dec_input, unit=(u.hour, u.degree))
    ra_coord = coord.ra.deg
    dec_coord = coord.dec.deg
    pm_ra_min = pmra - pm_r
    pm_ra_max = pmra + pm_r
    pm_dec_min = pmdec - pm_r
    pm_dec_max = pmdec + pm_r
    d_near_plx = (1 / d_near) * 1000
    d_far_plx = (1 / d_far) * 1000

    if dr != 3:
        gaia = 'gaiadr2.gaia_source'
    else:
        gaia = 'gaiaedr3.gaia_source'

    query = """SELECT parallax, parallax_error, pmra, pmra_error, pmdec, pmdec_error, bp_rp, phot_g_mean_mag, ra, dec \
                FROM {g} \
                WHERE CONTAINS(POINT('ICRS',{g}.ra,{g}.dec),CIRCLE('ICRS',{ra},{dec},{area}))=1 \
                AND bp_rp IS NOT NULL \
                AND parallax<{ig_near} AND parallax>{ig_far} \
                AND visibility_periods_used>8 \
                AND pmra<{ra_max} AND pmra>{ra_min} \
                AND pmdec<{dec_max} AND pmdec>{dec_min};""".format(g=gaia, ra=str(ra_coord), dec=str(dec_coord),
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
        print('RA:', ra_input, 'Dec:', dec_input, 'Rad:', radius)
        print('PM_RA:', pmra, 'PM_Dec:', pmdec, 'PM_Rad:', pm_r)
        print('Distance range: {:.0f} pc to {:.0f} pc'.format(d_near, d_far))


def search_name(name, radius=0.5, pm_r=2, dr=2):
    """Sends a cluster name, view radius and optional proper motion radius query to SIMBAD to download position,
    proper motion and distance information. Sends these to other functions to check for null values and so return
    defaults values. Once ready, sends values to search() function for Gaia query."""
    customsimbad = Simbad()
    customsimbad.add_votable_fields('pmra', 'pmdec', 'plx', 'fe_h')
    result_table = customsimbad.query_object(name)
    raval = result_table['RA']
    decval = result_table['DEC']
    pmraval = result_table['PMRA']
    pmdecval = result_table['PMDEC']
    plxval = result_table['PLX_VALUE']
    ra_coord = raval[0]
    dec_coord = decval[0]
    pmra = pmraval[0]
    pmdec = pmdecval[0]
    plx = plxval[0]
    d_near, d_far = distance_range(plx, radius)
    pmra, pmdec, pm_r = pm_range(pmra, pmdec, pm_r)
    search(ra_coord, dec_coord, radius, pmra, pmdec, pm_r, d_near, d_far, dr)


def distance_range(plx, radius):
    """If SIMBAD returns blank value for distance, calculates approximate distance from cluster diameter. Then
    produces distance range for search()."""
    if type(plx) != np.float64:
        est_dist = (10 / (2 * (np.tan(radius) / 2) * 0.0174533))
        d_near = (est_dist / 5)
        d_far = (est_dist * 3)
    else:
        est_dist = (1 / plx) * 1000
        d_near = est_dist - (est_dist / 2)
        d_far = est_dist + (est_dist / 2)
    return d_near, d_far


def pm_range(pmra, pmdec, pm_r):
    """If SIMBAD returns blank value for proper motion, calculates a default range for search()."""
    if type(pmra) != np.float64:
        pmra = 0
        pm_r = 30
    if type(pmdec) != np.float64:
        pmdec = 0
        pm_r = 30
    return pmra, pmdec, pm_r


class Refine:
    """Contains functions to aid more precise filtering of search()"""

    def load(self):
        """Loads initial output CSV from search() and adjusts/converts some values for plotting, etc."""
        # Load output table

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

        # Calculate cluster membership probability

        centroids_pm = k_means(data[['pmra', 'pmdec']])
        centroids_pos = k_means(data[['ra', 'dec']])

        sd_x = pd.DataFrame.std(data['ra'])
        sd_y = pd.DataFrame.std(data['dec'])
        sd_z = pd.DataFrame.std(data['distance'])
        sd_pmx = pd.DataFrame.std(data['pmra'])
        sd_pmy = pd.DataFrame.std(data['pmdec'])
        mean_d = data.loc[:, 'distance'].mean()
        data['probx'] = (data['ra'] - centroids_pos[:, 0])
        data['proby'] = (data['dec'] - centroids_pos[:, 1])
        data['probz'] = np.abs((data['distance']) - mean_d)
        data['probpmx'] = (data['pmra'] - centroids_pm[:, 0])
        data['probpmy'] = (data['pmdec'] - centroids_pm[:, 1])
        data['prob'] = (np.exp(-0.5 * ((data['probx'] / sd_x) ** 2 + (data['proby'] / sd_y) ** 2
                                       + (data['probpmx'] / sd_pmx) ** 2 + (data['probpmy'] / sd_pmy) ** 2
                                       + (data['probz'] / sd_z) ** 2))) * 100

        data.drop(columns=['probx', 'proby', 'probz', 'probpmx', 'probpmy'], inplace=True)

        # Save to CSV

        data.to_csv(download_dir + r'\complete.csv', index=None)

        return data

    def prob(self, pb):
        """Refine dataframe by membership probability"""
        df = self[self['prob'] > pb]
        df.to_csv(download_dir + r'\high_prob.csv', index=None)
        print('Number of stars with prob > {:.0f}%: {:.0f}'.format(pb, len(df)))

        return df

    def pm_plot(self):
        """Produces Pandas scatter plot of proper motion values."""
        self.plot.scatter(x='pmra', y='pmdec', s=1, alpha=0.25)

    def pm_kde(self):
        """Produces 2D Seaborn kernel distribution density plot of proper motion values."""
        sns.kdeplot(x=self['pmra'], y=self['pmdec'], fill=True, thresh=0, levels=100, cmap="mako")

    def map(self):
        """Produces Pandas scatter plot map of RA and Dec values. Size value is proportional to magnitude."""
        self.plot.scatter(x='ra', y='dec', s=(((1 / self['m_v_tycho']) * 20) ** 4), alpha=0.25, figsize=(8, 8))

    def d_plot(self):
        """Produces 1D Seaborn kernel distribution plot of distance values. Peak shows approx. cluster distance."""
        sns.kdeplot(self['distance'])


class Info:
    """Contains some astrophysical calculations."""

    def dist(self):
        """Uses gaussian kernel distribution of distance values and finds float value of peak as well as 5% and 95%
        quartiles."""
        cluster_d_p = dp(self['distance'], self.loc[:, 'distance'].quantile(0.05),
                         self.loc[:, 'distance'].quantile(0.95))
        print(r'Distance: {:.0f} pc'.format(cluster_d_p))
        print(r'5%: {:.0f} pc - 95%: {:.0f}'.format(self.loc[:, 'distance'].quantile(0.05),
                                                    self.loc[:, 'distance'].quantile(0.95)))

    def trgb(self, eb_v):
        """Estimate of distance via the Tip of the Reg Giant Branch standard candle (TRGB) from globular
        cluster plot. Currently just """
        df = self.nsmallest(5, 'm_v_tycho')
        mag = df['m_v_tycho'].mean()
        cluster_d_p = 10 ** ((mag + 3 + 5 - (3.1 * eb_v)) / 5)
        print(r'Distance: {:.0f} pc'.format(cluster_d_p))

    def radial_velocity(self):
        """Estimate of the radial velocity of the cluster center based on the ."""
        # Calculate radial velocity in km per sec
        r_v_mas = np.mean(np.sqrt(self['pmra'] ** 2 + self['pmdec'] ** 2))
        # Convert radial velocity from mas per year to rad per sec.
        r_v_rad = r_v_mas * 1.537289632e-16
        # Convert distance from pc to km
        cluster_d_p = dp(self['distance'], self.loc[:, 'distance'].quantile(0.05),
                         self.loc[:, 'distance'].quantile(0.95))
        d_km = cluster_d_p * 3.086e13
        # Perform final calculation
        r_v = d_km * (np.tan(r_v_rad / 2))
        # Calculate angle of rv
        centroids_pm = k_means(self[['pmra', 'pmdec']])
        centroids_pos = k_means(self[['ra', 'dec']])
        cpm = centroids_pm[0]
        cp = centroids_pos[0]
        r_v_theta = np.arctan(cpm[1] / cpm[0])
        r_v_theta = r_v_theta * 180 / np.pi
        if cpm[0] > 0:
            r_v_theta = 90 - r_v_theta
        else:
            r_v_theta = 270 - r_v_theta
        print(r'Radial Vel: %.1f km^s | Angle: %.0f deg' % (r_v, r_v_theta))


class Plotting:
    """Contains some further plotting functions."""

    def cmd(self):
        """Produces a Pandas scatter plot as a Colour Magnitude Diagram using original Gaia filters."""
        ax = self.plot.scatter(x='bp_rp', y='phot_g_mean_mag', s=3, alpha=0.25, figsize=(8, 8))
        ax.invert_yaxis()

    def cmd2(self):
        """Produces a Pandas scatter plot as a Colour Magnitude Diagram using values converted to more standardised
        filters (Tycho Vmag and Johnson B-V colour index)."""
        ax = self.plot.scatter(x='b_v', y='m_v_tycho', s=3, alpha=0.25, figsize=(8, 8))
        ax.invert_yaxis()
