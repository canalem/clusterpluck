from astropy import units as u
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia
from astroquery.simbad import Simbad
import os
import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt
from ipywidgets import interact, Layout, FloatSlider
import ipywidgets as widgets
import scipy.stats as stats
import pickle
from PyAstronomy import pyasl
from pylab import *
from matplotlib.pyplot import axes
from scipy.stats import gaussian_kde
from .cp_clustering import k_means, dp
from IPython.display import display, clear_output

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
                                 output_file=download_dir + r'\output.csv')
    r = job1.get_results()
    if len(r) == 0:
        print('No values returned')
    else:
        print('Number of stars:', len(r))
        print('RA:', ra_input, 'Dec:', dec_input, 'Rad:', radius)
        print('PM_RA:', pmra, 'PM_Dec:', pmdec, 'PM_Rad:', pm_r)
        print('Distance range: {:.0f} pc to {:.0f} pc'.format(d_near, d_far))

    stats = {'ra': ra_coord, 'dec': dec_coord, 'radius': radius, 'pmra': pmra, 'pmdec': pmdec,
             'pm_radius': pm_r, 'dr': dr}
    pickle_out = open("dict.pickle", "wb")
    pickle.dump(stats, pickle_out)
    pickle_out.close()


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
    fe_hvalue = result_table['Fe_H_Fe_H']
    ra_coord = raval[0]
    dec_coord = decval[0]
    pmra = pmraval[0]
    pmdec = pmdecval[0]
    plx = plxval[0]
    fe_h = fe_hvalue[0]
    d_near, d_far = distance_range(plx, radius)
    pmra, pmdec, pm_r = pm_range(pmra, pmdec, pm_r)

    # Initialise data and store values in pickle

    simbad_stats = {'name': name, 'ra': ra_coord, 'dec': dec_coord, 'radius': radius, 'pmra': pmra, 'pmdec': pmdec,
                    'pm_radius': pm_r, 'plx_value': plx, 'fe_h': fe_h, 'dr': dr}
    pickle_out = open("dict.pickle", "wb")
    pickle.dump(simbad_stats, pickle_out)
    pickle_out.close()

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


def load():
    """Loads initial output CSV from search() and adjusts/converts some values for plotting, etc."""
    # Load output table

    data = pd.read_csv(download_dir + r'\output.csv')
    data = data.sort_values(by=['phot_g_mean_mag']).reset_index(drop=True)
    pickle_in = open("dict.pickle", "rb")
    g_data = pickle.load(pickle_in)
    dr = g_data['dr']

    # Adjust parallax for median offset in data (Lindegren, et al. 2018, 2020)
    if dr == 2:
        data['parallax'] = data['parallax'] + 0.029
    else:
        data['parallax'] = data['parallax'] + 0.017

    # Adjust parallax error to 2 sigma confidence level

    data['parallax_error'] = data['parallax_error'] * 2

    # Find trivial distance and error.

    data['distance'] = 1 / (data['parallax'] / 1000)
    data['distance_error'] = data['distance'] - (1 / ((data['parallax'] + data['parallax_error']) / 1000))

    # Calculate M_v Tycho (Jordi, et al. 2010)

    if dr == 2:
        data['m_v_tycho'] = data['phot_g_mean_mag'] - (
                (0.02157 * (data['bp_rp'] ** 3)) - (0.2346 * (data['bp_rp'] ** 2)) - (
                0.06629 * data['bp_rp']) - 0.01842)
    else:
        data['m_v_tycho'] = data['phot_g_mean_mag'] - (
                (0.02342 * (data['bp_rp'] ** 3)) - (0.2387 * (data['bp_rp'] ** 2)) - (0.0682 * data['bp_rp']) - 0.01077)

    # Calculate (B - V) (Jordi, et al. 2010)

    data['b_v'] = 0.0043372 * (316.97 * np.sqrt(
        10046700000000 * data['bp_rp'] ** 2 - 43399787980000 * data['bp_rp'] + 701395108514151) + 1004670000 * data[
                                   'bp_rp'] - 2169989399) ** (1 / 3) - 17506 / (316.97 * np.sqrt(
        10046700000000 * data['bp_rp'] ** 2 - 43399787980000 * data['bp_rp'] + 701395108514151) + 1004670000 * data[
                                                                                    'bp_rp'] - 2169989399) ** (
                          1 / 3) + 1.4699

    # else: data['b_v'] = (5451 * np.sqrt(3) * np.sqrt( 891402030000000000 * data['bp_rp'] ** 2 - 2049490307788600000 *
    # data['bp_rp'] + 1494892818678820000) + 8914020300000 * data['bp_rp'] - 10247451538943) ** (1 / 3) / (5451 * 2
    # ** (2 / 3) * 5 ** (1/3)) - (41332604 * 2 ** (2 / 3) * 5 ** (1 / 3)) / (5451 * (5451 * np.sqrt(3)*np.sqrt(
    # 891402030000000000 * data['bp_rp'] ** 2 - 2049490307788600000 * data['bp_rp'] + 1494892818678820000) +
    # 8914020300000 * data['bp_rp'] - 10247451538943) ** (1 / 3)) + 6176 / 5451

    # Calculate absolute magnitudes, temperature and luminosity

    data['abs'] = (data['m_v_tycho'] + 5 + (5 * np.log10(data['parallax'] / 1000)))
    ball = pyasl.BallesterosBV_T()
    data['t_k'] = ball.bv2T(data['bp_rp'])
    data['lum_s'] = (3.85e26 * (10 ** ((4.77 - data['abs']) / 2.5))) / 3.38e26

    # Calculate some stats

    mu = np.mean(data['distance'])
    sigma = np.std(data['distance'])
    dist_med = np.median(data['distance'])
    conf = stats.norm.interval(0.95, loc=mu, scale=sigma)
    density = gaussian_kde(data['distance'])
    d_near_pc = data['distance'].min()
    d_far_pc = data['distance'].max()
    xs = np.linspace(d_near_pc, d_far_pc, 1000)
    ys = density(xs)
    dist_index = np.argmax(ys)
    cluster_d_p = xs[dist_index]

    # Calculate cluster membership probability

    centroids_pm = k_means(data[['pmra', 'pmdec']])
    centroids_pos = k_means(data[['ra', 'dec']])

    sd_x = pd.DataFrame.std(data['ra'])
    sd_y = pd.DataFrame.std(data['dec'])
    sd_z = pd.DataFrame.std(data['distance'])
    sd_pmx = pd.DataFrame.std(data['pmra'])
    sd_pmy = pd.DataFrame.std(data['pmdec'])
    data['probx'] = (data['ra'] - centroids_pos[:, 0])
    data['proby'] = (data['dec'] - centroids_pos[:, 1])
    data['probz'] = np.abs((data['distance']) - cluster_d_p)
    data['probpmx'] = (data['pmra'] - centroids_pm[:, 0])
    data['probpmy'] = (data['pmdec'] - centroids_pm[:, 1])
    data['prob'] = (np.exp(-0.5 * ((data['probx'] / sd_x) ** 2 + (data['proby'] / sd_y) ** 2
                                   + (data['probpmx'] / sd_pmx) ** 2 + (data['probpmy'] / sd_pmy) ** 2
                                   + (data['probz'] / sd_z) ** 2))) * 100

    data.drop(columns=['probx', 'proby', 'probz', 'probpmx', 'probpmy'], inplace=True)

    # SIMBAD!!! Build a SIMBAD database of objects in the FOV.

    pickle_in = open("dict.pickle", "rb")
    cluster_stats = pickle.load(pickle_in)
    ra_coord = cluster_stats['ra']
    dec_coord = cluster_stats['dec']
    radius = cluster_stats['radius']

    customSimbad = Simbad()
    customSimbad.add_votable_fields('otypes')

    result_table2 = customSimbad.query_region(SkyCoord(ra=ra_coord, dec=dec_coord,
                                                       unit=(u.deg, u.deg), frame='icrs'), radius * u.deg)
    raValue2 = result_table2['RA']
    decValue2 = result_table2['DEC']
    coord2 = SkyCoord(ra=raValue2, dec=decValue2, unit=(u.hour, u.degree), frame='icrs')
    ra_coord2 = coord2.ra.deg
    dec_coord2 = coord2.dec.deg
    dfid = pd.DataFrame(result_table2['MAIN_ID'], columns=['Name'], dtype="string")
    dfot = pd.DataFrame(result_table2['OTYPES'], columns=['Otypes'], dtype="string")
    dfra = pd.DataFrame(ra_coord2, columns=['ra'])
    dfdec = pd.DataFrame(dec_coord2, columns=['dec'])
    df_simbad = pd.concat([dfid, dfot, dfra, dfdec], axis=1, sort=False)

    data.insert(0, 'name', '', True)
    data.insert(1, 'otypes', '', True)

    for row in df_simbad.itertuples():
        for row2 in data.itertuples():
            if (row2.ra - 0.00015) < row.ra < (row2.ra + 0.00015) and (row2.dec - 0.00015) < row.dec < (row2.dec + 0.00015):
                data.at[row2[0], 'name'] = row.Name
                data.at[row2[0], 'otypes'] = row.Otypes

    # Initialise data and store in pickle

    hist_stats = {'mu': mu, 'sigma': sigma, 'dist_med': dist_med, 'conf': conf, 'cluster_d': cluster_d_p}
    pickle_out = open("dict.pickle", "wb")
    pickle.dump(hist_stats, pickle_out)
    pickle_out.close()

    # Save to CSV

    data.to_csv(download_dir + r'\complete.csv', index=None)

    return data


class Refine:
    """Contains functions to aid more precise filtering of search()"""

    def __init__(self):
        self.plot = None

    def prob(self, pb):
        """Refine dataframe by membership probability"""
        df = self[self['prob'] > pb]
        df.to_csv(download_dir + r'\high_prob.csv', index=None)
        print('Number of stars with prob > {:.0f}%: {:.0f}'.format(pb, len(df)))

        # Calculate some stats

        mu = np.mean(df['distance'])
        sigma = np.std(df['distance'])
        dist_med = np.median(df['distance'])
        conf = stats.norm.interval(0.95, loc=mu, scale=sigma)
        density = gaussian_kde(df['distance'])
        d_near_pc = df['distance'].min()
        d_far_pc = df['distance'].max()
        xs = np.linspace(d_near_pc, d_far_pc, 1000)
        ys = density(xs)
        dist_index = np.argmax(ys)
        cluster_d_p = xs[dist_index]

        # Initialise data and store in pickle

        hist_stats_prob = {'mu': mu, 'sigma': sigma, 'dist_med': dist_med, 'conf': conf, 'cluster_d': cluster_d_p}
        pickle_out = open("dict.pickle_prob", "wb")
        pickle.dump(hist_stats_prob, pickle_out)
        pickle_out.close()

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

    def d_hist(self):
        """Produces histogram of object distances derived from parallax."""
        fig = plt.figure()
        ax = fig.add_subplot(111)
        n, bins, patches = ax.hist(self['distance'], bins=50)
        nmax = np.max(n)
        arg_max = None
        for j, _n in enumerate(n):
            if _n == nmax:
                arg_max = j
                break


class Info:
    """Contains some astrophysical calculations."""

    def __init__(self):
        self.loc = None
        self.nsmallest = None

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
        cluster plot. Needs E(B-V) value - from 'The Globular Cluster Database'?. Currently just in beta."""
        df = self.nsmallest(5, 'm_v_tycho')
        mag = df['m_v_tycho'].mean()
        cluster_d_p = 10 ** ((mag + 3 + 5 - (3.1 * eb_v)) / 5)
        print(r'Distance: {:.0f} pc'.format(cluster_d_p))

    def eb_v(self, dist_x):
        """Estimate of extinction via the Tip of the Reg Giant Branch standard candle (TRGB) from globular
        cluster plot. Uses distance as measured from the parallax!!! Not accurate yet. Currently just in beta."""
        df = self.nsmallest(5, 'm_v_tycho')
        mag = df['m_v_tycho'].mean()
        eb_v = (- (5 * log10(dist_x) - (mag + 3 + 5))) / 3.1
        print(r'E(B - V): {:.2f}'.format(eb_v))

    def radial_velocity(self):
        """Estimate of the radial velocity of the cluster center."""
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
        cpm = centroids_pm[0]
        r_v_theta = np.arctan(cpm[1] / cpm[0])
        r_v_theta = r_v_theta * 180 / np.pi
        if cpm[0] > 0:
            r_v_theta = 90 - r_v_theta
        else:
            r_v_theta = 270 - r_v_theta
        print(r'Radial Vel: %.1f km^s | Angle: %.0f deg' % (r_v, r_v_theta))


class Plotting:
    """Contains some further plotting functions."""

    def __init__(self):
        self.plot = None

    def cmd(self):
        """Produces a Pandas scatter plot as a Colour Magnitude Diagram using original Gaia filters."""
        ax = self.plot.scatter(x='bp_rp', y='phot_g_mean_mag', s=3, alpha=0.25, figsize=(10, 10))
        ax.invert_yaxis()

    def cmd2(self):
        """Produces a Pandas scatter plot as a Colour Magnitude Diagram using values converted to more standardised
        filters (Tycho Vmag and Johnson B-V colour index)."""
        ax = self.plot.scatter(x='b_v', y='m_v_tycho', s=3, alpha=0.25, figsize=(10, 10))
        ax.invert_yaxis()

    def d_hist(self):
        """Produces histogram of object distances derived from parallax with information plotted such as mean, median
        and mode, etc."""
        pickle_in = open("dict.pickle", "rb")
        hist_stats = pickle.load(pickle_in)
        sigma = hist_stats['sigma']
        mu = hist_stats['mu']
        dist_med = hist_stats['dist_med']
        fig = plt.figure()
        ax = fig.add_subplot()
        n, bins, patches = ax.hist(self['distance'], bins=50)
        nmax = np.max(n)
        arg_max = None
        for j, _n in enumerate(n):
            if _n == nmax:
                arg_max = j
                break
        x_max = bins[arg_max]
        plt.close()
        fig, ax = plt.subplots(figsize=(12, 8))
        n, bins, patches = ax.hist(self['distance'], bins='auto', density=1, color='w', edgecolor='black')
        y = ((1 / (np.sqrt(2 * np.pi) * sigma)) * np.exp(-0.5 * (1 / sigma * (bins - mu)) ** 2))
        plt.axvline(mu, color='blue', linestyle='--', linewidth=1)
        plt.axvline(dist_med, color='green', linestyle=':', linewidth=1)
        plt.axvline(x_max, color='yellow', linestyle=':', linewidth=1)
        ax.plot(bins, y, color='blue', linestyle='--', linewidth=1)
        ax.set_xlabel('Distance / pc', fontsize=12)
        plt.yticks([])
        ax.set_title('Distance to Stars in Search Area', fontsize=18)
        fig.tight_layout()
        textstr = '\n'.join((
            r'$\mu=%.2f$' % (mu,),
            r'$\sigma=%.2f$' % (sigma,),
            r'median=%.2f' % (dist_med,),
            r'mode=%.2f' % (x_max,)))
        ax.text(0.85, 0.95, textstr, transform=ax.transAxes, fontsize=12, verticalalignment='top')

    def d_hist_ref(self):
        """Produces histogram of object distances derived from refined parallax with information plotted such as mean,
        median and mode, etc."""
        pickle_in = open("dict.pickle_prob", "rb")
        hist_stats_prob = pickle.load(pickle_in)
        sigma = hist_stats_prob['sigma']
        mu = hist_stats_prob['mu']
        dist_med = hist_stats_prob['dist_med']
        fig = plt.figure()
        ax = fig.add_subplot()
        n, bins, patches = ax.hist(self['distance'], bins=50)
        nmax = np.max(n)
        arg_max = None
        for j, _n in enumerate(n):
            if _n == nmax:
                arg_max = j
                break
        x_max = bins[arg_max]
        plt.close()
        fig, ax = plt.subplots(figsize=(12, 8))
        n, bins, patches = ax.hist(self['distance'], bins='auto', density=1, color='w', edgecolor='black')
        y = ((1 / (np.sqrt(2 * np.pi) * sigma)) * np.exp(-0.5 * (1 / sigma * (bins - mu)) ** 2))
        plt.axvline(mu, color='blue', linestyle='--', linewidth=1)
        plt.axvline(dist_med, color='green', linestyle=':', linewidth=1)
        plt.axvline(x_max, color='yellow', linestyle=':', linewidth=1)
        ax.plot(bins, y, color='blue', linestyle='--', linewidth=1)
        ax.set_xlabel('Distance / pc', fontsize=12)
        plt.yticks([])
        ax.set_title('Distance to Stars in Search Area', fontsize=18)
        fig.tight_layout()
        textstr = '\n'.join((
            r'$\mu=%.2f$' % (mu,),
            r'$\sigma=%.2f$' % (sigma,),
            r'median=%.2f' % (dist_med,),
            r'mode=%.2f' % (x_max,)))
        ax.text(0.85, 0.95, textstr, transform=ax.transAxes, fontsize=12, verticalalignment='top')

    def three_d(self):
        """Produces 3D plot of objects. If used with jupyter notebook, %matplotlib notebook is required for
        interactivity."""
        dot_size = ((np.log10((self['lum_s'] * 1000))) - 2) ** 4
        fig = plt.figure(figsize=(12, 12))
        ax = fig.add_subplot(1, 1, 1, projection="3d")
        x = self['ra']
        y = self['distance']
        z = self['dec']
        ax.scatter(x, y, z, marker='o', s=dot_size, c='blue')
        plt.show()

    def iso_match_alt(self):
        """Create adjustable CMD diagram to determine approximate cluster age and reddening using matplotlib.
        Requires %matplotlib notebook."""
        pickle_in = open("dict.pickle", "rb")
        iso_stats = pickle.load(pickle_in)
        fig9, ax_cmd = plt.subplots(figsize=(12, 10))
        iso = pd.read_csv(download_dir + r'\iso\iso_00.csv')
        plt.subplots_adjust(left=0.2, bottom=0.25)
        dm = (-5) + 5 * np.log10(iso_stats['cluster_d'])
        g = self['m_v_tycho']
        r = self['b_v']
        dminit = dm
        brinit = 0

        ax_g = axes([0.25, 0.05, 0.65, 0.03])
        ax_br = axes([0.25, 0.10, 0.65, 0.03])

        sdm = Slider(ax_g, 'Dist_Mod', 0, 20, valinit=dminit)
        sbr = Slider(ax_br, 'Reddening', 0, 1.0, valinit=brinit)
        plt.rcParams['lines.linewidth'] = 0.5
        plt.rcParams['lines.linestyle'] = '--'
        ax_cmd.set_ylim(-5, 9)
        ax_cmd.invert_yaxis()

        ax_cmd.plot(iso['B_V_100'], iso['V_100'], label=r'1x$10^{8}$ yr')
        ax_cmd.plot(iso['B_V_200'], iso['V_200'], c='Black', alpha=0.1)
        ax_cmd.plot(iso['B_V_300'], iso['V_300'], c='Black', alpha=0.1)
        ax_cmd.plot(iso['B_V_400'], iso['V_400'], c='Black', alpha=0.1)
        ax_cmd.plot(iso['B_V_500'], iso['V_500'], label=r'5x$10^{8}$ yr')
        ax_cmd.plot(iso['B_V_600'], iso['V_600'], c='Black', alpha=0.1)
        ax_cmd.plot(iso['B_V_700'], iso['V_700'], c='Black', alpha=0.1)
        ax_cmd.plot(iso['B_V_800'], iso['V_800'], c='Black', alpha=0.1)
        ax_cmd.plot(iso['B_V_900'], iso['V_900'], c='Black', alpha=0.1)
        ax_cmd.plot(iso['B_V_1000'], iso['V_1000'], label=r'1x$10^{9}$ yr')
        ax_cmd.scatter(r, g - dm, s=5, marker='s', c='w', edgecolor='black')
        leg = ax_cmd.legend()
        ax_cmd.set_ylabel(r'$M_{V}$')
        ax_cmd.set_xlabel(r'$B-V$')

        def update(val):
            ax_cmd.cla()
            ax_cmd.set_ylim(-5, 9)
            ax_cmd.invert_yaxis()
            ax_cmd.plot(iso['B_V_100'], iso['V_100'], label=r'1x$10^{8}$ yr')
            ax_cmd.plot(iso['B_V_200'], iso['V_200'], c='Black', alpha=0.1)
            ax_cmd.plot(iso['B_V_300'], iso['V_300'], c='Black', alpha=0.1)
            ax_cmd.plot(iso['B_V_400'], iso['V_400'], c='Black', alpha=0.1)
            ax_cmd.plot(iso['B_V_500'], iso['V_500'], label=r'5x$10^{8}$ yr')
            ax_cmd.plot(iso['B_V_600'], iso['V_600'], c='Black', alpha=0.1)
            ax_cmd.plot(iso['B_V_700'], iso['V_700'], c='Black', alpha=0.1)
            ax_cmd.plot(iso['B_V_800'], iso['V_800'], c='Black', alpha=0.1)
            ax_cmd.plot(iso['B_V_900'], iso['V_900'], c='Black', alpha=0.1)
            ax_cmd.plot(iso['B_V_1000'], iso['V_1000'], label=r'1x$10^{9}$ yr')
            ax_cmd.scatter(r - sbr.val, g - sdm.val, s=5, marker='s', c='w', edgecolor='black')
            ax_cmd.legend()
            ax_cmd.set_ylabel(r'$M_{V}$')
            ax_cmd.set_xlabel(r'$B-V$')
            draw()

        sdm.on_changed(update)
        sbr.on_changed(update)

        show()

    def iso_match(self):
        """Create adjustable CMD diagram to determine approximate cluster age and reddening using ipyWidgets. Can be
        used inline within a Jupyter notebook."""

        pickle_in = open("dict.pickle", "rb")
        iso_stats = pickle.load(pickle_in)
        fig9, ax_cmd = plt.subplots(figsize=(12, 10))
        iso = pd.read_csv(download_dir + r'\iso\iso_00.csv')
        plt.subplots_adjust(left=0.2, bottom=0.25)
        dm = (-5) + 5 * np.log10(iso_stats['cluster_d'])
        g = self['m_v_tycho']
        r = self['b_v']
        dminit = dm
        brinit = 0
        out = widgets.Output(layout=widgets.Layout())

        plt.rcParams['lines.linewidth'] = 0.5
        plt.rcParams['lines.linestyle'] = '--'
        ax_cmd.set_ylim(-5, 9)
        ax_cmd.invert_yaxis()
        plt.close()

        def update(sdm=dminit, sbr=brinit):
            ax_cmd.cla()
            ax_cmd.set_ylim(-5, 9)
            ax_cmd.invert_yaxis()
            ax_cmd.plot(iso['B_V_100'], iso['V_100'], label=r'1x$10^{8}$ yr')
            ax_cmd.plot(iso['B_V_200'], iso['V_200'], c='Black', alpha=0.1)
            ax_cmd.plot(iso['B_V_300'], iso['V_300'], c='Black', alpha=0.1)
            ax_cmd.plot(iso['B_V_400'], iso['V_400'], c='Black', alpha=0.1)
            ax_cmd.plot(iso['B_V_500'], iso['V_500'], label=r'5x$10^{8}$ yr')
            ax_cmd.plot(iso['B_V_600'], iso['V_600'], c='Black', alpha=0.1)
            ax_cmd.plot(iso['B_V_700'], iso['V_700'], c='Black', alpha=0.1)
            ax_cmd.plot(iso['B_V_800'], iso['V_800'], c='Black', alpha=0.1)
            ax_cmd.plot(iso['B_V_900'], iso['V_900'], c='Black', alpha=0.1)
            ax_cmd.plot(iso['B_V_1000'], iso['V_1000'], label=r'1x$10^{9}$ yr')
            ax_cmd.scatter(r - sbr, g - sdm, s=5, marker='s', c='w', edgecolor='black')
            ax_cmd.legend()
            ax_cmd.set_ylabel(r'$M_{V}$')
            ax_cmd.set_xlabel(r'$B-V$')
            out.clear_output(wait=True)
            display(fig9, out)

        interact(update, sdm=FloatSlider(value=dminit, min=0, max=20, step=0.1, layout=Layout(width='75%'),
                                         description='Dist_Mod:'),
                 sbr=FloatSlider(value=brinit, min=0, max=1, step=0.01, layout=Layout(width='75%'),
                                 description='Reddening:'))
