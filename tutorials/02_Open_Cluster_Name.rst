How to conduct a simple name search
===================================

After installing the clusterpluck package, this tutorial demonstrates
how to search for an open cluster using just its name.

First, import the module and classes.

.. code:: ipython3

    import clusterpluck as cp
    from clusterpluck.gaia import Refine, Info, Plotting


.. parsed-literal::

    Created TAP+ (v1.2.1) - Connection:
    	Host: gea.esac.esa.int
    	Use HTTPS: True
    	Port: 443
    	SSL Port: 443
    Created TAP+ (v1.2.1) - Connection:
    	Host: geadata.esac.esa.int
    	Use HTTPS: True
    	Port: 443
    	SSL Port: 443
    

Then, perform the search, downloading the basic cluster data from the
`SIMBAD Astronomical Database <http://simbad.u-strasbg.fr/simbad/>`__
and the individual star data from the `Gaia
archive <https://gea.esac.esa.int/archive/>`__.

The name of the cluster must be in a string and must also be in a
recognised format, i.e. M.. for Messier catalogue, NGC.. for New General
Catalogue, etc. It is possible that lesser known catalogues can be used
but they may cause an error if the format doesn’t match SIMBAD’s.

The results will be downloaded and stored as a CSV file.

.. code:: ipython3

    cp.search_name('M47')


.. parsed-literal::

    Number of stars: 994
    RA: 07 36 35 Dec: -14 29.0 Rad: 0.5
    PM_RA: -7.02 PM_Dec: 0.9592 PM_Rad: 2
    Distance range: 242 pc to 725 pc
    

Messier 47 is an open cluster in Cancer. The output of the search
contains the following information:

-  The number of stars downloaded from the Gaia database.
-  The right ascention (RA), declination (Dec) and search radius of the
   cluster. This is its position in the sky and it’s size.
-  The proper motion RA, Dec and search radius. This is the rate of
   cluster’s apparent movement across the sky. Cluster stars will all
   share approximately the same apparent movement and so will form a
   tight group when these data are plotted. In fact it is the main
   method used to identify cluster members.
-  Finally, the distance range of the search. This helps filter out lots
   of stars that are unrelated but can also cause cluster stars to be
   lost. In particular any objects further than 1 kpc (1000 pc) away can
   suffer from this.

Any or all of these can be amended by using the general ``search()``
method but that is for another tutorial.

Now let’s use the ``load()`` method to load the data into a Pandas
dataframe.

.. code:: ipython3

    t = cp.load()

We can check the dataframe using simple pandas commands.

.. code:: ipython3

    t.info()


.. parsed-literal::

    <class 'pandas.core.frame.DataFrame'>
    RangeIndex: 994 entries, 0 to 993
    Data columns (total 15 columns):
     #   Column           Non-Null Count  Dtype  
    ---  ------           --------------  -----  
     0   parallax         994 non-null    float64
     1   parallax_error   994 non-null    float64
     2   pmra             994 non-null    float64
     3   pmra_error       994 non-null    float64
     4   pmdec            994 non-null    float64
     5   pmdec_error      994 non-null    float64
     6   bp_rp            994 non-null    float64
     7   phot_g_mean_mag  994 non-null    float64
     8   ra               994 non-null    float64
     9   dec              994 non-null    float64
     10  distance         994 non-null    float64
     11  distance_error   994 non-null    float64
     12  m_v_tycho        994 non-null    float64
     13  b_v              994 non-null    float64
     14  prob             994 non-null    float64
    dtypes: float64(15)
    memory usage: 116.6 KB
    

The dataframe is in descending g (green) magnitude order.

.. code:: ipython3

    t




.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>parallax</th>
          <th>parallax_error</th>
          <th>pmra</th>
          <th>pmra_error</th>
          <th>pmdec</th>
          <th>pmdec_error</th>
          <th>bp_rp</th>
          <th>phot_g_mean_mag</th>
          <th>ra</th>
          <th>dec</th>
          <th>distance</th>
          <th>distance_error</th>
          <th>m_v_tycho</th>
          <th>b_v</th>
          <th>prob</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>2.077637</td>
          <td>0.207615</td>
          <td>-7.109151</td>
          <td>0.155182</td>
          <td>0.780030</td>
          <td>0.133099</td>
          <td>-0.118344</td>
          <td>5.627999</td>
          <td>114.016190</td>
          <td>-14.492770</td>
          <td>481.315959</td>
          <td>43.727498</td>
          <td>5.641895</td>
          <td>-0.150834</td>
          <td>80.099210</td>
        </tr>
        <tr>
          <th>1</th>
          <td>2.012069</td>
          <td>0.106906</td>
          <td>-7.138949</td>
          <td>0.080708</td>
          <td>1.145077</td>
          <td>0.083712</td>
          <td>0.081778</td>
          <td>6.473912</td>
          <td>114.411640</td>
          <td>-14.440757</td>
          <td>497.000887</td>
          <td>25.074589</td>
          <td>6.499310</td>
          <td>-0.011235</td>
          <td>38.982746</td>
        </tr>
        <tr>
          <th>2</th>
          <td>2.021288</td>
          <td>0.100697</td>
          <td>-7.310861</td>
          <td>0.086194</td>
          <td>0.785335</td>
          <td>0.074871</td>
          <td>-0.002212</td>
          <td>6.655011</td>
          <td>114.171852</td>
          <td>-14.443610</td>
          <td>494.734135</td>
          <td>23.477293</td>
          <td>6.673285</td>
          <td>-0.069918</td>
          <td>81.064037</td>
        </tr>
        <tr>
          <th>3</th>
          <td>1.743823</td>
          <td>0.067370</td>
          <td>-7.792587</td>
          <td>0.054587</td>
          <td>2.608254</td>
          <td>0.045973</td>
          <td>2.395475</td>
          <td>6.861200</td>
          <td>114.285469</td>
          <td>-14.324596</td>
          <td>573.452773</td>
          <td>21.330566</td>
          <td>8.088121</td>
          <td>1.639623</td>
          <td>0.560803</td>
        </tr>
        <tr>
          <th>4</th>
          <td>1.865338</td>
          <td>0.091669</td>
          <td>-7.078137</td>
          <td>0.070441</td>
          <td>1.084258</td>
          <td>0.063201</td>
          <td>-0.031709</td>
          <td>6.872384</td>
          <td>114.150463</td>
          <td>-14.484610</td>
          <td>536.095954</td>
          <td>25.111499</td>
          <td>6.888938</td>
          <td>-0.090495</td>
          <td>71.562474</td>
        </tr>
        <tr>
          <th>...</th>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
        </tr>
        <tr>
          <th>989</th>
          <td>2.101019</td>
          <td>2.697285</td>
          <td>-5.451604</td>
          <td>1.943831</td>
          <td>-0.055180</td>
          <td>2.318828</td>
          <td>1.512276</td>
          <td>20.853565</td>
          <td>114.108059</td>
          <td>-14.141460</td>
          <td>475.959477</td>
          <td>267.552516</td>
          <td>21.434158</td>
          <td>1.004453</td>
          <td>0.785147</td>
        </tr>
        <tr>
          <th>990</th>
          <td>3.006464</td>
          <td>3.076210</td>
          <td>-5.971694</td>
          <td>2.067749</td>
          <td>-0.149494</td>
          <td>2.225318</td>
          <td>1.086000</td>
          <td>20.883945</td>
          <td>114.298768</td>
          <td>-14.845272</td>
          <td>332.616656</td>
          <td>168.215278</td>
          <td>21.223415</td>
          <td>0.699234</td>
          <td>0.314587</td>
        </tr>
        <tr>
          <th>991</th>
          <td>3.919975</td>
          <td>4.900724</td>
          <td>-5.411766</td>
          <td>3.079569</td>
          <td>-0.184064</td>
          <td>4.386642</td>
          <td>1.416664</td>
          <td>20.885458</td>
          <td>114.524786</td>
          <td>-14.436349</td>
          <td>255.103698</td>
          <td>141.733991</td>
          <td>21.407289</td>
          <td>0.935868</td>
          <td>0.006988</td>
        </tr>
        <tr>
          <th>992</th>
          <td>4.023713</td>
          <td>3.740082</td>
          <td>-6.316390</td>
          <td>3.189567</td>
          <td>1.229143</td>
          <td>2.255890</td>
          <td>1.818579</td>
          <td>20.899426</td>
          <td>113.957239</td>
          <td>-14.413898</td>
          <td>248.526653</td>
          <td>119.723675</td>
          <td>21.684544</td>
          <td>1.224519</td>
          <td>0.507653</td>
        </tr>
        <tr>
          <th>993</th>
          <td>1.840975</td>
          <td>3.254297</td>
          <td>-5.187855</td>
          <td>2.186801</td>
          <td>2.466972</td>
          <td>2.374998</td>
          <td>3.004002</td>
          <td>20.947468</td>
          <td>114.296388</td>
          <td>-14.796842</td>
          <td>543.190356</td>
          <td>346.929973</td>
          <td>22.697336</td>
          <td>2.076619</td>
          <td>0.047302</td>
        </tr>
      </tbody>
    </table>
    <p>994 rows × 15 columns</p>
    </div>



Once the data are loaded to a variable we can check to see if we have a
cluster. Let’s see how the proper motion plot, ``pm_plot()``, looks.

.. code:: ipython3

    Refine.pm_plot(t)



.. image:: img/02_Open_Cluster_Name_11_0.png


We can see the cluster’s stars are forming a group right in the middle
of the plot! That means SIMBAD has given us good proper motion data and
the default proper motion radius is acceptable. Any closer and we would
lose relevent star data. Further away and there would be mnore chance of
unrelated stars included.

Now let’s look at a plot of the cluster as a star map.

.. code:: ipython3

    Refine.map(t)



.. image:: 02_Open_Cluster_Name_files%5C02_Open_Cluster_Name_13_0.png


This plot shows a star map of the search with the star size proportional
to their magnitude. It can help show us if our search radius is too wide
or narrow. This looks pretty good as the cluster appears fully
contained.

Now let’s see if the distance filter has correct values.

.. code:: ipython3

    Refine.d_plot(t)



.. image:: 02_Open_Cluster_Name_files%5C02_Open_Cluster_Name_15_0.png


There is a very clear, tall peak in the middle of our graph that tells
us the cluster stars are clearly outnumbering the unrelated stars. We
can also see roughly how far away the cluster is in parsecs just by
looking. The distance filter doesn’t need refining either. So that is
all the parameters involved with the search.

However let’s have a look at two of the features we can draw from this
data; a colour magnitude diagram and a more precise measurement of the
distance.

Using the ``Plotting`` class, we can call the ``cmd2()`` instance which
uses calculated values of apparent visual magnitude and the standardised
B-V colour index…

.. code:: ipython3

    Plotting.cmd2(t)



.. image:: 02_Open_Cluster_Name_files%5C02_Open_Cluster_Name_17_0.png


… and we have a beautiful CMD plot with the classic *main sequence* of
stars running from top left to bottom right. These stars are in the
middle of their lives, burning Helium in their cores in a relatively
stable way just like our sun. The thin line of stars above running
parallel to it are multiple star systems that have a slightly higher
luminosity.

Other features are the brightest stars at the top which appear to just
be ‘curling’ upwards. This is called the *main sequence turn off*. The
stars here are running low on core Helium and starting to evolve into
*red giants*. They’re not quite at that point but the position of the
turn off is a major method of ageing clusters. At the other end are the
red and white dwarfs.

Finally, use the ``Info`` class ``dist()`` instance to extract a precise
calculated distance from the parallax data including a 2-sigma range.
This distance isn’t to be used in a scientific context yet as it simply
uses an inverted parallax method to calculate.

.. code:: ipython3

    Info.dist(t)


.. parsed-literal::

    Distance: 474 pc
    5%: 341 pc - 95%: 633
    

We can also show the radial velocity of the cluster based on the
distance and proper motion.

.. code:: ipython3

    Info.radial_velocity(t)


.. parsed-literal::

    Radial Vel: 7.9 km^s | Angle: 278 deg
    

Next, we will look at how to refine the cluster data if name search
doesn’t quite give us what we need.
