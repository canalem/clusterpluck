from sklearn.cluster import KMeans
from scipy.stats import gaussian_kde
import numpy as np


def k_means(ax):
    dfk = ax
    dfk.columns = ['x', 'y']
    kmeans = KMeans(n_clusters=1).fit(dfk)
    centroids = kmeans.cluster_centers_
    return centroids


def dp(dist, lq, uq):
    density = gaussian_kde(dist)
    xs = np.linspace(lq, uq, 1000)
    ys = density(xs)
    dist_index = np.argmax(ys)
    return xs[dist_index]
