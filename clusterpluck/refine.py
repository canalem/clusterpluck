import os
import pandas as pd
import numpy as np
import seaborn as sns
from bokeh.plotting import figure, output_file, show
from bokeh.models import ColumnDataSource
from bokeh.models.tools import HoverTool
from bokeh.palettes import Spectral5
from bokeh.transform import factor_cmap

download_dir = str(os.getcwd() + r'\clusterpluck\data')


class Refine:
    None
