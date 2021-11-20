#%%
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.io import loadmat
from sklearn.decomposition import PCA



coord=SkyCoord(ra=280,dec=60, unit=(u.degree, u.degree), frame='icrs')
width=u.Quantity(0.1, u.deg)
height=u.Quantity(0.1, u.deg)
job = Gaia.launch_job_async("SELECT TOP 1000 gaia_source.paralla \
    x,gaia_source.phot_g_mean_mag,gaia_source.bp_rp \
        FROM gaiaedr3.gaia_source \
            WHERE CONTAINS(POINT('ICRS',gaiaedr3.gaia_source.ra,gaiaedr3.gaia_source.dec),CIRCLE('ICRS',0,0,180))=1")
print(job)
r=job.get_results()
r.pprint()
#%%
R= r.to_pandas()
#%%
X = R["bp_rp"]
Y = R["phot_g_mean_mag"] + 5*np.log10(R["parallax"]/100)
plt.scatter(X, Y)
plt.gca().invert_yaxis()
# %%
