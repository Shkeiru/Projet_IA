#%%
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.io import loadmat
from sklearn.decomposition import PCA
import seaborn as sns



coord=SkyCoord(ra=280,dec=60, unit=(u.degree, u.degree), frame='icrs')
width=u.Quantity(0.1, u.deg)
height=u.Quantity(0.1, u.deg)
job = Gaia.launch_job_async("SELECT ra, dec, parallax, phot_g_mean_mag, bp_rp, \
astrometric_excess_noise, \
phot_bp_rp_excess_factor \
FROM gaiaedr3.gaia_source \
WHERE parallax > 10 \
AND parallax_over_error > 10 \
AND phot_bp_mean_flux_over_error > 10 \
AND phot_rp_mean_flux_over_error > 10 \
AND astrometric_excess_noise < 1 ")
print(job)
r=job.get_results()
r.pprint()
#%%
R= r.to_pandas()
ext = R["astrometric_excess_noise"].max()-R["astrometric_excess_noise"].min()
#%%
fig = plt.figure(figsize=(12,8))
X = R["bp_rp"]
Y = R["phot_g_mean_mag"] #+ 5*np.log10(R["parallax"]/100)
sns.relplot(x=X, y=-Y, hue=R["astrometric_excess_noise"], size=0.1)
plt.savefig("./fig.png")
plt.show()

# %%
