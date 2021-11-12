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
r=Gaia.query_object_async(coordinate=coord,width=width, height=height)
r.pprint()
#%%
R= r.to_pandas()

# %%


R=R.drop(columns=['designation',
'astrometric_primary_flag','phot_variable_flag',
'datalink_url'])

R = R.fillna(0)
print(R)
#%%
R=R.iloc[:,0:92]
#%%

pca=PCA()
pca.fit(R)



# %%
