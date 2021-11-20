
#%%
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.io import loadmat
from sklearn.decomposition import PCA
from sklearn.manifold import Isomap


coord=SkyCoord(ra=280,dec=60, unit=(u.degree, u.degree), frame='icrs')
width=u.Quantity(0.1, u.deg)
height=u.Quantity(0.1, u.deg)
r=Gaia.query_object_async(coordinate=coord,width=width, height=height)
r.pprint()
#%%
X=r.to_pandas()
R= r.to_pandas()

# %%
X=R.drop(columns=['designation',
'astrometric_primary_flag','phot_variable_flag',
'datalink_url'])
R=R.drop(columns=['designation',
'astrometric_primary_flag','phot_variable_flag',
'datalink_url'])
X = X.fillna(0)
R = R.fillna(0)
print(R)

#%%
X=X.iloc[:,0:92]
R=R.iloc[:,0:92]
#%%
#process of the data with the PCA
pca=PCA()
pca.fit(R)
R=pca.fit_transform(R)
plt.plot(np.cumsum(pca.explained_variance_ratio_))
plt.xlabel('n components')
plt.ylabel('cumulative variance')


# %%
V=pca.components_
# %%
#plotting some of the 
plt.plot(X['dist'],X['phot_bp_mean_flux'],'.')
plt.plot(X['dist'],X['phot_g_mean_flux'],'+')
plt.plot(X['dist'],X['phot_rp_mean_flux'],'*')
# %%
#seeing if there are some correlations between the phot mean data
plt.plot(X['phot_g_mean_flux'],X['phot_bp_mean_flux'],'.')
plt.plot(X['phot_g_mean_flux'],X['phot_rp_mean_flux'],'+')
plt.plot(X['phot_rp_mean_flux'],X['phot_bp_mean_flux'],'*')
