
# coding: utf-8

# In[1]:

get_ipython().magic('cat 0Source_Citation.txt')


# In[2]:

get_ipython().magic('matplotlib inline')
# %matplotlib notebook # for interactive


# For high dpi displays.

# In[3]:

get_ipython().magic("config InlineBackend.figure_format = 'retina'")


# # 0. General note

# - This notebook shows an example of how to include anharmonic effects and/or electronic effects in the P-V-T EOS fitting.
# 
# - We use data on SiC from [Nisr et al. (2017, JGR-Planet)](http://onlinelibrary.wiley.com/doi/10.1002/2016JE005158/full).
# 
# - Note that the code here is for demonstration purpose.  In fact, we did not find any evidence for including anharmonic and electronic effects in fitting the data for SiC.  

# # 1. Global setup

# In[4]:

import numpy as np
import matplotlib.pyplot as plt
import uncertainties as uct
from uncertainties import unumpy as unp
import pytheos as eos


# # 2. Setup for fitting with different gold pressure scales

# The equations of state of gold from Fei et al. (2007, PNAS) and Dorogokupets and Dewaele (2007, HPR).  These equations are provided in `pytheos` as built-in classes.

# In[5]:

au_eos = {'Fei2007': eos.gold.Fei2007bm3(), 'Dorogokupets2007': eos.gold.Dorogokupets2007()}


# Because we use Birch-Murnaghan EOS version of Fei2007 and Dorogokupets2007 used Vinet EOS, we create a dictionary to provide different static compression EOSs for the different pressure scales used.

# In[6]:

st_model = {'Fei2007': eos.BM3Model, 'Dorogokupets2007': eos.VinetModel}


# Assign initial values for the EOS parameters.

# In[7]:

k0_3c = {'Fei2007': 241.2, 'Dorogokupets2007': 243.0}
k0p_3c = {'Fei2007': 2.84, 'Dorogokupets2007': 2.68}
k0_6h = {'Fei2007': 243.1, 'Dorogokupets2007': 245.5}
k0p_6h = {'Fei2007': 2.79, 'Dorogokupets2007': 2.59}


# Also for the thermal parameters.  In this example, we will use the constant $q$ equation for the thermal part of the EOS.

# In[8]:

gamma0 = 1.06
q = 1.
theta0 = 1200.


# We also setup for the physical constants of two different polymorphs of SiC.

# In[9]:

v0 = {'3C': 82.8042, '6H': 124.27}
n_3c = 2.; z_3c = 4.
n_6h = 2.; z_6h = 6.


# # 3. Data distribution (3C)

# The data set is provided in a `csv` file.

# In[10]:

data = np.recfromcsv('./data/3C-HiTEOS-final.csv', case_sensitive=True, deletechars='')


# Set up variables for the data.

# In[11]:

v_std = unp.uarray( data['V(Au)'], data['sV(Au)'])
temp = unp.uarray(data['T(3C)'], data['sT(3C)'])
v = unp.uarray(data['V(3C)'], data['sV(3C)'])


# Plot $P$-$V$-$T$ data in the $P$-$V$ and $P$-$T$ spaces.

# In[12]:

for key, value in au_eos.items():
    p = au_eos[key].cal_p(v_std, temp)
    eos.plot.thermal_data({'p': p, 'v': v, 'temp': temp}, title=key)


# # 4. Data fitting with constq equation (3C) with electronic effects

# Normally weight for each data point can be calculated from $\sigma(P)$.  In this case, using `uncertainties`, we can easily propagate the temperature and volume uncertainties to get the value.

# In[13]:

for key, value in au_eos.items():
    # calculate pressure
    p = au_eos[key].cal_p(v_std, temp)
    # add prefix to the parameters.  this is important to distinguish thermal and static parameters
    eos_st = st_model[key](prefix='st_') 
    eos_th = eos.ConstqModel(n_3c, z_3c, prefix='th_')
    eos_el = eos.ZharkovElecModel(n_3c, z_3c, prefix='el_')
    # define initial values for parameters
    params = eos_st.make_params(v0=v0['3C'], k0=k0_3c[key], k0p=k0p_3c[key])
    params += eos_th.make_params(v0=v0['3C'], gamma0=gamma0, q=q, theta0=theta0)
    params += eos_el.make_params(v0=v0['3C'], e0=0.1e-6, g=0.01)
    # construct PVT eos
    pvteos = eos_st + eos_th + eos_el
    # fix static parameters and some other well known parameters
    params['th_v0'].vary=False; params['th_gamma0'].vary=False; params['th_theta0'].vary=False
    params['th_q'].vary=False
    params['st_v0'].vary=False; params['st_k0'].vary=False; params['st_k0p'].vary=False
    params['el_v0'].vary=False#; params['el_e0'].vary=False#; params['el_g'].vary=False
    # calculate weights.  setting it None results in unweighted fitting
    weights = 1./unp.std_devs(p) #None
    fit_result = pvteos.fit(unp.nominal_values(p), params, v=unp.nominal_values(v), 
                                   temp=unp.nominal_values(temp), weights=weights)
    print('********'+key)
    print(fit_result.fit_report())
    # plot fitting results
    eos.plot.thermal_fit_result(fit_result, p_err=unp.std_devs(p), v_err=unp.std_devs(v), title=key)


# # 5. Data fitting with constq equation (3C) with anharmonic effects

# In[14]:

for key, value in au_eos.items():
    # calculate pressure
    p = au_eos[key].cal_p(v_std, temp)
    # add prefix to the parameters.  this is important to distinguish thermal and static parameters
    eos_st = st_model[key](prefix='st_') 
    eos_th = eos.ConstqModel(n_3c, z_3c, prefix='th_')
    eos_anh = eos.ZharkovAnhModel(n_3c, z_3c, prefix='anh_')
    # define initial values for parameters
    params = eos_st.make_params(v0=v0['3C'], k0=k0_3c[key], k0p=k0p_3c[key])
    params += eos_th.make_params(v0=v0['3C'], gamma0=gamma0, q=q, theta0=theta0)
    params += eos_anh.make_params(v0=v0['3C'], a0=0.1e-6, m=0.01)
    # construct PVT eos
    pvteos = eos_st + eos_th + eos_anh
    # fix static parameters and some other well known parameters
    params['th_v0'].vary=False; params['th_gamma0'].vary=False; params['th_theta0'].vary=False
    params['th_q'].vary=False
    params['st_v0'].vary=False; params['st_k0'].vary=False; params['st_k0p'].vary=False
    params['anh_v0'].vary=False#; params['el_e0'].vary=False#; params['el_g'].vary=False
    # calculate weights.  setting it None results in unweighted fitting
    weights = 1./unp.std_devs(p) #None
    fit_result = pvteos.fit(unp.nominal_values(p), params, v=unp.nominal_values(v), 
                                   temp=unp.nominal_values(temp), weights=weights)
    print('********'+key)
    print(fit_result.fit_report())
    # plot fitting results
    eos.plot.thermal_fit_result(fit_result, p_err=unp.std_devs(p), v_err=unp.std_devs(v), title=key)


# In[ ]:



