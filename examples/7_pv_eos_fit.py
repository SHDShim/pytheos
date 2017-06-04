
# coding: utf-8

# ## Source and citation
# 
# - This notebook is part of the `pytheos` package ([Github]()). 
# 
# - __[Citation]__ S.-H. Shim (2017) Pytheos - python equations of state tools. doi:

# In[3]:

get_ipython().magic('matplotlib inline')
# %matplotlib notebook # for interactive


# For high dpi displays.

# In[4]:

get_ipython().magic("config InlineBackend.figure_format = 'retina'")


# # 0. General note

# - Through high pressure experiments, we can measure the molar volume of minerals at high pressure.  The data can be then fit to static equations of state, in order to obtain bulk modulus $(K_0)$ and pressure derivative of bulk modulus $(K'_0)$.  
# 
# - `pytheos` uses `lmfit` to provide equation of state fitting for static compression data.  It has a few convenient features for the flexibility, such as fixing/varying particular parameters, propagating uncertainties, etc.  This notebook shows a few examples.

# # 1. Global setup

# In[5]:

import numpy as np
import matplotlib.pyplot as plt
import pytheos as eos


# # 2. Setup fitting models

# Using MgO, we create models for three equations for static EOS.  Volume can be either in unit-cell volume or molar volume.  But `v` and `v0` should be in the same unit.

# In[8]:

v0 = 11.244; k0 = 160.; k0p = 4.0


# Three different equations for static EOS of MgO.

# In[9]:

exp_bm3 = eos.BM3Model()
exp_vinet = eos.VinetModel()
exp_kunc = eos.KuncModel(order=5) # this way I can handle order


# Assign initial values to the parameters.

# In[10]:

params = exp_bm3.make_params(v0=v0, k0=k0, k0p=k0p)


# # 3. Synthesize data

# We make data with random noise.

# In[12]:

v_data = v0 * np.linspace(0.99,0.6,20)


# Calculate pressures from three different equations.

# In[13]:

p_bm3 = exp_bm3.eval(params, v=v_data)
p_vinet = exp_vinet.eval(params, v=v_data)
p_kunc = exp_kunc.eval(params, v=v_data)


# Create random noise (`noise`) and add it to pressure value to simulate the data.

# In[17]:

sp = np.random.normal(0.,2.,20)
p_data = p_bm3 + sp


# Plot the synthetic data together with true values.

# In[19]:

plt.plot(p_data, v_data, 'ko', label='data')
plt.plot(p_bm3, v_data, label='bm3')
plt.plot(p_vinet, v_data, label='vinet')
plt.plot(p_kunc, v_data, label='kunc')
plt.xlabel('Pressure (GPa)'); plt.ylabel('Molar volume (cm$^3$/mol)')
plt.legend();


# The cell below shows the systematic differences between the equations.

# In[21]:

plt.plot(p_bm3, p_vinet-p_bm3, label='Vinet - BM3')
plt.plot(p_bm3, p_kunc-p_bm3, label='Kunc - BM3')
plt.xlabel('Pressure (GPa)'); plt.ylabel('$\Delta P$')
plt.axhline(0, c='k', ls=':')
plt.legend();


# # 4. EOS fitting

# We fix $V_0$ as it is typically well known.  We can also fix any parameters in `lmfit`.

# In[22]:

params['v0'].vary = False
#params['k0'].vary = False


# It is common practice not to include weight for the $P$-$V$ fit because it will de-emphasize the high pressure data points.  But if you need, you can set weight in the following way.  `None` performs unweighted fitting

# In[27]:

weight = None # 1./sp#None


# The cell below performs fitting with BM3 equation.

# In[29]:

fitresult_bm3 = exp_bm3.fit(p_data, params, v=v_data, weights=weight)


# You may print out the fitting result.

# In[31]:

print(fitresult_bm3.fit_report())


# Fitting can be performed for other equations.

# In[32]:

fitresult_vinet = exp_vinet.fit(p_data, params, v=v_data, weights=weight)
print(fitresult_vinet.fit_report())


# In[33]:

fitresult_kunc = exp_kunc.fit(p_data, params, v=v_data, weights=weight)
print(fitresult_kunc.fit_report())


# # 5. Play with the results

# The example below shows how to get individual parameters from the fitting result.

# In[34]:

fitresult_vinet.params['k0p'].value


# You can also get estimated uncertainties.

# In[35]:

fitresult_vinet.params['k0p'].stderr


# The cells below show how to calculate fitted pressure values at the data points.  First, we get the fit parameters.

# In[37]:

v0_r = fitresult_vinet.params['v0']
k0_r = fitresult_vinet.params['k0']
k0p_r = fitresult_vinet.params['k0p']


# Then we get the function used for the fitting.

# In[38]:

f = fitresult_vinet.model


# `eval` conducts calculation for given volume value.

# In[39]:

f.eval(params, v=v_data)


# In[41]:

f.func(v_data, v0, k0, k0p)


# If you want to get array of (p, v) for plotting smooth fitting curve:

# In[42]:

v_fitline = np.linspace(v0,v_data.min(),100)


# In[43]:

p_fitline = f.func(v_fitline, v0_r, k0_r, k0p_r)


# In[45]:

plt.plot(p_fitline, v_fitline)
plt.plot(p_data_bm3, v_data, 'k.', label='data')
plt.xlabel('Pressure (GPa)'); plt.ylabel('Molar volume (cm$^3$/mol)');


# The cell below shows how to get the fitting results with uncertainties.

# In[46]:

fitresult_bm3.params


# The cell below shows how to get the uncertainties of the fit curve.

# In[53]:

from uncertainties import ufloat
p_err = f.func(v_fitline, ufloat(fitresult_bm3.params['v0'].value, 0.),
         ufloat(fitresult_bm3.params['k0'].value, fitresult_bm3.params['k0'].stderr),
         ufloat(fitresult_bm3.params['k0p'].value, fitresult_bm3.params['k0p'].stderr))


# The cell below shows the results in a table.  It only shows the first five data points, but `p_err` contains 100 points as we setup that way above.

# In[54]:

import pandas as pd
df = pd.DataFrame()
df['p'] = p_err
df.head()


# `pytheos` has internal plot script, which takes care of the plotting the fitting results with useful information.

# In[55]:

eos.plot.static_fit_result(fitresult_bm3, p_err=sp)


# How about the fit results for the other functions?

# In[57]:

eos.plot.static_fit_result(fitresult_vinet)


# In[58]:

eos.plot.static_fit_result(fitresult_kunc)


# Why the estimated uncertainty shown in the bottom parts of the figures are greater in Vinet and Kunc fittings?  It is perhaps due to the systematic differences among the equations.  Note that the synthetic data set was created from BM3 equation.
