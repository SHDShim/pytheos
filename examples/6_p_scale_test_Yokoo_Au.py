
# coding: utf-8

# ## Source and citation
# 
# - This notebook is part of the `pytheos` package ([Github]()). 
# 
# - __[Citation]__ S.-H. Shim (2017) Pytheos - python equations of state tools. doi:

# In[1]:

get_ipython().magic('matplotlib inline')
# %matplotlib notebook # for interactive


# For high dpi displays.

# In[2]:

get_ipython().magic("config InlineBackend.figure_format = 'retina'")


# # 0. General note

# This example compares pressure calculated from `pytheos` and original publication for the gold scale by Yokoo 2009.

# # 1. Global setup

# In[3]:

import matplotlib.pyplot as plt
import numpy as np
from uncertainties import unumpy as unp
import pytheos as eos


# # 3. Compare

# In[4]:

eta = np.linspace(1., 0.60, 21)
print(eta)


# In[5]:

yokoo_au = eos.gold.Yokoo2009()


# In[6]:

yokoo_au.print_equations()


# In[7]:

yokoo_au.print_equations()


# In[8]:

yokoo_au.print_parameters()


# In[9]:

v0 = 67.84742110765599


# In[10]:

yokoo_au.three_r


# In[11]:

v = v0 * (eta) 
temp = 3000.


# In[12]:

p = yokoo_au.cal_p(v, temp * np.ones_like(v))


# <img src='./tables/Yokoo_Au.png'>

# In[13]:

print('for T = ', temp)
for eta_i, p_i in zip(eta, p):
    print("{0: .3f} {1: .2f}".format(eta_i, p_i))


# In[18]:

v = yokoo_au.cal_v(p, temp * np.ones_like(p), min_strain=0.6)
print(1.-(v/v0))


# - I cannot quite reproduce the table values.  The mismatch is about 3 GPa at 3000 K and 380 GPa.
# 
# - This means his parameters may have been rounded.
# 
# - Therefore, I readjusted the eos parameters from Yokoo to match their table values better.  Users have a choice if they use the table values or the parameter values.  If `reproduce_table` sets to `True`, the difference reduces to 0.1 GPa.

# In[19]:

yokoo_au = eos.gold.Yokoo2009(reproduce_table=True)


# In[20]:

p = yokoo_au.cal_p(v, temp * np.ones_like(v))

