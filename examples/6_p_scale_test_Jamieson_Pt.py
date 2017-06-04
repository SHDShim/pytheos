
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

# This example compares pressure calculated from `pytheos` and original publication for the platinum scale by Jamieson 1983.

# # 1. Global setup

# In[3]:

import matplotlib.pyplot as plt
import numpy as np
from uncertainties import unumpy as unp
import pytheos as eos


# # 3. Compare

# In[4]:

eta = np.linspace(0., 0.18, 37)


# In[5]:

jamieson_pt = eos.platinum.Jamieson1982()


# In[6]:

jamieson_pt.print_equations()


# In[7]:

jamieson_pt.print_parameters()


# In[8]:

v0 = 60.421757141035926


# In[9]:

jamieson_pt.three_r


# In[10]:

v = v0 * (1.-eta) 
temp = 1500.


# In[11]:

p = jamieson_pt.cal_p(v, temp * np.ones_like(v))


# <img src='./tables/Jamieson_Pt_1.png'>

# In[12]:

print('for T = ', temp)
for eta_i, p_i in zip(eta, p):
    print("{0: .3f} {1: .2f}".format(eta_i, p_i))


# In[13]:

v = jamieson_pt.cal_v(p, temp * np.ones_like(p), min_strain=0.6)
print(1.-(v/v0))


# In[ ]:



