
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

# This example compares pressure calculated from `pytheos` and original publication for the gold scale by Jamieson 1983.

# # 1. Global setup

# In[3]:

import matplotlib.pyplot as plt
import numpy as np
from uncertainties import unumpy as unp
import pytheos as eos


# # 3. Compare

# In[4]:

eta = np.linspace(0., 0.225, 46)
print(eta)


# In[5]:

jamieson_aul = eos.gold.Jamieson1982L()
jamieson_auh = eos.gold.Jamieson1982H()


# In[6]:

jamieson_aul.print_equations()


# In[7]:

jamieson_auh.print_equations()


# In[8]:

jamieson_aul.print_parameters()


# In[9]:

jamieson_auh.print_parameters()


# In[10]:

v0 = 67.84747902176544


# In[11]:

jamieson_aul.three_r


# In[12]:

jamieson_auh.three_r


# In[13]:

v = v0 * (1.-eta) 
temp = 1500.


# In[14]:

p = jamieson_auh.cal_p(v, temp * np.ones_like(v))


# <img src='./tables/Jamieson_Au_1.png'>
# <img src='./tables/Jamieson_Au_2.png'>

# In[15]:

print('for T = ', temp)
for eta_i, p_i in zip(eta, p):
    print("{0: .3f} {1: .2f}".format(eta_i, p_i))


# In[16]:

v = jamieson_auh.cal_v(p, temp * np.ones_like(p), min_strain=0.6)
print(1.-(v/v0))

