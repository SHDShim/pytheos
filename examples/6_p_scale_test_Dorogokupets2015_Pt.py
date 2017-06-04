
# coding: utf-8

# ## Source and citation
# 
# - This notebook is a part of the `pytheos` package ([Github](http://github.com/SHDShim/pytheos)). 
# 
# - __[How to cite]__ S.-H. Shim (2017) Pytheos - a python tool set for equations of state. DOI:

# In[1]:

get_ipython().magic('matplotlib inline')
# %matplotlib notebook # for interactive


# For high dpi displays.

# In[2]:

get_ipython().magic("config InlineBackend.figure_format = 'retina'")


# # 0. General note

# This example compares pressure calculated from `pytheos` and original publication for the platinum scale by Dorogokupets 2015.

# # 1. Global setup

# In[3]:

import matplotlib.pyplot as plt
import numpy as np
from uncertainties import unumpy as unp
import pytheos as eos


# # 3. Compare

# In[4]:

eta = np.linspace(1., 0.70, 7)
print(eta)


# In[5]:

dorogokupets2015_pt = eos.platinum.Dorogokupets2015()


# In[6]:

help(eos.platinum.Dorogokupets2015)


# In[7]:

dorogokupets2015_pt.print_equations()


# In[8]:

dorogokupets2015_pt.print_equations()


# In[9]:

dorogokupets2015_pt.print_parameters()


# In[10]:

v0 = 60.37930856339099


# In[11]:

dorogokupets2015_pt.three_r


# In[12]:

v = v0 * (eta) 
temp = 3000.


# In[13]:

p = dorogokupets2015_pt.cal_p(v, temp * np.ones_like(v))


# In[14]:

print('for T = ', temp)
for eta_i, p_i in zip(eta, p):
    print("{0: .3f} {1: .2f}".format(eta_i, p_i))


# The table is not given in this publication.

# In[15]:

v = dorogokupets2015_pt.cal_v(p, temp * np.ones_like(p), min_strain=0.6)
print((v/v0))

