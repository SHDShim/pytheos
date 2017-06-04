
# coding: utf-8

# ## Source and citation
# 
# - This notebook is a part of the `pytheos` package ([Github](http://github.com/SHDShim/pytheos)). 
# 
# - __[Citation]__ S.-H. Shim (2017) Pytheos - python equations of state tools. doi:

# In[1]:

get_ipython().magic('matplotlib inline')
# %matplotlib notebook # for interactive


# For high dpi displays.

# In[2]:

get_ipython().magic("config InlineBackend.figure_format = 'retina'")


# # 0. General note

# This example compares pressure calculated from `pytheos` and original publication for the MgO scale by Dorogokupets 2015.

# # 1. Global setup

# In[3]:

import matplotlib.pyplot as plt
import numpy as np
from uncertainties import unumpy as unp
import pytheos as eos


# # 3. Compare

# In[4]:

eta = np.linspace(1., 0.6, 9)
print(eta)


# In[5]:

dorogokupets2015_mgo = eos.periclase.Dorogokupets2015()


# In[6]:

help(eos.periclase.Dorogokupets2015)


# In[7]:

dorogokupets2015_mgo.print_equations()


# In[8]:

dorogokupets2015_mgo.print_equations()


# In[9]:

dorogokupets2015_mgo.print_parameters()


# In[10]:

v0 = 74.698


# In[11]:

dorogokupets2015_mgo.three_r


# In[12]:

v = v0 * (eta) 
temp = 3000.


# In[13]:

p = dorogokupets2015_mgo.cal_p(v, temp * np.ones_like(v))


# Table is not given in this publication.

# In[15]:

print('for T = ', temp)
for eta_i, p_i in zip(eta, p):
    print("{0: .3f} {1: .2f} ".format(eta_i, p_i))


# In[17]:

v = dorogokupets2015_mgo.cal_v(p, temp * np.ones_like(p), min_strain=0.6)
print((v/v0))

