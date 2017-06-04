
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

# This example compares pressure calculated from `pytheos` and original publication for the MgO scale by Jamieson 1983.

# # 1. Global setup

# In[4]:

import matplotlib.pyplot as plt
import numpy as np
from uncertainties import unumpy as unp
import pytheos as eos


# # 3. Compare

# In[5]:

eta = np.linspace(0., 0.26, 53)
print(eta)


# In[6]:

jamieson_mgo = eos.periclase.Jamieson1982()


# In[7]:

jamieson_mgo.print_equations()


# In[8]:

jamieson_mgo.print_parameters()


# In[9]:

v0 = 74.67451012663052


# In[10]:

jamieson_mgo.three_r


# In[11]:

v = v0 * (1.-eta) 
temp = 1500.


# In[12]:

p = jamieson_mgo.cal_p(v, temp * np.ones_like(v))


# <img src='./tables/Jamieson_MgO_1.png'>
# <img src='./tables/Jamieson_MgO_2.png'>

# In[13]:

print('for T = ', temp)
for eta_i, p_i in zip(eta, p):
    print("{0: .3f} {1: .2f}".format(eta_i, p_i))


# In[14]:

v = jamieson_mgo.cal_v(p, temp * np.ones_like(p))
print(1.-(v/v0))

