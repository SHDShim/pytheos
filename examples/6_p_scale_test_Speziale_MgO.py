
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

# This example compares pressure calculated from `pytheos` and original publication for the gold scale by Speiale 2001.

# # 1. Global setup

# In[4]:

import matplotlib.pyplot as plt
import numpy as np
from uncertainties import unumpy as unp
import pytheos as eos


# # 3. Compare

# In[5]:

eta = np.linspace(1., 0.60, 21)
print(eta)


# In[6]:

speziale_mgo = eos.periclase.Speziale2001()


# In[7]:

speziale_mgo.print_equations()


# In[8]:

speziale_mgo.print_equations()


# In[9]:

speziale_mgo.print_parameters()


# In[10]:

v0 = 74.698


# In[11]:

speziale_mgo.three_r


# In[16]:

v = v0 * (eta) 
temp = 300.


# In[17]:

p = speziale_mgo.cal_p(v, temp * np.ones_like(v))


# In[18]:

print('for T = ', temp)
for eta_i, p_i in zip(eta, p):
    print("{0: .3f} {1: .2f}".format(eta_i, p_i))


# In[19]:

v = speziale_mgo.cal_v(p, temp * np.ones_like(p), min_strain=0.6)
print((v/v0))


# In[ ]:



