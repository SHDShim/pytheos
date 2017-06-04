
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

# This example compares pressure calculated from `pytheos` and original publication for the platinum scale by Dorogokupets 2007.

# # 1. Global setup

# In[4]:

import matplotlib.pyplot as plt
import numpy as np
from uncertainties import unumpy as unp
import pytheos as eos


# # 3. Compare

# In[5]:

eta = np.linspace(1., 0.70, 7)
print(eta)


# In[6]:

dorogokupets2007_pt = eos.platinum.Dorogokupets2007()


# In[7]:

help(dorogokupets2007_pt)


# In[8]:

dorogokupets2007_pt.print_equations()


# In[9]:

dorogokupets2007_pt.print_equations()


# In[10]:

dorogokupets2007_pt.print_parameters()


# In[11]:

v0 = 60.37930856339099


# In[12]:

dorogokupets2007_pt.three_r


# In[13]:

v = v0 * (eta) 
temp = 3000.


# In[14]:

p = dorogokupets2007_pt.cal_p(v, temp * np.ones_like(v))


# <img src='./tables/Dorogokupets2007_Pt.png'>

# In[15]:

print('for T = ', temp)
for eta_i, p_i in zip(eta, p):
    print("{0: .3f} {1: .2f}".format(eta_i, p_i))


# In[16]:

v = dorogokupets2007_pt.cal_v(p, temp * np.ones_like(p), min_strain=0.6)
print(1.-(v/v0))


# In[ ]:



