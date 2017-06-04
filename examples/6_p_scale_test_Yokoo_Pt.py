
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

# This example compares pressure calculated from `pytheos` and original publication for the platinum scale by Yokoo 2009.

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

yokoo_pt = eos.platinum.Yokoo2009()


# In[6]:

yokoo_pt.print_equations()


# In[7]:

yokoo_pt.print_equations()


# In[8]:

yokoo_pt.print_parameters()


# In[9]:

v0 = 60.37930856339099


# In[10]:

yokoo_pt.three_r


# In[11]:

v = v0 * (eta) 
temp = 3000.


# In[12]:

p = yokoo_pt.cal_p(v, temp * np.ones_like(v))


# <img src='./tables/Yokoo_Pt.png'>

# In[13]:

print('for T = ', temp)
for eta_i, p_i in zip(eta, p):
    print("{0: .3f} {1: .2f}".format(eta_i, p_i))


# It is alarming that even 300 K isotherm does not match with table value.  The difference is 1%.

# In[14]:

v = yokoo_pt.cal_v(p, temp * np.ones_like(p), min_strain=0.6)
print(1.-(v/v0))

