
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

# This example compares pressure calculated from `pytheos` and original publication for the gold scale by Jamieson 1983.

# # 1. Global setup

# In[4]:

import matplotlib.pyplot as plt
import numpy as np
from uncertainties import unumpy as unp
import pytheos as eos


# # 3. Compare

# In[5]:

eta = np.linspace(0., 0.225, 46)
print(eta)


# In[6]:

jamieson_aul = eos.gold.Jamieson1982L()
jamieson_auh = eos.gold.Jamieson1982H()


# In[7]:

jamieson_aul.print_equations()


# In[8]:

jamieson_auh.print_equations()


# In[9]:

jamieson_aul.print_parameters()


# In[10]:

jamieson_auh.print_parameters()


# In[11]:

v0 = 67.84747902176544


# In[12]:

jamieson_aul.three_r


# In[13]:

jamieson_auh.three_r


# In[14]:

v = v0 * (1.-eta) 
temp = 1500.


# In[15]:

p = jamieson_auh.cal_p(v, temp * np.ones_like(v))


# <img src='./tables/Jamieson_Au_1.png'>
# <img src='./tables/Jamieson_Au_2.png'>

# In[16]:

print('for T = ', temp)
for eta_i, p_i in zip(eta, p):
    print("{0: .3f} {1: .2f}".format(eta_i, p_i))


# In[17]:

v = jamieson_auh.cal_v(p, temp * np.ones_like(p), min_strain=0.6)
print(1.-(v/v0))

