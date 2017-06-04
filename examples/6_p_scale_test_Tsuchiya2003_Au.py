
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

# This example compares pressure calculated from `pytheos` and original publication for the gold scale by Tsuchiya 2003.

# # 1. Global setup

# In[4]:

import matplotlib.pyplot as plt
import numpy as np
from uncertainties import unumpy as unp
import pytheos as eos


# # 3. Compare

# In[5]:

eta = np.linspace(0., 0.34, 18)
print(eta)


# In[6]:

tsuchiya_au = eos.gold.Tsuchiya2003()


# In[7]:

help(tsuchiya_au)


# In[8]:

tsuchiya_au.print_equations()


# In[9]:

tsuchiya_au.print_equations()


# In[10]:

tsuchiya_au.print_parameters()


# In[11]:

v0 = 67.84742110765599


# In[12]:

tsuchiya_au.three_r


# In[13]:

v = v0 * (1.-eta) 
temp = 2500.


# In[14]:

p = tsuchiya_au.cal_p(v, temp * np.ones_like(v))


# <img src='./tables/Tsuchiya_Au.png'>

# In[15]:

print('for T = ', temp)
for eta_i, p_i in zip(eta, p):
    print("{0: .3f} {1: .2f}".format(eta_i, p_i))


# In[16]:

v = tsuchiya_au.cal_v(p, temp * np.ones_like(p), min_strain=0.6)
print(1.-(v/v0))


# - I cannot quite reproduce the table values.  The mismatch is about 1 GPa at 2500 K and 240 GPa.
# 
# - This means his parameters may have been rounded.
# 
# - Therefore, I readjusted the eos parameters from Tsuchiya to match their table values better.  Users have a choice if they use the table values or the parameter values.  If `reproduce_table` sets to `True`, the difference reduces to 0.1 GPa.

# In[17]:

tsuchiya_au = eos.gold.Tsuchiya2003(reproduce_table=True)


# In[18]:

p = tsuchiya_au.cal_p(v, temp * np.ones_like(v))


# In[19]:

print('for T = ', temp)
for eta_i, p_i in zip(eta, p):
    print("{0: .3f} {1: .2f}".format(eta_i, p_i))


# In[ ]:



