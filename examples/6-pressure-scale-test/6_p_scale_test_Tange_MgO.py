
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('cat', '0Source_Citation.txt')


# In[2]:


get_ipython().run_line_magic('matplotlib', 'inline')
# %matplotlib notebook # for interactive


# For high dpi displays.

# In[3]:


get_ipython().run_line_magic('config', "InlineBackend.figure_format = 'retina'")


# # 0. General note

# This example compares pressure calculated from `pytheos` and original publication for the gold scale by Tange 2008.

# # 1. Global setup

# In[4]:


import matplotlib.pyplot as plt
import numpy as np
from uncertainties import unumpy as unp
import pytheos as eos


# # 3. Compare

# In[5]:


eta = np.linspace(1., 0.65, 36)
print(eta)


# In[6]:


tange_mgo = eos.periclase.Tange2009()


# In[7]:


tange_mgo.print_equations()


# In[8]:


tange_mgo.print_equations()


# In[9]:


tange_mgo.print_parameters()


# In[10]:


v0 = 74.698


# In[11]:


tange_mgo.three_r


# In[12]:


v = v0 * (eta) 
temp = 3000.


# In[13]:


p = tange_mgo.cal_p(v, temp * np.ones_like(v))


# <img src='./tables/Tange_MgO.png'>

# In[14]:


print('for T = ', temp)
for eta_i, p_i in zip(eta, p):
    print("{0: .3f} {1: .2f}".format(eta_i, p_i))


# Make comparison with the Vinet column in the table.

# In[15]:


v = tange_mgo.cal_v(p, temp * np.ones_like(p), min_strain=0.6)
print(1.-(v/v0))

