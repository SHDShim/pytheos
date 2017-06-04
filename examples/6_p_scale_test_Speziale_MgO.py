
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

# This example compares pressure calculated from `pytheos` and original publication for the gold scale by Speiale 2001.

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


speziale_mgo = eos.periclase.Speziale2001()


# In[6]:


speziale_mgo.print_equations()


# In[7]:


speziale_mgo.print_equations()


# In[8]:


speziale_mgo.print_parameters()


# In[9]:


v0 = 74.698


# In[10]:


speziale_mgo.three_r


# In[11]:


v = v0 * (eta) 
temp = 3000.


# In[12]:


p = speziale_mgo.cal_p(v, temp * np.ones_like(v))


# In[13]:


print('for T = ', temp)
for eta_i, p_i in zip(eta, p):
    print("{0: .3f} {1: .2f}".format(eta_i, p_i))


# In[15]:


v = speziale_mgo.cal_v(p, temp * np.ones_like(p), min_strain=0.6)
print((v/v0))


# In[ ]:




