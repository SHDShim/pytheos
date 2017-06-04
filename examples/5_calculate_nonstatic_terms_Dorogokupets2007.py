
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

# - This notebook shows the magnitude of different non-static pressure terms in the EOS of platinum by Dorogokupets and Dewaele (2007, HPR).

# # 1. General setup

# In[4]:

import uncertainties as uct
import numpy as np
import matplotlib.pyplot as plt
from uncertainties import unumpy as unp
import pytheos as eos


# In[5]:

v0 = 3.9231**3
v = np.linspace(v0, v0 * 0.8, 20)


# # 2. Calculate thermal pressure

# In[6]:

p_th = eos.dorogokupets2007_pth(v, 2000., v0, 2.82, 1.83, 8.11, 220., 1, 4)


# # 3. Calculate pressure from anharmonicity

# In[7]:

help(eos.zharkov_panh)


# In[8]:

p_anh = eos.zharkov_panh(v, 2000., v0, -166.9e-6, 4.32, 1, 4)


# # 4. Calculate pressure from electronic effects

# In[9]:

help(eos.zharkov_pel)


# In[10]:

p_el = eos.zharkov_pel(v, 2000., v0, 260.e-6, 2.4, 1, 4)


# # 5. Plot with respect to volume

# In[11]:

plt.plot(v, p_th, label='$P_{th}$')
plt.plot(v, p_el, label='$P_{el}$')
plt.plot(v, p_anh, label='$P_{anh}$')
plt.legend();


# # 5. Plot with respect to pressure

# We call the built-in dorogokupets2007 scale in `pytheos`.

# In[12]:

dorogokupets2007_pt = eos.platinum.Dorogokupets2007()


# In[13]:

help(dorogokupets2007_pt)


# In[14]:

p = dorogokupets2007_pt.cal_p(v, 2000.)


# In[15]:

plt.plot(unp.nominal_values(p), p_th, label='$P_{th}$')
plt.plot(unp.nominal_values(p), p_el, label='$P_{el}$')
plt.plot(unp.nominal_values(p), p_anh, label='$P_{anh}$')
plt.legend();


# In[ ]:



