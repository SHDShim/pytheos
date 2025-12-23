#!/usr/bin/env python
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

# Quick pressure calculator using `pytheos`.

# # 1. Global setup

# In[4]:


import matplotlib.pyplot as plt
import numpy as np
from uncertainties import unumpy as unp
from uncertainties import ufloat
import pytheos as eos


# In[5]:


scales = {'MgO': [eos.periclase.Tange2009(), eos.periclase.Dorogokupets2015(),
              eos.periclase.Dorogokupets2015() , eos.periclase.Jamieson1982()],
          'Au': [eos.gold.Dorogokupets2007(), eos.gold.Dorogokupets2015(),
                   eos.gold.Fei2007bm3(), eos.gold.Fei2007vinet(), 
                   eos.gold.Heinz1984(), eos.gold.Shim2002(), 
                   eos.gold.Jamieson1982H(), eos.gold.Jamieson1982L()],
          'Pt': [eos.platinum.Dorogokupets2007(), eos.platinum.Dorogokupets2015(),
                 eos.platinum.Fei2007bm3(), eos.platinum.Fei2007vinet(), eos.platinum.Holmes1989(),
                 eos.platinum.Jamieson1982(), eos.platinum.Yokoo2009()] }


# # 2. Calculator

# In[6]:


# change for'Au', 'Pt', 'MgO'
my_standard = 'Pt' 

# See below about how to make inputs
my_volume = 60.0 
my_temperature = 2000. 


# _Examples for volume and temperature inputs_
# 
# For single data point with errors
# ```python
# my_volume = ufloat(60.0, 0.01)  
# my_temperature = ufloat(1500, 100) 
# ```
# 
# For array of data points without errors
# ```python
# my_volume = np.array([70, 60]) 
# my_temperature = np.array([1500, 1700]) 
# ```
# 
# For arrays with error bars
# ```python
# nominal = np.array([1.0, 2.0, 3.0])
# error = np.array([0.1, 0.2, 0.1])
# my_volume = unp.uarray(nominal, error)
# nominal = np.array([2100, 1210, 2100])
# error = np.array([200, 150, 160])
# my_temperature = unp.uarray(nominal, error)
# ```

# In[7]:


for scale in scales[my_standard]:
    eos_i = scale
    eos_i.print_reference()
    print(eos_i.cal_p(my_volume, my_temperature))


# In[ ]:





# In[ ]:




