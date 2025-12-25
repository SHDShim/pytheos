#!/usr/bin/env python
# coding: utf-8

# In[62]:


get_ipython().run_line_magic('cat', '0Source_Citation.txt')


# In[63]:


get_ipython().run_line_magic('matplotlib', 'inline')
# %matplotlib notebook # for interactive


# For high dpi displays.

# In[64]:


get_ipython().run_line_magic('config', "InlineBackend.figure_format = 'retina'")


# # 0. General note

# This example compares pressure calculated from `pytheos` and original publication for the NaCl-B2 scales.

# # 1. Global setup

# In[65]:


import matplotlib.pyplot as plt
import numpy as np
from uncertainties import unumpy as unp
import pytheos as eos


# # 2. Compare

# In[66]:


dorogokupets2007 = eos.sodium_chloride_b2.Dorogokupets2007()
fei2007v = eos.sodium_chloride_b2.Fei2007vinet()
fei2007bm = eos.sodium_chloride_b2.Fei2007bm3()


# In[67]:


help(eos.sodium_chloride_b2.Dorogokupets2007)


# In[68]:


dorogokupets2007.print_equations()


# In[69]:


dorogokupets2007.print_parameters()


# In[70]:


fei2007bm.print_parameters()


# In[71]:


v0_mol_dorogokupets2007 = 24.53 # value given in the paper


# In[72]:


from scipy.constants import N_A
v0_dorogokupets2007 = v0_mol_dorogokupets2007 * 1.e24 / N_A
v0_dorogokupets2007


# In[73]:


v0_fei2007 = 41.35
v0_mol_fei2007 = v0_fei2007 / 1.e24 * N_A
v0_mol_fei2007


# In[74]:


dorogokupets2007.three_r


# <img src='./tables/Dorogokupets2007_NaCl-B2.png' width=500>

# In[75]:


temp = np.ones(7) * 2500.


# In[76]:


vmol = np.linspace(17,11,7)
v = vmol/v0_mol_dorogokupets2007 * v0_dorogokupets2007
p = dorogokupets2007.cal_p(v, temp)
print('for T = ', temp)
for eta_i, p_i in zip(eta, p):
    print("{0: .3f} {1: .2f} ".format(eta_i, p_i))


# In[77]:


vmol = np.linspace(17,11,7)
v = vmol/v0_mol_fei2007 * v0_fei2007
p = fei2007bm.cal_p(v, temp)
print('for T = ', temp)
for eta_i, p_i in zip(eta, p):
    print("{0: .3f} {1: .2f} ".format(eta_i, p_i))


# I cannot reproduce Dorogokupets 2007 table for Fei and that is likely due to use of different v0 between these two studies.  For Dorogokupets2007 scale, I was able to reproduce their table in the paper.  
# 
# Fei2007 did not provide tables for their EOSs.  However, I digitized their figures and compared with my calculation.  They match well.
