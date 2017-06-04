
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

# - There exist different formulations for the Gruneisen parameter and thermal pressure.  In this notebook we will compare the different formulations for MgO.

# # 1. General setup

# In[4]:

import uncertainties as uct
import numpy as np
from uncertainties import unumpy as unp
import matplotlib.pyplot as plt
import pytheos as eos


# # 2. Calculate Gruneisen parameter

# We test three different ways to get Gruneisen parameter at high pressure-temperature.

# In[5]:

help(eos.constq_grun)


# In[6]:

help(eos.speziale_grun)


# In[7]:

help(eos.altshuler_grun)


# We will test for MgO.

# In[8]:

v0 = 74.698
v = np.linspace(v0, v0*0.8, 10)


# - From Dewaele et al. (2000, JGR):
# 
# $$\dfrac{\gamma}{\gamma_0} = \left(\dfrac{V}{V_0}\right)^q$$
# 
# where $\gamma_0 = 1.45$ and $q = 0.8 \pm 0.5$.
# 
# - From Speziale et al. (2001, JGR):
# 
# $$\gamma = \gamma_0 \exp\left\{ \dfrac{q_0}{q_1} \left[ \left(\dfrac{V}{V_0}\right)^{q_1} -1 \right] \right\}$$
# 
# where $\gamma_0 = 1.49 \pm 0.03$, $q_0 = 1.65 \pm 0.4$, and $q_1 = 11.8 \pm 0.2$.
# 
# - From Dorogokupets and Dewaele (2007, HPR):
# 
# $$q = \beta x^\beta \dfrac{\gamma_0 - \gamma_\inf}{\gamma} $$
# 
# $$\gamma = \gamma_0 x^\beta$$
# 
# where $\gamma_0 = 1.50$, $\gamma_\inf = 0.75$, and $\beta = 2.96$.
# 
# 
# 

# In[9]:

gamma_de = eos.constq_grun(v, v0, 1.45, 0.8)
gamma_sp = eos.speziale_grun(v, v0, uct.ufloat(1.49, 0.03), uct.ufloat(1.65, 0.4), uct.ufloat(11.8, 0.2))
gamma_do = eos.altshuler_grun(v, v0, 1.50, 0.75, 2.96) 


# In[10]:

plt.plot(v, gamma_de, label='Dewaele2000')
plt.errorbar(v, unp.nominal_values(gamma_sp), yerr=unp.std_devs(gamma_sp), label='Speziale2001')
plt.plot(v, gamma_do, label='Dorogokupets2007')
plt.xlabel('Unit-cell volume ($\mathrm{\AA}^3$)')
plt.ylabel(r"$\mathrm{Gr{\"u}neisen parameter}$")
plt.legend();


# # 3. Calculate Debye temperature

# In[11]:

help(eos.constq_debyetemp)


# In[12]:

help(eos.speziale_debyetemp)


# In[13]:

help(eos.altshuler_debyetemp)


# In[14]:

theta_de = eos.constq_debyetemp(v, v0, 1.45, 0.8, 800)
theta_sp = eos.speziale_debyetemp(v, v0, uct.ufloat(1.49, 0.03), uct.ufloat(1.65, 0.4), uct.ufloat(11.8, 0.2), 773)
theta_do = eos.altshuler_debyetemp(v, v0, 1.50, 0.75, 2.96, 760) 


# In[15]:

plt.plot(v, theta_de, label='Dewaele2000')
plt.errorbar(v, unp.nominal_values(theta_sp), yerr=unp.std_devs(theta_sp), label='Speziale2001')
plt.plot(v, theta_do, label='Dorogokupets2007')
plt.xlabel('Unit-cell volume ($\mathrm{\AA}^3$)')
plt.ylabel("Debye Temperature (K)")
plt.legend();


# # 4. Calculate thermal pressure

# In[16]:

help(eos.constq_pth)


# In[17]:

help(eos.speziale_pth)


# In[18]:

help(eos.dorogokupets2007_pth)


# In[19]:

n = 2; z = 4
temp = np.ones_like(v) * 2000.


# In[20]:

pth_de = eos.constq_pth(v, temp, v0, 1.45, 0.8, 800, n, z)
pth_sp = eos.speziale_pth(v, temp, v0, uct.ufloat(1.49, 0.03), uct.ufloat(1.65, 0.4), uct.ufloat(11.8, 0.2), 773, n, z)
pth_do = eos.dorogokupets2007_pth(v, temp, v0, 1.50, 0.75, 2.96, 760, n, z) 


# In[21]:

plt.plot(v, pth_de, label='Dewaele2000')
plt.errorbar(v, unp.nominal_values(pth_sp), yerr=unp.std_devs(pth_sp), label='Speziale2001')
plt.plot(v, pth_do, label='Dorogokupets2007')
plt.xlabel('Unit-cell volume ($\mathrm{\AA}^3$)')
plt.ylabel("Thermal pressure (GPa)")
plt.legend();

