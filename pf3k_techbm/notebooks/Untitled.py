
# coding: utf-8

# In[2]:

import numpy as np


# In[3]:

a = np.array([0, 1, 4, 6, 12])
b = np.array([-1, 1, 2, 5, 9, 12, 22])
print(np.searchsorted(a, b))
print(np.searchsorted(a, b, side='right'))
print(np.searchsorted(a, b) % 5)


# In[4]:

b - a[np.searchsorted(a, b)]


# In[5]:

b - a[np.searchsorted(a, b)-1]


# In[6]:

abs(b - a[np.searchsorted(a, b) % len(a)])


# In[7]:

abs(b - a[np.searchsorted(a, b)-1])


# In[30]:

np.minimum(abs(b - a[np.searchsorted(a, b) % len(a)]), abs(b - a[np.searchsorted(a, b)-1]))


# In[8]:

b - a[np.searchsorted(a, b, side='right')-1]


# In[8]:

b = np.array([0, 1, 4, 6, 12])
a = np.array([-1, 1, 2, 5, 9, 12, 22])


# In[17]:

def find_nearest(a=a, b=b):
    nearest_before = abs(a - b[np.searchsorted(b, a)-1])
    print(nearest_before)
    nearest_after = abs(a - b[np.searchsorted(b, a) % len(b)])
    print(nearest_after)
    return(np.minimum(nearest_before, nearest_after))
    


# In[18]:

print(a)
print(b)
print()
find_nearest()


# In[ ]:



