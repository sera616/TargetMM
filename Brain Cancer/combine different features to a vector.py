
# coding: utf-8

# ## pssm+ss+pdo:20+3+1=24

# In[1]:


import pandas as pd 

##------------------  a
df1 = pd.read_csv('./data/APAACN.csv', encoding='gbk') 
df2 = pd.read_csv('./data/PAACN.csv', encoding='gbk') 
df3 = pd.read_csv('./data/WAPSSM_bc60.csv', encoding='gbk') 
df3=df3.fillna(0)
df4 = pd.read_csv('./data/WASS_bc9.csv', encoding='gbk')
df4=df4.fillna(0)
df5 = pd.read_csv('./data/WAPDO_bc3.csv', encoding='gbk')
df5=df5.fillna(0)

dfc2 = df1.merge(df2.loc[:,'proteinID':],left_on='proteinID', right_on='proteinID')
dfc3 = dfc2.merge(df3.loc[:,'proteinID':],left_on='proteinID', right_on='proteinID')
dfc4 = dfc3.merge(df4.loc[:,'proteinID':],left_on='proteinID', right_on='proteinID')
dfc5 = dfc4.merge(df5.loc[:,'proteinID':],left_on='proteinID', right_on='proteinID')

dfc5.to_csv('./data/bc202f.csv', encoding='gbk',sep=',',index=False,header=True) 


# In[ ]:


## k=11_nonWA
import pandas as pd 
##------------------  a
df1 = pd.read_csv('./data/result/PSSM/w11_pssm.csv', encoding='gbk') 
df2 = pd.read_csv('./data/result/SS/w11_ss.csv', encoding='gbk') 
df3 = pd.read_csv('./data/result/PDO/w11_pdo.csv', encoding='gbk') 

df1=df1.fillna(0)
df2=df2.fillna(0)
df3=df3.fillna(0)

dfc2 = df1.merge(df2.loc[:,'proteinID':],left_on='proteinID', right_on='proteinID')
dfc3 = dfc2.merge(df3.loc[:,'proteinID':],left_on='proteinID', right_on='proteinID')

dfc3.to_csv('./data/w11_nonWA.csv', encoding='gbk',sep=',',index=False,header=True) 


# In[5]:


## k=11_WA
import pandas as pd 
##------------------  a
df1 = pd.read_csv('./data/result/PSSM/WAPSSM_bc60.csv', encoding='gbk') 
df2 = pd.read_csv('./data/result/SS/WASS_bc9.csv', encoding='gbk') 
df3 = pd.read_csv('./data/result/PDO/WAPDO_bc3.csv', encoding='gbk') 

df1=df1.fillna(0)
df2=df2.fillna(0)
df3=df3.fillna(0)

dfc2 = df1.merge(df2.loc[:,'proteinID':],left_on='proteinID', right_on='proteinID')
dfc3 = dfc2.merge(df3.loc[:,'proteinID':],left_on='proteinID', right_on='proteinID')

dfc3.to_csv('./data/w11_WA.csv', encoding='gbk',sep=',',index=False,header=True) 


# In[1]:


## Global feature
import pandas as pd 
df1 = pd.read_csv('./data/APAACN.csv', encoding='gbk') 
df2 = pd.read_csv('./data/PAACN.csv', encoding='gbk') 
df1=df1.fillna(0)
df2=df2.fillna(0)

dfc2 = df1.merge(df2.loc[:,'proteinID':],left_on='proteinID', right_on='proteinID')

dfc2.to_csv('./data/bc_global.csv', encoding='gbk',sep=',',index=False,header=True)


