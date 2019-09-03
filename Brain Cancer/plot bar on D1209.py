
# coding: utf-8

# ## 通过调调整Std,来适当调整每个柱状图上面数字的显示位置

# ## 10-flod-cross-validation

# In[3]:


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from collections import namedtuple
get_ipython().run_line_magic('matplotlib', 'inline')


n_groups = 6

##------------------- SN -------------------------##
SN = (0.2595,0.3419,0.1900,0.4819,0.4605,0.1421)
SNStd = (0.017, 0.017, 0.015, 0.022,0.021, 0.013)


##------------------- Sp -------------------------##
SP=(0.9943,0.9558,0.9900,1.0000,0.9638,0.9917)
SPStd = (0.015, 0.015, 0.015,0.015, 0.015, 0.015)



##------------------- ACC -------------------------##
ACC=(0.8289,0.8178,0.8101,0.8918,0.8587,0.8145)
ACCStd = (0.020, 0.020, 0.018, 0.025, 0.015, 0.025)



##------------------- MCC -------------------------##
MCC=(0.4329,0.3946,0.3380,0.6487,0.5249,0.1912)
MCCStd = (0.021, 0.021, 0.028, 0.021,0.040, 0.017)


font1 = {'family': 'Times New Roman',
         'weight': 'normal',
         'size': 8,
         }

plt.rc('font',family='Times New Roman')###-----------------

fig, ax = plt.subplots()
# index = np.arange(n_groups)
index = np.arange(1,18,3)
bar_width = 0.5

opacity = 0.7
error_config = {'ecolor': '0.5'}

rects1 = ax.bar(index, SN, bar_width,
                alpha=opacity, color='mediumpurple',
                yerr=SNStd, error_kw=error_config,
                label='SN')

rects2 = ax.bar(index + bar_width, SP, bar_width,
                alpha=opacity, color='olive',
                yerr=SPStd, error_kw=error_config,
                label='SP')

rects6 = ax.bar(index + 2*bar_width, ACC, bar_width,
                alpha=opacity, color='darkturquoise',
                yerr=ACCStd, error_kw=error_config,
                label='ACC')

rects7 = ax.bar(index + 3*bar_width, MCC, bar_width,
                alpha=opacity, color='seagreen',
                yerr=MCCStd, error_kw=error_config,
                label='MCC')


# ax.set_xlabel('Local      Global      Local+Global',fontsize=10)
ax.set_ylabel('Evaluation Values',fontsize=10)
ax.set_title('Feature Evalution using Three classifiers on D1209',fontsize=10)
ax.set_xticks(index + bar_width*1.5 )
ax.set_xticklabels(('RF', 'SVM', 'GPC', 'RF', 'SVM', 'GPC'),fontsize=10)

# ax.legend(fontsize=8,loc='best',bbox_to_anchor=(0.85,0.97))
ax.legend(fontsize=8,loc='best',bbox_to_anchor=(0.85,0.67),framealpha=0.5)

for rect in rects1:
    height = rect.get_height()
    ax.text(rect.get_x() + rect.get_width()/2., 1.06*height,
            '%0.4f' % (height),ha='center', va='bottom',fontsize=8)
    
for rect in rects2:
    height = rect.get_height()
    ax.text(rect.get_x() + rect.get_width()/2., 1.015*height,
            '%0.4f' % (height),ha='center', va='bottom',fontsize=8)    



for rect in rects6:
    height = rect.get_height()
    ax.text(rect.get_x() + rect.get_width()/2., 1.03*height,
            '%0.4f' % (height),ha='center', va='bottom',fontsize=8)         

for rect in rects7:
    height = rect.get_height()
    ax.text(rect.get_x() + rect.get_width()/2., 1.08*height,
            '%0.4f' % (height),ha='center', va='bottom',fontsize=8)  
    


plt.rcParams['figure.figsize'] = (16, 14) # 设置figure_size尺寸
plt.rcParams['savefig.dpi'] = 300 #图片像素
plt.rcParams['figure.dpi'] = 300 #分辨率
# 默认的像素：[6.0,4.0]，分辨率为100，图片尺寸为 600&400
# 指定dpi=200，图片尺寸为 1200*800
# 指定dpi=300，图片尺寸为 1800*1200
# 设置figsize可以在不改变分辨率情况下改变比例
plt.savefig('./data/BC_10fold-cross-validation.tif', dpi=300,bbox_inches ='tight') #指定分辨率保存

# fig.tight_layout()

plt.show()


# In[5]:


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from collections import namedtuple
get_ipython().run_line_magic('matplotlib', 'inline')


n_groups = 6

##------------------- SN -------------------------##
SN = (0.3529,0.3922,0.3333,0.5455,0.5000,0.3182)
SNStd = (0.017, 0.017, 0.015, 0.022,0.021, 0.013)


##------------------- Sp -------------------------##
SP=(0.9881,0.9484,0.9881,0.9916,0.9620,0.9873)
SPStd = (0.015, 0.015, 0.015,0.015, 0.015, 0.015)



##------------------- ACC -------------------------##
ACC=(0.8812,0.8548,0.8779,0.8944,0.8614,0.8416)
ACCStd = (0.020, 0.020, 0.018, 0.025, 0.015, 0.025)



##------------------- MCC -------------------------##
MCC=(0.5024,0.4090,0.4844,0.6693,0.5519,0.4670)
MCCStd = (0.021, 0.021, 0.028, 0.021,0.040, 0.017)


font1 = {'family': 'Times New Roman',
         'weight': 'normal',
         'size': 8,
         }

plt.rc('font',family='Times New Roman')###-----------------

fig, ax = plt.subplots()
# index = np.arange(n_groups)
index = np.arange(1,18,3)
bar_width = 0.5

opacity = 0.7
error_config = {'ecolor': '0.5'}

rects1 = ax.bar(index, SN, bar_width,
                alpha=opacity, color='mediumpurple',
                yerr=SNStd, error_kw=error_config,
                label='SN')

rects2 = ax.bar(index + bar_width, SP, bar_width,
                alpha=opacity, color='olive',
                yerr=SPStd, error_kw=error_config,
                label='SP')

rects6 = ax.bar(index + 2*bar_width, ACC, bar_width,
                alpha=opacity, color='darkturquoise',
                yerr=ACCStd, error_kw=error_config,
                label='ACC')

rects7 = ax.bar(index + 3*bar_width, MCC, bar_width,
                alpha=opacity, color='seagreen',
                yerr=MCCStd, error_kw=error_config,
                label='MCC')


# ax.set_xlabel('Local      Global      Local+Global',fontsize=10)
ax.set_ylabel('Evaluation Values',fontsize=10)
ax.set_title('Feature Evalution using Three classifiers on D1209',fontsize=10)
ax.set_xticks(index + bar_width*1.5 )
ax.set_xticklabels(('RF', 'SVM', 'GPC', 'RF', 'SVM', 'GPC'),fontsize=10)

# ax.legend(fontsize=8,loc='best',bbox_to_anchor=(0.85,0.97))
ax.legend(fontsize=8,loc='best',bbox_to_anchor=(0.85,0.67),framealpha=0.5)

for rect in rects1:
    height = rect.get_height()
    ax.text(rect.get_x() + rect.get_width()/2., 1.06*height,
            '%0.4f' % (height),ha='center', va='bottom',fontsize=8)
    
for rect in rects2:
    height = rect.get_height()
    ax.text(rect.get_x() + rect.get_width()/2., 1.015*height,
            '%0.4f' % (height),ha='center', va='bottom',fontsize=8)    



for rect in rects6:
    height = rect.get_height()
    ax.text(rect.get_x() + rect.get_width()/2., 1.03*height,
            '%0.4f' % (height),ha='center', va='bottom',fontsize=8)         

for rect in rects7:
    height = rect.get_height()
    ax.text(rect.get_x() + rect.get_width()/2., 1.08*height,
            '%0.4f' % (height),ha='center', va='bottom',fontsize=8)  
    


plt.rcParams['figure.figsize'] = (16, 14) # 设置figure_size尺寸
plt.rcParams['savefig.dpi'] = 300 #图片像素
plt.rcParams['figure.dpi'] = 300 #分辨率
# 默认的像素：[6.0,4.0]，分辨率为100，图片尺寸为 600&400
# 指定dpi=200，图片尺寸为 1200*800
# 指定dpi=300，图片尺寸为 1800*1200
# 设置figsize可以在不改变分辨率情况下改变比例
plt.savefig('./data/BC_Independent Test.tif', dpi=300,bbox_inches ='tight') #指定分辨率保存

# fig.tight_layout()

plt.show()

