
# coding: utf-8

# ## 通过调调整Std,来适当调整每个柱状图上面数字的显示位置

# ## 10-flod-cross-validation

# In[58]:


# In[ ]:


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from collections import namedtuple
get_ipython().run_line_magic('matplotlib', 'inline')


n_groups = 9

##------------------- SN -------------------------##
SN = (0.1662,0.2792,0.0355,0.7832,0.6753,0.6428)
SNStd = (0.017, 0.017, 0.012,  0.050,0.035, 0.035)


##------------------- Sp -------------------------##
SP=(0.9060,0.8035,0.9795,0.7418,0.5912,0.5927)
SPStd = (0.0085, 0.0085, 0.0085, 0.0085, 0.0085,0.0085)



##------------------- ACC -------------------------##
ACC=(0.6084,0.5926,0.5997,0.7626,0.6334,0.6178)
ACCStd = (0.033, 0.033, 0.033,0.042, 0.039, 0.036)



##------------------- MCC -------------------------##
MCC=(0.1079,0.0973,0.0436,0.5258,0.2675,0.2360)
MCCStd = (0.015, 0.017, 0.008, 0.042,0.022, 0.019)


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
ax.set_title('Feature Evalution using Three classifiers on D5223',fontsize=10)
ax.set_xticks(index + bar_width*1.5 )
ax.set_xticklabels(('RF', 'SVM', 'GPC', 'RF', 'SVM', 'GPC'),fontsize=10)

# ax.legend(fontsize=8,loc='best',bbox_to_anchor=(0.85,0.97))
ax.legend(fontsize=8,loc='best',bbox_to_anchor=(0.83,0.75),framealpha=0.5)

for rect in rects1:
    height = rect.get_height()
    ax.text(rect.get_x() + rect.get_width()/2., 1.08*height,
            '%0.4f' % (height),ha='center', va='bottom',fontsize=8)
    
for rect in rects2:
    height = rect.get_height()
    ax.text(rect.get_x() + rect.get_width()/2., 1.*height,
            '%0.4f' % (height),ha='center', va='bottom',fontsize=8)    



for rect in rects6:
    height = rect.get_height()
    ax.text(rect.get_x() + rect.get_width()/2., 1.065*height,
            '%0.4f' % (height),ha='center', va='bottom',fontsize=8)         

for rect in rects7:
    height = rect.get_height()
    ax.text(rect.get_x() + rect.get_width()/2., 1.08*height,
            '%0.4f' % (height),ha='center', va='bottom',fontsize=8)  
    

# plt.ylim((0, 1.05))
plt.rcParams['figure.figsize'] = (16, 14) # 设置figure_size尺寸
plt.rcParams['savefig.dpi'] = 300 #图片像素
plt.rcParams['figure.dpi'] = 300 #分辨率
# 默认的像素：[6.0,4.0]，分辨率为100，图片尺寸为 600&400
# 指定dpi=200，图片尺寸为 1200*800
# 指定dpi=300，图片尺寸为 1800*1200
# 设置figsize可以在不改变分辨率情况下改变比例
plt.savefig('./data/CV_10fold-cross-validation.tif', dpi=300,bbox_inches ='tight') #指定分辨率保存

# fig.tight_layout()

plt.show()


# In[56]:


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from collections import namedtuple
get_ipython().run_line_magic('matplotlib', 'inline')


n_groups = 9

##------------------- SN -------------------------##
SN = (0.1768,0.2683,0.0549,0.7760,0.6314,0.5894)
SNStd = (0.0, 0.0, 0.0, 0.0,0.0, 0.0)


##------------------- Sp -------------------------##
SP=(0.9029,0.7961,0.9705,0.7828,0.6259,0.6350)
SPStd = (0.025, 0.028, 0.030, 0.025,0.018, 0.020)



##------------------- ACC -------------------------##
ACC=(0.6294,0.5972,0.6256,0.7795,0.6286,0.6126)
ACCStd = (0.078, 0.078, 0.075, 0.098,0.075, 0.072)



##------------------- MCC -------------------------##
MCC=(0.1161,0.0743,0.0635,0.5589,0.2573,0.2247)
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
ax.set_title('Feature Evalution using Three classifiers on D5223',fontsize=10)
ax.set_xticks(index + bar_width*1.5 )
ax.set_xticklabels(('RF', 'SVM', 'GPC', 'RF', 'SVM', 'GPC'),fontsize=10)

# ax.legend(fontsize=8,loc='best',bbox_to_anchor=(0.85,0.97))
ax.legend(fontsize=8,loc='best',bbox_to_anchor=(0.83,0.75),framealpha=0.5)

for rect in rects1:
    height = rect.get_height()
    ax.text(rect.get_x() + rect.get_width()/2., 0.97*height,
            '%0.4f' % (height),ha='center', va='bottom',fontsize=8)
    
for rect in rects2:
    height = rect.get_height()
    ax.text(rect.get_x() + rect.get_width()/2., 1.03*height,
            '%0.4f' % (height),ha='center', va='bottom',fontsize=8)    



for rect in rects6:
    height = rect.get_height()
    ax.text(rect.get_x() + rect.get_width()/2., 1.13*height,
            '%0.4f' % (height),ha='center', va='bottom',fontsize=8)         

for rect in rects7:
    height = rect.get_height()
    ax.text(rect.get_x() + rect.get_width()/2., 1.04*height,
            '%0.4f' % (height),ha='center', va='bottom',fontsize=8)  
    

plt.ylim((0, 1.05))
plt.rcParams['figure.figsize'] = (16, 14) # 设置figure_size尺寸
plt.rcParams['savefig.dpi'] = 300 #图片像素
plt.rcParams['figure.dpi'] = 300 #分辨率
# 默认的像素：[6.0,4.0]，分辨率为100，图片尺寸为 600&400
# 指定dpi=200，图片尺寸为 1200*800
# 指定dpi=300，图片尺寸为 1800*1200
# 设置figsize可以在不改变分辨率情况下改变比例
plt.savefig('./data/CV_Independent Test.tif', dpi=300,bbox_inches ='tight') #指定分辨率保存

# fig.tight_layout()

plt.show()

