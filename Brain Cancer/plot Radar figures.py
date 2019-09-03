
"""
author : sera
date: 2019
"""

# coding: utf-8

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
plt.rc('font',family='Times New Roman')
get_ipython().run_line_magic('matplotlib', 'inline')
from IPython.core.pylabtools import figsize # import figsize

def result_pic(result):
    """
    雷达图的绘制
    :param result: 分类数据
    :return: 雷达图
    """
    labels = ['SN', 'SP', 'Precision','NPV','F1', 'ACC', 'MCC']
    kinds = list(result.iloc[:, 0])
    
    result = pd.concat([result, result[['SN']]], axis=1)
    centers = np.array(result.iloc[:, 1:])
    
    # circle
    n = len(labels)
    angle = np.linspace(0, 2 * np.pi, n, endpoint=False)
    angle = np.concatenate((angle, [angle[0]]))

    # plot
    fig = plt.figure()
    ax = fig.add_subplot(111, polar=True)    
    ax.spines['polar'].set_visible(False) 
    ax.set_rlim(0,1) 

    # plot line
    for i in range(len(kinds)):
        ax.plot(angle, centers[i], linewidth=2, label=kinds[i])

    ax.set_thetagrids(angle * 180 / np.pi, labels)
    

    plt.rcParams['savefig.dpi'] = 300 #
    plt.rcParams['figure.dpi'] = 300 
    
    plt.title('Evaluation Values of Different Classifiers')
    ax.set_xlabel('D1209')
    plt.legend(loc=(1.05, 0.05))
    plt.savefig('./data/D1209-CV.jpg', dpi=300) # save figure 
    plt.show()
   

if __name__ == '__main__':

    result = pd.read_csv('./data/D1209-CV.csv', sep=',')
    result_pic(result)
