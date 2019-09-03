
# coding: utf-8

# In[11]:


def comE(y_true,y_pred):
    from sklearn.metrics import confusion_matrix
    #####------------   定义：tn，fp, fn, tp  ------------###
    def tn(y_true, y_pred): 
        return metrics.confusion_matrix(y_true, y_pred)[0, 0]
    def fp(y_true, y_pred): 
        return metrics.confusion_matrix(y_true, y_pred)[0, 1]
    def fn(y_true, y_pred): 
        return metrics.confusion_matrix(y_true, y_pred)[1, 0]
    def tp(y_true, y_pred): 
        return metrics.confusion_matrix(y_true, y_pred)[1, 1]

    TN = tn(y_true, y_pred)
    FP = fp(y_true, y_pred)
    TP = tp(y_true, y_pred)
    FN = fn(y_true, y_pred)

    #sensitivity, recall, hit rate, true positive rate ：TPR = TP / (TP + FN)
    SN = TP*1.0/(TP + FN)*1.0 
    #specificity, true negative rate:TNR = TN / (TN + FP)
    SP = TN*1.0 / (TN + FP)*1.0  
    #precision, prositive predictive value:PPV = TP / (TP + FP)
    precision = TP*1.0 / (TP + FP)*1.0
    #negative predictive value:NPV = TN / (TN + FN)
    NPV = TN*1.0 / (TN + FN)*1.0
    # F1 score is the harmonic mean of precision and sensitivity
    F1= 2*TP / (2*TP + FP+FN)*1.0
    print("TP, FP, FN, TN:\t", TP, FP, FN, TN)
    
    return (SN, SP,precision, NPV,F1)
    


# In[12]:

import pandas as pd 
data = pd.read_csv('./data/bc202f.csv')
y = data['label']
X = data.loc[:,'APAAC1':]

from sklearn.model_selection import train_test_split
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.25, random_state=12)

w11train = pd.DataFrame()
X_train.insert(0,'label',y_train)
w11train=X_train
w11train.to_csv('./data/bc202ftrain.csv')

w11test = pd.DataFrame()
X_test.insert(0,'label',y_test)
w11test=X_test
w11test.to_csv('./data/bc202ftest.csv')


data = pd.read_csv('./data/bc202ftrain.csv')
y = data['label']
X = data.loc[:,'APAAC1':]

data = pd.read_csv('./data/bc202ftest.csv')
yy = data['label']
XX = data.loc[:,'APAAC1':]
print('数据加载完毕！')


# In[21]:


from sklearn.model_selection import cross_val_score
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import VotingClassifier
from sklearn.linear_model.stochastic_gradient import SGDClassifier
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import RBF
from sklearn.tree import DecisionTreeClassifier
from sklearn import metrics
from sklearn.ensemble import AdaBoostClassifier

clf1 =SVC(probability=False, C=9, gamma=0.15, kernel='rbf')
clf2 = RandomForestClassifier( criterion='gini', n_estimators=34 ,random_state=12)
clf3 = GaussianProcessClassifier(kernel = 1.0 * RBF(1.0),n_restarts_optimizer=1,
                              max_iter_predict=50,random_state=2)

eclf = VotingClassifier(estimators=[('svc', clf1), ('rf', clf2), ('gpc', clf3)], voting='hard',
                        weights = [2,5,2])
# for clf, label in zip([clf1, clf2, clf3, eclf], ['SVC', 'RF', 'GPC', 'Ensemble']):
#     scores = cross_val_score(clf, X, y, cv=10, scoring='accuracy')
#     print("Accuracy: %0.2f (+/- %0.2f) [%s]" % (scores.mean(), scores.std(), label))
eclf4 = eclf.fit(X, y)

   
    
print('\n---------Independent test set ----------\n ')
yy_true,yy_pred = yy,eclf4.predict(XX)
print("\tmatthews_corrcoef: %.4f\n" % metrics.matthews_corrcoef(yy_true, yy_pred))
MCC=metrics.matthews_corrcoef(yy_true, yy_pred)
print("\taccuracy_score: %.4f\n" % metrics.accuracy_score(yy_true, yy_pred))
ACC=metrics.accuracy_score(yy_true, yy_pred)
from sklearn.metrics import confusion_matrix
tn, fp, fn, tp = confusion_matrix(yy_true, yy_pred).ravel()
tn, fp, fn, tp
SN = tp/(tp+fn)
SP = tn/(tn+fp)
precision = tp/(tp+fp)
NPV = tn/(tn+fn)
F1 = 2*tp / (2*tp + fp+fn)*1.0
print("SN : %.4f \n"% SN)
print("SP :%.4f\n"% SP)
print("precision : %.4f \n"% precision)
print("NPV :%.4f \n"% NPV)
print("F1: %.4f \n"% F1)
    