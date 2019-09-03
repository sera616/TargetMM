
# coding: utf-8

# ## DT

# In[27]:

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
    SN = TP*1.0/(TP + FN)*1.0 ## 也就是：SN
    #specificity, true negative rate:TNR = TN / (TN + FP)
    SP = TN / (TN + FP)*1.0  ## 也就是：SP
    #precision, prositive predictive value:PPV = TP / (TP + FP)
    precision = TP / (TP + FP)*1.0
    #negative predictive value:NPV = TN / (TN + FN)
    NPV = TN / (TN + FN)*1.0
    # F1 score is the harmonic mean of precision and sensitivity
    F1= 2*TP / (2*TP + FP+FN)*1.0
    
    return SN, SP,precision, NPV,F1


# In[43]:

####----------------- find out the optimal parameters for classifiers------------------##
import pandas as pd 
data = pd.read_csv('./data/result/PSSM/cv202f.csv')
y = data['label']
X = data.loc[:,'APAAC1':]

from sklearn.model_selection import train_test_split
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.25, random_state=12)

## --  split train and test data  --------------------##
cv202ftrain = pd.DataFrame()

X_train.insert(0,'label',y_train)
cv202ftrain=X_train
cv202ftrain.to_csv('./data/result/PSSM/cv202ftrain.csv')

cv202ftest = pd.DataFrame()
X_test.insert(0,'label',y_test)
cv202ftest=X_test
cv202ftest.to_csv('./data/result/PSSM/cv202ftest.csv')


# ## SVM
# In[59]:
data = pd.read_csv('./data/result/PSSM/cv202ftrain.csv')
y = data['label']
X = data.loc[:,'APAAC1':]


from sklearn.model_selection import KFold,StratifiedKFold
from sklearn.svm import SVC

C_range = [9,9.4,9.8,10,10.5,10.8]
gamma_range = [0.05,0.08,0.11,0.12,0.15]
kernel_range = ['rbf','linear','sigmoid','poly']
    
cv = StratifiedKFold(n_splits=10) 
 
MCCmean=[]
ACCmean=[]
SNmean=[]
SPmean=[]
precisionmean=[]
NPVmean=[]
F1mean=[]

paramArray1=[]
paramArray2=[]
paramArray3=[]
aa=open('./data/result/PSSM/cv202f_evalution','a')
aa.write("\n\n\n\n-------------------SVM------------------------\n")
aa.write("\n-------------------10-fold corss-validation result--------------------\n")
for i in C_range: 
    for g in gamma_range:
        for j in kernel_range:
#             print("\tthe param is :%s" % i)
#             print("\tthe param is :%s" % j)
            MCC=[]
            ACC=[]
            SN =[]
            SP=[]
            precision=[]
            NPV=[]
            F1=[]

            k=1  
            for train, test in cv.split(X, y):
                svc =SVC(probability=False, C=i, gamma=g, kernel=j)
                y_true,y_pred =y[test], svc.fit(X.iloc[train], y[train]).predict(X.iloc[test])
    #             print("\tmatthews_corrcoef: %1.3f" % metrics.matthews_corrcoef(y_true, y_pred))
                MCCv=metrics.matthews_corrcoef(y_true, y_pred)
                MCC.append(MCCv)
        #         print("\taccuracy_score: %1.3f\n" % metrics.accuracy_score(y_true, y_pred))
                ACCv=metrics.accuracy_score(y_true, y_pred)
                ACC.append(ACCv)
                SNv,SPv,precisionv,NPVv,F1v = comE(y_true, y_pred)
                SN.append(SNv)
                SP.append(SPv)
                precision.append(precisionv)
                NPV.append(NPVv)
                F1.append(F1v)
                
                aa.write("参数组合为:(C,gamma,kernel)\t"+str(i)+'\t'+str(g)+'\t'+str(j)+'\t')
                aa.write("第"+str(k)+"折交叉验证"+'\t')
                aa.write("\tMCC: \t"+str(MCCv))
                aa.write("\tACC: \t"+str(ACCv))
                aa.write("\tSN: \t"+str(SNv))
                aa.write("\tSP: \t"+str(SPv))
                aa.write("\tprecision: \t"+str(precisionv))
                aa.write("\tNPV: \t"+str(NPVv))
                aa.write("\tF1: \t"+str(F1v))
                aa.write('\n')

                k=k+1  

#             print("MCC mean: %.3f +/- %.3f :"% (np.mean(MCC),np.std(MCC)))
#             print("ACC mean: %.3f +/- %.3f :"% (np.mean(ACC),np.std(ACC)))
#             print("SN mean: %.3f +/- %.3f :"% (np.mean(SN),np.std(SN)))
#             print("SP mean: %.3f +/- %.3f :"% (np.mean(SP),np.std(SP)))
#             print("precision mean: %.3f +/- %.3f :"% (np.mean(precision),np.std(precision)))
#             print("NPV mean: %.3f +/- %.3f :"% (np.mean(NPV),np.std(NPV)))
#             print("F1 mean: %.3f +/- %.3f :"% (np.mean(F1),np.std(F1)))
           
            MCCmean.append(np.mean(MCC))
            ACCmean.append(np.mean(ACC))
            SNmean.append(np.mean(SN))
            SPmean.append(np.mean(SP))
            precisionmean.append(np.mean(precision))
            NPVmean.append(np.mean(NPV))
            F1mean.append(np.mean(F1))

            paramArray1.append(i)
            paramArray2.append(g)
            paramArray3.append(j)
            
            ####--------------------------------
            aa.write('该组参数下的评价指标平均值：\n')
            aa.write("\tMCCmean: \t"+str(np.mean(MCC)))
            aa.write("\tACCmean: \t"+str(np.mean(ACC)))
            aa.write("\tSNmean: \t"+str(np.mean(SN)))
            aa.write("\tSPvmean: \t"+str(np.mean(SP)))
            aa.write("\tprecisionmean: \t"+str(np.mean(precision)))
            aa.write("\tNPVmean: \t"+str(np.mean(NPV)))
            aa.write("\tF1mean: \t"+str(np.mean(F1)))
            aa.write('\n\n\n')
        

import numpy as np
maxMCC=np.max(MCCmean) 
index=MCCmean.index(maxMCC)
# print(maxMCC)
# print(paramArray1[index])
# print(paramArray2[index])
# print(paramArray3[index])


param1=paramArray1[index]
param2=paramArray2[index]
param3=paramArray3[index]

aa.write("\n-------------------------------   SVM  --------------------------\n")
aa.write("\n-----The maximum MCC_mean on the cross-validation set with the corresponding parameters-----\n")
aa.write("MCCmean：\t"+str(maxMCC))
aa.write("\tACCmean：\t"+str(ACCmean[index]))
aa.write("\tSNmean：\t"+str(SNmean[index]))
aa.write("\tSPmean：\t"+str(SPmean[index]))
aa.write("\tprecisionmean：\t"+str(precisionmean[index]))
aa.write("\tNPVmean：\t"+str(NPVmean[index]))
aa.write("\tF1mean：\t"+str(F1mean[index]))
aa.write("\n对应的参数组合为：\t"+"C=:\t"+str(param1)+'\t'+"gamma=:\t"+str(param2)+'\t'+"kernel=:\t"+str(param3))


svc =SVC(probability=False, C=param1, gamma=param2, kernel=param3)
svc.fit(X,y)


data = pd.read_csv('./data/result/PSSM/cv202ftest.csv')
yy = data['label']
XX = data.loc[:,'APAAC1':]

yy_true,yy_pred=yy,svc.predict(XX)
# print("\tmatthews_corrcoef: %1.3f" % metrics.matthews_corrcoef(yy_true, yy_pred))
MCCT=metrics.matthews_corrcoef(yy_true, yy_pred)
# print("\taccuracy_score: %1.3f" % metrics.accuracy_score(yy_true, yy_pred))
ACCT= metrics.accuracy_score(yy_true, yy_pred)

SNT,SPT,precisionT,NPVT,F1T = comE(yy_true, yy_pred) 
# print("\t SN : %1.3f" % SNT) 
# print("\t SP : %1.3f" % SPT) 
# print("\t Precision : %1.3f" % precisionT) 
# print("\t NPV: %1.3f" % NPVT) 
# print("\t F1: %1.3f" % F1T) 

aa.write("\n\n\n-------------------independent test result-------------\n")
aa.write("\tMCC: \t"+str(MCCT))
aa.write("\tACC: \t"+str(ACCT))
aa.write("\tSN: \t"+str(SNT))
aa.write("\tSP: \t"+str(SPT))
aa.write("\tprecision: \t"+str(precisionT))
aa.write("\tNPV: \t"+str(NPVT))
aa.write("\tF1: \t"+str(F1T))
aa.close()




# ## RF 
# In[53]:
## 第二步

data = pd.read_csv('./data/result/PSSM/cv202ftrain.csv')
y = data['label']
X = data.loc[:,'APAAC1':]

from sklearn.model_selection import KFold,StratifiedKFold
from sklearn.ensemble import RandomForestClassifier
from sklearn import metrics
import numpy as np

criterion_range = ['gini','entropy']  
n_estimators_range=[ 10,  12,  14,  16,  18,  20,  22,  24,  26,  28,  30,  32,  34,
            36,  38,  40,  42,  44,  46,  48,  50,  52,  54,  56,  58,  60,
            62,  64,  66,  68,  70,  72,  74,  76,  78,  80,  82,  84,  86,
            88,  90,  92,  94,  96,  98, 100]

cv = StratifiedKFold(n_splits=10) 
 
MCCmean=[]
ACCmean=[]
SNmean=[]
SPmean=[]
precisionmean=[]
NPVmean=[]
F1mean=[]

paramArray1=[]
paramArray2=[]
aa=open('./data/result/PSSM/cv202f_evalution','a')
aa.write("\n-------------------RF-------------\n")
aa.write("\n-------------------10-fold corss-validation result--------------------\n")
for i in criterion_range: 
    for j in n_estimators_range:
#         print("\tthe param is :%s" % i)
#         print("\tthe param is :%s" % j)
        MCC=[]
        ACC=[]
        SN =[]
        SP=[]
        precision=[]
        NPV=[]
        F1=[]
        
        k=1  
        for train, test in cv.split(X, y):
            rfc = RandomForestClassifier( criterion=i, n_estimators=j )
            y_true,y_pred =y[test], rfc.fit(X.iloc[train], y[train]).predict(X.iloc[test])
#             print("\tmatthews_corrcoef: %1.3f" % metrics.matthews_corrcoef(y_true, y_pred))
            MCCv=metrics.matthews_corrcoef(y_true, y_pred)
            MCC.append(MCCv)
    #         print("\taccuracy_score: %1.3f\n" % metrics.accuracy_score(y_true, y_pred))
            ACCv=metrics.accuracy_score(y_true, y_pred)
            ACC.append(ACCv)
            
            SNv,SPv,precisionv,NPVv,F1v = comE(y_true, y_pred)
            ## y_true, y_pred 
            SN.append(SNv)
            SP.append(SPv)
            precision.append(precisionv)
            NPV.append(NPVv)
            F1.append(F1v)
            
            aa.write("参数组合为:(criterion,n_estimators)\t"+str(i)+'\t'+str(j)+'\t')
            aa.write("第"+str(k)+"折交叉验证"+'\t')
            aa.write("\tMCC: \t"+str(MCCv))
            aa.write("\tACC: \t"+str(ACCv))
            aa.write("\tSN: \t"+str(SNv))
            aa.write("\tSP: \t"+str(SPv))
            aa.write("\tprecision: \t"+str(precisionv))
            aa.write("\tNPV: \t"+str(NPVv))
            aa.write("\tF1: \t"+str(F1v))
            aa.write('\n')
            
            k=k+1  
            
#         print("MCC mean: %.3f +/- %.3f :"% (np.mean(MCC),np.std(MCC)))
#         print("ACC mean: %.3f +/- %.3f :"% (np.mean(ACC),np.std(ACC)))
#         print("SN mean: %.3f +/- %.3f :"% (np.mean(SN),np.std(SN)))
#         print("SP mean: %.3f +/- %.3f :"% (np.mean(SP),np.std(SP)))
#         print("precision mean: %.3f +/- %.3f :"% (np.mean(precision),np.std(precision)))
#         print("NPV mean: %.3f +/- %.3f :"% (np.mean(NPV),np.std(NPV)))
#         print("F1 mean: %.3f +/- %.3f :"% (np.mean(F1),np.std(F1)))
        
        MCCmean.append(np.mean(MCC))
        ACCmean.append(np.mean(ACC))
        SNmean.append(np.mean(SN))
        SPmean.append(np.mean(SP))
        precisionmean.append(np.mean(precision))
        NPVmean.append(np.mean(NPV))
        F1mean.append(np.mean(F1))
        
        paramArray1.append(i)
        paramArray2.append(j)
        
        aa.write('该组参数下的评价指标平均值：\n')
        aa.write("\tMCCmean: \t"+str(np.mean(MCC)))
        aa.write("\tACCmean: \t"+str(np.mean(ACC)))
        aa.write("\tSNmean: \t"+str(np.mean(SN)))
        aa.write("\tSPvmean: \t"+str(np.mean(SP)))
        aa.write("\tprecisionmean: \t"+str(np.mean(precision)))
        aa.write("\tNPVmean: \t"+str(np.mean(NPV)))
        aa.write("\tF1mean: \t"+str(np.mean(F1)))
        aa.write('\n\n\n')
        


maxMCC=np.max(MCCmean) 
index=MCCmean.index(maxMCC)
# print(maxMCC)
# print(paramArray1[index])
# print(paramArray2[index])

param1=paramArray1[index]
param2=paramArray2[index]

aa.write("\n-------------------------------   RF  --------------------------\n")
aa.write("\n-----The maximum MCC_mean on the cross-validation set with the corresponding parameters-----\n")
aa.write("MCCmean：\t"+str(maxMCC))
aa.write("\tACCmean：\t"+str(ACCmean[index]))
aa.write("\tSNmean：\t"+str(SNmean[index]))
aa.write("\tSPmean：\t"+str(SPmean[index]))
aa.write("\tprecisionmean：\t"+str(precisionmean[index]))
aa.write("\tNPVmean：\t"+str(NPVmean[index]))
aa.write("\tF1mean：\t"+str(F1mean[index]))
aa.write("\n对应的参数组合为：\t"+"criterion:\t"+str(param1)+'\t\t'+"n_estimators:\t"+str(param2))


rfc = RandomForestClassifier(criterion=param1,n_estimators=param2)
rfc.fit(X,y)

data = pd.read_csv('./data/result/PSSM/cv202ftest.csv')
yy = data['label']
XX = data.loc[:,'APAAC1':]

yy_true,yy_pred=yy,rfc.predict(XX)
# print("\tmatthews_corrcoef: %1.3f" % metrics.matthews_corrcoef(yy_true, yy_pred))
MCCT=metrics.matthews_corrcoef(yy_true, yy_pred)
# print("\taccuracy_score: %1.3f" % metrics.accuracy_score(yy_true, yy_pred))
ACCT= metrics.accuracy_score(yy_true, yy_pred)

SNT,SPT,precisionT,NPVT,F1T = comE(yy_true, yy_pred) 
# print("\t SN : %1.3f" % SNT) 
# print("\t SP : %1.3f" % SPT) 
# print("\t Precision : %1.3f" % precisionT) 
# print("\t NPV: %1.3f" % NPVT) 
# print("\t F1: %1.3f" % F1T) 


aa.write("\n\n\n------------------- independent test result -------------\n")
aa.write("\tMCC: \t"+str(MCCT))
aa.write("\tACC: \t"+str(ACCT))
aa.write("\tSN: \t"+str(SNT))
aa.write("\tSP: \t"+str(SPT))
aa.write("\tprecision: \t"+str(precisionT))
aa.write("\tNPV: \t"+str(NPVT))
aa.write("\tF1: \t"+str(F1T))
aa.close()



# ## GPC
# In[58]:
## 第二步
data = pd.read_csv('./data/result/PSSM/cv202ftrain.csv')
y = data['label']
X = data.loc[:,'APAAC1':]


from sklearn.model_selection import KFold,StratifiedKFold
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import RBF
from sklearn import metrics
import numpy as np

n_restarts_optimizerR = [0,1]
max_iter_predictR =[50,100]

cv = StratifiedKFold(n_splits=10) 
 
MCCmean=[]
ACCmean=[]
SNmean=[]
SPmean=[]
precisionmean=[]
NPVmean=[]
F1mean=[]

paramArray1=[]
paramArray2=[]
aa=open('./data/result/PSSM/cv202f_evalution','a')
aa.write("\n-------------------GPC------------\n")
aa.write("\n-------------------10-fold corss-validation result--------------------\n")
for i in n_restarts_optimizerR: 
    for j in max_iter_predictR:
#         print("\tthe param is :%s" % i)
#         print("\tthe param is :%s" % j)
        MCC=[]
        ACC=[]
        SN =[]
        SP=[]
        precision=[]
        NPV=[]
        F1=[]
        
        k=1  
        for train, test in cv.split(X, y):
            gpc=GaussianProcessClassifier(kernel = 1.0 * RBF(1.0),n_restarts_optimizer=i,max_iter_predict=j)
            y_true,y_pred =y[test], gpc.fit(X.iloc[train], y[train]).predict(X.iloc[test])
#             print("\tmatthews_corrcoef: %1.3f" % metrics.matthews_corrcoef(y_true, y_pred))
            MCCv=metrics.matthews_corrcoef(y_true, y_pred)
            MCC.append(MCCv)
    #         print("\taccuracy_score: %1.3f\n" % metrics.accuracy_score(y_true, y_pred))
            ACCv=metrics.accuracy_score(y_true, y_pred)
            ACC.append(ACCv)
            
            SNv,SPv,precisionv,NPVv,F1v = comE(y_true, y_pred)
            ## y_true, y_pred 
            SN.append(SNv)
            SP.append(SPv)
            precision.append(precisionv)
            NPV.append(NPVv)
            F1.append(F1v)
            aa.write("参数组合为(n_restarts_optimizer,max_iter_predict):\t"+str(i)+'\t'+str(j)+'\t')
            aa.write("第"+str(k)+"折交叉验证"+'\t')
            aa.write("\tMCC: \t"+str(MCCv))
            aa.write("\tACC: \t"+str(ACCv))
            aa.write("\tSN: \t"+str(SNv))
            aa.write("\tSP: \t"+str(SPv))
            aa.write("\tprecision: \t"+str(precisionv))
            aa.write("\tNPV: \t"+str(NPVv))
            aa.write("\tF1: \t"+str(F1v))
            aa.write('\n')
            
            k=k+1  
            
#         print("MCC mean: %.3f +/- %.3f :"% (np.mean(MCC),np.std(MCC)))
#         print("ACC mean: %.3f +/- %.3f :"% (np.mean(ACC),np.std(ACC)))
#         print("SN mean: %.3f +/- %.3f :"% (np.mean(SN),np.std(SN)))
#         print("SP mean: %.3f +/- %.3f :"% (np.mean(SP),np.std(SP)))
#         print("precision mean: %.3f +/- %.3f :"% (np.mean(precision),np.std(precision)))
#         print("NPV mean: %.3f +/- %.3f :"% (np.mean(NPV),np.std(NPV)))
#         print("F1 mean: %.3f +/- %.3f :"% (np.mean(F1),np.std(F1)))
        MCCmean.append(np.mean(MCC))
        ACCmean.append(np.mean(ACC))
        SNmean.append(np.mean(SN))
        SPmean.append(np.mean(SP))
        precisionmean.append(np.mean(precision))
        NPVmean.append(np.mean(NPV))
        F1mean.append(np.mean(F1))
        
        paramArray1.append(i)
        paramArray2.append(j)
       
        aa.write('该组参数下的评价指标平均值：\n')
        aa.write("\tMCCmean: \t"+str(np.mean(MCC)))
        aa.write("\tACCmean: \t"+str(np.mean(ACC)))
        aa.write("\tSNmean: \t"+str(np.mean(SN)))
        aa.write("\tSPvmean: \t"+str(np.mean(SP)))
        aa.write("\tprecisionmean: \t"+str(np.mean(precision)))
        aa.write("\tNPVmean: \t"+str(np.mean(NPV)))
        aa.write("\tF1mean: \t"+str(np.mean(F1)))
        aa.write('\n\n\n')
        


maxMCC=np.max(MCCmean) 
index=MCCmean.index(maxMCC)
# print(maxMCC)
# print(paramArray1[index])
# print(paramArray2[index])

param1=paramArray1[index]
param2=paramArray2[index]

aa.write("\n-------------------------------  GPC  --------------------------\n")
aa.write("\n-----The maximum MCC_mean on the cross-validation set with the corresponding parameters-----\n")
aa.write("MCCmean：\t"+str(maxMCC))
aa.write("\tACCmean：\t"+str(ACCmean[index]))
aa.write("\tSNmean：\t"+str(SNmean[index]))
aa.write("\tSPmean：\t"+str(SPmean[index]))
aa.write("\tprecisionmean：\t"+str(precisionmean[index]))
aa.write("\tNPVmean：\t"+str(NPVmean[index]))
aa.write("\tF1mean：\t"+str(F1mean[index]))
aa.write("\n对应的参数组合为：\t"+"n_restarts_optimizer:\t"+str(param1)+'\t\t'+"max_iter_predict:\t"+str(param2))

gpc=GaussianProcessClassifier(kernel = 1.0 * RBF(1.0),n_restarts_optimizer=param1,max_iter_predict=param2)
gpc.fit(X,y)


data = pd.read_csv('./data/result/PSSM/cv202ftest.csv')
yy = data['label']
XX = data.loc[:,'APAAC1':]

yy_true,yy_pred=yy,gpc.predict(XX)
# print("\tmatthews_corrcoef: %1.3f" % metrics.matthews_corrcoef(yy_true, yy_pred))
MCCT=metrics.matthews_corrcoef(yy_true, yy_pred)
# print("\taccuracy_score: %1.3f" % metrics.accuracy_score(yy_true, yy_pred))
ACCT= metrics.accuracy_score(yy_true, yy_pred)

SNT,SPT,precisionT,NPVT,F1T = comE(yy_true, yy_pred) 
# print("\t SN : %1.3f" % SNT) 
# print("\t SP : %1.3f" % SPT) 
# print("\t Precision : %1.3f" % precisionT) 
# print("\t NPV: %1.3f" % NPVT) 
# print("\t F1: %1.3f" % F1T) 

aa.write("\n\n\n-------------------indepedent test result------------\n")
aa.write("\tMCC: \t"+str(MCCT))
aa.write("\tACC: \t"+str(ACCT))
aa.write("\tSN: \t"+str(SNT))
aa.write("\tSP: \t"+str(SPT))
aa.write("\tprecision: \t"+str(precisionT))
aa.write("\tNPV: \t"+str(NPVT))
aa.write("\tF1: \t"+str(F1T))
aa.close()

