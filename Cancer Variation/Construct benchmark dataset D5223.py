
# coding: utf-8

# ## 提取突变部分，并把突变前后的AA和num分开

# ## Disease-associated  / SNP
# 
# ## 将fasta 文件，根据突变个数，改写成多个sequences
# 
# 
# ##  总共的序列数：  32540
# ##  总共的蛋白质数：  21558
# 
# ##    --------------  最终取了如下数据：
# ##  总共的序列数：  2608
# ##  总共的蛋白质数：  2247
# 
# 
# ## 每个蛋白最多留了2个突变，预防过拟合

# In[1]:


cancerV=open('./data/original data/Ensembl79_homo_cancer_variation_protein.txt','r')
CV_Mut=open('./data/temp/CV_Mut.txt','w')
CV_pssm=open('./data/temp/CV_pssm.txt','w')##做pssm时需要这样的格式，仅有>,和一个字符串 ##用于存放1029个蛋白质序列：257（bc）+772=1029
Forilearn=open('./data/temp/Forilearn.txt','w')  ##存放ilearn格式的蛋白质序列，并设置75%的training, 25%的testing
CV_MutFN=open('./data/temp/CV_MutFN_5.txt','w')  ##存放ilearn格式的蛋白质序列，并设置75%的training, 25%的testing
ForilearnAA=open('./data/temp/ForilearnAA.txt','w')  ##用于存放11个AA序列，并以ilearn格式存入
proteinIDPOS=open('./data/temp/proteinIDPOS.txt','w')##把proteinID写入文件，用于后面的多个特征文件合并：proteinID作为主键-
PSA=open('./data/temp/PSA.txt','w')  ##求PSA时，要求长度小于1000，这里取了突变点前后各5个AA，用于求PSA
proteinIDPOSAA=open('./data/temp/proteinIDPOSAA.txt','w')##把proteinID,突变位置，突变前后AA，写入文件，用以进行突变预测分析
proteinIDPOSAAfathmm=open('./data/temp/proteinIDPOSAAfathmm.txt','w')##为：fathmm使用的格式，把proteinID,突变位置，突变前后AA，写入文件
proteinID=open('./data/temp/proteinID.txt','w')##把proteinID写入文件，用于提取PSSM,SS,PDO文件


##############################-----------------------处理序列突变数，把突变后为*的突变去掉----------########################
cv_lines=cancerV.readlines()
total_Seqnumcv=0##用于记录总共的序列数
total_proteincv = 0 #用于蛋白质条数
train_number=0 ##用于记录训练集的蛋白质条数
train_number_AA =0 ##用于记录训练集的蛋白质条数:11个AA
ProteinID_CV=[] ##用于存放cancer_variation 的蛋白质ID，用以判断后面的SNP与这个不重复

# for i in range(len(cv_lines)): ## 用于判断有多少个蛋白质和多少个序列。决定使用欧冠多少数据进行试验后，用下面的语句
for i in range(2540):  ##  取2540个序列
        flag = 0 #用于判断该蛋白质有没有写入文件
     ##----------------fasta头部前两项---------------#
        Tou=cv_lines[i].split('\t')[:-2]  ## mut之前的fasta名称,除了突变信息
        number=Tou[0]
        name=Tou[1]
        
    ##----------------sequence的后面*去掉---------------#   
        b=cv_lines[i].split()
        seq=b[-1]
        if (seq[-1]=='*'):
            seq=seq[:-1]
        else:
            seq=seq[:]
        
        temp=cv_lines[i].split('\t')[:-1]##fasta名称中，最后一个\t第表示突变名称和突变信息的
        mut=temp[-1].split(';')#每个突变分开
        k=1##用于记录是第几个突变序列
        mutationN =0 ## 记录该蛋白质中的突变个数，每个蛋白质只取有效的3个突变
        for i in range(len(mut)):
            if(len(mut) < 3): #####################-------如果突变个数小于3，直接写--
                mut2=mut[i].split(':')##突变名臣与突变信息分开
                qian=mut2[1][0]
                hou=mut2[1][-1]
                num=mut2[1][1:-1]
              ###--------- 筛选条件：---------------------------------------------------------##
              ## ---------1.突变前后的AA都不是*，并且(int(num)+6)< len(seq)，即没超过序列长度
              ## ---------2.int(num)-5)> 0，即突变点前面至少有5个AA，用于求微环境好用的
             ## 为了用iLearn计算APAAC, PACC,需要每个蛋白质的长度大于31
                if (hou != '*' and qian !='*'and (int(num)+6)< len(seq) and (int(num)-5)> 0 and len(seq)>31):## 把突变后为*的序列去掉
                    ProteinID_CV.append(name) ###用于存放cancer_variation 的蛋白质ID，用以判断后面的SNP与这个不重复
                    
                    ##------------------------------ 突变后的AA替换掉突变前的AA :然后写入到下面的所有文件中---------------#
                    n2=int(num)-1   ## int(num)-1 :才是突变的真是下标
                    newseq = seq[:n2] + hou + seq[n2+1:] 
                    ##------------------------------ 突变位置前后5个位置的AA  ----保存，用于计算feature邻居的特征 --------------#
                    mutfn=newseq[n2-5: n2+6]  ## 以突变点为中心，提取前后各5个AA---------
                    CV_MutFN.write('>'+str(number)+'\t'+str(name)+'_'+str(k)+'\t')
                    CV_MutFN.write(mutfn+'\n')

                    #---------------------写突变文件-----------------------------##
                    CV_Mut.write('>'+str(number)+'\t'+str(name)+'_'+str(k)+'\t')
                    CV_Mut.write(str(mut2)+'\t'+str(qian)+'\t'+str(num)+'\t'+str(hou)+'\t'+str(newseq)+'\n')
                    #--------------------写PSSM格式需要的文件-------------------##
                    CV_pssm.write('>'+str(name)+'_'+str(k)+'\n')
                    CV_pssm.write(str(newseq)+'\n')
                    
                    #------------------------  将前后5个AA（共11个AA）按照ilearn的格式写入文件 -------#
                    #--------------------for ilearn 写入文件-------------------##
                    ## training = 2608*0.75 = 1956---------
                    if(train_number_AA < 1956):
                        ForilearnAA.write('>'+str(name)+'_'+str(k)+'|'+'1'+'|'+'training'+'\n'+str(mutfn)+'\n')
                    else:
                        ForilearnAA.write('>'+str(name)+'_'+str(k)+'|'+'1'+'|'+'testing'+'\n'+str(mutfn)+'\n')  
                    train_number_AA = train_number_AA+1

                    #--------------------for ilearn 写入文件-------------------##
                     ## training = 2608*0.75 = 1956---------
                    if(train_number <1956):
                        Forilearn.write('>'+str(name)+'_'+str(k)+'|'+'1'+'|'+'training'+'\n'+str(newseq)+'\n')
                    else:
                        Forilearn.write('>'+str(name)+'_'+str(k)+'|'+'1'+'|'+'testing'+'\n'+str(newseq)+'\n')
                    train_number = train_number+1

                    #----------------   将proteinID 写入文件，用于后面的多个特征文件合并：proteinID作为主键--------#
                    proteinIDPOS.write(str(name)+'_'+str(k)+'\t'+str(num)+'\t'+str(len(newseq))+'\n')
                    
                    proteinID.write(str(name)+'_'+str(k)+'\t'+str(num)+'\t'+'1'+'\n')

                    ## -----------------##把proteinID,突变位置，突变前后AA，写入文件，用以进行突变预测分析---------#
                    proteinIDPOSAA.write(str(name)+'_'+str(k)+'\t'+str(num)+'\t'+qian+'\t'+hou+'\n')

                    ##-----------------  ##为：fathmm使用的格式，把proteinID,突变位置，突变前后AA，写入文件 ------#
                    proteinIDPOSAAfathmm.write(str(name)+'_'+str(k)+'\t'+qian+str(num)+hou+'\n')

                    ##-------------##求PSA时，要求长度小于1000，这里取了突变点前后各5个AA，用于求PSA  --------------------##  
                    PSA.write('>'+str(name)+'_'+str(k))
                    PSA.write('\n'+str(mutfn)+'\n')


                    total_Seqnumcv=total_Seqnumcv+1
                    k=k+1 ##用于记录是第几个突变序列
                    flag = 1 #用于判断该蛋白质有没有写入文件：1表示该蛋白质写入了文件
                else:###如果突变为*，则跳过，即不写入
                    continue    

            else:#####################-------如果突变个数大于3，只写前面3个---
                    if(mutationN < 2): ##突变个数小于3，就写
                        mut2=mut[i].split(':')##突变名臣与突变信息分开
                        qian=mut2[1][0]
                        hou=mut2[1][-1]
                        num=mut2[1][1:-1]
                          ###--------- 筛选条件：---------------------------------------------------------##
                          ## ---------1.突变前后的AA都不是*，并且(int(num)+6)< len(seq)，即没超过序列长度
                          ## ---------2.int(num)-5)> 0，即突变点前面至少有5个AA，用于求微环境好用的
                         ## 为了用iLearn计算APAAC, PACC,需要每个蛋白质的长度大于31
                        if (hou != '*' and qian !='*'and (int(num)+6)< len(seq) and (int(num)-5)> 0 and len(seq)>31):## 把突变后为*的序列去掉
                            ProteinID_CV.append(name) ###用于存放cancer_variation 的蛋白质ID，用以判断后面的SNP与这个不重复
                            mutationN =  mutationN +1 ## 每个蛋白质只取有效的3个突变
                            ##------------------------------ 突变后的AA替换掉突变前的AA :然后写入到下面的所有文件中---------------#
                            n2=int(num)-1   ## int(num)-1 :才是突变的真是下标
                            newseq = seq[:n2] + hou + seq[n2+1:] 
                            ##------------------------------ 突变位置前后5个位置的AA  ----保存，用于计算feature邻居的特征 --------------#
                            mutfn=newseq[n2-5: n2+6]  ## 以突变点为中心，提取前后各5个AA---------
                            CV_MutFN.write('>'+str(number)+'\t'+str(name)+'_'+str(k)+'\t')
                            CV_MutFN.write(mutfn+'\n')

                            #---------------------写突变文件-----------------------------##
                            CV_Mut.write('>'+str(number)+'\t'+str(name)+'_'+str(k)+'\t')
                            CV_Mut.write(str(mut2)+'\t'+str(qian)+'\t'+str(num)+'\t'+str(hou)+'\t'+str(newseq)+'\n')

                            #--------------------写PSSM格式需要的文件-------------------##
                            CV_pssm.write('>'+str(name)+'_'+str(k)+'\n')
                            CV_pssm.write(str(newseq)+'\n')
                               

                            #------------------------  将前后5个AA（共11个AA）按照ilearn的格式写入文件 -------#
                            #--------------------for ilearn 写入文件-------------------##
                            ## training = 2608*0.75 = 1956---------
                            if(train_number_AA < 1956):
                                ForilearnAA.write('>'+str(name)+'_'+str(k)+'|'+'1'+'|'+'training'+'\n'+str(mutfn)+'\n')
                            else:
                                ForilearnAA.write('>'+str(name)+'_'+str(k)+'|'+'1'+'|'+'testing'+'\n'+str(mutfn)+'\n')  
                            train_number_AA = train_number_AA+1

                            #--------------------for ilearn 写入文件-------------------##
                             ## training = 2608*0.75 = 1956---------
                            if(train_number <1956):
                                Forilearn.write('>'+str(name)+'_'+str(k)+'|'+'1'+'|'+'training'+'\n'+str(newseq)+'\n')
                            else:
                                Forilearn.write('>'+str(name)+'_'+str(k)+'|'+'1'+'|'+'testing'+'\n'+str(newseq)+'\n')
                            train_number = train_number+1

                            #----------------   将proteinID 写入文件，用于后面的多个特征文件合并：proteinID作为主键--------#
                            proteinIDPOS.write(str(name)+'_'+str(k)+'\t'+str(num)+'\t'+str(len(newseq))+'\n')
                            
                            proteinID.write(str(name)+'_'+str(k)+'\t'+str(num)+'\t'+'1'+'\n')

                            ## -----------------##把proteinID,突变位置，突变前后AA，写入文件，用以进行突变预测分析---------#
                            proteinIDPOSAA.write(str(name)+'_'+str(k)+'\t'+str(num)+'\t'+qian+'\t'+hou+'\n')

                            ##-----------------  ##为：fathmm使用的格式，把proteinID,突变位置，突变前后AA，写入文件 ------#
                            proteinIDPOSAAfathmm.write(str(name)+'_'+str(k)+'\t'+qian+str(num)+hou+'\n')

                            ##-------------##求PSA时，要求长度小于1000，这里取了突变点前后各5个AA，用于求PSA  --------------------##  
                            PSA.write('>'+str(name)+'_'+str(k))
                            PSA.write('\n'+str(mutfn)+'\n')

                            mutationN =mutationN +1
                            total_Seqnumcv=total_Seqnumcv+1
                            k=k+1 ##用于记录是第几个突变序列
                            flag = 1 #用于判断该蛋白质有没有写入文件：1表示该蛋白质写入了文件
                        else:###如果突变为*，则跳过，即不写入
                            continue  
                    else:##突变个数大于3，直接退出这个蛋白质，继续下个蛋白质
                        break
        if(flag == 1): ## 该蛋白质写入了文件
            total_proteincv = total_proteincv+1                
                                
CV_Mut.close()
cancerV.close()
CV_pssm.close()
Forilearn.close()
CV_MutFN.close()
ForilearnAA.close()
proteinIDPOS.close()
PSA.close()
proteinIDPOSAA.close()
proteinIDPOSAAfathmm.close()
print("总共的序列数： ",total_Seqnumcv)
print("总共的蛋白质数： ",total_proteincv)
proteinID.close()


# In[2]:


len(ProteinID_CV)


# ## Ensembl79_homo_dbSNP_variation_protein.txt
# 
# ## 总共的序列数：  210379
# ## 总共的蛋白质数：  76337
#  
# ##    --------------  最终取了如下数据：
# ##  总共的序列数： 2615
# ##  总共的蛋白质数：1377
# 
# 
# 
# ## 每个蛋白最多留了2个突变，预防过拟合

# In[3]:


SNP=open('./data/original data/Ensembl79_homo_dbSNP_variation_protein.txt','r')
SNP_Mut=open('./data/temp/SNP_Mutall.txt','w')
SNP_pssm=open('./data/temp/SNP_pssmall.txt','w')##做pssm时需要这样的格式，仅有>,和一个字符串
CV_pssm=open('./data/temp/CV_pssm.txt','a')##用于存放1029个蛋白质序列：257（bc）+772=1029
Forilearn=open('./data/temp/Forilearn.txt','a') ##存放ilearn格式的蛋白质序列，并设置75%的training, 25%的testing
SNP_MutFN=open('./data/temp/SNP_MutFN_5.txt','w')  ##存放ilearn格式的蛋白质序列，并设置75%的training, 25%的testing
ForilearnAA=open('./data/temp/ForilearnAA.txt','a')  ##用于存放11个AA序列，并以ilearn格式存入
proteinIDPOS=open('./data/temp/proteinIDPOS.txt','a')##把proteinID写入文件，用于后面的多个特征文件合并：proteinID作为主键
PSA=open('./data/temp/PSA.txt','a')##求PSA时，要求长度小于1000，这里取了突变点前后各5个AA，用于求PSA
proteinIDPOSAA=open('./data/temp/proteinIDPOSAA.txt','a')##把proteinID,突变位置，突变前后AA，写入文件，用以进行突变预测分析
proteinIDPOSAAfathmm=open('./data/temp/proteinIDPOSAAfathmm.txt','a')##为：fathmm使用的格式，把proteinID,突变位置，突变前后AA，写入文件
proteinID=open('./data/temp/proteinID.txt','a')##把proteinID写入文件，用于提取PSSM,SS,PDO文件
##############################-----------------------处理序列突变数，把突变后为*的突变去掉----------########################
bc_lines=SNP.readlines()
total_proteinSNP = 0 #用于蛋白质条数
total_SeqnumSNP = 0 #用于记录总共的序列数
train_numberSNP=0 ##用于记录训练集的蛋白质条数
train_numberSNPAA =0 ##用于记录训练集的蛋白质条数:11个AA
count=0 ##只取前773个序列，所有count是用于记录是否到了773个
shuzi=['0','1','2','3','4','5','6','7','8','9'] ##用于判断突变中间数字是否都是数字？
# for i in range(len(bc_lines)):## 用于判断有多少个蛋白质和多少个序列。决定使用欧冠多少数据进行试验后，用下面的语句
for i in range( 1798):  ##  取1798个序列
        flag = 0 #用于判断该蛋白质有没有写入文件
     ##----------------fasta头部前两项---------------#
        Tou=bc_lines[i].split('\t')[:-2]  ## mut之前的fasta名称,除了突变信息
#         print(Tou)
        number=Tou[0]
        name=Tou[1]
        
    ##----------------sequence的后面*去掉---------------#   
        b=bc_lines[i].split()
        seq=b[-1]
        if (seq[-1]=='*'):
            seq=seq[:-1]
        else:
            seq=seq[:]
        
        temp=bc_lines[i].split('\t')[:-1]##fasta名称中，最后一个\t第表示突变名称和突变信息的
        mut=temp[-1].split(';')#每个突变分开
        k=1##用于记录是第几个突变序列
        mutationN =0 ## 记录该蛋白质中的突变个数，每个蛋白质只取有效的3个突变
        for i in range(len(mut)):
            if(len(mut) < 3): #####################-------如果突变个数小于3，直接写--
                youzimu = 0
                mut2=mut[i].split(':')##突变名臣与突变信息分开
                qian=mut2[1][0]
                hou=mut2[1][-1]
                num=mut2[1][1:-1]
                for j in range(len(num)):
                    if(num[j] not in shuzi):#如果突变数字中有字母，则把youzimu设置为1，用于后面判断是否写入文件
                        youzimu=1
                if( youzimu == 0): #如果突变数字中只有数字，则写入文件
                      ###--------- 筛选条件：---------------------------------------------------------##
                      ## ---------1.突变前后的AA都不是*，并且(int(num)+6)< len(seq)，即没超过序列长度
                      ## ---------2.int(num)-5)> 0，即突变点前面至少有5个AA，用于求微环境好用的
                     ## 为了用iLearn计算APAAC, PACC,需要每个蛋白质的长度大于31  
                     ## 保证所有的SNP蛋白质不在CV中出现，即保证蛋白质不同时出现在正负样本中
                      
                    if (hou != '*' and qian !='*'and (int(num)+6)< len(seq) and (int(num)-5)>0 and len(seq)>31 and name not in ProteinID_CV):## 把突变后为*的序列去掉
                            ##------------------------------ 突变后的AA替换掉突变前的AA :然后写入到下面的所有文件中-----------#
                            n2=int(num)-1   ## int(num)-1 :才是突变的真是下标
                            newseq = seq[:n2] + hou + seq[n2+1:] 
                            ##------------------------------ 突变位置前后5个位置的AA  ----保存，用于计算feature邻居的特征 -----#
                            mutfn=newseq[n2-5: n2+6]  ## 以突变点为中心，提取前后各5个AA---------
                            SNP_MutFN.write('>'+str(number)+'\t'+str(name)+'_'+str(k)+'\t')
                            SNP_MutFN.write(mutfn+'\n')                        

                        ###------------------------------------   all    ------------------------##    
                            #---------------------写突变文件-----------------------------##
                            SNP_Mut.write('>'+str(number)+'\t'+str(name)+'_'+str(k)+'\t')
                            SNP_Mut.write(str(mut2)+'\t'+str(qian)+'\t'+str(num)+'\t'+str(hou)+'\t'+str(newseq)+'\n')

                          #--------------------写PSSM格式需要的文件----all---------------##
                            SNP_pssm.write('>'+str(name)+'_'+str(k)+'\t')
                            SNP_pssm.write(str(newseq)+'\n')

                            CV_pssm.write('>'+str(name)+'_'+str(k)+'S'+'\n')
                            CV_pssm.write(str(newseq)+'\n')
                                                
                            #----------------   将proteinID 写入文件，用于后面的多个特征文件合并：proteinID作为主键--------#
                            proteinIDPOS.write(str(name)+'_'+str(k)+'\t'+str(num)+'\t'+str(len(newseq))+'\n')
                            
                            proteinID.write(str(name)+'_'+str(k)+'\t'+str(num)+'\t'+'0'+'\n')
                            
                           ## -----------------##把proteinID,突变位置，突变前后AA，写入文件，用以进行突变预测分析---------#
                            proteinIDPOSAA.write(str(name)+'_'+str(k)+'\t'+str(num)+'\t'+qian+'\t'+hou+'\n')
                            
                            ##-----------------  ##为：fathmm使用的格式，把proteinID,突变位置，突变前后AA，写入文件 ------#
                            proteinIDPOSAAfathmm.write(str(name)+'_'+str(k)+'\t'+qian+str(num)+hou+'\n')
                
                
                           ##-------------##求PSA时，要求长度小于1000，这里取了突变点前后各5个AA，用于求PSA  --------------------##  
                            PSA.write('>'+str(name)+'_'+str(k))
                            PSA.write('\n'+str(mutfn)+'\n')

                            #------------------------  将前后5个AA（共11个AA）按照ilearn的格式写入文件 -------#
                            #--------------------for ilearn 写入文件-------------------##
                            ##------------------------training = 2615*0.75  =1961--------------
                            if (train_numberSNP < 1961):
                                ForilearnAA.write('>'+str(name)+'_'+str(k)+'|'+'0'+'|'+'training'+'\n'+str(mutfn)+'\n')
                            else:
                                ForilearnAA.write('>'+str(name)+'_'+str(k)+'|'+'0'+'|'+'testing'+'\n'+str(mutfn)+'\n')
                            train_numberSNPAA = train_numberSNPAA+1


                            ###---------------------  ilearn 格式的写入文件------------------
                            ##------------------------training = 2615*0.75  =1961-------------
                            if (train_numberSNP < 1961):
                                Forilearn.write('>'+str(name)+'_'+str(k)+'|'+'0'+'|'+'training'+'\n'+str(newseq)+'\n')
                            else:
                                Forilearn.write('>'+str(name)+'_'+str(k)+'|'+'0'+'|'+'testing'+'\n'+str(newseq)+'\n')
                            train_numberSNP = train_numberSNP+1

                            total_SeqnumSNP=total_SeqnumSNP+1
                            k=k+1 ##用于记录是第几个突变序列
                            flag = 1 #用于判断该蛋白质有没有写入文件：1表示该蛋白质写入了文件
                    else:###如果突变为*，则跳过，即不写入
                        continue                         
                else: # 如果突变数字中包含字母，则不写入文件，调到下一个突变信息的判断中
                    continue     
                
            else:
                mut2=mut[i].split(':')##突变名臣与突变信息分开
                qian=mut2[1][0]
                hou=mut2[1][-1]
                num=mut2[1][1:-1]
                youzimu = 0
                if(mutationN < 2): ##突变个数小于3，就写
                    for j in range(len(num)):
                        if(num[j] not in shuzi):#如果突变数字中有字母，则把youzimu设置为1，用于后面判断是否写入文件
                            youzimu=1
                    if(youzimu == 0): #如果突变数字中只有数字，则写入文件
                      ###--------- 筛选条件：---------------------------------------------------------##
                      ## ---------1.突变前后的AA都不是*，并且(int(num)+6)< len(seq)，即没超过序列长度
                      ## ---------2.int(num)-5)> 0，即突变点前面至少有5个AA，用于求微环境好用的
                        if (hou != '*' and qian !='*'and (int(num)+6)< len(seq) and (int(num)-5)> 0 and len(seq)>31 and name not in ProteinID_CV):## 把突变后为*的序列去掉
                            ##------------------------------ 突变后的AA替换掉突变前的AA :然后写入到下面的所有文件中-----------#
                            n2=int(num)-1   ## int(num)-1 :才是突变的真是下标
                            newseq = seq[:n2] + hou + seq[n2+1:] 
                            ##------------------------------ 突变位置前后5个位置的AA  ----保存，用于计算feature邻居的特征 -----#
                            mutfn=newseq[n2-5: n2+6]  ## 以突变点为中心，提取前后各5个AA---------
                            SNP_MutFN.write('>'+str(number)+'\t'+str(name)+'_'+str(k)+'\t')
                            SNP_MutFN.write(mutfn+'\n')                        

                            ###------------------------------------   all    ------------------------##    
                            #---------------------写突变文件-----------------------------##
                            SNP_Mut.write('>'+str(number)+'\t'+str(name)+'_'+str(k)+'\t')
                            SNP_Mut.write(str(mut2)+'\t'+str(qian)+'\t'+str(num)+'\t'+str(hou)+'\t'+str(newseq)+'\n')

                            #--------------------写PSSM格式需要的文件----all---------------##
                            SNP_pssm.write('>'+str(name)+'_'+str(k)+'\t')
                            SNP_pssm.write(str(newseq)+'\n')
                            
                            CV_pssm.write('>'+str(name)+'_'+str(k)+'S'+'\n')
                            CV_pssm.write(str(newseq)+'\n')
                                                
                            #----------------   将proteinID 写入文件，用于后面的多个特征文件合并：proteinID作为主键--------#
                            proteinIDPOS.write(str(name)+'_'+str(k)+'\t'+str(num)+'\t'+str(len(newseq))+'\n')
                        
                            proteinID.write(str(name)+'_'+str(k)+'\t'+str(num)+'\t'+'0'+'\n')
                            
                           ## -----------------##把proteinID,突变位置，突变前后AA，写入文件，用以进行突变预测分析---------#
                            proteinIDPOSAA.write(str(name)+'_'+str(k)+'\t'+str(num)+'\t'+qian+'\t'+hou+'\n')
                            
                            ##-----------------  ##为：fathmm使用的格式，把proteinID,突变位置，突变前后AA，写入文件 ------#
                            proteinIDPOSAAfathmm.write(str(name)+'_'+str(k)+'\t'+qian+str(num)+hou+'\n')
                
                
                           ##-------------##求PSA时，要求长度小于1000，这里取了突变点前后各5个AA，用于求PSA  --------------------##  
                            PSA.write('>'+str(name)+'_'+str(k))
                            PSA.write('\n'+str(mutfn)+'\n')

                            #------------------------  将前后5个AA（共11个AA）按照ilearn的格式写入文件 -------#
                            #--------------------for ilearn 写入文件-------------------##
                            ##------------------------training = 2615*0.75  =1961-------------
                            if (train_numberSNP < 1961):
                                ForilearnAA.write('>'+str(name)+'_'+str(k)+'|'+'0'+'|'+'training'+'\n'+str(mutfn)+'\n')
                            else:
                                ForilearnAA.write('>'+str(name)+'_'+str(k)+'|'+'0'+'|'+'testing'+'\n'+str(mutfn)+'\n')
                            train_numberSNPAA = train_numberSNPAA+1


                            ###---------------------  ilearn 格式的写入文件------------------
                            ##------------------------training = 2615*0.75  =1961--------------
                            if (train_numberSNP < 1961):
                                Forilearn.write('>'+str(name)+'_'+str(k)+'|'+'0'+'|'+'training'+'\n'+str(newseq)+'\n')
                            else:
                                Forilearn.write('>'+str(name)+'_'+str(k)+'|'+'0'+'|'+'testing'+'\n'+str(newseq)+'\n')
                            train_numberSNP = train_numberSNP+1
                            
                            mutationN =mutationN +1
                            total_SeqnumSNP=total_SeqnumSNP+1
                            k=k+1 ##用于记录是第几个突变序列
                            flag = 1 #用于判断该蛋白质有没有写入文件：1表示该蛋白质写入了文件
                        else:###如果突变为*，则跳过，即不写入
                            continue                         
                    else: # 如果突变数字中包含字母，则不写入文件，调到下一个突变信息的判断中
                        continue     
                else:##突变个数大于3，直接退出这个蛋白质，继续下个蛋白质
                     break  
        if(flag == 1): ## 该蛋白质写入了文件
            total_proteinSNP= total_proteinSNP+1                
                                
            
SNP_Mut.close()
SNP_pssm.close()
SNP_MutFN.close()
Forilearn.close()
SNP.close()
ForilearnAA.close()
proteinIDPOS.close()
CV_pssm.close()
PSA.close()
proteinIDPOSAA.close()
proteinIDPOSAAfathmm.close()
print("总共的序列数： ",total_SeqnumSNP)
print("总共的蛋白质数： ",total_proteinSNP)
proteinID.close()

