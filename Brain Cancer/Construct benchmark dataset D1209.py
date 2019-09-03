
# coding: utf-8

# ## 提取突变部分，并把突变前后的AA和num分开
# ## brain cancer
# 
# ## 将fasta 文件，根据突变个数，改写成多个sequences
# 
# ##    --------------  最终取了如下数据：
# ##  总共的序列数：  255
# ##  总共的蛋白质数： 177
# 
# ## 每个蛋白最多留了2个突变，预防过拟合

# In[1]:

bcancer=open('./data/original data/brain cancer.txt','r')
BC_Mut=open('./data/temp/BC_Mut.txt','w')
BC_pssm=open('./data/temp/BC_pssm.txt','w')##做pssm时需要这样的格式，仅有>,和一个字符串 ##用于存放1029个蛋白质序列：257（bc）+772=1029
Forilearn=open('./data/temp/Forilearn.txt','w')  ##存放ilearn格式的蛋白质序列，并设置75%的training, 25%的testing
BC_MutFN=open('./data/temp/BC_MutFN_5.txt','w')  ##存放ilearn格式的蛋白质序列，并设置75%的training, 25%的testing
ForilearnAA=open('./data/temp/ForilearnAA.txt','w')  ##用于存放11个AA序列，并以ilearn格式存入
proteinIDPOS=open('./data/temp/proteinIDPOS.txt','w')##把proteinID写入文件，用于后面的多个特征文件合并：proteinID作为主键-
PSA=open('./data/temp/PSA1030.txt','w')  ##求PSA时，要求长度小于1000，这里取了突变点前后各5个AA，用于求PSA
proteinIDPOSAA=open('./data/temp/proteinIDPOSAA.txt','w')##把proteinID,突变位置，突变前后AA，写入文件，用以进行突变预测分析
proteinIDPOSAAfathmm=open('./data/temp/proteinIDPOSAAfathmm.txt','w')##为：fathmm使用的格式，把proteinID,突变位置，突变前后AA，写入文件

##############################-----------------------处理序列突变数，把突变后为*的突变去掉----------########################
bc_lines=bcancer.readlines()
total_Seqnumbc=0##用于记录总共的序列数
total_proteinbc = 0 #用于蛋白质条数
train_number=0 ##用于记录训练集的蛋白质条数
train_number_AA =0 ##用于记录训练集的蛋白质条数:11个AA
ProteinID_BC=[] ##用于存放brain cancer的蛋白质ID，用以判断后面的SNP与这个不重复

for i in range(len(bc_lines)):
        flag = 0 #用于判断该蛋白质有没有写入文件
     ##----------------fasta头部前两项---------------#
        Tou=bc_lines[i].split('\t')[:-2]  ## mut之前的fasta名称,除了突变信息
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
                mut2=mut[i].split(':')##突变名臣与突变信息分开
                qian=mut2[1][0]
                hou=mut2[1][-1]
                num=mut2[1][1:-1]
              ###--------- 筛选条件：---------------------------------------------------------##
              ## ---------1.突变前后的AA都不是*，并且(int(num)+6)< len(seq)，即没超过序列长度
              ## ---------2.int(num)-5)> 0，即突变点前面至少有5个AA，用于求微环境好用的
             ## 为了用iLearn计算APAAC, PACC,需要每个蛋白质的长度大于31
                if (hou != '*' and qian !='*'and (int(num)+6)< len(seq) and (int(num)-5)> 0 and len(seq)>31):## 把突变后为*的序列去掉
                    ProteinID_BC.append(name) ###用于存放cancer_variation 的蛋白质ID，用以判断后面的SNP与这个不重复
                    ##------------------------------ 突变后的AA替换掉突变前的AA :然后写入到下面的所有文件中---------------#
                    n2=int(num)-1   ## int(num)-1 :才是突变的真是下标
                    newseq = seq[:n2] + hou + seq[n2+1:] 
                    ##------------------------------ 突变位置前后5个位置的AA  ----保存，用于计算feature邻居的特征 --------------#
                    mutfn=newseq[n2-5: n2+6]  ## 以突变点为中心，提取前后各5个AA---------
                    BC_MutFN.write('>'+str(number)+'\t'+str(name)+'_'+str(k)+'\t')
                    BC_MutFN.write(mutfn+'\n')
                    
                    #---------------------写突变文件-----------------------------##
                    BC_Mut.write('>'+str(number)+'\t'+str(name)+'_'+str(k)+'\t')
                    BC_Mut.write(str(mut2)+'\t'+str(qian)+'\t'+str(num)+'\t'+str(hou)+'\t'+str(newseq)+'\n')

                    #--------------------写PSSM格式需要的文件-------------------##
                    BC_pssm.write('>'+str(name)+'_'+str(k)+'C'+'\n')
                    BC_pssm.write(str(newseq)+'\n')

                    #------------------------  将前后5个AA（共11个AA）按照ilearn的格式写入文件 -------#
                    #--------------------for ilearn 写入文件-------------------##
                    ## training = 255*0.75 = 191  ---------
                    if(train_number_AA < 191):
                        ForilearnAA.write('>'+str(name)+'_'+str(k)+'|'+'1'+'|'+'training'+'\n'+str(mutfn)+'\n')
                    else:
                        ForilearnAA.write('>'+str(name)+'_'+str(k)+'|'+'1'+'|'+'testing'+'\n'+str(mutfn)+'\n')  
                    train_number_AA = train_number_AA+1
 
                    #--------------------for ilearn 写入文件-------------------##
                    ## training = 257*0.75 = 192  ---------
                    if(train_number < 192):
                        Forilearn.write('>'+str(name)+'_'+str(k)+'|'+'1'+'|'+'training'+'\n'+str(newseq)+'\n')
                    else:
                        Forilearn.write('>'+str(name)+'_'+str(k)+'|'+'1'+'|'+'testing'+'\n'+str(newseq)+'\n')                    
                    

                    #----------------   将proteinID 写入文件，用于后面的多个特征文件合并：proteinID作为主键--------#
                    proteinIDPOS.write(str(name)+'_'+str(k)+'C'+'\t'+str(num)+'\t'+'1'+'\n')

#                     ## -----------------##把proteinID,突变位置，突变前后AA，写入文件，用以进行突变预测分析---------#
#                     proteinIDPOSAA.write(str(name)+'\t'+str(num)+'\t'+qian+'\t'+hou+'\n')

#                     ##-----------------  ##为：fathmm使用的格式，把proteinID,突变位置，突变前后AA，写入文件 ------#
#                     proteinIDPOSAAfathmm.write(str(name)+'\t'+qian+str(num)+hou+'\n')
                      ## -----------------##把proteinID,突变位置，突变前后AA，写入文件，用以进行突变预测分析---------#
                    proteinIDPOSAA.write(str(name)+'_'+str(k)+'\t'+str(num)+'\t'+qian+'\t'+hou+'\n')

                      ##-----------------  ##为：fathmm使用的格式，把proteinID,突变位置，突变前后AA，写入文件 ------#
                    proteinIDPOSAAfathmm.write(str(name)+'_'+str(k)+'\t'+qian+str(num)+hou+'\n')

                    ##-------------##求PSA时，要求长度小于1000，这里取了突变点前后各5个AA，用于求PSA  --------------------##  
                    PSA.write('>'+str(name)+'_'+str(k))
                    PSA.write('\n'+str(mutfn)+'\n') 

                    total_Seqnumbc=total_Seqnumbc+1
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
                            ProteinID_BC.append(name) ###用于存放cancer_variation 的蛋白质ID，用以判断后面的SNP与这个不重复
                            ##------------------------------ 突变后的AA替换掉突变前的AA :然后写入到下面的所有文件中---------------#
                            n2=int(num)-1   ## int(num)-1 :才是突变的真是下标
                            newseq = seq[:n2] + hou + seq[n2+1:] 
                            ##------------------------------ 突变位置前后5个位置的AA  ----保存，用于计算feature邻居的特征 --------------#
                            mutfn=newseq[n2-5: n2+6]  ## 以突变点为中心，提取前后各5个AA---------
                            BC_MutFN.write('>'+str(number)+'\t'+str(name)+'_'+str(k)+'\t')
                            BC_MutFN.write(mutfn+'\n')

                            #---------------------写突变文件-----------------------------##
                            BC_Mut.write('>'+str(number)+'\t'+str(name)+'_'+str(k)+'\t')
                            BC_Mut.write(str(mut2)+'\t'+str(qian)+'\t'+str(num)+'\t'+str(hou)+'\t'+str(newseq)+'\n')

                            #--------------------写PSSM格式需要的文件-------------------##
                            BC_pssm.write('>'+str(name)+'_'+str(k)+'C'+'\n')
                            BC_pssm.write(str(newseq)+'\n')

                            #------------------------  将前后5个AA（共11个AA）按照ilearn的格式写入文件 -------#
                            #--------------------for ilearn 写入文件-------------------##
                            ## training = 255*0.75 = 191  ---------
                            if(train_number_AA < 191):
                                ForilearnAA.write('>'+str(name)+'_'+str(k)+'|'+'1'+'|'+'training'+'\n'+str(mutfn)+'\n')
                            else:
                                ForilearnAA.write('>'+str(name)+'_'+str(k)+'|'+'1'+'|'+'testing'+'\n'+str(mutfn)+'\n')  
                            train_number_AA = train_number_AA+1
                            
                            #--------------------for ilearn 写入文件-------------------##
                            ## training = 255*0.75 = 191  ---------
                            if(train_number < 191):
                                Forilearn.write('>'+str(name)+'_'+str(k)+'|'+'1'+'|'+'training'+'\n'+str(newseq)+'\n')
                            else:
                                Forilearn.write('>'+str(name)+'_'+str(k)+'|'+'1'+'|'+'testing'+'\n'+str(newseq)+'\n')  

                            #----------------   将proteinID 写入文件，用于后面的多个特征文件合并：proteinID作为主键--------#
                            proteinIDPOS.write(str(name)+'_'+str(k)+'C'+'\t'+str(num)+'\t'+'1'+'\n')

#                             ## -----------------##把proteinID,突变位置，突变前后AA，写入文件，用以进行突变预测分析---------#
#                             proteinIDPOSAA.write(str(name)+'\t'+str(num)+'\t'+qian+'\t'+hou+'\n')

#                             ##-----------------  ##为：fathmm使用的格式，把proteinID,突变位置，突变前后AA，写入文件 ------#
#                             proteinIDPOSAAfathmm.write(str(name)+'\t'+qian+str(num)+hou+'\n')
                            ## -----------------##把proteinID,突变位置，突变前后AA，写入文件，用以进行突变预测分析---------#
                            proteinIDPOSAA.write(str(name)+'_'+str(k)+'\t'+str(num)+'\t'+qian+'\t'+hou+'\n')

                            ##-----------------  ##为：fathmm使用的格式，把proteinID,突变位置，突变前后AA，写入文件 ------#
                            proteinIDPOSAAfathmm.write(str(name)+'_'+str(k)+'\t'+qian+str(num)+hou+'\n')

                            ##-------------##求PSA时，要求长度小于1000，这里取了突变点前后各5个AA，用于求PSA  --------------------##  
                            PSA.write('>'+str(name)+'_'+str(k))
                            PSA.write('\n'+str(mutfn)+'\n')

                            total_Seqnumbc=total_Seqnumbc+1
                            k=k+1 ##用于记录是第几个突变序列
                            flag = 1 #用于判断该蛋白质有没有写入文件：1表示该蛋白质写入了文件
                        else:###如果突变为*，则跳过，即不写入
                            continue               
                    else:##突变个数大于2，直接退出这个蛋白质，继续下个蛋白质
                        break
        if(flag == 1): ## 该蛋白质写入了文件
            total_proteinbc = total_proteinbc+1                
                                
BC_Mut.close()
bcancer.close()
BC_pssm.close()
Forilearn.close()
BC_MutFN.close()
ForilearnAA.close()
proteinIDPOS.close()
PSA.close()
proteinIDPOSAA.close()
proteinIDPOSAAfathmm.close()
print("总共的序列数： ",total_Seqnumbc)
print("总共的蛋白质数： ",total_proteinbc)


# ## Ensembl79_homo_dbSNP_variation_protein.txt
# 
# ## Ensembl79_homo_dbSNP_variation_protein.txt
# 
# ## 总共的序列数：  210379
# ## 总共的蛋白质数：  76337
#  
# ##    --------------  最终取了如下数据：
# ##  总共的序列数： 954
# ##  总共的蛋白质数： 495
# 
# 
# ## 每个蛋白最多留了2个突变，预防过拟合

# In[2]:


SNP=open('./data/original data/Ensembl79_homo_dbSNP_variation_protein.txt','r')
SNP_Mut=open('./data/temp/SNP_Mutall.txt','w')
SNP_pssm=open('./data/temp/SNP_pssmall.txt','w')##做pssm时需要这样的格式，仅有>,和一个字符串
BC_pssm=open('./data/temp/BC_pssm.txt','a')##用于存放1029个蛋白质序列：257（bc）+772=1029
Forilearn=open('./data/temp/Forilearn.txt','a') ##存放ilearn格式的蛋白质序列，并设置75%的training, 25%的testing
SNP_MutFN=open('./data/temp/SNP_MutFN_5.txt','w')  ##存放ilearn格式的蛋白质序列，并设置75%的training, 25%的testing
ForilearnAA=open('./data/temp/ForilearnAA.txt','a')  ##用于存放11个AA序列，并以ilearn格式存入
proteinIDPOS=open('./data/temp/proteinIDPOS.txt','a')##把proteinID写入文件，用于后面的多个特征文件合并：proteinID作为主键
PSA=open('./data/temp/PSA.txt','a')##求PSA时，要求长度小于1000，这里取了突变点前后各5个AA，用于求PSA
proteinIDPOSAA=open('./data/temp/proteinIDPOSAA.txt','a')##把proteinID,突变位置，突变前后AA，写入文件，用以进行突变预测分析
proteinIDPOSAAfathmm=open('./data/temp/proteinIDPOSAAfathmm.txt','a')##为：fathmm使用的格式，把proteinID,突变位置，突变前后AA，写入文件

##############################-----------------------处理序列突变数，把突变后为*的突变去掉----------########################
bc_lines=SNP.readlines()
total_proteinSNP = 0 #用于蛋白质条数
total_SeqnumSNP = 0 #用于记录总共的序列数
train_numberSNP=0 ##用于记录训练集的蛋白质条数
train_numberSNPAA =0 ##用于记录训练集的蛋白质条数:11个AA
count=0 ##只取前773个序列，所有count是用于记录是否到了773个
shuzi=['0','1','2','3','4','5','6','7','8','9'] ##用于判断突变中间数字是否都是数字？
# for i in range(len(bc_lines)):## 用于判断有多少个蛋白质和多少个序列。决定使用欧冠多少数据进行试验后，用下面的语句
for i in range( 500):  ##  取1800个序列
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
                      ## 判断即将写入的蛋白质不在BC中，用是否在proteinID_BC这个list来判断
                     ## 为了用iLearn计算APAAC, PACC,需要每个蛋白质的长度大于31
                    if (hou != '*' and qian !='*'and (int(num)+6)< len(seq) and (int(num)-5)> 0 and len(seq)>31 and name not in ProteinID_BC):
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

                            BC_pssm.write('>'+str(name)+'_'+str(k)+'S'+'\n')
                            BC_pssm.write(str(newseq)+'\n')
                                                
                            #----------------   将proteinID 写入文件，用于后面的多个特征文件合并：proteinID作为主键--------#
                            proteinIDPOS.write(str(name)+'_'+str(k)+'S'+'\t'+str(num)+'\t'+'0'+'\n')
                        
                           ## -----------------##把proteinID,突变位置，突变前后AA，写入文件，用以进行突变预测分析---------#
#                             proteinIDPOSAA.write(str(name)+'\t'+str(num)+'\t'+qian+'\t'+hou+'\n')
                            
#                             ##-----------------  ##为：fathmm使用的格式，把proteinID,突变位置，突变前后AA，写入文件 ------#
#                             proteinIDPOSAAfathmm.write(str(name)+'\t'+qian+str(num)+hou+'\n')
                           # -----------------##把proteinID,突变位置，突变前后AA，写入文件，用以进行突变预测分析---------#
                            proteinIDPOSAA.write(str(name)+'_'+str(k)+'\t'+str(num)+'\t'+qian+'\t'+hou+'\n')
                            
                            ##-----------------  ##为：fathmm使用的格式，把proteinID,突变位置，突变前后AA，写入文件 ------#
                            proteinIDPOSAAfathmm.write(str(name)+'_'+str(k)+'\t'+qian+str(num)+hou+'\n')                
                
                           ##-------------##求PSA时，要求长度小于1000，这里取了突变点前后各5个AA，用于求PSA  --------------------##  
                            PSA.write('>'+str(name)+'_'+str(k)+'S')
                            PSA.write('\n'+str(mutfn)+'\n')

                            #------------------------  将前后5个AA（共11个AA）按照ilearn的格式写入文件 -------#
                            #--------------------for ilearn 写入文件-------------------##
                            ##------------------------training = 954*0.75  =716--------------
                            if (train_numberSNP < 716):
                                ForilearnAA.write('>'+str(name)+'_'+str(k)+'S'+'|'+'0'+'|'+'training'+'\n'+str(mutfn)+'\n')
                            else:
                                ForilearnAA.write('>'+str(name)+'_'+str(k)+'S'+'|'+'0'+'|'+'testing'+'\n'+str(mutfn)+'\n')
                            train_numberSNPAA = train_numberSNPAA+1
                            ###---------------------  ilearn 格式的写入文件------------------
                            if (train_numberSNP < 716):
                                Forilearn.write('>'+str(name)+'_'+str(k)+'S'+'|'+'0'+'|'+'training'+'\n'+str(newseq)+'\n')
                            else:
                                Forilearn.write('>'+str(name)+'_'+str(k)+'S'+'|'+'0'+'|'+'testing'+'\n'+str(newseq)+'\n')
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
                      ## 判断即将写入的蛋白质不在BC中，用是否在proteinID_BC这个list来判断
                     ## 为了用iLearn计算APAAC, PACC,需要每个蛋白质的长度大于31
                        if (hou != '*' and qian !='*'and (int(num)+6)< len(seq) and (int(num)-5)> 0 and len(seq)>31 and name not in ProteinID_BC):

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
                            SNP_pssm.write('>'+str(name)+'_'+str(k)+'S'+'\t')
                            SNP_pssm.write(str(newseq)+'\n')
                            
                            BC_pssm.write('>'+str(name)+'_'+str(k)+'S'+'\n')
                            BC_pssm.write(str(newseq)+'\n')
                                                
                            #----------------   将proteinID 写入文件，用于后面的多个特征文件合并：proteinID作为主键--------#
                            proteinIDPOS.write(str(name)+'_'+str(k)+'S'+'\t'+str(num)+'\t'+'0'+'\n')
                        
#                            ## -----------------##把proteinID,突变位置，突变前后AA，写入文件，用以进行突变预测分析---------#
#                             proteinIDPOSAA.write(str(name)+'\t'+str(num)+'\t'+qian+'\t'+hou+'\n')
                            
#                             ##-----------------  ##为：fathmm使用的格式，把proteinID,突变位置，突变前后AA，写入文件 ------#
#                             proteinIDPOSAAfathmm.write(str(name)+'\t'+qian+str(num)+hou+'\n')
                           # -----------------##把proteinID,突变位置，突变前后AA，写入文件，用以进行突变预测分析---------#
                            proteinIDPOSAA.write(str(name)+'_'+str(k)+'\t'+str(num)+'\t'+qian+'\t'+hou+'\n')
                            
                            ##-----------------  ##为：fathmm使用的格式，把proteinID,突变位置，突变前后AA，写入文件 ------#
                            proteinIDPOSAAfathmm.write(str(name)+'_'+str(k)+'\t'+qian+str(num)+hou+'\n')                
                
                
                           ##-------------##求PSA时，要求长度小于1000，这里取了突变点前后各5个AA，用于求PSA  --------------------##  
                            PSA.write('>'+str(name)+'_'+str(k)+'S')
                            PSA.write('\n'+str(mutfn)+'\n')

                            #------------------------  将前后5个AA（共11个AA）按照ilearn的格式写入文件 -------#
                            #--------------------for ilearn 写入文件-------------------##
                            ##------------------------training = 954*0.75  =716--------------
                            if (train_numberSNP < 716):
                                ForilearnAA.write('>'+str(name)+'_'+str(k)+'S'+'|'+'0'+'|'+'training'+'\n'+str(mutfn)+'\n')
                            else:
                                ForilearnAA.write('>'+str(name)+'_'+str(k)+'S'+'|'+'0'+'|'+'testing'+'\n'+str(mutfn)+'\n')
                            train_numberSNPAA = train_numberSNPAA+1
                            ###---------------------  ilearn 格式的写入文件------------------
                            if (train_numberSNP < 716):
                                Forilearn.write('>'+str(name)+'_'+str(k)+'S'+'|'+'0'+'|'+'training'+'\n'+str(newseq)+'\n')
                            else:
                                Forilearn.write('>'+str(name)+'_'+str(k)+'S'+'|'+'0'+'|'+'testing'+'\n'+str(newseq)+'\n')
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
BC_pssm.close()
PSA.close()
proteinIDPOSAA.close()
proteinIDPOSAAfathmm.close()
print("总共的序列数： ",total_SeqnumSNP)
print("总共的蛋白质数： ",total_proteinSNP)

