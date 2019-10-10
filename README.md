# TargetMM

TargetMM:
Missense mutation is one kind of Single Nucleotide Variations, which might likely cause diseade or even cancer.TargetMM utilize evaluationary information, such as PSSM and predicted Secondary Structure information, and Disorder region of proteins. In addition, we also consider protein global features, such as APAAC and PAAC. On the training data, we trained and validated SVM, RF,and GPC models. Then,intergate them together to form an ensemble classifier model. We compare TargetMM with four existing prediction methods,such as, SIFT, PROVEAN, FATHMM, and PolyPhen-2. The comparision results showed that TargetMM has much higher MCC than that of all these four methods.


## CanProVar 2.0 protein missense mutation database
Address: http://canprovar2.zhang-lab.org/datadownload.php
CanProVar provides the download of human protein database (Ensembl v79) in the fasta format, in which variation information is recorded in the header line of each sequence. 
MS-CanProVar (version 2.0) is a protein sequence database that includes variation information to facilitate peptide variant detection in shotgun proteomics. In the .fasta file, each variant peptide is included as an independent entry; variations are annotated in the header line; variations are labeled as "rs" for SNPs and "cs" for cancer-related mutations. Please refer to A bioinformatics workflow for variant peptide detection in shotgun proteomics. Li et al., MCP, 2011 for details about the MS-CanProVar database. The current version of MS-CanProVar is based on Ensembl V79. 


## Extracting features
### 1. Position Specific Score Matrix (PSSM)
You could use PSI-BLAST software searing swiss-port database. Software address :https://www.ebi.ac.uk/seqdb/confluence/display/THD/PSI-BLAST.

### 2.Predicted Secondary Structure (PSS)
PSIPRED is one kind of tools to predict secondary structure, for more datails at https://hpc.nih.gov/apps/PSIPRED.html

### 3.Disorder
Dynamically disordered regions appear to be relatively abundant in eukaryotic proteomes. The DISOPRED server allows users to submit a protein sequence, and returns a probability estimate of each residue in the sequence being disordered. The results are sent in both plain text and graphical formats, and the server can also supply predictions of secondary structure to provide further structural information.

See more information from this paper below:
Jonathan J. Ward, Liam J. McGuffin, Kevin Bryson, Bernard F. Buxton, David T. Jones, The DISOPRED server for the prediction of protein disorder, Bioinformatics, Volume 20, Issue 13, 1 September 2004, Pages 2138–2139, https://doi.org/10.1093/bioinformatics/bth195

### 4. Amphiphilic Pseudo-Amino Acid Composition and Pseudo-Amino Acid Composition (APAAC and PAAC)
We gather APAAC and PAAC from iLearn (http://ilearn.erc.monash.edu.au/), which is an integrated platform and meta-learner for feature engineering, machine-learning analysis and modeling of DNA, RNA and protein sequence data.

If you use feature extracting tools,please cite the following paper:
Chen Z, Zhao P, Li F, et al. iLearn: an integrated platform and meta-learner for feature engineering, machine-learning analysis and modeling of DNA, RNA and protein sequence data[J]. Briefings in bioinformatics, 2019.


## Machine Learning Technology
In this work, we utilize three classifiers: Random Forest, Support Vector Machine, Gauss Process Classification from scikit-learn, which is a "Simple and efficient tools for data mining and data analysis", "Accessible to everybody, and reusable in various contexts",  "Built on NumPy, SciPy, and matplotlib", and "Open source, commercially usable - BSD license". See more datails from https://scikit-learn.org/stable/index.html

## Contact 
If you are interested in our work,OR, if you have any suggestiongs and questions about our source code, PLEASE contact with us. 
E-mail address : gfang0616@njust.edu.cn

## Cite
paper: Fang Ge, Yi-heng Zhu, Wen-wen Kan, Xu Chen, Dong-jun Yu. Accurate Prediction of Missense Mutation by Combing Local and Global Sequence Information with Classifier Ensemble. Under review.



