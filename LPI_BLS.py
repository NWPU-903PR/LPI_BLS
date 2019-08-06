# -*- coding: utf-8 -*-
"""

"""

import numpy as np
#import GetFeatures
from sklearn import preprocessing
from sklearn.model_selection import StratifiedKFold
from sklearn.linear_model import LogisticRegression
import BroadLearningSystem
from sklearn import metrics
import matplotlib.pyplot as plt 
from sklearn.metrics import auc
import argparse
import os


# ######function readInteractions()
def readInteractions(datafile):
    RPI_Names = {'Protein':[],
                 'lncRNA':[]}
    with open(datafile, 'r') as file_to_read:
        while True:
            lines = file_to_read.readline()
            if not lines:
                break
            lines = lines.strip()
            lines = lines.split('\t')
            RPI_Names['Protein'].append(lines[1])
            RPI_Names['lncRNA'].append(lines[0])
    return RPI_Names

######function GetKmer()
def get_4_trid_RNA():
    nucle_com = []
    chars = ['A', 'C', 'G', 'U']
    base = len(chars)
    end = len(chars)**4
    for i in range(end):
        n = i
        ch0 = chars[int(n%base)]
        n = n/base
        ch1=chars[int(n%base)]
        n = n/base
        ch2 = chars[int(n%base)]
        n = n/base
        ch3 = chars[int(n%base)]
        nucle_com.append(ch0 + ch1 + ch2 + ch3)
    return  nucle_com

def TransDict_from_list():
    tar_list = ['0', '1', '2', '3', '4', '5', '6']
    result = {}
    index = 0
    groups = ['AGV', 'ILFP', 'YMTS', 'HNQW', 'RK', 'DE', 'C']
    for group in groups:
        g_members = sorted(group)
        for c in g_members:
            result[c] = str(tar_list[index]) 
        index = index + 1
    return result

def translate_sequence(seq, group_dict):
    from_list = []
    to_list = []
    for k,v in group_dict.items():
        from_list.append(k)
        to_list.append(v)
    TRANS_seq = seq.translate(str.maketrans(str(from_list), str(to_list)))
    return TRANS_seq

def get_3_protein_trids():
    nucle_com = []
    chars = ['0', '1', '2', '3', '4', '5', '6']
    base = len(chars)
    end = len(chars)**3
    for i in range(0,end):
        n = i
        ch0 = chars[int(n%base)]
        n = n/base
        ch1 = chars[int(n%base)]
        n = n/base
        ch2 = chars[int(n%base)]
        nucle_com.append(ch0 + ch1 + ch2)
    return  nucle_com


def get_protein_trids(seq, group_dict):
    tran_seq = translate_sequence (seq, group_dict)
    return tran_seq

def GetKmer_featrue(RNA_file, Protein_file, RNA_4_trid, group_dict, protein_3_tris):
    Kmer_Features = {}
    RNA = list()
    with open(RNA_file, 'r') as Rfile_to_read:
        while True:
            lines = Rfile_to_read.readline()
            if not lines:
                break
            if lines.startswith('>'):
                name = lines.replace('>','').split()[0]
                Kmer_Features[name] = ''
                RNA.append(name)
            else:
                Rseq = lines.replace('\n','').strip()
                seq_len = len(Rseq)
                Rkmer = []
                for val in RNA_4_trid:
                    num = Rseq.count(val)
                    Rkmer.append(float(num)/seq_len)
                Rkmer = np.array(Rkmer)
                Kmer_Features[name] = Rkmer.reshape(1,len(Rkmer))

    protein = list()
    with open(Protein_file, 'r') as Pfile_to_read:
        while True:
            lines = Pfile_to_read.readline()
            if not lines:
                break
            if lines.startswith('>'):
                name = lines.replace('>','').split()[0]
                Kmer_Features[name] = ''
                protein.append(name)
            else:
                Pseq = lines.replace('\n','').strip()
                Pseq = translate_sequence(Pseq, group_dict)
                seq_len = len(Pseq)
                Pkmer = []
                for val in protein_3_tris:
                    num = Pseq.count(val)
                    Pkmer.append(float(num)/seq_len)
                Pkmer = np.array(Pkmer)
                Kmer_Features[name] = Pkmer.reshape(1,len(Pkmer))
    return Kmer_Features, RNA, protein

######function Get_pse_featrue()
def Get_pse_featrue(PseProtein_file, PseRNA_file, RNA, protein):
    Pse_Features = {}
    with open(PseProtein_file, 'r') as Pfile_to_read:
        k = 0
        while True:
            lines = Pfile_to_read.readline()
            if not lines:
                break
            ff = lines.replace('\n','').strip()
            ff = [float(i) for i in ff.split()]
            ff = np.array(ff)
            Pse_Features[protein[k]] = ff.reshape(1,len(ff))
            k += 1
    
    with open(PseRNA_file, 'r') as Rfile_to_read:
        i = 0
        while True:
            lines = Rfile_to_read.readline()
            if not lines:
                break
            ff = lines.replace('\n','').strip()
            ff = [float(i) for i in ff.split()]
            ff = np.array(ff)
            Pse_Features[RNA[i]] = ff.reshape(1,len(ff))
            i += 1
    return Pse_Features


def Get_kmer_data(RPI_Names, RNA_file, Protein_file):
    RNA_4_trid = get_4_trid_RNA()
    group_dict = TransDict_from_list()
    protein_3_tris = get_3_protein_trids()
    Kmer_Features, RNA, protein = GetKmer_featrue(RNA_file, Protein_file, RNA_4_trid, group_dict, protein_3_tris)
    
    Kmer_Matrix = np.empty(shape=[0, 599])
    num = len(RPI_Names['lncRNA'])
    for i in range(num):
        proteinName = RPI_Names['Protein'][i]
        if proteinName in Kmer_Features:
            protein_features = Kmer_Features[proteinName]
        else:
            print('please input all protein sequence provided in lncRNA-protein pairs file, \n No', proteinName, 'in fasta file of protein sequence')
            break
        
        RNAName = RPI_Names['lncRNA'][i]
        if RNAName in Kmer_Features:
            RNA_features = Kmer_Features[RNAName]
        else:
            print('please input all lncRNA sequence provided in lncRNA-protein pairs file, \n No', RNAName, 'in fasta file of lncRNA sequence')
            break
        ff = np.concatenate((RNA_features, protein_features),axis=1)
        Kmer_Matrix = np.concatenate((Kmer_Matrix, ff),axis=0)
    return Kmer_Matrix, RNA, protein

def Get_pse_data(RPI_Names, PseProtein_file, PseRNA_file, RNA, protein):
    Pse_Features = Get_pse_featrue(PseProtein_file, PseRNA_file, RNA, protein)
    
    Pse_Matrix = np.empty(shape=[0, 51])
    num = len(RPI_Names['lncRNA'])
    for i in range(num):
        proteinName = RPI_Names['Protein'][i]
        if proteinName in Pse_Features:
            protein_features = Pse_Features[proteinName]
        else:
            print('cannot find correct pse feature file of protein')
            break
        
        RNAName = RPI_Names['lncRNA'][i]
        if RNAName in Pse_Features:
            RNA_features = Pse_Features[RNAName]
        else:
            print('cannot find correct pse feature file of lncRNA')
            break
        ff = np.concatenate((RNA_features, protein_features),axis=1)
        Pse_Matrix = np.concatenate((Pse_Matrix, ff),axis=0)
   
    return Pse_Matrix


def generate_features_for_new_pairs(pairs_file, lncRNA_fa, protein_fa, lncRNA_pse, protein_pse):
    RPI_Names = readInteractions(pairs_file)

    print('---------get kmer feature---------------')
    Data_kmer, RNA, protein = Get_kmer_data(RPI_Names, lncRNA_fa, protein_fa)
    print('---------get ktuple feature---------------')
    Data_ktuple = Get_pse_data(RPI_Names, protein_pse, lncRNA_pse, RNA, protein)
    return RPI_Names, Data_kmer, Data_ktuple

def generate_RPI488_features():
    Data_kmer = np.load('./data/Data_kmer_RPI488.npy')
    Data_ktuple = np.load('./data/Data_ktuple_RPI488.npy')
    labels = np.load('./data/labels_RPI488.npy')
    N1 = 3
    N2 = 100
    N3 = 100
    skf = StratifiedKFold(n_splits = 5, random_state = 512374, shuffle = True)
    y = labels[:,0]
    return skf, Data_kmer, Data_ktuple, y, N1, N2, N3


def generate_RPI7317_features():
    Data_kmer = np.load('./data/Data_kmer_RPI7317.npy')
    Data_ktuple = np.load('./data/Data_ktuple_RPI7317.npy')
    labels = np.load('./data/labels_RPI7317.npy')
    N1 = 3
    N2 = 60
    N3 = 900
    skf = StratifiedKFold(n_splits = 5, random_state = 1, shuffle = True)
    y = labels[:,0]
    return skf, Data_kmer, Data_ktuple, y, N1, N2, N3
 

def show_performance(prob, true_y):
    TP = 0
    FP = 0
    FN = 0
    TN = 0
    for i in range(len(true_y)):
        if true_y[i] == 1:
            if prob[i] >= 0.5:
                TP += 1
            else:
                FN += 1
        if true_y[i] == 0:
            if prob[i] >= 0.5:
                FP += 1
            else:
                TN += 1
    ACC = (TP + TN)/len(true_y)
    Sn = TP/(TP + FN)
    Sp = TN/(FP + TN)
    MCC = ((TP * TN) - (FP * FN))/np.sqrt((TP+FP) * (TP+FN) * (TN+FP) * (TN+FN))
    return ACC, Sn, Sp, MCC


def run_stacking(x_train_cv1, x_train_cv2, x_train_cv3, x_train_cv4, x_train_cv5, y_train_cv, 
                 x_test_cv1, x_test_cv2, x_test_cv3, x_test_cv4, x_test_cv5, y_test_cv, N1, N2, N3):
    
    
    skf0 = StratifiedKFold(n_splits = 10, random_state = 1, shuffle = True)

    blend_train = np.zeros((x_train_cv1.shape[0], 5))
    blend_test = np.zeros((x_test_cv1.shape[0], 5))
    blend_test_1 = np.zeros((x_test_cv1.shape[0], 10))
    blend_test_2 = np.zeros((x_test_cv1.shape[0], 10))
    blend_test_3 = np.zeros((x_test_cv1.shape[0], 10))
    blend_test_4 = np.zeros((x_test_cv1.shape[0], 10))
    blend_test_5 = np.zeros((x_test_cv1.shape[0], 10))
    num_cv = 0
    
    for train_index_cv, test_index_cv in skf0.split(x_train_cv1, y_train_cv): 
        xx_train1, xx_test1 = x_train_cv1[train_index_cv,:], x_train_cv1[test_index_cv,:]
        xx_train2, xx_test2 = x_train_cv2[train_index_cv,:], x_train_cv2[test_index_cv,:]
        xx_train3, xx_test3 = x_train_cv3[train_index_cv,:], x_train_cv3[test_index_cv,:]
        xx_train4, xx_test4 = x_train_cv4[train_index_cv,:], x_train_cv4[test_index_cv,:]
        xx_train5, xx_test5 = x_train_cv5[train_index_cv,:], x_train_cv5[test_index_cv,:]
        yy_train, yy_test = y_train_cv[train_index_cv], y_train_cv[test_index_cv]
        train_prob_cv_1 = BroadLearningSystem.BLS(xx_train1, yy_train, xx_test1, yy_test, 1, 10**-10, N1, N2, N3)
        train_prob_cv_2 = BroadLearningSystem.BLS(xx_train2, yy_train, xx_test2, yy_test, 1, 10**-10, N1, N2, N3)
        train_prob_cv_3 = BroadLearningSystem.BLS(xx_train3, yy_train, xx_test3, yy_test, 1, 10**-10, N1, N2, N3)
        train_prob_cv_4 = BroadLearningSystem.BLS(xx_train4, yy_train, xx_test4, yy_test, 1, 10**-10, N1, N2, N3)
        train_prob_cv_5 = BroadLearningSystem.BLS(xx_train5, yy_train, xx_test5, yy_test, 1, 10**-10, N1, N2, N3)
        blend_train[test_index_cv, 0] = train_prob_cv_1
        blend_train[test_index_cv, 1] = train_prob_cv_2
        blend_train[test_index_cv, 2] = train_prob_cv_3
        blend_train[test_index_cv, 3] = train_prob_cv_4
        blend_train[test_index_cv, 4] = train_prob_cv_5
        test_prob_cv_1 = BroadLearningSystem.BLS(xx_train1, yy_train, x_test_cv1, y_test_cv, 1, 10**-10, N1, N2, N3)
        test_prob_cv_2 = BroadLearningSystem.BLS(xx_train2, yy_train, x_test_cv2, y_test_cv, 1, 10**-10, N1, N2, N3)
        test_prob_cv_3 = BroadLearningSystem.BLS(xx_train3, yy_train, x_test_cv3, y_test_cv, 1, 10**-10, N1, N2, N3)
        test_prob_cv_4 = BroadLearningSystem.BLS(xx_train4, yy_train, x_test_cv4, y_test_cv, 1, 10**-10, N1, N2, N3)
        test_prob_cv_5 = BroadLearningSystem.BLS(xx_train5, yy_train, x_test_cv5, y_test_cv, 1, 10**-10, N1, N2, N3)
        blend_test_1[:,num_cv] = test_prob_cv_1
        blend_test_2[:,num_cv] = test_prob_cv_2
        blend_test_3[:,num_cv] = test_prob_cv_3
        blend_test_4[:,num_cv] = test_prob_cv_4
        blend_test_5[:,num_cv] = test_prob_cv_5
        num_cv += 1
    blend_test[:, 0] = blend_test_1.mean(axis = 1)
    blend_test[:, 1] = blend_test_2.mean(axis = 1)
    blend_test[:, 2] = blend_test_3.mean(axis = 1)
    blend_test[:, 3] = blend_test_4.mean(axis = 1)
    blend_test[:, 4] = blend_test_5.mean(axis = 1)
    
#    blend_train = preprocessing.scale(blend_train)
#    blend_test = preprocessing.scale(blend_test)

    bclf = LogisticRegression(solver='liblinear', penalty = 'l1')
    bclf.fit(blend_train, y_train_cv)
    test_predict = bclf.predict_proba(blend_test)
    
    return test_predict[:,1]


def LPI_BL_new(pairs_file, lncRNA_fa, protein_fa, lncRNA_pse, protein_pse):
    skf, Data_kmer, Data_ktuple, y, N1, N2, N3 = generate_RPI488_features()
    x_train_1 = Data_kmer
    x_train_2 = Data_ktuple
    x_train_3 = np.hstack([Data_kmer[:, 0:256], Data_ktuple[:, 22:51]])
    x_train_4 = np.hstack([Data_ktuple[:, 0:22], Data_kmer[:, 256:599]])
    x_train_5 = np.hstack([Data_kmer, Data_ktuple])
    y_train = y
        
    RPI_Names, test_kmer, test_ktuple = generate_features_for_new_pairs(pairs_file, lncRNA_fa, protein_fa, lncRNA_pse, protein_pse)
    x_test_1 = test_kmer
    x_test_2 = test_ktuple
    x_test_3 = np.hstack([test_kmer[:, 0:256], test_ktuple[:, 22:51]])
    x_test_4 = np.hstack([test_ktuple[:, 0:22], test_kmer[:, 256:599]])
    x_test_5 = np.hstack([test_kmer, test_ktuple])
    y_test = 0
    print('---predicting---')
    test_prob_6 = run_stacking(x_train_1, x_train_2, x_train_3, x_train_4, x_train_5, y_train, 
                               x_test_1, x_test_2, x_test_3, x_test_4, x_test_5, y_test, N1, N2, N3)
    np.savetxt('./results/predicted_probs.txt', test_prob_6)


def LPI_BL_cv(dataset):
    if dataset == "RPI488":
        skf, Data_kmer, Data_ktuple, y, N1, N2, N3 = generate_RPI488_features()
    if dataset == "RPI7317":
        skf, Data_kmer, Data_ktuple, y, N1, N2, N3 = generate_RPI7317_features()
        
    
    profermance = np.zeros((5,6,4))
    n_cv = 0
    for train_index, test_index in skf.split(Data_kmer, y):
        x_train_1, x_test_1 = Data_kmer[train_index,:], Data_kmer[test_index,:]
        x_train_2, x_test_2 = Data_ktuple[train_index,:], Data_ktuple[test_index,:]
        
        x_train_3 = np.hstack([Data_kmer[train_index, 0:256], Data_ktuple[train_index, 22:51]])
        x_test_3 = np.hstack([Data_kmer[test_index, 0:256], Data_ktuple[test_index, 22:51]])
        
        x_train_4 = np.hstack([Data_ktuple[train_index, 0:22], Data_kmer[train_index, 256:599]])
        x_test_4 = np.hstack([Data_ktuple[test_index, 0:22], Data_kmer[test_index, 256:599]])
        
        x_train_5 = np.hstack([Data_kmer[train_index,:], Data_ktuple[train_index,:]])
        x_test_5 = np.hstack([Data_kmer[test_index,:], Data_ktuple[test_index,:]])
        
        y_train, y_test = y[train_index], y[test_index]
        
        print('---1 run BLS1---')
        test_prob_1 = BroadLearningSystem.BLS(x_train_1, y_train, x_test_1, y_test, 1, 10**-10, N1, N2, N3)
        acc1, sn1, sp1, mcc1 = show_performance(test_prob_1, y_test)
        profermance[n_cv,0,:] = np.array((acc1, sn1, sp1, mcc1))
        print(acc1, sn1, sp1, mcc1 )
        
        print('---2 run BLS2---')
        test_prob_2 = BroadLearningSystem.BLS(x_train_2, y_train, x_test_2, y_test, 1, 10**-10, N1, N2, N3)
        acc2, sn2, sp2, mcc2 = show_performance(test_prob_2, y_test)
        profermance[n_cv,1,:] = np.array((acc2, sn2, sp2, mcc2))
        print(acc2, sn2, sp2, mcc2)
        
        print('---3 run BLS3---')
        test_prob_3 = BroadLearningSystem.BLS(x_train_3, y_train, x_test_3, y_test, 1, 10**-10, N1, N2, N3)
        acc3, sn3, sp3, mcc3 = show_performance(test_prob_3, y_test)
        profermance[n_cv,2,:] = np.array((acc3, sn3, sp3, mcc3))
        print(acc3, sn3, sp3, mcc3)
        
        print('---4 run BLS4---')
        test_prob_4 = BroadLearningSystem.BLS(x_train_4, y_train, x_test_4, y_test, 1, 10**-10, N1, N2, N3)
        acc4, sn4, sp4, mcc4 = show_performance(test_prob_4, y_test)
        profermance[n_cv,3,:] = np.array((acc4, sn4, sp4, mcc4))
        print(acc4, sn4, sp4, mcc4)
        
        print('---5 run BLS5---')
        test_prob_5 = BroadLearningSystem.BLS(x_train_5, y_train, x_test_5, y_test, 1, 10**-10, N1, N2, N3)
        acc5, sn5, sp5, mcc5 = show_performance(test_prob_5, y_test)
        profermance[n_cv,4,:] = np.array((acc5, sn5, sp5, mcc5))  
        print(acc5, sn5, sp5, mcc5)
        
        print('---6 run stacked ensemble ---')
        test_prob_6 = run_stacking(x_train_1, x_train_2, x_train_3, x_train_4, x_train_5, y_train, 
                                   x_test_1, x_test_2, x_test_3, x_test_4, x_test_5, y_test, N1, N2, N3)
        acc6, sn6, sp6, mcc6 = show_performance(test_prob_6, y_test)
        profermance[n_cv,5,:] = np.array((acc6, sn6, sp6, mcc6))
        print(acc6, sn6, sp6, mcc6)
        
        n_cv += 1
        
    print('BLS1---', 'acc---', np.mean(profermance[:,0,0]), 'sn---', np.mean(profermance[:,0,1]), 
          'sp---', np.mean(profermance[:,0,2]), 'mcc---', np.mean(profermance[:,0,3]))
    print('BLS2---', 'acc---', np.mean(profermance[:,1,0]), 'sn---', np.mean(profermance[:,1,1]), 
          'sp---', np.mean(profermance[:,1,2]), 'mcc---', np.mean(profermance[:,1,3]))
    print('BLS3---', 'acc---', np.mean(profermance[:,2,0]), 'sn---', np.mean(profermance[:,2,1]), 
          'sp---', np.mean(profermance[:,2,2]), 'mcc---', np.mean(profermance[:,2,3]))        
    print('BLS4---', 'acc---', np.mean(profermance[:,3,0]), 'sn---', np.mean(profermance[:,3,1]), 
          'sp---', np.mean(profermance[:,3,2]), 'mcc---', np.mean(profermance[:,3,3]))
    print('BLS5---', 'acc---', np.mean(profermance[:,4,0]), 'sn---', np.mean(profermance[:,4,1]), 
          'sp---', np.mean(profermance[:,4,2]), 'mcc---', np.mean(profermance[:,4,3]))        
    print('stacked ensemble ---', 'acc---', np.mean(profermance[:,5,0]), 'sn---', np.mean(profermance[:,5,1]), 
          'sp---', np.mean(profermance[:,5,2]), 'mcc---', np.mean(profermance[:,5,3]))   

    

parser = argparse.ArgumentParser(description="LPI-BLS: Predicting lncRNA-protein interactions with a broad learning system-based stacked ensemble classifier")
parser.add_argument('-dataset', type=str, help='dataset for 5-fold cross validation, (RPI488 or RPI7317)')
parser.add_argument('-pair', type=str, help='txt file for lncRNA-protein pairs which you want to be predicted, \n see ./data/example_pairs.txt')
parser.add_argument('-rf', type=str, help='fasta file of lncRNA sequence, \n see ./data/example_lncRNA.fa')
parser.add_argument('-pf', type=str, help='fasta file of protein sequence, \n see ./data/example_protein.fa')
parser.add_argument('-rP', type=str, help='file of lncRNA pse feature, \n see ./data/example_lncRNA_pse_feature.fa')
parser.add_argument('-pP', type=str, help='file of protein pse feature, \n see ./data/example_protein_pse_feature.fa')


args = parser.parse_args()
dataset = args.dataset
if dataset is not None:
    print('LPI-BLS performs on', dataset, 'in 5-fold cross validation')
    LPI_BL_cv(dataset)
else:
    pairs_file = args.pair
    lncRNA_fa = args.rf
    protein_fa = args.pf
    lncRNA_pse = args.rP
    protein_pse = args.pP
    if lncRNA_fa is None or protein_fa is None:
        print('please input both lncRNA file and protein file')
    try:
        f = open(protein_pse)
        f.close()
    except IOError:
        print("protein pse feature file is not accessible.")
        print('please make sure that you have run the pse-in-one first')
    try:
        f = open(lncRNA_pse)
        f.close()
    except IOError:
        print("lncRNA pse feature file is not accessible.")
        print('please make sure that you have run the pse-in-one first') 
    LPI_BL_new(pairs_file, lncRNA_fa, protein_fa, lncRNA_pse, protein_pse)
