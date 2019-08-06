# -*- coding: utf-8 -*-
"""
Created on Mon Sep 17 15:09:38 2018

@author: HAN_RUIZHI yb77447@umac.mo OR  501248792@qq.com

This code is the first version of BLS Python. 
If you have any questions about the code or find any bugs
   or errors during use, please feel free to contact me.
If you have any questions about the original paper, 
   please contact the authors of related paper.
"""

import numpy as np
from sklearn import preprocessing
from numpy import random
from scipy import linalg as LA 


'''
激活函数
'''
def tansig(x):
    return (2/(1+np.exp(-2*x)))-1

def sigmoid(data):
    return 1.0/(1+np.exp(-data))
    
def linear(data):
    return data
    
def tanh(data):
    return (np.exp(data)-np.exp(-data))/(np.exp(data)+np.exp(-data))
    
def relu(data):
    return np.maximum(data,0)

def pinv(A,reg):
    return np.mat(reg*np.eye(A.shape[1])+A.T.dot(A)).I.dot(A.T)
'''
参数压缩
'''
def shrinkage(a,b):
    z = np.maximum(a - b, 0) - np.maximum( -a - b, 0)
    return z
'''
参数稀疏化
'''
def sparse_bls(A,b):
    lam = 0.001
    itrs = 50
    AA = A.T.dot(A)   
    m = A.shape[1]
    n = b.shape[1]
    x1 = np.zeros([m,n])
    wk = x1
    ok = x1
    uk = x1
    L1 = np.mat(AA + np.eye(m)).I
    L2 = (L1.dot(A.T)).dot(b)
    for i in range(itrs):
        ck = L2 + np.dot(L1,(ok - uk))
        ok = shrinkage(ck + uk, lam)
        uk = uk + ck - ok
        wk = ok   
    return wk


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def BLS(train_x,train_y,test_x,test_y,s,c,N1,N2,N3):

    train_x = preprocessing.scale(train_x,axis = 1)# ,with_mean = '0') #处理数据 
    FeatureOfInputDataWithBias = np.hstack([train_x, 0.1 * np.ones((train_x.shape[0],1))])
    OutputOfFeatureMappingLayer = np.zeros([train_x.shape[0],N2*N1])
    Beta1OfEachWindow = []

    distOfMaxAndMin = []
    minOfEachWindow = []
    ymin = 0
    ymax = 1

    for i in range(N2):
        random.seed(i)
        weightOfEachWindow = 2 * random.randn(train_x.shape[1]+1,N1)-1; #生成每个窗口的权重系数，最后一行为偏差
#        WeightOfEachWindow([],[],i) = weightOfEachWindow; #存储每个窗口的权重系数
        FeatureOfEachWindow = np.dot(FeatureOfInputDataWithBias,weightOfEachWindow) #生成每个窗口的特征
        #压缩每个窗口特征到[-1，1]
        scaler1 = preprocessing.MinMaxScaler(feature_range=(0, 1)).fit(FeatureOfEachWindow)
        FeatureOfEachWindowAfterPreprocess = scaler1.transform(FeatureOfEachWindow)
        #通过稀疏化计算映射层每个窗口内的最终权重
        betaOfEachWindow  =  sparse_bls(FeatureOfEachWindowAfterPreprocess,FeatureOfInputDataWithBias).T
        #存储每个窗口的系数化权重
        Beta1OfEachWindow.append(betaOfEachWindow)
        #每个窗口的输出 T1
        outputOfEachWindow = np.dot(FeatureOfInputDataWithBias,betaOfEachWindow)
#        print('Feature nodes in window: max:',np.max(outputOfEachWindow),'min:',np.min(outputOfEachWindow))
        distOfMaxAndMin.append(np.max(outputOfEachWindow,axis =0) - np.min(outputOfEachWindow,axis=0))
        minOfEachWindow.append(np.min(outputOfEachWindow,axis = 0))
        outputOfEachWindow = (outputOfEachWindow-minOfEachWindow[i])/distOfMaxAndMin[i]
        OutputOfFeatureMappingLayer[:,N1*i:N1*(i+1)] = outputOfEachWindow
        del outputOfEachWindow 
        del FeatureOfEachWindow 
        del weightOfEachWindow 
        
    #生成强化层
    #以下为映射层输出加偏置（强化层输入）
    InputOfEnhanceLayerWithBias = np.hstack([OutputOfFeatureMappingLayer, 0.1 * np.ones((OutputOfFeatureMappingLayer.shape[0],1))])
    #生成强化层权重
    if N1*N2>=N3:
        random.seed(67797325)
#        dim = N1*N2+1
#        temp_matric = stats.ortho_group(dim)
#        weightOfEnhanceLayer = temp_matric[:,0:N3]
        weightOfEnhanceLayer = LA.orth(2 * random.randn(N2*N1+1,N3))-1
    else:
        random.seed(67797325)
        weightOfEnhanceLayer = LA.orth(2 * random.randn(N2*N1+1,N3).T-1).T
    
    tempOfOutputOfEnhanceLayer = np.dot(InputOfEnhanceLayerWithBias,weightOfEnhanceLayer)
#    print('Enhance nodes: max:',np.max(tempOfOutputOfEnhanceLayer),'min:',np.min(tempOfOutputOfEnhanceLayer))

    parameterOfShrink = s/np.max(tempOfOutputOfEnhanceLayer)

    OutputOfEnhanceLayer = sigmoid(tempOfOutputOfEnhanceLayer * parameterOfShrink)
    
    #生成最终输入
    InputOfOutputLayer = np.hstack([OutputOfFeatureMappingLayer,OutputOfEnhanceLayer])
    pinvOfInput = pinv(InputOfOutputLayer,c)
    OutputWeight = np.dot(pinvOfInput,train_y) #全局违逆
    OutputWeight = OutputWeight.T


    
    #训练输出
#    OutputOfTrain = np.dot(InputOfOutputLayer,OutputWeight)
#    print('Training time is ',trainTime,'s')
    #测试过程
    test_x = preprocessing.scale(test_x,axis = 1)#,with_mean = True,with_std = True) #处理数据 x = (x-mean(x))/std(x) x属于[-1，1]
    FeatureOfInputDataWithBiasTest = np.hstack([test_x, 0.1 * np.ones((test_x.shape[0],1))])
    OutputOfFeatureMappingLayerTest = np.zeros([test_x.shape[0],N2*N1])
    
#  映射层
    for i in range(N2):
        outputOfEachWindowTest = np.dot(FeatureOfInputDataWithBiasTest,Beta1OfEachWindow[i])
        OutputOfFeatureMappingLayerTest[:,N1*i:N1*(i+1)] =(ymax-ymin)*(outputOfEachWindowTest-minOfEachWindow[i])/distOfMaxAndMin[i]-ymin
#  强化层
    InputOfEnhanceLayerWithBiasTest = np.hstack([OutputOfFeatureMappingLayerTest, 0.1 * np.ones((OutputOfFeatureMappingLayerTest.shape[0],1))])
    tempOfOutputOfEnhanceLayerTest = np.dot(InputOfEnhanceLayerWithBiasTest,weightOfEnhanceLayer)
#  强化层输出
    OutputOfEnhanceLayerTest = sigmoid(tempOfOutputOfEnhanceLayerTest * parameterOfShrink)    
#  最终层输入
    InputOfOutputLayerTest = np.hstack([OutputOfFeatureMappingLayerTest,OutputOfEnhanceLayerTest])
#  最终测试输出   
    OutputOfTest = np.dot(InputOfOutputLayerTest,OutputWeight)
    OutProbOfTest = np.zeros(len(OutputOfTest))
    for i in range(len(OutputOfTest)):
        OutProbOfTest[i] = float(OutputOfTest[i])
    return OutProbOfTest

