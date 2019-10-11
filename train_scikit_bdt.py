#!/usr/bin/python
# this is the example script to use xgboost to train
import numpy as np
import pandas as pd

from ROOT import TFile, TTree
from array import array


from sklearn.ensemble import GradientBoostingClassifier

def getGradBDT(options={}):
    """ Standard BDT classifier based on GradienBoosting"""
    bdt = GradientBoostingClassifier(n_estimators=20,
                                     learning_rate=0.08,
                                     max_depth=6,
                                     min_weight_fraction_leaf=0.08,
                                     random_state=0,
                                     verbose=4)
    bdt.set_params(**options)
    return bdt

def scale_linear_bycolumn(rawpoints, high=100.0, low=0.0):
    mins = np.min(rawpoints, axis=0)
    maxs = np.max(rawpoints, axis=0)
    rng = maxs - mins
    return high - (((high - low) * (maxs - rawpoints)) / rng)

# path to where the data lies
dpath_train = 'data/training_1l_10ji3bi.csv'
dpath_test = 'data/test_1l_10ji3bi.csv'
modelfile = 'weights/sm4top_1l_10ji3bi_0'
output_file_name = 'output_1l_10ji3bi_0.root'


variable_size = 16

# load in training data, directly use numpy
dtrain = np.loadtxt( dpath_train, delimiter=',', skiprows=1, converters={variable_size+2: lambda x:int(x=='s'.encode('utf-8')) } )
dtest = np.loadtxt( dpath_test, delimiter=',', skiprows=1, converters={variable_size+2: lambda x:int(x=='s'.encode('utf-8')) } )
print ('finish loading from csv ')

label_train  = dtrain[:,variable_size+2]
label_test  = dtest[:,variable_size+2]
data_train   = dtrain[:,1:variable_size+1]
data_test   = dtest[:,1:variable_size+1]
weight_train = dtrain[:,variable_size+1] * float(len(label_test)) / float(len(label_train))
weight_test = dtest[:,variable_size+1] * float(len(label_test)) / float(len(label_train))

sum_wpos = sum( weight_train[i] for i in range(len(label_train)) if label_train[i] == 1.0  )
sum_wneg = sum( weight_train[i] for i in range(len(label_train)) if label_train[i] == 0.0  )

for i in range(len(weight_train)):
    if (label_train[i]==1.0):
        weight_train[i] = weight_train[i] * (0.5/sum_wpos)
    else:
        weight_train[i] = weight_train[i] * (0.5/sum_wneg)

data_train_scaled = scale_linear_bycolumn(data_train,-1.0,1.0)
data_test_scaled = scale_linear_bycolumn(data_test,-1.0,1.0)

noe = len(label_train)
weight_train = noe*weight_train
weight_test = noe*weight_test
        
sum_wpos_check = sum( weight_train[i] for i in range(len(label_train)) if label_train[i] == 1.0  )
sum_wneg_check = sum( weight_train[i] for i in range(len(label_train)) if label_train[i] == 0.0  )

# print weight statistics
print ('weight statistics: wpos=%g, wneg=%g, ratio=%g' % ( sum_wpos_check, sum_wneg_check, sum_wneg_check/sum_wpos_check ))

GradBDT = getGradBDT()
print (GradBDT)
GradBDT.fit(data_train_scaled,label_train,sample_weight = weight_train)

# evaluate test data
ypred_train = GradBDT.predict_proba( data_train_scaled )
ypred_test = GradBDT.predict_proba( data_test_scaled )

print ( type(ypred_train), ypred_train )

maxn = len(data_train[0])
print (maxn)

class_id = array('i', [ 0 ] )
variables = array('f', maxn*[ 0. ] )
weight = array('f', [ 0. ] )
score = array('f', [ 0. ] )

f = TFile(output_file_name,'recreate')
t1 = TTree('TestTree','Tree for Test Data')
t2 = TTree('TrainTree','Tree for Training Data')

t1.Branch('class_id',class_id,'class_id/I')
t1.Branch('variables',variables,'variables[]/F')
t1.Branch('weight',weight,'weight/F')
t1.Branch('score',score,'score/F')

t2.Branch('class_id',class_id,'class_id/I')
t2.Branch('variables',variables,'variables[]/F')
t2.Branch('weight',weight,'weight/F')
t2.Branch('score',score,'score/F')


for i in range(len(label_test)):
    class_id[0] = int(label_test[i])
    for j in range(maxn):
        variables[j] = data_test[i][j]
    weight[0] = weight_test[i]
    score[0] = ypred_test[i]
    t1.Fill()

for i in range(len(label_train)):
    class_id[0] = int(label_train[i])
    for j in range(maxn):
        variables[j] = data_train[i][j]
    weight[0] = weight_train[i]
    score[0] = ypred_train[i]
    t2.Fill()

t1.Write()
t2.Write()

f.Write()
f.Close()

print (ypred_train)
print (ypred_test)

print ('Number of Test Events: ', len(label_test))
print ('Number of Train Events: ', len(label_train))
print ('finish training')
