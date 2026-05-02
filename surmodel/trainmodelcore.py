'''
Train PWR 3D core on Boron and cycle length
'''

'''
Model 
'''
import os
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"
os.environ["AUTOGRAPH_VERBOSITY"] = "0"
os.environ["TF_AUTOGRAPH_DISABLE"] = "1"
os.environ["TMPDIR"] = "/tmp"
import sys
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from deepxde.backend import tf
import joblib
import pickle
import shutil

import logging
tf.get_logger().setLevel("ERROR")
tf.config.optimizer.set_jit(False)
import os
import gc
import deepxde as dde
from deepxde.nn.tensorflow_compat_v1.mionet import MIONetCartesianProd_custom7
from deepxde.data.quadruple import QuadrupleCartesianProd
import midas_data 

# dde.config.set_default_float("float32")


BatchSampler = dde.data.sampler.BatchSampler
Data = dde.data.Data
## create new class
class PWRCustomCartesianProd(Data):
    """Cartesian Product input data format for MIONet custom architecture for the Rx core test.

    Args:
        X_train: A tuple of five NumPy arrays.
        The first element has the shape (`N1`,`dim1`),
        the second element has the shape (`N1`, `dim2`),
        ... to the forth

        and the last
            element has the shape (`N2`, `dim3`).
        y_train: A NumPy array of shape (`N1`, `N2`).
    """

    def __init__(self, X_train, y_train, X_test, y_test):
        # ## need to check for general release esle no
        # if (
        #     len(X_train[0]) * len(X_train[2]) != y_train.size
        #     or len(X_train[1]) * len(X_train[2]) != y_train.size
        #     or len(X_train[0]) != len(X_train[1])
        # ):
        #     raise ValueError(
        #         "The training dataset does not have the format of Cartesian product."
        #     )
        # if (
        #     len(X_test[0]) * len(X_test[2]) != y_test.size
        #     or len(X_test[1]) * len(X_test[2]) != y_test.size
        #     or len(X_test[0]) != len(X_test[1])
        # ):
        #     raise ValueError(
        #         "The testing dataset does not have the format of Cartesian product."
        #     )

        self.train_x, self.train_y = X_train, y_train
        self.test_x, self.test_y = X_test, y_test

        self.branch_sampler = BatchSampler(len(X_train[0]), shuffle=True)
        self.trunk_sampler = BatchSampler(len(X_train[-1]), shuffle=True) # take the last array as it is the trunk data

    def losses(self, targets, outputs, loss_fn, inputs, model, aux=None):
        return loss_fn(targets, outputs)

    def train_next_batch(self, batch_size=None):
        if batch_size is None:
            return self.train_x, self.train_y
        if not isinstance(batch_size, (tuple, list)):
            indices = self.branch_sampler.get_next(batch_size)
            return (
                self.train_x[0][indices],
                self.train_x[1][indices],
                self.train_x[2][indices],
                self.train_x[3][indices],
                self.train_x[4][indices],
                self.train_x[5][indices],
                self.train_x[6][indices],
                self.train_x[7],
            ), self.train_y[indices]
        indices_branch = self.branch_sampler.get_next(batch_size[0])
        indices_trunk = self.trunk_sampler.get_next(batch_size[1])
        return (
            self.train_x[0][indices_branch],
            self.train_x[1][indices_branch],
            self.train_x[2][indices_branch],
            self.train_x[3][indices_branch],
            self.train_x[4][indices_branch],
            self.train_x[5][indices_branch],
            self.train_x[6][indices_branch],
            self.train_x[7][indices_trunk],
        ), self.train_y[indices_branch][:, indices_trunk]

    def test(self):
        return self.test_x, self.test_y

def minmaxscale(X, globalmax, globalmin):
    '''
    Scaling X (need to be flatten to 2D) based on global min/max
    '''
    X = X.reshape(X.shape[0],-1)
    temp = (X-globalmin)/(globalmax-globalmin + 1e-20) ## prevent zero division

    return temp

def minmaxReverse(Xscaled, globalmax, globalmin):
    '''
    Scale back to original array
    '''
    temp = Xscaled * (globalmax-globalmin+1e-20) + globalmin

    return temp

def create_scaler(all_data,r=0.5,argmin=1.0, argmax=1.0):
    train_min = all_data.min(axis=0)
    train_max = all_data.max(axis=0)
    range_extension = r * (train_max - train_min)
    global_min = argmin*train_min - range_extension
    global_max = argmax*train_max + range_extension
    return global_min, global_max

def traincoremodel_update():
    tf.keras.backend.clear_session()
    tf.compat.v1.reset_default_graph()
    gc.collect()
    batch_size = 32
    seed = 199478
    tf.set_random_seed(seed) # for tf 2.0 .compat.v1.set_random_seed(seed)
    datapath = midas_data.__path_new_data__[1]
    datapath_new = midas_data.__path_to_store_retrain_data__
    X1 = np.load(datapath_new+'branchXS_1.npy').astype(np.float32)
    X2 = np.load(datapath_new+'branchXS_2.npy').astype(np.float32)
    X3 = np.load(datapath_new+'branchXS_3.npy').astype(np.float32)
    X4 = np.load(datapath_new+'branchXS_4.npy').astype(np.float32)
    X5 = np.load(datapath_new+'branchXS_5.npy').astype(np.float32)
    X6 = np.load(datapath_new+'branchXS_6.npy').astype(np.float32)
    X7 = np.load(datapath_new+'branchXS_7.npy').astype(np.float32)
    Y = np.load(datapath_new+'outputcore.npy').astype(np.float32)
    # trunk = np.load(datapath+'trunk.npy').astype(np.float32)
    
    ## reshape
    X1 = X1.reshape(X1.shape[0], -1)
    X2 = X2.reshape(X2.shape[0], -1)
    X3 = X3.reshape(X3.shape[0], -1)
    X4 = X4.reshape(X4.shape[0], -1)
    X5 = X5.reshape(X5.shape[0], -1)
    X6 = X6.reshape(X6.shape[0], -1)
    X7 = X7.reshape(X7.shape[0], -1)
    Y  = Y.reshape(Y.shape[0], -1)
    ## create scaler 


    # ## set negative value to 0??
    # X1[X1<0]=0
    # X2[X2<0]=0
    # X3[X3<0]=0
    # X4[X4<0]=0
    # X5[X5<0]=0
    # X6[X6<0]=0
    # X7[X7<0]=0
    ## load global max min value
    # #scalerparam = pickle.load(open(datapath+'scaler50.pkl','rb'))
    # b1min, b1max = create_scaler(X1)
    # b2min, b2max = create_scaler(X2)
    # b3min, b3max = create_scaler(X3)
    # b4min, b4max = create_scaler(X4)
    # b5min, b5max = create_scaler(X5)
    # b6min, b6max = create_scaler(X6)
    # b7min, b7max = create_scaler(X7)
    # omin, omax = create_scaler(Y,r=0.5,argmin=1.0, argmax=1.0)
    # scalerparam ={"b1min":b1min,
    #               "b1max":b1max,
    #               "b2min":b2min,
    #               "b2max":b2max,
    #               "b3min":b3min,
    #               "b3max":b3max,
    #               "b4min":b4min,
    #               "b4max":b4max,
    #               "b5min":b5min,
    #               "b5max":b5max,
    #               "b6min":b6min,
    #               "b6max":b6max,
    #               "b7min":b7min,
    #               "b7max":b7max,
    #               "omin":omin,
    #               "omax":omax, 
    #               }
    ## dumping scaler
    #joblib.dump(scalerparam, datapath+'scaler.joblib')
    scalerparam = joblib.load(datapath+'scaler.joblib')
    X1s = minmaxscale(X1,scalerparam['b1max'], scalerparam['b1min'])
    X2s = minmaxscale(X2,scalerparam['b2max'], scalerparam['b2min'])
    X3s = minmaxscale(X3,scalerparam['b3max'], scalerparam['b3min'])
    X4s = minmaxscale(X4,scalerparam['b4max'], scalerparam['b4min'])
    X5s = minmaxscale(X5,scalerparam['b5max'], scalerparam['b5min'])
    X6s = minmaxscale(X6,scalerparam['b6max'], scalerparam['b6min'])
    X7s = minmaxscale(X7,scalerparam['b7max'], scalerparam['b7min'])
    Ys  = minmaxscale(Y,scalerparam['omax'], scalerparam['omin'])
    ## splitting
    # Then split
    (X1_train, X1_test,
     X2_train, X2_test,
     X3_train, X3_test,
     X4_train, X4_test,
     X5_train, X5_test,
     X6_train, X6_test,
     X7_train, X7_test,
     y_train, y_test) = train_test_split(
        X1s, X2s, X3s, X4s,
        X5s, X6s, X7s, Ys,
        test_size=0.3,
        random_state=42,
        shuffle=True
    ) ## no scale on Y
    trunks = np.array([1,2]).reshape(2,1)

    X_train = (X1_train,
               X2_train,
               X3_train,
               X4_train,
               X5_train,
               X6_train,
               X7_train, trunks)
    # y_train = y_train.reshape(y_train.shape[0],-1) # no scaling
    y_train = y_train # scaling
    X_test =  (X1_test,
               X2_test,
               X3_test,
               X4_test,
               X5_test,
               X6_test,
               X7_test, trunks)
    # y_test = y_test.reshape(y_test.shape[0],-1) # no scaling
    y_test = y_test

    data = PWRCustomCartesianProd(X_train, y_train, X_test, y_test)
    m = 1296
    net = dde.maps.mionet.MIONetCartesianProd_custom7(
    	[m,800, 600,400,200], [m,800, 600,400,200],
    	[m,800, 600,400,200], [m,800, 600,400,200],
    	[m,800, 600,400,200], [m,800, 600,400,200],
    	[m,800, 600,400,200],
        [1, 500,300,200],
         {"branch1":  "relu", "branch2": "relu", "branch3": "relu",
          "branch4":  "relu", "branch5": "relu", "branch6": "relu",
          "branch7": "relu", "trunk":  "relu"},
           "Glorot normal"
    )
    model = dde.Model(data, net)
    # early_stopping_callback = dde.callbacks.EarlyStopping(
    #     monitor="loss_train",
    #     min_delta=1e-5,
    #     patience=5000,
    # )

    model.compile(
    	"adam",
    	lr=5e-4,
    	#decay=("inverse time", 1, 1e-4),
    	#metrics=["l2 relative error"],
    	metrics=["mae"],
    	# metrics=["accuracy"],
    )
    iter = midas_data.__training_epochs_model_2__
    # # # # # Compile and Train
    model.train(epochs=iter,batch_size=batch_size,
                model_restore_path=midas_data.__path_to_restore_model__[1],
                #    callbacks=[early_stopping_callback],
                model_save_path=midas_data.__path_to_save_model__[1])
    
    ## update path 
    temp = midas_data.__path_to_save_model__[1]+"-"+str(iter)+".ckpt"
    midas_data.__path_to_restore_model__[1] = temp 
    # return
    
