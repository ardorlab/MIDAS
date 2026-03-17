'''
Load model and create some functions
Model 6 final
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
import deepxde as dde
import gc
import time as TT
import midas_data

import logging
# tf.get_logger().setLevel("ERROR")
# tf.config.optimizer.set_jit(False)
seed = 199478
tf.keras.backend.clear_session()
# tf.keras.utils.set_random_seed(seed) # fro tf 1.0
tf.set_random_seed(seed) # for tf 2.0 .compat.v1.set_random_seed(seed)
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


st = TT.time()
batch_size = 32
datapath = midas_data.__path_base_data__[0]



## randomly create test and train data
X1_train = np.random.rand(1,1296)
X2_train = np.random.rand(1,1296)
X3_train = np.random.rand(1,1296)
X4_train = np.random.rand(1,1296)
X5_train = np.random.rand(1,1296)
X6_train = np.random.rand(1,1296)
X7_train = np.random.rand(1,1296)
trunk = np.load(datapath+'trunk.npy', mmap_mode='r').astype(np.float32)
y_train = np.random.rand(1,1296)
### test data 
X1_test = np.random.rand(1,1296)
X2_test = np.random.rand(1,1296)
X3_test = np.random.rand(1,1296)
X4_test = np.random.rand(1,1296)
X5_test = np.random.rand(1,1296)
X6_test = np.random.rand(1,1296)
X7_test = np.random.rand(1,1296)
y_test  = np.random.rand(1,1296)
## load global max min value 
scalerparam = joblib.load(datapath+'scaler.joblib', mmap_mode='r') ## load joblib
# scalerparam = pickle.load(open(datapath+'scaler.pkl','rb'))
##
trunks = minmaxscale(trunk,scalerparam['tmax'], scalerparam['tmin'])

X_train = (X1_train,
           X2_train,
           X3_train, 
           X4_train, 
           X5_train, 
           X6_train, 
           X7_train, trunks)
y_train = y_train
X_test =  (X1_test,
           X2_test,
           X3_test, 
           X4_test, 
           X5_test, 
           X6_test, 
           X7_test, trunks)
y_test = y_test

from deepxde.nn.tensorflow_compat_v1.mionet import MIONetCartesianProd_custom7
from deepxde.data.quadruple import QuadrupleCartesianProd
data = PWRCustomCartesianProd(X_train, y_train, X_test, y_test)
m = 1296
net = dde.maps.mionet.MIONetCartesianProd_custom7(
	[m,800, 600,400,150], [m,800, 600,400,150],
	[m,800, 600,400,150], [m,800, 600,400,150],
	[m,800, 600,400,150], [m,800, 600,400,150],
	[m,800, 600,400,150],
    [3, 500,300,150],
     {"branch1":  "relu", "branch2": "relu", "branch3": "relu",
      "branch4":  "relu", "branch5": "relu", "branch6": "relu",
      "branch7": "relu", "trunk":  "relu"},
       "Glorot normal"
)
model = dde.Model(data, net)
model.compile(loss='customMSE', ## old is MSE
	optimizer="adam",
	lr=5e-4,
	# decay=("inverse time", 1, 1e-4),
	#metrics=["l2 relative error"],
	metrics=["MSE"],
	# metrics=["accuracy"],
)
# model.restore(r'./PWR-model03/MIONet_PWR3D_03-200000.ckpt') 
model.restore(midas_data.__path_to_save_model__[0]+"-"+str(midas_data.__training_epochs_model_1__)+".ckpt") 
print('Load model rpf time: ', TT.time()-st)