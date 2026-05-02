## train pin with deeper unet (No Scaling + Difference Plot)
import h5py
import numpy as np
import joblib
import os
import matplotlib.pyplot as plt

# --- DeepXDE Backend Imports ---
import deepxde as dde
from deepxde.backend import tf
# ---

# --- START FIX (Eager Execution & Imports) ---
tf.enable_eager_execution()

Model = tf.keras.models.Model
Input = tf.keras.layers.Input
Conv2D = tf.keras.layers.Conv2D
Dropout = tf.keras.layers.Dropout
MaxPooling2D = tf.keras.layers.MaxPooling2D
UpSampling2D = tf.keras.layers.UpSampling2D
concatenate = tf.keras.layers.concatenate
Cropping2D = tf.keras.layers.Cropping2D
EarlyStopping = tf.keras.callbacks.EarlyStopping
ModelCheckpoint = tf.keras.callbacks.ModelCheckpoint
Adam = tf.keras.optimizers.Adam
custom_object_scope = tf.keras.utils.custom_object_scope
# --- END FIX ---

from sklearn.model_selection import train_test_split
# MinMaxScaler removed
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score

# --- 0. Verify TF Version ---
print(f"DeepXDE Backend: {dde.backend.backend_name}")
print(f"TensorFlow Version: {tf.__version__}")
if not tf.__version__.startswith("1."):
    print("Warning: This script is intended for TensorFlow 1.x.")

# --- 1. Configuration ---
# --- 1. Configuration ---
H5_FILE = 'pin_data_157case05.h5'
INPUT_DATASET = 'calpin'
OUTPUT_DATASET = 'parcpin'

# --- 2. Load and Prepare Data ---
print(f"Loading data from {H5_FILE}...")

with h5py.File(H5_FILE, 'r') as hf:
    # Load as (num_samples, 153, 153)
    X_data = hf[INPUT_DATASET][:]
    y_data = hf[OUTPUT_DATASET][:]

print(f"Original X shape: {X_data.shape}")
print(f"Original y shape: {y_data.shape}")
tetetetetet = y_data.reshape(37,16,153,153,1)
tetetetetet=tetetetetet[:,::-1,:,:,:]
print(np.max(tetetetetet, axis=(1, 2, 3, 4)))   # shape: (M,)
print(np.unravel_index(np.argmax(tetetetetet), tetetetetet.shape))
