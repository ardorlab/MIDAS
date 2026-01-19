import os
import tensorflow as tf
from tensorflow.keras.layers import Input, Conv2D, MaxPooling2D, UpSampling2D, concatenate, Cropping2D
from tensorflow.keras.models import Model

# Architecture config
IMG_DIM = 153
INPUT_SHAPE = (IMG_DIM, IMG_DIM, 1)
SAVED_MODEL_NAME = 'pin_power_unet_deep_noscale.h5'
pathmodel = '/home/khnguy22/Deeponet-midas/MIDAS/surmodel/pinmodel/'

def build_unet_manually(input_shape):
    inputs = Input(input_shape)
    # Encoder
    c1 = Conv2D(16, (3, 3), activation='relu', padding='same')(inputs)
    c1 = Conv2D(16, (3, 3), activation='relu', padding='same')(c1)
    p1 = MaxPooling2D((2, 2), padding='same')(c1)
    
    c2 = Conv2D(32, (3, 3), activation='relu', padding='same')(p1)
    c2 = Conv2D(32, (3, 3), activation='relu', padding='same')(c2)
    p2 = MaxPooling2D((2, 2), padding='same')(c2)
    
    c3 = Conv2D(64, (3, 3), activation='relu', padding='same')(p2)
    c3 = Conv2D(64, (3, 3), activation='relu', padding='same')(c3)
    p3 = MaxPooling2D((2, 2), padding='same')(c3)
    
    # Bottleneck
    b = Conv2D(128, (3, 3), activation='relu', padding='same')(p3)
    b = Conv2D(128, (3, 3), activation='relu', padding='same')(b)
    
    # Decoder
    u3 = UpSampling2D((2, 2))(b)
    u3_cropped = Cropping2D(((1, 0), (1, 0)))(u3)
    u3 = concatenate([u3_cropped, c3])
    c4 = Conv2D(64, (3, 3), activation='relu', padding='same')(u3)
    c4 = Conv2D(64, (3, 3), activation='relu', padding='same')(c4)
    
    u2 = UpSampling2D((2, 2))(c4)
    u2_cropped = Cropping2D(((1, 0), (1, 0)))(u2)
    u2 = concatenate([u2_cropped, c2])
    c5 = Conv2D(32, (3, 3), activation='relu', padding='same')(u2)
    c5 = Conv2D(32, (3, 3), activation='relu', padding='same')(c5)
    
    u1 = UpSampling2D((2, 2))(c5)
    u1_cropped = Cropping2D(((1, 0), (1, 0)))(u1)
    u1 = concatenate([u1_cropped, c1])
    c6 = Conv2D(16, (3, 3), activation='relu', padding='same')(u1)
    c6 = Conv2D(16, (3, 3), activation='relu', padding='same')(c6)
    
    outputs = Conv2D(1, (1, 1), activation='linear')(c6)
    return Model(inputs=[inputs], outputs=[outputs])

def load_pin_model():
    full_path = os.path.join(pathmodel, SAVED_MODEL_NAME)
    model = build_unet_manually(INPUT_SHAPE)
    # Loading weights only avoids the "dtype" initializer error
    model.load_weights(full_path)
    return model

# We will call this inside the worker to prevent global hanging