'''
Training script for pin power correction with UNET-CNN model 
Update log: Convert to TF 2.15
'''
import h5py
import numpy as np
import os
import matplotlib.pyplot as plt
import time

import deepxde as dde # for compability when integrated with DeepONet surrogate models

import tensorflow.compat.v1 as tf
tf.disable_eager_execution()          
tf.disable_v2_behavior()              

from tensorflow.compat.v1.keras.models import Model
from tensorflow.compat.v1.keras.layers import (
    Input, Conv2D, MaxPooling2D, UpSampling2D,
    concatenate, Cropping2D, Dropout
)
from tensorflow.compat.v1.keras.callbacks import EarlyStopping, ModelCheckpoint
from tensorflow.compat.v1.keras.optimizers import Adam
import tensorflow.compat.v1.keras.backend as K

from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score

print(f"DeepXDE Backend : {dde.backend.backend_name}")
print(f"TensorFlow Version: {tf.__version__}")
print(f"Eager execution enabled: {tf.executing_eagerly()}")  

# --- Configuration, chnage base on user's preference ---
H5_FILE           = 'pin_data_ne512.h5'
INPUT_DATASET     = 'calpin'
OUTPUT_DATASET    = 'parcpin'
SAVED_MODEL_NAME  = 'pin_power_unet_deep_noscale_updated_ne512.h5'
OLD_MODEL_NAME    = '/home/khnguy22/Deeponet-midas/MIDAS/surmodel/pinmodel/pin_power_unet_deep_noscale.h5'
PLOT_RESULTS_FILE = 'model_predictions_unet_deep_diff_ne512.png'

IMG_DIM          = 153
INPUT_SHAPE      = (IMG_DIM, IMG_DIM, 1)
BATCH_SIZE       = 16
EPOCHS           = 10
VALIDATION_SPLIT = 0.2

print(f"Loading data from {H5_FILE}...")
with h5py.File(H5_FILE, 'r') as hf:
    X_data = hf[INPUT_DATASET][:]
    y_data = hf[OUTPUT_DATASET][:]

print(f"Original X shape: {X_data.shape}, Original y shape: {y_data.shape}")

sampled_arrX, _, sampled_arrY, _ = train_test_split(
    X_data, y_data, train_size=0.2, random_state=42
)
# only take 20% of data for CNN model training
print(f"Sample X shape: {sampled_arrX.shape}, Sample y shape: {sampled_arrY.shape}")

X_train, X_test, y_train, y_test = train_test_split(
    sampled_arrX, sampled_arrY, test_size=VALIDATION_SPLIT, random_state=42
)

print("Reshaping data...")
X_train_cnn = X_train.reshape(-1, IMG_DIM, IMG_DIM, 1)
y_train_cnn = y_train.reshape(-1, IMG_DIM, IMG_DIM, 1)
X_test_cnn  = X_test.reshape(-1, IMG_DIM, IMG_DIM, 1)
y_test_cnn  = y_test.reshape(-1, IMG_DIM, IMG_DIM, 1)

print(f"Training samples: {X_train_cnn.shape[0]}, Testing samples: {X_test_cnn.shape[0]}")


def create_composite_loss(alpha=0.1):
    '''Custom Loss function here for CNN model'''
    def masked_mse_loss(y_true, y_pred):
        mask = tf.cast(tf.not_equal(y_true, 0), tf.float32)
        masked_sq_err = tf.square(y_true - y_pred) * mask
        return tf.reduce_sum(masked_sq_err) / (tf.reduce_sum(mask) + 1e-7)

    def peak_loss(y_true, y_pred):
        mask = tf.cast(tf.not_equal(y_true, 0), tf.float32)
        max_true = tf.reduce_max(y_true,        axis=[1, 2, 3])
        max_pred = tf.reduce_max(y_pred * mask, axis=[1, 2, 3])
        return tf.reduce_mean(tf.square(max_true - max_pred))

    def composite_loss(y_true, y_pred):
        return (1.0 - alpha) * masked_mse_loss(y_true, y_pred) \
                     + alpha * peak_loss(y_true, y_pred)

    composite_loss.__name__ = 'composite_loss'
    return composite_loss

composite_loss_fn = create_composite_loss(alpha=0.1)

def build_deep_unet(input_shape):
    inputs = Input(input_shape)
    c1 = Conv2D(16,  (3,3), activation='relu', kernel_initializer='he_normal', padding='same')(inputs)
    c1 = Conv2D(16,  (3,3), activation='relu', kernel_initializer='he_normal', padding='same')(c1)
    p1 = MaxPooling2D((2,2), padding='same')(c1)
    c2 = Conv2D(32,  (3,3), activation='relu', kernel_initializer='he_normal', padding='same')(p1)
    c2 = Conv2D(32,  (3,3), activation='relu', kernel_initializer='he_normal', padding='same')(c2)
    p2 = MaxPooling2D((2,2), padding='same')(c2)
    c3 = Conv2D(64,  (3,3), activation='relu', kernel_initializer='he_normal', padding='same')(p2)
    c3 = Conv2D(64,  (3,3), activation='relu', kernel_initializer='he_normal', padding='same')(c3)
    p3 = MaxPooling2D((2,2), padding='same')(c3)
    b  = Conv2D(128, (3,3), activation='relu', kernel_initializer='he_normal', padding='same')(p3)
    b  = Conv2D(128, (3,3), activation='relu', kernel_initializer='he_normal', padding='same')(b)
    u3 = UpSampling2D((2,2))(b)
    u3 = concatenate([Cropping2D(((1,0),(1,0)))(u3), c3])
    c4 = Conv2D(64,  (3,3), activation='relu', kernel_initializer='he_normal', padding='same')(u3)
    c4 = Conv2D(64,  (3,3), activation='relu', kernel_initializer='he_normal', padding='same')(c4)
    u2 = UpSampling2D((2,2))(c4)
    u2 = concatenate([Cropping2D(((1,0),(1,0)))(u2), c2])
    c5 = Conv2D(32,  (3,3), activation='relu', kernel_initializer='he_normal', padding='same')(u2)
    c5 = Conv2D(32,  (3,3), activation='relu', kernel_initializer='he_normal', padding='same')(c5)
    u1 = UpSampling2D((2,2))(c5)
    u1 = concatenate([Cropping2D(((1,0),(1,0)))(u1), c1])
    c6 = Conv2D(16,  (3,3), activation='relu', kernel_initializer='he_normal', padding='same')(u1)
    c6 = Conv2D(16,  (3,3), activation='relu', kernel_initializer='he_normal', padding='same')(c6)
    outputs = Conv2D(1, (1,1), activation='linear')(c6)
    return Model(inputs=[inputs], outputs=[outputs])

custom_objects = {'composite_loss': composite_loss_fn}

if os.path.exists(OLD_MODEL_NAME):
    print(f"Loading existing model from {OLD_MODEL_NAME} ...")
    model = tf.keras.models.load_model(
        OLD_MODEL_NAME,
        custom_objects=custom_objects,
        compile=False
    )
    print("Old model loaded successfully.")
else:
    print(f"'{OLD_MODEL_NAME}' not found — building fresh U-Net.")
    model = build_deep_unet(INPUT_SHAPE)

model.compile(
    optimizer=Adam(learning_rate=0.0005),
    loss=composite_loss_fn,
    metrics=['mean_absolute_error']
)
model.summary()

print("Training model...")
early_stopper    = EarlyStopping(patience=5, monitor='val_loss', mode='min', verbose=1)
model_checkpoint = ModelCheckpoint(
    SAVED_MODEL_NAME, monitor='val_loss', mode='min',
    save_best_only=True, verbose=1
)

st = time.time()
history = model.fit(
    X_train_cnn, y_train_cnn,
    batch_size=BATCH_SIZE,
    epochs=EPOCHS,
    validation_data=(X_test_cnn, y_test_cnn),
    callbacks=[early_stopper, model_checkpoint],
    verbose=1
)
print(f"Training complete. Best model saved to {SAVED_MODEL_NAME} "
      f"({time.time() - st:.1f}s)")

print("Evaluating model...")
best_model = tf.keras.models.load_model(
    SAVED_MODEL_NAME,
    custom_objects=custom_objects,
    compile=False
)

y_pred_raw  = best_model.predict(X_test_cnn)
keep_mask   = (y_test[0] != 0).reshape(1, IMG_DIM, IMG_DIM, 1)
y_pred_corr = y_pred_raw * keep_mask

mse = mean_squared_error(y_test.flatten(), y_pred_corr.flatten())
mae = mean_absolute_error(y_test.flatten(), y_pred_corr.flatten())
r2  = r2_score(y_test.flatten(), y_pred_corr.flatten())

print("\n--- Evaluation Metrics (Post-Processed, No Scale) ---")
print(f"  MSE : {mse:.6e}")
print(f"  MAE : {mae:.6f}")
print(f"  R²  : {r2:.6f}")
print("-----------------------------------------------------\n")

print("Generating prediction plots...")
y_pred_img = y_pred_corr.reshape(-1, IMG_DIM, IMG_DIM)
X_test_img = X_test.reshape(-1, IMG_DIM, IMG_DIM)

num_examples = 3
indices = np.random.choice(len(X_test_cnn), num_examples, replace=False)

fig, axes = plt.subplots(num_examples, 4, figsize=(20, 5 * num_examples))
fig.suptitle("DEEP U-Net Prediction (No Scaling) vs. Ground Truth", fontsize=16)

for i, idx in enumerate(indices):
    input_img = X_test_img[idx]
    true_img  = y_test[idx]
    pred_img  = y_pred_img[idx]
    diff_img  = pred_img - true_img

    print(f"Sample {idx} — Max True: {true_img.max():.4f}, "
          f"Max Pred: {pred_img.max():.4f}, Max Input: {input_img.max():.4f}")

    ax = axes[i, 0]
    im0 = ax.imshow(input_img, cmap='viridis', origin='lower')
    ax.set_title(f"Input (Sample {idx})")
    fig.colorbar(im0, ax=ax, label='Pin Power')

    ax = axes[i, 1]
    im1 = ax.imshow(true_img, cmap='viridis', origin='lower')
    ax.set_title("True")
    fig.colorbar(im1, ax=ax, label='Pin Power')

    ax = axes[i, 2]
    im2 = ax.imshow(pred_img, cmap='viridis', origin='lower')
    ax.set_title("Predicted (Corrected)")
    fig.colorbar(im2, ax=ax, label='Pin Power')

    ax = axes[i, 3]
    diff_limit = max(abs(diff_img.min()), abs(diff_img.max()))
    im3 = ax.imshow(diff_img, cmap='bwr', origin='lower',
                    vmin=-diff_limit, vmax=diff_limit)
    ax.set_title("Difference (Pred - True)")
    fig.colorbar(im3, ax=ax, label='Error')

plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.savefig(PLOT_RESULTS_FILE)
print(f"Saved plots to {PLOT_RESULTS_FILE}")