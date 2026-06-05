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
H5_FILE = 'pin_data_ne512.h5'
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

### load input/ouput 
# INPUT_DATASET = np.load('./testdata/calculatedcoremap.npy')
# OUTPUT_DATASET =np.load('./testdata/PARCScoremap.npy')
# N=INPUT_DATASET.shape[0]
# ## 
# INPUT_DATASET = INPUT_DATASET.reshape(100,34,16,153,153)
# INPUT_DATASET = INPUT_DATASET[:,:,::-1,:,:]
# INPUT_DATASET = INPUT_DATASET.reshape(N,153,153)



SAVED_MODEL_NAME = 'pin_power_unet_deep_noscale_updated_ne512.h5' 
OLD_MODEL_NAME = '/home/khnguy22/Deeponet-midas/MIDAS/surmodel/pinmodel/pin_power_unet_deep_noscale.h5' 
PLOT_RESULTS_FILE = 'model_predictions_unet_deep_diff_ne512.png' 

IMG_DIM = 153 # 153x153
INPUT_SHAPE = (IMG_DIM, IMG_DIM, 1)
BATCH_SIZE = 16
EPOCHS = 10
VALIDATION_SPLIT = 0.2

# --- 2. Load and Prepare Data ---
# print(f"Loading data from {H5_FILE}...")
# with h5py.File(H5_FILE, 'r') as hf:
#     X_data = hf[INPUT_DATASET][:]
#     y_data = hf[OUTPUT_DATASET][:]
# X_data=INPUT_DATASET
# y_data = OUTPUT_DATASET
# stop
# N = X_data.shape[0] 
# X_data = X_data.reshape(500,34,16,153,153)
# X_data = X_data[:,:,::-1,:,:]
# X_data = X_data.reshape(N,153,153)



print(f"Original X shape: {X_data.shape}, Original y shape: {y_data.shape}")

### take a portion of it 
sampled_arrX, restX, sampled_arrY, restY = train_test_split(
    X_data, y_data, train_size=0.2, random_state=42
)
print(f"Sample X shape: {sampled_arrX.shape}, Sample y shape: {sampled_arrY.shape}")

X_train, X_test, y_train, y_test = train_test_split(
    sampled_arrX, sampled_arrY, test_size=VALIDATION_SPLIT, random_state=42
)

# --- 3. Reshape Data (NO SCALING) ---
print("Reshaping data (No Scaling)...")

# Directly reshape to (N, 153, 153, 1)
X_train_cnn = X_train.reshape(-1, IMG_DIM, IMG_DIM, 1)
y_train_cnn = y_train.reshape(-1, IMG_DIM, IMG_DIM, 1)
X_test_cnn = X_test.reshape(-1, IMG_DIM, IMG_DIM, 1)
y_test_cnn = y_test.reshape(-1, IMG_DIM, IMG_DIM, 1)

print(f"Training samples: {X_train_cnn.shape[0]}, Testing samples: {X_test_cnn.shape[0]}")

# --- 4. Build the DEEPER U-Net Model ---
print("Building DEEPER U-Net model...")

def build_deep_unet(input_shape):
    inputs = Input(input_shape)
    # Encoder
    c1 = Conv2D(16, (3, 3), activation='relu', kernel_initializer='he_normal', padding='same')(inputs)
    c1 = Conv2D(16, (3, 3), activation='relu', kernel_initializer='he_normal', padding='same')(c1)
    p1 = MaxPooling2D((2, 2), padding='same')(c1)
    c2 = Conv2D(32, (3, 3), activation='relu', kernel_initializer='he_normal', padding='same')(p1)
    c2 = Conv2D(32, (3, 3), activation='relu', kernel_initializer='he_normal', padding='same')(c2)
    p2 = MaxPooling2D((2, 2), padding='same')(c2)
    c3 = Conv2D(64, (3, 3), activation='relu', kernel_initializer='he_normal', padding='same')(p2)
    c3 = Conv2D(64, (3, 3), activation='relu', kernel_initializer='he_normal', padding='same')(c3)
    p3 = MaxPooling2D((2, 2), padding='same')(c3)
    # Bottleneck
    b = Conv2D(128, (3, 3), activation='relu', kernel_initializer='he_normal', padding='same')(p3)
    b = Conv2D(128, (3, 3), activation='relu', kernel_initializer='he_normal', padding='same')(b)
    # Decoder
    u3 = UpSampling2D((2, 2))(b)
    u3_cropped = Cropping2D(((1, 0), (1, 0)))(u3)
    u3 = concatenate([u3_cropped, c3])
    c4 = Conv2D(64, (3, 3), activation='relu', kernel_initializer='he_normal', padding='same')(u3)
    c4 = Conv2D(64, (3, 3), activation='relu', kernel_initializer='he_normal', padding='same')(c4)
    u2 = UpSampling2D((2, 2))(c4)
    u2_cropped = Cropping2D(((1, 0), (1, 0)))(u2)
    u2 = concatenate([u2_cropped, c2])
    c5 = Conv2D(32, (3, 3), activation='relu', kernel_initializer='he_normal', padding='same')(u2)
    c5 = Conv2D(32, (3, 3), activation='relu', kernel_initializer='he_normal', padding='same')(c5)
    u1 = UpSampling2D((2, 2))(c5)
    u1_cropped = Cropping2D(((1, 0), (1, 0)))(u1)
    u1 = concatenate([u1_cropped, c1])
    c6 = Conv2D(16, (3, 3), activation='relu', kernel_initializer='he_normal', padding='same')(u1)
    c6 = Conv2D(16, (3, 3), activation='relu', kernel_initializer='he_normal', padding='same')(c6)
    outputs = Conv2D(1, (1, 1), activation='linear')(c6)
    model = Model(inputs=[inputs], outputs=[outputs])
    return model

# --- Custom Composite Loss Function ---
def create_composite_loss(alpha=0.5):
    def masked_mse_loss(y_true, y_pred):
        mask = tf.cast(tf.not_equal(y_true, 0), tf.float32)
        squared_error = tf.square(y_true - y_pred)
        masked_squared_error = squared_error * mask
        # Add epsilon for numerical stability
        loss = tf.reduce_sum(masked_squared_error) / (tf.reduce_sum(mask) + 1e-7)
        return loss

    def peak_loss(y_true, y_pred):
        mask = tf.cast(tf.not_equal(y_true, 0), tf.float32)
        max_true = tf.reduce_max(y_true, axis=[1, 2, 3])
        max_pred = tf.reduce_max(y_pred * mask, axis=[1, 2, 3])
        peak_loss = tf.reduce_mean(tf.square(max_true - max_pred))
        return peak_loss

    def composite_loss(y_true, y_pred):
        mse = masked_mse_loss(y_true, y_pred)
        peak = peak_loss(y_true, y_pred)
        return (1.0 - alpha) * mse + alpha * peak
    
    return composite_loss

# NOTE: When data is NOT scaled, loss values will be larger (since values > 1).
# You might need to adjust alpha or learning rate if training is unstable.
# Here we keep alpha=0.3
composite_loss_fn = create_composite_loss(alpha=0.1)
model = build_deep_unet(INPUT_SHAPE)
optimizer = Adam(lr=0.0005) 
model.compile(optimizer=optimizer, loss=composite_loss_fn, metrics=['mean_absolute_error'])
model.summary()
## load model first
# Load the best model
with custom_object_scope({'composite_loss': composite_loss_fn}):
    model = tf.keras.models.load_model(OLD_MODEL_NAME)
    print('loaded old model completed!')
# --- 5. Train the Model ---
print("Training model...")
early_stopper = EarlyStopping(patience=5, monitor='val_loss', mode='min', verbose=1)
model_checkpoint = ModelCheckpoint(SAVED_MODEL_NAME, monitor='val_loss', mode='min', save_best_only=True, verbose=1)
import time 
st = time.time()
history = model.fit(
   X_train_cnn, y_train_cnn, 
   batch_size=BATCH_SIZE,
   epochs=EPOCHS,
   validation_data=(X_test_cnn, y_test_cnn),
   callbacks=[early_stopper, model_checkpoint],
   verbose=1
)
print(f"Model training complete. Best model saved to {SAVED_MODEL_NAME}, time completed {time.time()-st}")

# --- 6. Evaluate Metrics (with Post-Processing Correction) ---
print("Evaluating model performance with post-processing correction...")

# Load the best model
with custom_object_scope({'composite_loss': composite_loss_fn}):
    best_model = tf.keras.models.load_model(SAVED_MODEL_NAME)

# Get predictions (These are already in original scale)
y_pred_raw_cnn = best_model.predict(X_test_cnn)

# --- CORRECTION STEP ---
# Create a mask of "pixels to keep" (non-zero) from y_test (which is already unscaled)
keep_mask = (y_test[0] != 0) 
keep_mask_4d = keep_mask.reshape(1, IMG_DIM, IMG_DIM, 1) 

# Apply the mask to the predictions
y_pred_corrected_cnn = y_pred_raw_cnn * keep_mask_4d

# Flatten for metrics
y_test_flat = y_test.flatten()
y_pred_corrected_flat = y_pred_corrected_cnn.flatten()

# Calculate metrics
mse = mean_squared_error(y_test_flat, y_pred_corrected_flat)
mae = mean_absolute_error(y_test_flat, y_pred_corrected_flat)
r2 = r2_score(y_test_flat, y_pred_corrected_flat)

print("\n--- Model Evaluation Metrics (Post-Processed, No Scale) ---")
print(f"  Mean Squared Error (MSE): {mse:.6e}")
print(f"  Mean Absolute Error (MAE): {mae:.6f}")
print(f"  R-squared (R²) Score:     {r2:.6f}")
print("-------------------------------------------------\n")

# --- 7. Visualize Results (with Difference Plot) ---
print("Generating and saving prediction plots...")

# Reshape for plotting
y_pred_corrected_img = y_pred_corrected_cnn.reshape(-1, IMG_DIM, IMG_DIM)
X_test_img = X_test.reshape(-1, IMG_DIM, IMG_DIM)

# Plot 3 random examples with 4 columns
num_examples = 3
indices = np.random.choice(range(len(X_test_cnn)), num_examples, replace=False)

fig, axes = plt.subplots(num_examples, 4, figsize=(20, 5 * num_examples))
fig.suptitle("DEEP U-Net Prediction (No Scaling) vs. Ground Truth", fontsize=16)

for i, idx in enumerate(indices):
    input_img = X_test_img[idx]
    true_img = y_test[idx]
    pred_img = y_pred_corrected_img[idx]
    
    # Calculate Difference (Prediction - Truth)
    diff_img = pred_img - true_img

    # Determine consistent min/max for the main plots
    vmin = min(true_img.min(), pred_img.min())
    vmax = max(true_img.max(), pred_img.max())
    
    print(f'Sample {idx} - Max True: {true_img.max():.4f}, Max Pred: {pred_img.max():.4f}, Max old {input_img.max():.4f}')

    # 1. Plot Input
    ax = axes[i, 0]
    im0 = ax.imshow(input_img, cmap='viridis', origin='lower')
    ax.set_title(f"Input (Sample {idx})")
    fig.colorbar(im0, ax=ax, label='Pin Power')

    # 2. Plot Ground Truth
    ax = axes[i, 1]
    im1 = ax.imshow(true_img, cmap='viridis', origin='lower')#, vmin=vmin, vmax=vmax)
    ax.set_title(f"True")
    fig.colorbar(im1, ax=ax, label='Pin Power')

    # 3. Plot Predicted
    ax = axes[i, 2]
    im2 = ax.imshow(pred_img, cmap='viridis', origin='lower')#, vmin=vmin, vmax=vmax)
    ax.set_title(f"Predicted (Corrected)")
    fig.colorbar(im2, ax=ax, label='Pin Power')
    
    # 4. Plot Difference (Pred - True)
    ax = axes[i, 3]
    # Center the colormap at 0 to show + (red) and - (blue) errors
    diff_limit = max(abs(diff_img.min()), abs(diff_img.max()))
    im3 = ax.imshow(diff_img, cmap='bwr', origin='lower', vmin=-diff_limit, vmax=diff_limit)
    ax.set_title(f"Difference (Pred - True)")
    fig.colorbar(im3, ax=ax, label='Error')

plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.savefig(PLOT_RESULTS_FILE)
print(f"Saved prediction plots to {PLOT_RESULTS_FILE}")

print("\n✅ Process complete.")
