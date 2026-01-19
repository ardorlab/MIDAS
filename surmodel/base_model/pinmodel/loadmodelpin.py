
## train pin with deeper unet 
import h5py
import numpy as np
import joblib
import os
import matplotlib.pyplot as plt
import random
# --- DeepXDE Backend Imports ---
from deepxde.backend import tf

import deepxde as dde
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
from sklearn.preprocessing import MinMaxScaler
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score

IMG_DIM = 153 # 153x153
INPUT_SHAPE = (IMG_DIM, IMG_DIM, 1)
BATCH_SIZE = 16
EPOCHS = 50 
VALIDATION_SPLIT = 0.2
## models and scaler 
SAVED_MODEL_NAME = 'pin_power_unet_deep_noscale.h5' 
X_SCALER_FILE = 'x_scaler.gz'
Y_SCALER_FILE = 'y_scaler.gz'
PLOT_RESULTS_FILE ='testpred.png'
pathmodel ='/home/nguykhan/DeepOnet/PARCS-data/pinrescontruction/'

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
    """
    This function returns a loss function that is a weighted sum of
    masked MSE (for overall agreement) and a peak error (for max value).
    
    alpha = 0.0 -> 100% MSE
    alpha = 1.0 -> 100% Peak Error
    alpha = 0.5 -> 50% MSE, 50% Peak Error
    """
    
    def masked_mse_loss(y_true, y_pred):
        # 1. Masked Mean Squared Error (for overall agreement)
        mask = tf.cast(tf.not_equal(y_true, 0), tf.float32)
        squared_error = tf.square(y_true - y_pred)
        masked_squared_error = squared_error * mask
        # Add epsilon for numerical stability
        loss = tf.reduce_sum(masked_squared_error) / (tf.reduce_sum(mask) + 1e-7)
        return loss

    def peak_loss(y_true, y_pred):
        # 2. Peak Value Error (for max value)
        mask = tf.cast(tf.not_equal(y_true, 0), tf.float32)
        
        # Find the max of the true data (shape: [batch_size])
        max_true = tf.reduce_max(y_true, axis=[1, 2, 3])
        
        # Find the max of the predicted data *after* masking out padded areas
        max_pred = tf.reduce_max(y_pred * mask, axis=[1, 2, 3])
        
        # Calculate the mean squared error of the peak values
        peak_loss = tf.reduce_mean(tf.square(max_true - max_pred))
        return peak_loss

    def composite_loss(y_true, y_pred):
        # Combine the two losses
        mse = masked_mse_loss(y_true, y_pred)
        peak = peak_loss(y_true, y_pred)
        
        return (1.0 - alpha) * mse + alpha * peak
    
    return composite_loss

composite_loss_fn = create_composite_loss(alpha=0.1)
model = build_deep_unet(INPUT_SHAPE)
optimizer = Adam(learning_rate=0.0005) 
model.compile(optimizer=optimizer, loss=composite_loss_fn, metrics=['mean_absolute_error'])
#model.summary()


# Load the best model
with custom_object_scope({'composite_loss': composite_loss_fn}):
    best_model = tf.keras.models.load_model(pathmodel+SAVED_MODEL_NAME)



# # Get predictions (These are already in original scale)
# y_pred_raw_cnn = best_model.predict(X_test_cnn)

# # --- CORRECTION STEP ---
# # Create a mask of "pixels to keep" (non-zero) from y_test (which is already unscaled)
# keep_mask = (coremap[0] != 0) 
# keep_mask_4d = keep_mask.reshape(1, IMG_DIM, IMG_DIM, 1) 

# # Apply the mask to the predictions
# y_pred_corrected_cnn = y_pred_raw_cnn * keep_mask_4d

# # Flatten for metrics
# y_test_flat = coremap.flatten()
# y_pred_corrected_flat = y_pred_corrected_cnn.flatten()

# # Calculate metrics
# mse = mean_squared_error(y_test_flat, y_pred_corrected_flat)
# mae = mean_absolute_error(y_test_flat, y_pred_corrected_flat)
# r2 = r2_score(y_test_flat, y_pred_corrected_flat)

# print("\n--- Model Evaluation Metrics (Post-Processed, No Scale) ---")
# print(f"  Mean Squared Error (MSE): {mse:.6e}")
# print(f"  Mean Absolute Error (MAE): {mae:.6f}")
# print(f"  R-squared (R²) Score:     {r2:.6f}")
# print("-------------------------------------------------\n")

# # --- 7. Visualize Results (with Difference Plot) ---
# print("Generating and saving prediction plots...")

# # Reshape for plotting
# y_pred_corrected_img = y_pred_corrected_cnn.reshape(-1, IMG_DIM, IMG_DIM)
# X_test_img = nofluxmap.reshape(-1, IMG_DIM, IMG_DIM)

# ## print Fq
# print(f'Fq - True: {coremap.max():.4f}, Pred: {y_pred_corrected_img.max():.4f}, Old {X_test_img.max():.4f}')

# # Plot 3 random examples with 4 columns
# num_examples = 3
# indices = np.random.choice(range(len(X_test_cnn)), num_examples, replace=False)

# fig, axes = plt.subplots(num_examples, 4, figsize=(20, 5 * num_examples))
# fig.suptitle("DEEP U-Net Prediction (No Scaling) vs. Ground Truth", fontsize=16)

# for i, idx in enumerate(indices):
#     input_img = X_test_img[idx]
#     true_img = coremap[idx]
#     pred_img = y_pred_corrected_img[idx]
    
#     # Calculate Difference (Prediction - Truth)
#     diff_img = pred_img - true_img

#     # Determine consistent min/max for the main plots
#     vmin = min(true_img.min(), pred_img.min())
#     vmax = max(true_img.max(), pred_img.max())
    
#     print(f'Sample {idx} - Max True: {true_img.max():.4f}, Max Pred: {pred_img.max():.4f}, Max old {input_img.max():.4f}')

#     # 1. Plot Input
#     ax = axes[i, 0]
#     im0 = ax.imshow(input_img, cmap='viridis', origin='lower')
#     ax.set_title(f"Input (Sample {idx})")
#     fig.colorbar(im0, ax=ax, label='Pin Power')

#     # 2. Plot Ground Truth
#     ax = axes[i, 1]
#     im1 = ax.imshow(true_img, cmap='viridis', origin='lower')#, vmin=vmin, vmax=vmax)
#     ax.set_title(f"True")
#     fig.colorbar(im1, ax=ax, label='Pin Power')

#     # 3. Plot Predicted
#     ax = axes[i, 2]
#     im2 = ax.imshow(pred_img, cmap='viridis', origin='lower')#, vmin=vmin, vmax=vmax)
#     ax.set_title(f"Predicted (Corrected)")
#     fig.colorbar(im2, ax=ax, label='Pin Power')
    
#     # 4. Plot Difference (Pred - True)
#     ax = axes[i, 3]
#     # Center the colormap at 0 to show + (red) and - (blue) errors
#     diff_limit = max(abs(diff_img.min()), abs(diff_img.max()))
#     im3 = ax.imshow(diff_img, cmap='bwr', origin='lower', vmin=-diff_limit, vmax=diff_limit)
#     ax.set_title(f"Difference (Pred - True)")
#     fig.colorbar(im3, ax=ax, label='Error')

# plt.tight_layout(rect=[0, 0.03, 1, 0.95])
# plt.savefig(PLOT_RESULTS_FILE)
# print(f"Saved prediction plots to {PLOT_RESULTS_FILE}")
