import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import mean_squared_error, r2_score, mean_absolute_error

def plot_single_comparison(a_col, b_col, error_margin=0.05, column_name="Column"):
    """
    Plots a single 1D array of predictions (a_col) against ground truth (b_col).
    
    Parameters:
    a_col (numpy.ndarray): 1D array of predicted data, shape (N,).
    b_col (numpy.ndarray): 1D array of ground truth data, shape (N,).
    error_margin (float): The percentage for the error bands (default is 0.05 for 5%).
    column_name (str): Label to use in the title for clarity.
    """
    # Ensure inputs are flat 1D arrays
    a_col = np.asarray(a_col).squeeze()
    b_col = np.asarray(b_col).squeeze()
    
    if a_col.shape != b_col.shape:
        raise ValueError(f"Shape mismatch: a_col is {a_col.shape} but b_col is {b_col.shape}.")

    # --- Calculate Metrics ---
    # 1. RMSE & Normalized RMSE (by range)
    rmse = np.sqrt(mean_squared_error(b_col, a_col))
    range_b = np.max(b_col) - np.min(b_col)
    nrmse_range = (rmse / range_b) * 100
    
    # 2. Mean Relative Error (MRE) in %
    mre = np.mean(np.abs((b_col - a_col) / (b_col + 1e-10))) * 100
    
    # 3. R-squared Score
    r2 = r2_score(b_col, a_col)

    # --- Setup Plot ---
    fig, ax = plt.subplots(figsize=(7, 6))

    # --- Scatter Plot ---
    ax.scatter(b_col, a_col, alpha=0.2, s=3, color='red')

    # --- Error Bands and Perfect Match Line ---
    min_val = np.min(b_col)
    max_val = np.max(b_col)
    
    # Use linspace to guarantee perfectly straight lines regardless of data order
    line_x = np.linspace(min_val, max_val, 100)

    ax.plot(line_x, line_x, color='blue', linestyle='-', label='0% error')
    ax.plot(line_x, line_x * (1 + error_margin), color='black', linestyle=':', label=f'±{error_margin*100:.0f}% Error')
    ax.plot(line_x, line_x * (1 - error_margin), color='black', linestyle=':')

    # --- Formatting ---
    ax.set_xlabel('True Values')
    ax.set_ylabel('Predicted Values')
    ax.set_title(f'Evaluation on {column_name}', fontsize=14)

    # --- Text Box for Metrics ---
    textstr = f'NRMSE: {nrmse_range:.2f}%\nMRE: {mre:.2f}%\n$R^2$: {r2:.3f}'
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    
    ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=12,
            verticalalignment='top', bbox=props)
    
    # Place legend and grid
    ax.legend(loc='lower right')
    ax.grid(True, linestyle='--', alpha=0.6)
    
    plt.tight_layout()
    
    # Save the figure and cleanly close it so it doesn't eat up memory in loops
    plt.savefig(column_name + '.png', dpi=300)
    plt.close(fig)


pred = np.load('preddata_51.npy')
truth = np.load('truedata_51.npy')
colname = ['Cycle-Length_193rl_51','CBC_193rl_51','Fq_193rl_51','FdelH_193rl_51']
print(pred.shape)
N = pred.shape[0]
import random
idx = random.sample(range(N), 1000)
for col in range(pred.shape[1]):
    plot_single_comparison(pred[idx,col], truth[idx,col], error_margin=0.05, column_name=colname[col])
    