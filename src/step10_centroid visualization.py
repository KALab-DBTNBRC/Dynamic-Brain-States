import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat
import os

# Load centroids from MATLAB file
data = loadmat('/content/C patterns Rep set.mat')
C = data['C']  # shape = (4, 64)

print("Centroids shape:", C.shape)

# Create folder to store plots
os.makedirs("/content/centroid_plots_Rep", exist_ok=True)

# Plot and save
for i in range(C.shape[0]):
    fig, ax = plt.subplots(figsize=(12, 4))
    ax.bar(np.arange(1, C.shape[1] + 1), C[i, :])
    ax.set_title(f"Brain State {i+1} - ROI Activation Pattern")
    ax.set_xlabel("ROI Index")
    ax.set_ylabel("Mean Activation")

    # Save the figure (from the figure object, not plt)
    filename = f"/content/centroid_plots_Rep/centroid_state_{i+1}.png"
    fig.savefig(filename, dpi=300, bbox_inches='tight')
    plt.show()     # show it in Colab
    plt.close(fig) # close AFTER showing

print(" Saved all centroid plots successfully!")