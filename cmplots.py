import scipy.io as sio
import numpy as np
import matplotlib.pyplot as plt

data = sio.loadmat('total_cm.mat')
acc = sio.loadmat('total_acc.mat')

print(data)

# Create a figure and a grid of subplots
fig, axs = plt.subplots(nrows=3, ncols=6, figsize=(30, 15))

acc_list = [v for k, v in acc.items() if not k.startswith('__')]

# 
titles = ["MAV-VAR: VFvPinch", "MAV-VAR: VFvFlex", "MAV-VAR: VFvRest", "MAV-VAR: PinchvFlex", "MAV-VAR: PinchvRest", "MAV-VAR: FlexvRest", "MAV: VFvPinch", "MAV: VFvFlex", "MAV: VFvRest", "MAV: PinchvFlex", "MAV: PinchvRest", "MAV: FlexvRest", "VAR: VFvPinch", "VAR: VFvFlex", "VAR: VFvRest", "VAR: PinchvFlex", "VAR: PinchvRest", "VAR: FlexvRest"]

# Iterate over the keys in the data dictionary
for idx, key in enumerate([k for k in data.keys() if not k.startswith('__')]):
    # Extract the 2x2 matrix, round the numbers, and convert to integer
    matrix = np.round(data[key]).astype(int)
    matrix = np.flip(matrix)

    # Determine the row and column in the grid for the current subplot
    row = idx // 6
    col = idx % 6

    # Plot the matrix in the current subplot
    cax = axs[row, col].matshow(matrix, cmap='viridis')

    # Add numbers in each block with bigger font size
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            axs[row, col].text(j, i, matrix[i, j], ha='center', va='center', color='w', fontsize=20)

    # Add a title to the current subplot
    #axs[row, col].set_title(key)
    #axs[row, col].set_title(f"{key} (Accuracy: {acc[key]:.2f})")
    accuracy = acc_list[idx][0][0]

    title = titles[idx]


    # Add a title to the current subplot, including the accuracy data
    axs[row, col].set_title(f"{title} (Accuracy: {accuracy:.2f})", fontsize = 15)

# Add a colorbar to the figure
#fig.colorbar(cax, ax=axs.ravel().tolist(), orientation='horizontal')

# Save the figure as a PNG file
plt.savefig("all_matrices2.png")
# plt.close()
plt.show()