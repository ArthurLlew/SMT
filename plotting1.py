# Importing the required modules
import numpy as np
import matplotlib.pyplot as plt
  
# Generating data for the heat map
data = np.genfromtxt('results/task_v1_res.csv', delimiter=';')

N = data.shape[0]
M = data.shape[1]

# Init figure
fig = plt.figure(figsize=(7, 7))
    
# Figure title
fig.suptitle('Heat Map', fontsize=20, y=0.92)

# Subplot
subplot = fig.add_subplot()
  
# Ox and Oy axis points
subplot.set_xticks(np.linspace(0, M-1, 11, endpoint=True), ['{:.2f}'.format(i) for i in np.linspace(-3, 3, 11, endpoint=True)], rotation=45)
for label in subplot.get_xticklabels():
    label.set_horizontalalignment('right')
subplot.set_yticks(np.linspace(0, N-1, 11, endpoint=True), ['{:.2f}'.format(i) for i in np.flip(np.linspace(0, 4, 11, endpoint=True))])

# Axes labels
subplot.set_xlabel('X', labelpad=4)
subplot.set_ylabel('Y', rotation=0, labelpad=8)

# Function to show the heat map
axes = subplot.imshow(data, cmap='magma')

# Use colorbar with 2 values
cbar = fig.colorbar(axes, shrink=0.5, ticks=np.linspace(data.min(), data.max(), 3, endpoint=True))
# Color bar tick labels size
cbar.ax.tick_params(axis='both', which='major', labelsize=20)

plt.show()