import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import imread
from mpl_toolkits.mplot3d import Axes3D
import scipy.ndimage as ndimage
from plotly.subplots import make_subplots
import plotly.graph_objs as go

def process_image_data(image_data):
    # Convert the image data from base64 string to NumPy array
    nparr = np.fromstring(image_data, np.uint8)
    
    # Read the image from the NumPy array
    mat = plt.imread(nparr)
    
    # Convert to grayscale (assuming the image is RGB)
    if len(mat.shape) > 2:
        mat = mat[:, :, 0]  # get the first channel
    
    # Get dimensions of the image
    rows, cols = mat.shape
    
    # Create meshgrid for plotting
    xv, yv = np.meshgrid(range(cols), range(rows)[::-1])

    # Apply Gaussian filter to the image
    blurred = ndimage.gaussian_filter(mat, sigma=(5, 5), order=0)

    # Create a figure for visualization
    fig = plt.figure(figsize=(12, 8))

    # Plot original image
    ax1 = fig.add_subplot(221)
    ax1.imshow(mat, cmap='gray')
    ax1.set_title('Original Image')
    ax1.set_xlabel('Columns')
    ax1.set_ylabel('Rows')

    # Plot 3D surface of original image
    ax2 = fig.add_subplot(222, projection='3d')
    ax2.elev = 75
    ax2.plot_surface(xv, yv, mat, cmap='viridis')
    ax2.set_title('3D Plot of Original Image')
    ax2.set_xlabel('Columns')
    ax2.set_ylabel('Rows')
    ax2.set_zlabel('Intensity')

    # Plot blurred image
    ax3 = fig.add_subplot(223)
    ax3.imshow(blurred, cmap='gray')
    ax3.set_title('Blurred Image')
    ax3.set_xlabel('Columns')
    ax3.set_ylabel('Rows')

    # Plot 3D surface of blurred image
    ax4 = fig.add_subplot(224, projection='3d')
    ax4.elev = 75
    ax4.plot_surface(xv, yv, blurred, cmap='viridis')
    ax4.set_title('3D Plot of Blurred Image')
    ax4.set_xlabel('Columns')
    ax4.set_ylabel('Rows')
    ax4.set_zlabel('Intensity')

    # Display the plot
    plt.tight_layout()
    plt.show()

    # Optionally, return the figure for embedding in Flask app
    return fig
