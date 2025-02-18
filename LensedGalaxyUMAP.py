# JUST LENSED GALAXY
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import umap
from astropy.io import fits
from sklearn.preprocessing import StandardScaler
from matplotlib.colors import LogNorm, PowerNorm
# List of FITS files
fits_files = [
    'LENSED_f480m_cut.fits', 'LENSED_f444w_cut_reproj.fits', 'LENSED_f356w_cut_reproj.fits',
    'LENSED_f277w_cut_reproj.fits', 'LENSED_f200w_cut_reproj.fits', 'LENSED_f150w_cut_reproj.fits'
]

# File path directory
base_path = "/Users/alr/Desktop/LEGGOS/3DUMAP/justlens/"

# Initialize list to store the 3D cube (after stacking) and the flux intensities
cube_data = []

# Process each FITS file
for i, file in enumerate(fits_files):
    file_path = base_path + file

    with fits.open(file_path) as hdul:
        print(f"Processing {file}...")
        fdata = hdul['SCI'].data  # Extract image data
    
    # Replace NaNs and standardize
    fdata_fixed = np.nan_to_num(fdata)
    #min_flux = np.nanmax(fdata)  # Use np.nanmin to ignore NaNs
    flux_threshold = -1
    filtered_data = np.where(fdata_fixed > flux_threshold, fdata_fixed, -1)  #mask to remove low flux
    cube_data.append(filtered_data)  # Stack this 2D array into a 3D cube later
   
# Stack all 6 FITS images into a single 3D cube (shape: height x width x 6)
cube_data = np.stack(cube_data, axis=-1)

# Flatten the 3D cube into a 2D array (pixels x filters)
n_pixels = cube_data.shape[0] * cube_data.shape[1]  # Total number of pixels
flattened_data = cube_data.reshape(n_pixels, -1)  # Flatten into 2D array: (pixels x filters)
flux_intensities = flattened_data.mean(axis=1) 

# sample only a subset of pixels to speed up UMAP
sample_size = int(.1 * n_pixels)  # Take 10% of the pixels
indices = np.random.choice(n_pixels, sample_size, replace=False)
flattened_data_sampled = flattened_data[indices]
flux_intensities_sampled = flux_intensities[indices]
# Standardize the data
#scaled_data = StandardScaler().fit_transform(flattened_data_sampled)

# Apply UMAP to reduce dimensions from 6 to 3 (for 3D plotting)
reducer = umap.UMAP(n_neighbors=30, min_dist=1, n_components=3, random_state=42, metric='euclidean', n_jobs=-1)  # n_jobs=-1 uses all available CPU cores
embedding_3d = reducer.fit_transform(flattened_data_sampled) #replace with scaled_data

# ========== 3D UMAP Projection ==========
fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')

# Plot UMAP in 3D
scatter = ax.scatter(
    embedding_3d[:, 0], embedding_3d[:, 1], embedding_3d[:, 2],
    c=flux_intensities_sampled, cmap="plasma", s=10, norm=PowerNorm(gamma=.5)
)
ax.set_title("3D UMAP Projection of Combined FITS Cube")
ax.set_xlabel("Dim 1")
ax.set_ylabel("Dim 2")
ax.set_zlabel("Dim 3")
#plt.colorbar(flux_intensities_sampled, title = "Mean Flux Intensity")
plt.show()

fits_file = base_path + fits_files[0]

with fits.open(fits_file) as hdul:
    fdata = hdul['SCI'].data

# Reshape to match the original image shape
n_y, n_x = fdata.shape

# Normalize UMAP data to fit the image grid
#umap_img = np.zeros((n_y, n_x))

# Map UMAP points to image pixels
for i, idx in enumerate(indices):
    y, x = np.unravel_index(idx, (n_y, n_x))  # Get corresponding pixel coordinates
    umap_img[y, x] = flux_intensities_sampled[i]  # Assign flux intensities to corresponding pixels

# Plot the original FITS image with the UMAP overlay
fig, ax = plt.subplots(figsize=(10, 8))

# Display the FITS image (gray scale)
ax.imshow(fdata, cmap='gray', origin='lower', interpolation='nearest', alpha=.5)

# Overlay UMAP as a heatmap
ax.imshow(umap_img, cmap='plasma', origin='lower', interpolation='nearest', alpha=1)

# Add a colorbar to show flux intensity
#cbar = plt.colorbar(ax.imshow(umap_img, cmap='viridis', origin='lower', interpolation='nearest', alpha=0.5))
#cbar.set_label("Flux Intensity (from UMAP)")

plt.title("F480 FITS Image with UMAP Heatmap Overlay")
plt.show()
