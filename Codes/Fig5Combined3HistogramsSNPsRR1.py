import matplotlib.pyplot as plt
import matplotlib.image as mpimg

# 1. Define the file names of the previously saved plot images
file_overall = "/content/drive/MyDrive/mean_abs_effect_histogram_top_vs_random_cap20.png"
file_positive = "/content/drive/MyDrive/pos_snps_absolute_effect_histogram.png"
file_negative = "/content/drive/MyDrive/neg_snps_absolute_effect_histogram.png"

# 2. Setup the layout (3 rows, 1 column) for vertical orientation
# We use a taller figsize (e.g., width=10, height=24) to accommodate the stack
fig, axes = plt.subplots(3, 1, figsize=(10, 24))

# 3. Load and display each image in the respective subplot
image_paths = [file_overall, file_positive, file_negative]

for i, path in enumerate(image_paths):
    try:
        img = mpimg.imread(path)
        axes[i].imshow(img)
        axes[i].axis('off') # Turn off subplot axes to preserve original plot appearance
    except FileNotFoundError:
        print(f"Warning: {path} not found.")
        axes[i].text(0.5, 0.5, f"Missing Image:\n{path}", ha='center', va='center')
        axes[i].axis('off')

# 4. Final adjustments and saving
plt.subplots_adjust(hspace=0.05) # Minimize vertical space between panels
plt.tight_layout()
plt.savefig("/content/drive/MyDrive/histogram_3combined_vertical_layout.png", dpi=300, bbox_inches='tight')
plt.show()
