# Script for plotting HYSPLIT trajectories and sea ice extent for the Arctic
# Shaddy Ahmed
# 01/02/2022

# Load modules
import pandas as pd
import numpy as np
import matplotlib.colors
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
from grid_area import area_grid
import glob
import cartopy.io.shapereader as shpreader
from cartopy.feature import ShapelyFeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

# Get path to csv trajectory files and sort in sequential order
path = "./*.csv"
all_files = glob.glob(path)
all_files.sort(reverse=True)

# Read monthly averaged sea ice extent shape files
seaice_june = shpreader.Reader("./extent_N_201806_polygon_v3.0/extent_N_201806_polygon_v3.0.shp")
seaice_july = shpreader.Reader("./extent_N_201807_polygon_v3.0/extent_N_201807_polygon_v3.0.shp")
seaice_august = shpreader.Reader("./extent_N_201808_polygon_v3.0/extent_N_201808_polygon_v3.0.shp")

seaice_files = [seaice_june, seaice_july, seaice_august]

# File count (Used only to help with plotting the correct file)
count=0
# Names for the titles of the plots
months = ['June 2018', 'July 2018', 'August 2018']

# Set figure information
fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(15, 10),subplot_kw={'projection': ccrs.NorthPolarStereo()})
axs=axs.flatten()

# Loop over trajectory files
for fname in all_files:
    print("Reading "+fname[2:] )
    trajectory_file = pd.read_csv(fname, delimiter='\t')
# Get lat, lon, trajectory height, and PBL height columns in the trajectory file
    file_lat = trajectory_file.lat
    file_lon = trajectory_file.long
    file_height = trajectory_file.height
    file_mixdepth = trajectory_file.MIXDEPTH
# Create a 2D lat/lon grid
    density = np.zeros((360, 180))
# Create a 1D array grid for lat and lon
    latitude = np.arange(-90, 90, 1)
    longitude = np.arange(-180, 180, 1)
# Calculate the grid area in km2
    area_km2 = area_grid(latitude, longitude)
# Loop over each line in the trajectory file
    for i in range(len(trajectory_file)):
# Filter for trajectories below the PBL height
        if file_height[i] <= file_mixdepth[i]:
# Get lat and lon of each line in trajectory file
            traj_lat = file_lat[i]
            traj_lon = file_lon[i]
# Find the closest lat and lon indices in the 1D arrays
# Calculate the difference array to find the closest matching lat and lon
            lat_difference_array = np.absolute(latitude - traj_lat)
            lon_difference_array = np.absolute(longitude - traj_lon)
# find the index of closest lat and lon from the array
            lat_index = lat_difference_array.argmin()
            lon_index = lon_difference_array.argmin()
# Add one to the count at these indices in the density array
            density[lon_index, lat_index] = density[lon_index, lat_index] + 1
# Calculate the weighted density (in h/km2)
    density = np.transpose(density)
    weighted_density = density / area_km2

# -------------------------- PLOTTING --------------------------
    print("Plotting file "+fname[2:])
# Index for the trajectory plots
    count_traj=3+count
# Zoom into the map
    axs[count].set_extent([-180, 180, 90, 66], ccrs.PlateCarree())
    axs[count_traj].set_extent([-180, 180, 90, 66], ccrs.PlateCarree())
# Compute a circle in axes coordinates
    theta = np.linspace(0, 2 * np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(verts * radius + center)
    axs[count].set_boundary(circle, transform=axs[count].transAxes)
    axs[count_traj].set_boundary(circle, transform=axs[count_traj].transAxes)
# Add contour plot
    cp = axs[count_traj].pcolormesh(longitude, latitude, weighted_density, cmap='inferno_r', shading='auto',
                        norm=matplotlib.colors.LogNorm(), rasterized=True, transform=ccrs.PlateCarree())
    cp.set_clim(1e-3, 1e-1)
# Add station markers
    zeppelin = axs[count].plot(11.89, 78.9, marker="o", color='lime', transform=ccrs.PlateCarree())
    villum = axs[count].plot(-16.4, 81.36, marker="o", color='lime', transform=ccrs.PlateCarree())
    alert = axs[count].plot(-62.3, 82.5, marker="o", color='lime', transform=ccrs.PlateCarree())
    zeppelin = axs[count_traj].plot(11.89, 78.9, marker="o", color='lime', transform=ccrs.PlateCarree())
    villum = axs[count_traj].plot(-16.4, 81.36, marker="o", color='lime', transform=ccrs.PlateCarree())
    alert = axs[count_traj].plot(-62.3, 82.5, marker="o", color='lime', transform=ccrs.PlateCarree())
# Add subplot title
    axs[count].set_title(months[count], fontsize=16)
# Add sea ice extent
    shape_feature = ShapelyFeature(seaice_files[count].geometries(), ccrs.epsg(3411), facecolor="darkgray", edgecolor='blue', lw=2)
    axs[count].add_feature(shape_feature)
    shape_feature_traj = ShapelyFeature(seaice_files[count].geometries(), ccrs.epsg(3411), facecolor="None",
                                   edgecolor='blue', lw=2)
    axs[count_traj].add_feature(shape_feature_traj)
# Add map features
    gl=axs[count].gridlines(color='gray', alpha=0.25)
    axs[count].coastlines(resolution='50m')
    gl2=axs[count_traj].gridlines(color='gray', alpha=0.25)
    axs[count_traj].coastlines(resolution='50m')
# Add lat lines and labels
    gl.ylocator = mticker.FixedLocator([65, 70, 75, 80, 85])
    gl.yformatter = LATITUDE_FORMATTER
    gl2.ylocator = mticker.FixedLocator([65, 70, 75, 80, 85])
    gl2.yformatter = LATITUDE_FORMATTER
    axs[count].text(0.5, 0.605, r'85$^\circ$N', fontsize=8, fontweight='bold', transform=axs[count].transAxes)
    axs[count].text(0.5, 0.71, r'80$^\circ$N', fontsize=8, fontweight='bold', transform=axs[count].transAxes)
    axs[count].text(0.5, 0.815, r'75$^\circ$N', fontsize=8, fontweight='bold', transform=axs[count].transAxes)
    axs[count].text(0.5, 0.935, r'70$^\circ$N', fontsize=8, fontweight='bold', transform=axs[count].transAxes)
# Increase file count
    count=count+1
# Clear variables from memory
    del trajectory_file
    del density
    del latitude
    del longitude
    del area_km2
    del shape_feature

# Adjust the location of the subplots on the page to make room for the colorbar
fig.subplots_adjust(bottom=0.2, top=0.9, left=0.1, right=0.9,
                    wspace=0.02, hspace=0.02)

# Add a colorbar axis at the bottom of the graph (x pos, y pos, length, thickness)
cbar_ax = fig.add_axes([0.2, 0.1, 0.6, 0.04])

# Draw the colorbar
cbar = fig.colorbar(cp, cax=cbar_ax,orientation='horizontal')

# Add colorbar units and set label size
cbar.ax.set_title(label='h/km$^2$', fontsize=16, x=1.1, y=-0.2)
cbar.ax.tick_params(labelsize=14)

# Add subheading text for each row of plots
seaice_text = fig.text(0.027,0.735, 'Sea Ice', fontsize=16)
seaice2_text = fig.text(0.03,0.705, 'Extent', fontsize=16)
traj_text = fig.text(0.028,0.38, 'HYSPLIT', fontsize=16)
traj2_text = fig.text(0.02,0.35, 'trajectories', fontsize=16)

#plt.show()

# Save figure
plt.savefig('./trajectory_plots.pdf', dpi=300)
