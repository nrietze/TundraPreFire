import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import contextily as ctx
import geopandas as gpd
from glob import glob
import os
import datetime
import pandas as pd
import numpy as np 

# Set up function
def animate_stacked_polygons_with_basemap(df_viirs_filtered, 
                                          test_id: str, 
                                          save_path: str = None):
    """
    Animate polygons from filepaths, stacking them frame by frame on a dark basemap.

    Parameters:
    - df_viirs_filtered: DataFrame with a 'filepath' column containing polygon files
    - test_id: ID used to filter polygons for animation
    - save_path: Optional path to save the animation (e.g., 'stacked_polygons.mp4')
    """
    # Initialize figure
    fig, ax = plt.subplots(figsize=(10, 8))

    # Container to hold cumulative polygons
    polygons = []

    # Load and filter polygons for the given TEST_ID
    for filepath in df_viirs_filtered['filepath']:
        sd_perimeters = gpd.read_file(filepath)
        sd_perimeter_filtered = sd_perimeters[sd_perimeters.fireid == test_id]
        if not sd_perimeter_filtered.empty:
            polygons.append(sd_perimeter_filtered)

    # Ensure we have polygons to animate
    if not polygons:
        print("No polygons found for the given test ID.")
        return

    # Combine all polygons to compute global extent
    combined_gdf = gpd.GeoDataFrame(pd.concat(polygons, ignore_index=True),
                                    crs=polygons[0].crs)
    bounds = combined_gdf.total_bounds  # [minx, miny, maxx, maxy]

    # Color map setup
    cmap = plt.get_cmap('magma')
    num_frames = len(polygons)

    def update(frame):
        ax.clear()
        
        # Set map extent and plot basemap
        ax.set_xlim(bounds[0], bounds[2])
        ax.set_ylim(bounds[1], bounds[3])
        
        # Add dark background map from contextily
        ctx.add_basemap(ax, source=ctx.providers.CartoDB.Positron,
                        crs=polygons[0].crs.to_string())
        
        # Plot polygons up to the current frame with cumulative stacking
        for i in range(frame + 1):
            polygons[i].plot(ax=ax, edgecolor=cmap(i / num_frames), facecolor='none', linewidth=1.5)

        ax.set_title(f"Stacking Polygons - Frame: {frame + 1}")

    # Create the animation
    anim = FuncAnimation(fig, update, frames=num_frames, repeat=False)

    # Save or display the animation
    if save_path:
        anim.save(save_path, writer='ffmpeg', fps=1)
    else:
        plt.show()
        
def extract_datetime(PATH):
    FN = os.path.basename(PATH)
    datestr = FN.split(".")[0]
    return datetime.datetime.strptime(datestr, "%Y%m%d%p")
# %%

test_id = 14211

# Load list of sub-daily perimeters and assign datetime
PATH_VIIRS_PERIMETERS = "../data/feature_layers/fire_atlas/"

# Load VIIRS perimeters in Siberian tundra
FN_VIIRS_CAVM_PERIMETERS = os.path.join(PATH_VIIRS_PERIMETERS,"viirs_perimeters_in_cavm_e113.gpkg")

if os.path.exists(FN_VIIRS_CAVM_PERIMETERS):
    print("CAVM Fire perimeter file exists, loading.")
    merged_fire_perimeters = gpd.read_file(FN_VIIRS_CAVM_PERIMETERS)
else:
    print("Please extract VIIRS fire perimeters in CAVM extent first.")
    pass

FLIST_VIIRS_SUBDAILY = glob(os.path.join(PATH_VIIRS_PERIMETERS,
                                         "sub_daily/*/Snapshot/*M.gpkg"))
df_viirs_sub_daily = pd.DataFrame(data={"filepath":FLIST_VIIRS_SUBDAILY})
df_viirs_sub_daily['datetime'] = [extract_datetime(FN) for FN in FLIST_VIIRS_SUBDAILY]

df_viirs_sub_daily = df_viirs_sub_daily.sort_values("datetime")

perimeter = merged_fire_perimeters.loc[merged_fire_perimeters.fireid==test_id]
start_date = datetime.datetime(perimeter.tst_year.item(),
                               perimeter.tst_month.item(),
                               perimeter.tst_day.item())

end_date = datetime.datetime(perimeter.ted_year.item(),
                             perimeter.ted_month.item(),
                             perimeter.ted_day.item())

# filter sub-daily perimeters that are within the fire's start and end dates
cond = np.logical_and(df_viirs_sub_daily.datetime >= start_date,
                      df_viirs_sub_daily.datetime <= end_date)
df_viirs_filtered = df_viirs_sub_daily[cond]

animate_stacked_polygons_with_basemap(df_viirs_filtered,
                                      test_id,
                                      save_path='figures/polygon_animation.mp4')