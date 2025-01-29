# Blocking_movement

This repository provides the necessary scripts to reproduce the results and figures from paper the paper on blocking movement. It contains the following folders and content:

a. Python Scripts 

    1. BlockingIndex.py : Calculates the blocking index and blocking intensity from Z500 data

    2. Celltrack_figures.py: Uses the output of the 2D celltracking algorithm and converts the .txt data to usable .nc data. It also filters out the blocks lasting shorter dan 4 days or smaller than our blocking requirements. It results in Figures 4, 5, 6, 9 and S6, S8, S12, S13

    3. Ensemble_figures.py: Makes all ensemble figures. These are Figures 2, S2, S3, S4 and S7

    4. Hovmuller_ERA.py: Results in Figure 3b

    5. Hovmuller.py: Results in Figure S5

    6. Temperature.py: Creates temperature quadrants for each block and saves them into NetCDF files. It results in Figures 7, 8, S1, S9, S10, S11.
    
    7. z500_block.py: Creates and saves the z500 field around blocks to NetCDF files
