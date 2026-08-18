#!/usr/bin/env python3
"""
plot_ampar_1d_profiles.py

Description:
This script reads the 1D profile data (electrostatics and ion density) generated 
by `calc_ampar_1d_profiles.sh` and creates a publication-quality two-panel figure.

It plots:
  - Top Panel: 1D Electrostatic Potential (kT/e) along the Z-axis.
  - Bottom Panel: Ion Density (K+ and Cl-) along the Z-axis.

Outputs:
A high-resolution PNG/PDF figure of the 1D pore profiles.
"""

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

def main():
    # 1. Setup Argument Parsing
    parser = argparse.ArgumentParser(description="Plot AMPAR 1D Electrostatic and Ion Density Profiles")
    parser.add_argument("--files", nargs='+', required=True, 
                        help="List of profile text files to plot (e.g., run8_profile.txt)")
    parser.add_argument("--labels", nargs='+', required=True, 
                        help="List of labels for the legend (must match the number of files)")
    parser.add_argument("--out_file", default="1D_Profiles_Combined.png", 
                        help="Output filename for the plot (default: 1D_Profiles_Combined.png)")
    args = parser.parse_args()

    if len(args.files) != len(args.labels):
        print("Error: The number of --files must exactly match the number of --labels.")
        return

    print(f"Preparing to plot 1D profiles for {len(args.files)} systems...")

    # Set standard publication-quality aesthetics
    sns.set_theme(style="ticks", context="paper", font_scale=1.2)
    
    # Create a 2-panel plot sharing the X (Z-coordinate) axis
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8), sharex=True)
    colors = sns.color_palette("husl", len(args.files))

    data_found = False

    # 2. Read and Plot Data
    for idx, (input_file, label) in enumerate(zip(args.files, args.labels)):
        if not os.path.exists(input_file):
            print(f"Warning: Could not find {input_file}. Skipping...")
            continue
            
        print(f"Loading {input_file}...")
        
        try:
            # Read the data using sep=r'\s+' to avoid future pandas warnings
            df = pd.read_csv(input_file, sep=r'\s+', comment='#', 
                             names=['Z_mid', 'Potential', 'POT_density', 'CLA_density'])
            
            # --- Panel 1: Electrostatic Potential ---
            ax1.plot(df['Z_mid'], df['Potential'], label=label, color=colors[idx], linewidth=2.0)
            
            # --- Panel 2: Ion Densities ---
            # Plot Potassium (Solid line)
            ax2.plot(df['Z_mid'], df['POT_density'], label=f"{label} (K+)", 
                     color=colors[idx], linestyle='-', linewidth=1.5)
            # Plot Chloride (Dashed line)
            ax2.plot(df['Z_mid'], df['CLA_density'], label=f"{label} (Cl-)", 
                     color=colors[idx], linestyle='--', linewidth=1.5, alpha=0.7)
            
            data_found = True
            
        except Exception as e:
            print(f"Error processing {input_file}: {e}")

    # 3. Finalize and Save the Plot
    if data_found:
        # Format Top Panel
        ax1.set_title("1D Electrostatic Potential", weight='bold')
        ax1.set_ylabel("Potential (kT/e)", weight='bold')
        ax1.axhline(0, color='black', linewidth=0.8, linestyle='--') # 0 reference line
        ax1.legend(frameon=False, loc="best")
        
        # Format Bottom Panel
        ax2.set_title("Ion Density Distribution", weight='bold')
        ax2.set_xlabel(r"Z-Axis Coordinate ($\AA$)", weight='bold')
        ax2.set_ylabel("Density (counts/frame)", weight='bold')
        ax2.legend(frameon=False, loc="best", ncol=2) # 2 columns to fit both ions cleanly

        # General formatting
        sns.despine()
        plt.tight_layout()

        # Save figure
        print(f"Saving high-resolution plot to {args.out_file}...")
        plt.savefig(args.out_file, dpi=300, bbox_inches="tight")
        print("Plot generation complete!")
    else:
        print("Error: No valid data files were found to plot.")

if __name__ == "__main__":
    main()