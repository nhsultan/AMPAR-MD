#!/usr/bin/env python3
"""
calc_ampar_pore_radius.py

Description:
Calculates the pore radius of AMPA receptors using MDAnalysis and the HOLE2 program.
It analyzes the trajectory frame by frame, calculating the radius profile along the Z-axis.

The script automatically centers the initial search at the center of geometry of 
residue 586 (the Q/R site and main physical bottleneck). It uses the alpha carbons (CA) 
to ensure the script works universally for the Wild Type (Q586) and all mutants (E, G, R).

Outputs:
A single, easy-to-read CSV file containing the Frame, Z-coordinate, and Radius.
"""

import argparse
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import hole2
from timeit import default_timer as timer
import warnings

# Suppress minor warnings from MDAnalysis for cleaner terminal output
warnings.filterwarnings('ignore')

def main():
    start_time = timer()

    # 1. Setup Argument Parsing
    parser = argparse.ArgumentParser(description="Calculate AMPAR pore radius using MDAnalysis and HOLE2")
    parser.add_argument("--traj", required=True, help="Path to the trajectory file (.nc, .xtc)")
    parser.add_argument("--top", required=True, help="Path to the topology file (.psf, .pdb)")
    parser.add_argument("--out_prefix", required=True, help="Prefix for output files (e.g., Q586)")
    parser.add_argument("--executable", default="hole", help="Path to the HOLE executable binary (default: 'hole')")
    args = parser.parse_args()

    print(f"Loading topology: {args.top}")
    print(f"Loading trajectory: {args.traj}")
    
    # 2. Load the Universe
    u = mda.Universe(args.top, args.traj)

    # 3. Define the Pore Center
    # We use the alpha carbons of residue 586 (the selectivity filter) to define 
    # the starting point (cpoint) for the HOLE search algorithm. 
    # Using 'resid 586 and name CA' makes this robust for all Q586 mutations.
    pore_center_sel = u.select_atoms("resid 586 and name CA")
    
    # Calculate the initial center of geometry to seed the search.
    initial_cpoint = pore_center_sel.center_of_geometry()

    print(f"Initial search seeded at center of geometry for residue 586: {initial_cpoint}")
    print("Running HOLE analysis. This may take some time depending on trajectory length...")
    
    # 4. Initialize and Run HOLE analysis
    # cvect=[0,0,1] ensures it searches strictly along the Z-axis
    ha = hole2.HoleAnalysis(u, 
                            executable=args.executable,
                            cpoint=initial_cpoint,
                            cvect=[0, 0, 1])
    
    ha.run()

    # 5. Extract and Save Data
    # Aligning the output filename with the 'IonTracking.PREFIX.txt' convention
    out_file = f"PoreRadius.{args.out_prefix}.csv"
    print(f"Saving results to {out_file}...")
    
    with open(out_file, 'w') as f:
        f.write("Frame,Z_Coordinate,Radius\n")
        
        # ha.results.profiles is a list of record arrays containing the data per frame
        for frame_idx, profile in enumerate(ha.results.profiles):
            z_arr = profile['rxncoord'] # Z-axis coordinates
            rad_arr = profile['radius'] # Corresponding pore radii
            
            for z, r in zip(z_arr, rad_arr):
                f.write(f"{frame_idx},{z:.3f},{r:.3f}\n")

    end_time = timer()
    print(f"Pore radius calculation complete in {end_time - start_time:.2f} seconds!")

if __name__ == "__main__":
    main()
