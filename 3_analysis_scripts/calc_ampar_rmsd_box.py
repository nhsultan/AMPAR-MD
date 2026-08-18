#!/usr/bin/env python3
"""
calc_ampar_rmsd_box.py

Description:
This script calculates the average box dimensions and the Root Mean Square Deviation 
(RMSD) of the AMPA receptor backbone over the course of a molecular dynamics trajectory.

Key Features:
1. Extracts the X, Y, and Z unit cell lengths and computes their trajectory averages.
2. Calculates the RMSD of the protein backbone atoms relative to the first frame.
3. Converts all native nanometer (nm) measurements from MDTraj to Angstroms (Å).

Outputs:
1. BoxSummary.[PREFIX].txt - Contains average box dimensions and average RMSD.
2. RMSD.[PREFIX].txt - A time-series text file containing the RMSD value at every frame.
"""

import argparse
import numpy as np
import mdtraj as md
from timeit import default_timer as timer

def main():
    start_time = timer()

    # 1. Setup Argument Parsing
    parser = argparse.ArgumentParser(description="Calculate Box Dimensions and AMPAR Backbone RMSD")
    parser.add_argument("--traj", required=True, help="Path to the trajectory file (.nc, .xtc)")
    parser.add_argument("--top", required=True, help="Path to the topology file (.psf, .pdb)")
    parser.add_argument("--out_prefix", required=True, help="Prefix for output files (e.g., Q586)")
    args = parser.parse_args()

    # Align output filenames with the project's capitalization and dot-separated prefix standard
    summary_output = f"BoxSummary.{args.out_prefix}.txt"
    rmsd_output = f"RMSD.{args.out_prefix}.txt"

    print(f"Loading trajectory: {args.traj}")
    print(f"Loading topology: {args.top}")
    
    try:
        traj = md.load(args.traj, top=args.top)
    except Exception as e:
        print(f"Error loading trajectory or topology: {e}")
        return

    print(f"Loaded {traj.n_frames} frames.")

    # MDTraj uses nanometers; define conversion factor to Angstroms
    CONVERSION_FACTOR = 10.0

    # ==========================================
    # 2. Box Dimension Analysis
    # ==========================================
    print("Calculating box dimensions...")
    
    # traj.unitcell_lengths returns an array of shape (n_frames, 3) in nm
    if traj.unitcell_lengths is not None:
        box_lengths_nm = traj.unitcell_lengths
        avg_box_lengths_nm = np.mean(box_lengths_nm, axis=0)
        
        # Convert to Angstroms
        avg_x_A = avg_box_lengths_nm[0] * CONVERSION_FACTOR
        avg_y_A = avg_box_lengths_nm[1] * CONVERSION_FACTOR
        avg_z_A = avg_box_lengths_nm[2] * CONVERSION_FACTOR
        overall_avg_length_A = np.mean([avg_x_A, avg_y_A, avg_z_A])
    else:
        print("Warning: No unitcell (box) information found in this trajectory.")
        avg_x_A = avg_y_A = avg_z_A = overall_avg_length_A = 0.0

    # ==========================================
    # 3. Protein Backbone RMSD Analysis
    # ==========================================
    print("Calculating protein backbone RMSD relative to Frame 0...")
    
    # Focusing on the backbone removes noise from highly flexible sidechains
    backbone_selection = traj.topology.select('protein and backbone')

    if len(backbone_selection) == 0:
        print("Warning: No atoms matched the 'protein and backbone' selection. Skipping RMSD.")
        avg_rmsd_A = 0.0
    else:
        # md.rmsd automatically superposes the selected atoms before calculating
        rmsd_per_frame_nm = md.rmsd(traj, traj, 0, atom_indices=backbone_selection)
        rmsd_per_frame_A = rmsd_per_frame_nm * CONVERSION_FACTOR
        avg_rmsd_A = np.mean(rmsd_per_frame_A)

        # Save RMSD Time Series
        frame_indices = np.arange(traj.n_frames)
        rmsd_data_to_save = np.column_stack((frame_indices, rmsd_per_frame_A))
        
        header = "Frame_Index    RMSD_Backbone_(Angstroms)"
        np.savetxt(rmsd_output, rmsd_data_to_save, fmt='%d %.6f', header=header, comments='')
        print(f"RMSD time-series successfully written to {rmsd_output}")

    # ==========================================
    # 4. Save Summary Output
    # ==========================================
    with open(summary_output, 'w') as f:
        f.write("--- AMPAR Molecular Dynamics Analysis Summary ---\n")
        f.write(f"Topology File: {args.top}\n")
        f.write(f"Trajectory File: {args.traj}\n")
        f.write(f"Total Frames Processed: {traj.n_frames}\n\n")
        
        f.write("--- Average Box Lengths (Å) ---\n")
        f.write(f"X Dimension: {avg_x_A:.4f}\n")
        f.write(f"Y Dimension: {avg_y_A:.4f}\n")
        f.write(f"Z Dimension: {avg_z_A:.4f}\n")
        f.write(f"Overall Average: {overall_avg_length_A:.4f}\n\n")

        f.write("--- Protein RMSD Analysis ---\n")
        f.write(f"Reference: Frame 0\n")
        f.write(f"Selection: Protein Backbone\n")
        f.write(f"Average Backbone RMSD: {avg_rmsd_A:.4f} Å\n")

    print(f"Summary results successfully written to {summary_output}")

    end_time = timer()
    print(f"RMSD and Box analysis completed in {end_time - start_time:.2f} seconds.")

if __name__ == "__main__":
    main()
