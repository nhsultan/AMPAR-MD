#!/usr/bin/env python3
"""
calc_ampar_ion_flux.py

Description:
This script calculates the ion flux (permeation events) through the pore of 
AMPA receptors. It processes a molecular dynamics trajectory and tracks the 
movement of specific ions (POT, CLA). 

It defines the channel pore geometry using user-defined residue selections.
By default, it uses:
  - Bottom boundary: ASP 590
  - Top boundary: GLN 586
  - Z-axis reference: Midpoint between Top and Bottom boundaries

The script categorizes permeation events based on crossing the MidZ point:
  - TB: Top to Bottom permeation (Count: 1)
  - BB: Bottom to Bottom permeation (Count: 0)
  - BT: Bottom to Top permeation (Count: -1)
  - TT: Top to Top permeation (Count: 0)
  - MB/MT: Mid to Bottom/Top (Half events at trajectory start)
  - BM/TM: Bottom/Top to Mid (Half events at trajectory end)

Outputs:
Saves three text files containing the tracked ion events, pathway types, and ion identities.
"""

import argparse
import sys
import numpy as np
import mdtraj as md
from timeit import default_timer as timer

def main():
    start_time = timer()

    # 1. Setup Argument Parsing
    parser = argparse.ArgumentParser(description="Calculate AMPA receptor ion flux from MD trajectory")
    parser.add_argument("--traj", required=True, help="Path to the trajectory file (e.g., .nc, .xtc)")
    parser.add_argument("--top", required=True, help="Path to the topology file (e.g., .psf, .pdb)")
    parser.add_argument("--out_prefix", required=True, help="Prefix/path for output files (e.g., Q586)")
    parser.add_argument("--top_res", default="resname GLN and resid 586", help="MDTraj selection string for the top boundary")
    parser.add_argument("--bot_res", default="resname ASP and resid 590", help="MDTraj selection string for the bottom boundary")
    args = parser.parse_args()

    # Define output filenames based on prefix
    event_data_file = f"IonTracking.{args.out_prefix}.txt"
    pathway_file = f"IonPathway.{args.out_prefix}.txt"
    which_residue_file = f"Residue.{args.out_prefix}.txt"

    print(f"Loading trajectory: {args.traj}")
    print(f"Loading topology: {args.top}")
    traj = md.load(args.traj, top=args.top)
    top = traj.topology

    # 2. System Selection and Geometric Boundaries
    print(f"Calculating pore boundaries using:\n  Top: '{args.top_res}'\n  Bottom: '{args.bot_res}'")
    
    protein_idx = top.select('protein')
    avg_x = np.mean(traj.xyz[:, protein_idx, 0])
    avg_y = np.mean(traj.xyz[:, protein_idx, 1])

    bot_residue_idx = top.select(args.bot_res)
    if len(bot_residue_idx) == 0:
        print(f"Error: The bottom selection '{args.bot_res}' matched 0 atoms. Check your topology.", file=sys.stderr)
        sys.exit(1)
    bot_z_avg = np.mean(traj.xyz[:, bot_residue_idx, 2])

    top_residue_idx = top.select(args.top_res)
    if len(top_residue_idx) == 0:
        print(f"Error: The top selection '{args.top_res}' matched 0 atoms. Check your topology.", file=sys.stderr)
        sys.exit(1)
    top_z_avg = np.mean(traj.xyz[:, top_residue_idx, 2])

    # Z-reference point halfway between top and bottom residues
    mid_z = bot_z_avg + 0.5 * (top_z_avg - bot_z_avg)

    # 3. Track Ions
    target_ions = ['POT', 'CLA']
    
    total_events = 0
    pathway_list = []
    ion_number_list = []
    event_start_list = []
    event_end_list = []
    events_count_list = []

    print(f"Tracking ions: {target_ions}...")

    for ion_name in target_ions:
        md_ion_indices = top.select(f'resname {ion_name}')

        for ion_idx in md_ion_indices:
            coords = traj.xyz[:, ion_idx, :]
            
            # Vectorized calculation for radius and Z-axis mapping
            x_norm = coords[:, 0] - avg_x
            y_norm = coords[:, 1] - avg_y
            r_array = np.sqrt(x_norm**2 + y_norm**2)
            z_array = coords[:, 2]

            # Mask for when the ion is strictly inside the defined pore cylinder
            in_pore_mask = (z_array > bot_z_avg) & (z_array < top_z_avg) & (r_array < 2.0)
            
            radius_tracked = np.where(in_pore_mask, r_array, None)
            z_tracked = np.where(in_pore_mask, z_array, None)

            # 4. Identify Continuous Permeation Events
            event_starts = []
            event_ends = []
            
            for i in range(len(z_tracked)):
                if radius_tracked[i] is not None:
                    # Check for event start/end at the very boundaries of the trajectory
                    if i == 0 and radius_tracked[i+1] is not None:
                        event_starts.append(i)
                    if i == len(z_tracked) - 1 and radius_tracked[i-1] is not None:
                        event_ends.append(i)

                    # Check for event start/end during the trajectory
                    if 0 < i < len(z_tracked) - 1:
                        if radius_tracked[i-1] is None and radius_tracked[i+1] is not None:
                            event_starts.append(i)
                        if radius_tracked[i+1] is None and radius_tracked[i-1] is not None:
                            event_ends.append(i)

            # 5. Classify the Type of Event (Positive and Negative Event Tracking)
            event_type = []
            event_start_time = []
            event_end_time = []
            events_count = []

            for start, end in zip(event_starts, event_ends):
                z_start = z_tracked[start]
                z_end = z_tracked[end]
                
                e_type = None
                e_count = 0.0

                # Check for full events within the simulation timeframe
                if start != 0 and end != len(z_tracked) - 1:
                    if z_end < mid_z and z_start > mid_z:
                        e_type, e_count = 'TB', 1.0
                    elif z_end < mid_z and z_start < mid_z:
                        e_type, e_count = 'BB', 0.0
                    elif z_end > mid_z and z_start < mid_z:
                        e_type, e_count = 'BT', -1.0
                    elif z_end > mid_z and z_start > mid_z:
                        e_type, e_count = 'TT', 0.0
                
                # Check for half events - the ion is inside at the start
                elif start == 0:
                    if z_end < mid_z:
                        e_type, e_count = 'MB', 0.5
                    elif z_end > mid_z:
                        e_type, e_count = 'MT', -0.5
                
                # Check for half events - the ion is inside at the end
                elif end == len(z_tracked) - 1:
                    if z_start < mid_z:
                        e_type, e_count = 'BM', -0.5
                    elif z_start > mid_z:
                        e_type, e_count = 'TM', 0.5

                # Append verified events
                if e_type is not None:
                    event_type.append(e_type)
                    events_count.append(e_count)
                    event_start_time.append(start)
                    event_end_time.append(end)
                    total_events += 1

            if event_type:
                pathway_list.append(event_type)
                ion_number_list.append(top.atom(ion_idx))
                event_start_list.append(event_start_time)
                event_end_list.append(event_end_time)
                events_count_list.append(events_count)

    # 6. Format and Save Output Data
    print(f"Total pathway events detected: {total_events}")
    
    if total_events > 0:
        event_data_combined = np.zeros((total_events, 3))
        pathway_type_out = []
        ion_residue_out = []
        row_idx = 0

        for i in range(len(pathway_list)):
            for j in range(len(pathway_list[i])):
                event_data_combined[row_idx, 0] = events_count_list[i][j]
                event_data_combined[row_idx, 1] = event_start_list[i][j]
                event_data_combined[row_idx, 2] = event_end_list[i][j]
                ion_residue_out.append(str(ion_number_list[i].residue))
                pathway_type_out.append(pathway_list[i][j])
                row_idx += 1

        # Save arrays to text files
        np.savetxt(event_data_file, event_data_combined, fmt='%.1f')

        with open(which_residue_file, 'w') as fp:
            fp.write('\n'.join(ion_residue_out))

        with open(pathway_file, 'w') as fp:
            fp.write('\n'.join(pathway_type_out))
            
        print(f"Data successfully written to *.{args.out_prefix}.txt files.")
    else:
        print("No events detected. Output files were not generated.")

    end_time = timer()
    print(f"Analysis completed in {end_time - start_time:.2f} seconds.")

if __name__ == "__main__":
    main()