#!/usr/bin/env python3
"""
calc_ampar_conductance.py

Description:
This script calculates the ion channel conductance from ion flux and box dimension 
results generated in the AMPAR analysis pipeline. It processes multiple simulation 
runs to compute the total and per-ion currents and conductances based on the 
applied electric field and simulation time.

Outputs:
A summary text file containing overall and per-ion statistics (Current and Conductance) 
averaged across all processed simulation runs.
"""

import argparse
import re
import os
import numpy as np

def main():
    # 1. Setup Argument Parsing
    parser = argparse.ArgumentParser(description="Calculate AMPAR Ion Channel Conductance")
    parser.add_argument("--data_dir", default="RawData", 
                        help="Directory containing input files (default: RawData)")
    parser.add_argument("--out_file", default="AMPAR_Conductance_Analysis_Summary.txt", 
                        help="Output filename (default: AMPAR_Conductance_Analysis_Summary.txt)")
    parser.add_argument("--prefixes", nargs='+', required=True, 
                        help="List of system prefixes to process (e.g., Q586E_0 Q586E_1 wild_type)")
    parser.add_argument("--efield", type=float, default=0.11, 
                        help="Applied electric field in kcal/(mol*A*e) (default: 0.11)")
    parser.add_argument("--dt", type=float, default=100.0, 
                        help="Frame step in picoseconds (default: 100.0)")
    args = parser.parse_args()

    # 2. Initialize Variables
    # Ion Valences (Matched to standard AMPAR tracking: POT, CLA, SOD, CAL)
    valence_map = {
        'CAL': 2.0,  # Calcium (sometimes CAL or CAM depending on forcefield)
        'CAM': 2.0,  
        'POT': 1.0,  # Potassium
        'SOD': 1.0,  # Sodium
        'CLA': -1.0  # Chloride
    }

    # Dictionaries to collect data for statistics across runs
    results_current = {'Total': []}
    results_conductance = {'Total': []}
    for ion in valence_map.keys():
        results_current[ion] = []
        results_conductance[ion] = []

    num_runs = len(args.prefixes)
    print(f"--- AMPAR Ion Channel Conductance Analysis ({num_runs} Runs) ---")
    print(f"Results will be saved to: {args.out_file}\n")

    # 3. Open output file and Process Runs
    with open(args.out_file, 'w') as out_file:
        
        # Write Header
        out_file.write(f"--- AMPAR Ion Channel Conductance Analysis ({num_runs} Runs) ---\n\n")

        for run_prefix in args.prefixes:
            
            # Construct filenames dynamically based on your AMPAR script outputs
            file_box = os.path.join(args.data_dir, f'BoxSummary.{run_prefix}.txt')
            file_tracking = os.path.join(args.data_dir, f'IonTracking.{run_prefix}.txt')
            file_ions = os.path.join(args.data_dir, f'Residue.{run_prefix}.txt')

            header_msg = f"=== Processing {run_prefix} ==="
            print(header_msg)
            out_file.write(header_msg + "\n")

            # 4. Parse Box Info (Z-Dimension and Time)
            try:
                with open(file_box, 'r') as f:
                    content = f.read()
                    
                    z_match = re.search(r'Z Dimension:\s*([\d\.]+)', content)
                    if z_match:
                        z_length_angstrom = float(z_match.group(1))
                    else:
                        raise ValueError("Could not find 'Z Dimension' in BoxSummary file.")

                    frames_match = re.search(r'Total Frames Processed:\s*(\d+)', content)
                    if frames_match:
                        total_frames = int(frames_match.group(1))
                    else:
                        raise ValueError("Could not find 'Total Frames Processed' in BoxSummary file.")
            except FileNotFoundError:
                err_msg = f"Error: {file_box} not found. Skipping {run_prefix}.\n"
                print(err_msg)
                out_file.write(err_msg)
                continue
            except Exception as e:
                err_msg = f"Error reading {file_box}: {e}\n"
                print(err_msg)
                out_file.write(err_msg)
                continue

            # Calculate Time and Voltage
            total_time_ns = total_frames * (args.dt / 1000.0)
            total_time_s = total_time_ns * 1e-9

            voltage_energy = args.efield * z_length_angstrom 
            # Convert Energy to Volts: (kcal/mol * 4184) / F
            voltage_volts = (voltage_energy * 4184.0) / 96485.332

            # 5. Parse Ion Events
            try:
                with open(file_tracking, 'r') as ft, open(file_ions, 'r') as fi:
                    track_lines = ft.readlines()
                    ion_lines = fi.readlines()
            except FileNotFoundError:
                err_msg = f"Error: Tracking/Residue files not found for {run_prefix}. Skipping.\n"
                print(err_msg)
                out_file.write(err_msg)
                continue

            count_valid = 0
            total_charge_transported = 0.0
            charge_by_ion = {ion: 0.0 for ion in valence_map.keys()}
            ion_counts = {ion: 0.0 for ion in valence_map.keys()}

            num_events = min(len(track_lines), len(ion_lines))
            
            for i in range(num_events):
                track_line = track_lines[i].strip()
                ion_string = ion_lines[i].strip().upper()

                if not track_line or not ion_string:
                    continue

                parts = track_line.split()
                event_count = float(parts[0])

                matched_ion = None
                for key in valence_map.keys():
                    if key in ion_string:
                        matched_ion = key
                        break
                
                if matched_ion:
                    valence = valence_map[matched_ion]
                    charge_moved = valence * event_count
                    
                    total_charge_transported += charge_moved
                    charge_by_ion[matched_ion] += charge_moved
                    
                    ion_counts[matched_ion] += abs(event_count)
                    count_valid += 1

            # 6. Calculate Results for this run
            e_coulombs = 1.60217663e-19
            
            current_amps_total = (total_charge_transported * e_coulombs) / total_time_s
            current_pA_total = current_amps_total * 1e12
            
            conductance_pS_total = 0.0
            if voltage_volts != 0:
                conductance_siemens_total = current_amps_total / voltage_volts
                conductance_pS_total = conductance_siemens_total * 1e12

            results_current['Total'].append(current_pA_total)
            results_conductance['Total'].append(conductance_pS_total)

            run_report_lines = [
                f"  Box Z-Length: {z_length_angstrom:.2f} Angstroms",
                f"  Voltage:      {voltage_volts*1000:.2f} mV",
                f"  Time:         {total_time_ns:.2f} ns",
                f"  Events Tracker: POT: {ion_counts['POT']}, CLA: {ion_counts['CLA']}, SOD: {ion_counts['SOD']}",
                f"  --- Overall ---",
                f"  Total Net Charge:  {total_charge_transported:.4f} e",
                f"  Total Current:     {current_pA_total:.4f} pA",
                f"  Total Conductance: {conductance_pS_total:.4f} pS",
                f"  --- By Ion Type ---"
            ]

            for ion in valence_map.keys():
                charge = charge_by_ion[ion]
                amps = (charge * e_coulombs) / total_time_s
                pA = amps * 1e12
                
                pS = 0.0
                if voltage_volts != 0:
                    pS = (amps / voltage_volts) * 1e12
                    
                results_current[ion].append(pA)
                results_conductance[ion].append(pS)
                
                if ion_counts[ion] > 0:
                    run_report_lines.append(f"    {ion} -> Charge: {charge:>6.2f} e | Current: {pA:>8.4f} pA | Conductance: {pS:>8.4f} pS")

            report = "\n".join(run_report_lines) + "\n"
            print(report)
            out_file.write(report + "\n")

        # 7. Final Statistics
        stats_header = f"=== Final Statistics ({len(results_conductance['Total'])} Runs Processed) ==="
        print(stats_header)
        out_file.write(stats_header + "\n")

        if len(results_conductance['Total']) > 0:
            
            def format_stats(current_list, conductance_list, label):
                if not current_list or not conductance_list:
                    return ""
                
                avg_curr = np.mean(current_list)
                sem_curr = np.std(current_list, ddof=1) / np.sqrt(len(current_list)) if len(current_list) > 1 else 0.0
                
                avg_cond = np.mean(conductance_list)
                sem_cond = np.std(conductance_list, ddof=1) / np.sqrt(len(conductance_list)) if len(conductance_list) > 1 else 0.0
                
                return f"  {label:<6} | Current: {avg_curr:>8.4f} ± {sem_curr:>6.4f} pA  |  Conductance: {avg_cond:>8.4f} ± {sem_cond:>6.4f} pS\n"

            final_stats = "\n" + format_stats(results_current['Total'], results_conductance['Total'], 'Total')
            final_stats += "-" * 75 + "\n"
            
            for ion in valence_map.keys():
                if np.any(results_conductance[ion]):
                    final_stats += format_stats(results_current[ion], results_conductance[ion], ion)

            print(final_stats)
            out_file.write(final_stats)
        else:
            print("No valid runs processed.")
            out_file.write("No valid runs processed.\n")

if __name__ == "__main__":
    main()