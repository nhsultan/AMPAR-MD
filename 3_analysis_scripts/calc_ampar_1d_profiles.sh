#!/bin/bash
# ==============================================================================
# calc_ampar_1d_profiles.sh
# 
# Description:
# Runs VMD to calculate the 3D electrostatic potential (PMEPot) and extracts
# the Z-coordinates of permeating ions within a defined pore cylinder.
# An embedded Python script then bins these into a 1D profile along the Z-axis.
# ==============================================================================

# --- Argument Parsing & Configuration ---
if [ "$#" -lt 3 ]; then
    echo "Error: Missing required arguments."
    echo "Usage: $0 <SYSTEM_PREFIX> <RUN_START> <RUN_END> [CYL_RADIUS]"
    echo "Example: $0 5WEO.WildType.charmm 2 8"
    echo "         $0 5WEO.Q586E.charmm 1 5 5.0"
    exit 1
fi

SYSTEM_PREFIX=$1
RUN_START=$2
RUN_END=$3
CYL_RADIUS=${4:-5.0} # Defaults to 5.0 Angstroms if a 4th argument isn't provided

CYL_RADIUS_SQ=$(echo "$CYL_RADIUS * $CYL_RADIUS" | bc -l)

# --- Cleanup Trap ---
# Ensures temp files are deleted even if the user hits Ctrl+C
trap 'rm -f process_1d_profile.py temp.dx temp_ions.txt temp_frames.txt run*_analysis.tcl; exit' INT TERM EXIT

# ==============================================================================
# 1. Create the Python Processing Script
# ==============================================================================
cat << 'EOF' > process_1d_profile.py
import sys

run_idx = sys.argv[1]
n_frames = int(sys.argv[2])
r2_cyl = float(sys.argv[3])

dx_file = "temp.dx"
ion_file = "temp_ions.txt"

# Binning parameters
z_min = -150.0
z_max = 150.0
dz = 1.0
n_bins = int((z_max - z_min) / dz)

pot_counts = [0] * n_bins
cla_counts = [0] * n_bins

# Step A: Read ion Z-coordinates from VMD and bin them
try:
    with open(ion_file, "r") as f:
        for line in f:
            parts = line.split()
            if len(parts) == 2:
                res, z = parts[0], float(parts[1])
                bin_idx = int((z - z_min) / dz)
                if 0 <= bin_idx < n_bins:
                    if res == "POT":
                        pot_counts[bin_idx] += 1
                    elif res == "CLA":
                        cla_counts[bin_idx] += 1
except FileNotFoundError:
    print(f"Warning: No ion file found for run {run_idx}")

# Step B: Parse the PMEPot DX file and calculate 1D average within the cylinder
pot_sum = [0.0] * n_bins
pot_pts = [0] * n_bins

ox, oy, oz = 0.0, 0.0, 0.0
dx_step, dy_step, dz_step = 0.0, 0.0, 0.0
nx, ny, nz = 0, 0, 0

try:
    with open(dx_file, "r") as f:
        lines = f.readlines()
        
    data_start = 0
    for i, line in enumerate(lines):
        parts = line.strip().split()
        if not parts: continue
        if parts[0] == "object" and parts[1] == "1":
            nx, ny, nz = int(parts[-3]), int(parts[-2]), int(parts[-1])
        elif parts[0] == "origin":
            ox, oy, oz = float(parts[1]), float(parts[2]), float(parts[3])
        elif parts[0] == "delta":
            if dx_step == 0.0: dx_step = float(parts[1])
            elif dy_step == 0.0: dy_step = float(parts[2])
            elif dz_step == 0.0: dz_step = float(parts[3])
        elif "data follows" in line:
            data_start = i + 1
            break

    # DX format is written with Z varying fastest
    data_idx = 0
    for i in range(data_start, len(lines)):
        line = lines[i]
        if line.startswith("attribute") or line.startswith("object"):
            break
        for v in line.strip().split():
            val = float(v)
            iz = data_idx % nz
            iy = (data_idx // nz) % ny
            ix = data_idx // (ny * nz)
            
            x = ox + ix * dx_step
            y = oy + iy * dy_step
            z = oz + iz * dz_step
            
            # Filter for grid points inside the pore cylinder
            if (x*x + y*y) <= r2_cyl:
                bin_idx = int((z - z_min) / dz)
                if 0 <= bin_idx < n_bins:
                    pot_sum[bin_idx] += val
                    pot_pts[bin_idx] += 1
            data_idx += 1
except FileNotFoundError:
    print(f"Warning: No DX file found for run {run_idx}")

# Step C: Write final unified output
out_name = f"run{run_idx}_profile.txt"
with open(out_name, "w") as f:
    f.write("# Z_mid Potential_avg(kT/e) POT_density(count/frame) CLA_density(count/frame)\n")
    for i in range(n_bins):
        z_mid = z_min + (i + 0.5) * dz
        avg_pot = pot_sum[i] / pot_pts[i] if pot_pts[i] > 0 else 0.0
        
        # Avoid division by zero if n_frames is 0
        if n_frames > 0:
            avg_pot_count = pot_counts[i] / n_frames
            avg_cla_count = cla_counts[i] / n_frames
        else:
            avg_pot_count, avg_cla_count = 0.0, 0.0
        
        # Only print rows that map inside the system's Z-bounds
        if pot_pts[i] > 0 or avg_pot_count > 0 or avg_cla_count > 0:
            f.write(f"{z_mid:.2f} {avg_pot:.6f} {avg_pot_count:.6f} {avg_cla_count:.6f}\n")

print(f"Successfully wrote {out_name}")
EOF

# ==============================================================================
# 2. Main Loop Over Simulation Runs
# ==============================================================================
for i in $(seq $RUN_START $RUN_END); do
    echo "=================================================="
    echo "Processing Run $i..."
    echo "=================================================="

    PSF_FILE="${SYSTEM_PREFIX}${i}.EFz.psf"
    NC_FILE="${SYSTEM_PREFIX}${i}.EFz.nc"
    TCL_SCRIPT="run${i}_analysis.tcl"

    # Generate VMD TCL Script
    cat <<EOF > $TCL_SCRIPT
mol new {${PSF_FILE}} type psf
mol addfile {${NC_FILE}} type netcdf step 100 waitfor all

set num_frames [molinfo top get numframes]
if {\$num_frames == 0} {
    puts "Error: No frames loaded for run $i."
    quit
}

# Export frame count for the Python script
set f_out [open "temp_frames.txt" w]
puts \$f_out \$num_frames
close \$f_out

puts "Centering trajectory X and Y axes on GLN 586..."
set sel_all [atomselect top all]
set sel_ref [atomselect top "resname GLN and resid 586"]

for {set f 0} {\$f < \$num_frames} {incr f} {
    \$sel_all frame \$f
    \$sel_ref frame \$f

    set com [measure center \$sel_ref weight mass]
    set cx [lindex \$com 0]
    set cy [lindex \$com 1]
    
    # Move ONLY X and Y to origin (0, 0). Z remains unchanged.
    set vec [list [expr {-1.0 * \$cx}] [expr {-1.0 * \$cy}] 0.0]
    \$sel_all moveby \$vec
}

puts "Extracting Ion Z-coordinates inside cylinder..."
set out [open "temp_ions.txt" w]
for {set f 0} {\$f < \$num_frames} {incr f} {
    # Re-select per frame to capture dynamic ion movement in/out of cylinder
    set pot [atomselect top "resname POT and (x*x + y*y) <= ${CYL_RADIUS_SQ}" frame \$f]
    foreach z [\$pot get z] { puts \$out "POT \$z" }
    \$pot delete

    set cla [atomselect top "resname CLA and (x*x + y*y) <= ${CYL_RADIUS_SQ}" frame \$f]
    foreach z [\$cla get z] { puts \$out "CLA \$z" }
    \$cla delete
}
close \$out

puts "Calculating PMEPot for electrostatic profile..."
package require pmepot
pmepot -mol top -grid 1.0 -frames all -dxfile "temp.dx"

quit
EOF

    # Execute VMD in text mode
    vmd -dispdev text -e $TCL_SCRIPT
    
    # Read the frame count output from VMD
    if [ -f "temp_frames.txt" ]; then
        FRAMES=$(cat temp_frames.txt)
    else
        FRAMES=0
    fi

    # Execute Python aggregator
    python3 process_1d_profile.py "$i" "$FRAMES" "$CYL_RADIUS_SQ"
    
    # Clean up intermediate files for this specific run
    rm -f $TCL_SCRIPT temp.dx temp_ions.txt temp_frames.txt

done

echo "=================================================="
echo "Analysis Complete! Profiles saved as run[${RUN_START}-${RUN_END}]_profile.txt"
echo "=================================================="