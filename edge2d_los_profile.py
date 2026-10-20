"""
author: R. Cicioni
last updated: 31/07/2026

This code reconstruct the LoS signal from the 2D profile map of a JT-60SA edge2d simulation.
This code can be run on both JDC and EFGW servers. The run directory is found
automatically from the run name entered when the script starts.
1) To change the quantities of interest, change the values in the 'PROFILES' list.
2) To change the line of sight, simply change the coordinates of 'point1' and 'point2' variables.
"""
import sys
import os


def ask_run_name(prompt):
    """Ask for a non-empty run name."""
    while True:
        run_name = input(prompt).strip()
        if run_name:
            return run_name
        print("The run name cannot be empty.")


def get_run_name():
    """Read the run name from the command line, or ask interactively."""
    if len(sys.argv) > 2:
        raise ValueError("Usage: python jintrac_los_profile.py [run_name]")
    if len(sys.argv) == 2:
        return sys.argv[1].strip()
    return ask_run_name("Run name: ")


def find_run_directory(run_dir):
    """Find a run directory in the standard JDC and EFGW locations."""
    potential_paths = [
        f"/common/cmg/jnv7243/jetto/runs/{run_dir}",
        f"/pfs/work/g2rcicio/jetto/runs/{run_dir}",
        f"/common/cmg/jnv7243/edge2d/runs/{run_dir}",
        f"/pfs/work/g2rcicio/edge2d/runs/{run_dir}",
    ]

    for path in potential_paths:
        if os.path.exists(path):
            print(f"Run found at: {path}")
            return path

    raise FileNotFoundError(
        f"Run '{run_dir}' not found in any of the expected locations:\n"
        + "\n".join(potential_paths)
    )


run_name = get_run_name()
RUN_DIR = find_run_directory(run_name)

import matplotlib
matplotlib.use("TkAgg") 
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

sys.path.append('/home/lgarzot/Work/Python/JET') #jintrac python tools on JDC server
sys.path.append('/afs/eufus.eu/user/g/g2fjc/jintrac/37.0.0/python') #jintrac python tools on EFGW server
sys.path.append('/afs/eufus.eu/user/g/g2fjc/jintrac/37.0.0/libs') #eproc tools on EFGW server

import numpy as np
import h5py
import jetto_binary_tools
import eproc as ep
import matplotlib.pyplot as plt

TRAN_FILE = os.path.join(RUN_DIR, "tran")
PROFILES = ['TEVE', 'DENEL'] #signal names to process (add other TRAN variables here)
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
HDF5_OUTPUT = os.path.join(SCRIPT_DIR, f"{run_name}_los_profiles.h5")
# Define line of sight endpoints
point1 = np.array([4.5, 0.0]) #line of sight endpoint 1 (can be changed to any other coordinates)
point2 = np.array([2.5, 2.5]) #line of sight endpoint 2 (can be changed to any other coordinates)

rmesh=ep.data(TRAN_FILE,'RMESH') #centrid of the mesh in R direction
rmesh=rmesh.data[0:rmesh.nPts]
zmesh=ep.data(TRAN_FILE,'ZMESH') #centrid of the mesh in Z direction
zmesh=-zmesh.data[0:zmesh.nPts]

rvertp=ep.data(TRAN_FILE,'RVERTP') #vertices of the mesh in R direction (quintuple structure: 4 vertices + 1 centroid)
rvertp=rvertp.data[0:rvertp.nPts]
zvertp=ep.data(TRAN_FILE,'ZVERTP') #vertices of the mesh in Z direction (quintuple structure: 4 vertices + 1 centroid)
zvertp=-zvertp.data[0:zvertp.nPts]

vessel_r = ep.data(TRAN_FILE, 'RVESM2')
vessel_z = ep.data(TRAN_FILE, 'ZVESM2')
vessel_r = vessel_r.data[0:vessel_r.nPts]
vessel_z = -vessel_z.data[0:vessel_z.nPts]

sep_r = ep.data(TRAN_FILE, 'RSEPX')
sep_z = ep.data(TRAN_FILE, 'ZSEPX')
sep_r = sep_r.data[0:sep_r.nPts]
sep_z = -sep_z.data[0:sep_z.nPts]
sep_mask = np.isfinite(sep_r) & np.isfinite(sep_z)
sep_r = sep_r[sep_mask]
sep_z = sep_z[sep_mask]

# ============================================================================
# Extract nodes and quads from rvertp/zvertp quintuple structure
# ============================================================================
# Data structure: quintuple organization
# - Every 5 consecutive elements form one quad: 4 vertices + 1 centroid
# - Vertices are at indices [0,1,2,3], centroid at index [4]

n_total = len(rvertp)
n_quads = n_total // 5  # Number of quads

print(f"Extraxting mesh structure from rvertp/zvertp...")
print(f"Total number of points (rvertp,zvertp): {n_total}")
print(f"Number of quads: {n_quads}")

# Initialize nodes array (4 vertices per quad)
nodes = np.zeros((4 * n_quads, 2))
quads = np.zeros((n_quads, 4), dtype=int)

# Extract vertices (skip centroids which are at indices 5k+4)
node_idx = 0
for quad_idx in range(n_quads):
    # Extract the 4 vertices from this quintuple
    for vertex_in_quad in range(4):
        global_idx = 5 * quad_idx + vertex_in_quad
        nodes[node_idx, 0] = rvertp[global_idx]  # R coordinate
        nodes[node_idx, 1] = zvertp[global_idx]  # Z coordinate
        quads[quad_idx, vertex_in_quad] = node_idx
        node_idx += 1

# Printing info about the extracted mesh structure
print(f"\nExtracted mesh structure:")
print(f"  nodes shape: {nodes.shape}")
print(f"  quads shape: {quads.shape}")

# ============================================================================
# Map centroids to mesh points (rmesh, zmesh)
# ============================================================================
from scipy.spatial import cKDTree

# Extract centroids coordinates from rvertp/zvertp
centroid_indices = np.arange(4, len(rvertp), 5)
centroids = np.column_stack([rvertp[centroid_indices], zvertp[centroid_indices]])

# Build KDTree from mesh points for fast nearest-neighbor search
mesh_points = np.column_stack([rmesh, zmesh])
tree = cKDTree(mesh_points)

distances, indices = tree.query(centroids)
print(f"Centroid-to-mesh mapping:")
print(f"  Number of centroids: {len(centroids)}")
print(f"  Number of mesh points (rmesh): {len(rmesh)}")
print(f"  Max distance to nearest mesh point: {distances.max():.2e}")
print(f"{'='*70}\n")

# ============================================================================
# LoS intersection algorithm (bsed on Sutherland-Hodgman algorithm)
# ============================================================================

def cross2d(a, b):
    #Scalar cross product of 2-vectors.
    return a[0] * b[1] - a[1] * b[0]


def line_segment_intersect_t(origin, direction, p1, p2):
    """
    Intersect the infinite line X = origin + t*direction
    with the finite segment [p1, p2].
    Returns t and s parameters, or (None, None) if parallel.
    """
    d = direction
    v = p2 - p1
    denom = cross2d(d, v)
    if abs(denom) < 1e-12:
        return None, Noneh5py
    w = p1 - origin
    t = cross2d(w, v) / denom
    s = cross2d(w, d) / denom
    return t, s


def ray_quad_interval(origin, direction, quad_nodes):
    """
    Find interval [t_enter, t_exit] where line is inside quadrilateral.
    Returns (t_enter, t_exit) or None if line misses quad.
    """
    n = len(quad_nodes)
    ts = []
    for k in range(n):
        p1 = quad_nodes[k]
        p2 = quad_nodes[(k + 1) % n]
        t, s = line_segment_intersect_t(origin, direction, p1, p2)
        if t is not None and -1e-9 <= s <= 1 + 1e-9:
            ts.append(t)

    if len(ts) < 2:
        return None

    ts = sorted(ts)
    t_enter, t_exit = ts[0], ts[-1]
    if t_exit - t_enter < 1e-10:
        return None
    return t_enter, t_exit


def extract_profile(nodes, quads, values, origin, direction):
    """
    For each quad, compute the segment of the line inside it.
    Returns sorted list of (t_enter, t_exit, value).
    """
    segments = []
    for qi, quad in enumerate(quads):
        quad_nodes = nodes[quad]
        result = ray_quad_interval(origin, direction, quad_nodes)
        if result is not None:
            t_enter, t_exit = result
            segments.append((t_enter, t_exit, values[qi]))

    segments.sort(key=lambda s: s[0])
    return segments


def load_profile_data(profile_code):
    """
    Read one TRAN variable and map it from mesh points to cell centroids.
    """
    profile_obj = ep.data(TRAN_FILE, profile_code)
    profile_name = str(profile_obj.desc) if hasattr(profile_obj, 'desc') else profile_code
    profile_units = str(profile_obj.units) if hasattr(profile_obj, 'units') else 'N/A'
    profile_time = getattr(profile_obj, 'time', None)
    profile = profile_obj.data[0:profile_obj.nPts]
    if profile_obj.nPts == 0 or len(profile) == 0:
        print(f"WARNING: {profile_code} not found in TRAN file. Skipping this variable.")
        return None
    if len(profile) <= indices.max():
        print(
            f"WARNING: {profile_code} has {len(profile)} points, "
            f"but mesh mapping requires index {indices.max()}. Skipping this variable."
        )
        return None
    values = profile[indices]

    print(f"{profile_code} values mapping:")
    print(f"  Description: {profile_name}")
    print(f"  Units: {profile_units}")
    if profile_time is not None:
        print(f"  TRAN time: {profile_time:.3f} s")
    print(f"  Values extracted: {len(values)}")
    print(f"  Range: {values.min():.2f} - {values.max():.2f}")
    print(f"\nFirst 10 centroids and their {profile_code} values:")
    print("-" * 72)
    for i in range(min(10, len(centroids))):
        r_cent, z_cent = centroids[i]
        r_mesh, z_mesh = mesh_points[indices[i]]
        print(f"{i:>6}  {r_cent:>10.4f}  {z_cent:>10.4f}  {r_mesh:>10.4f}  {z_mesh:>10.4f}  {values[i]:>12.4f}  {distances[i]:>10.2e}")

    print(f"\n Successfully created {profile_code} values array.")
    print(f"{'='*70}\n")

    return {
        "code": profile_code,
        "description": profile_name,
        "units": profile_units,
        "time": profile_time,
        "values": values,
    }


def segment_arrays(segments, origin, direction):
    """
    Convert LoS segments to arrays useful for plotting and HDF5 export.
    """
    t_enter = np.array([s[0] for s in segments])
    t_exit = np.array([s[1] for s in segments])
    t_mid = 0.5 * (t_enter + t_exit)
    r_mid = origin[0] + t_mid * direction[0]
    z_mid = origin[1] + t_mid * direction[1]
    los_length = np.abs(t_exit - t_enter) * np.sqrt(direction[0]**2 + direction[1]**2)
    values = np.array([s[2] for s in segments])
    return t_enter, t_exit, t_mid, r_mid, z_mid, los_length, values


def export_los_profiles_hdf5(output_file, run_name, tran_file, point1, point2, profile_results):
    """
    Export LoS profiles to a clean HDF5 file for IDL post-processing.
    Physical variables are written as top-level datasets named after the TRAN variable.
    """
    with h5py.File(output_file, "w") as h5:
        h5.attrs["description"] = f"Line of sight profiles for run {run_name}"
        h5.attrs["run_name"] = run_name
        h5.attrs["tran_file"] = tran_file

        dset = h5.create_dataset("point1", data=point1)
        dset.attrs["units"] = "m"
        dset.attrs["description"] = "Line of sight start point [R, Z]"

        dset = h5.create_dataset("point2", data=point2)
        dset.attrs["units"] = "m"
        dset.attrs["description"] = "Line of sight end point [R, Z]"

        if profile_results:
            first = profile_results[0]
            for name, data, units, description in [
                ("los_t_enter", first["t_enter"], "dimensionless", "Entry parameter t for each intersected cell"),
                ("los_t_exit", first["t_exit"], "dimensionless", "Exit parameter t for each intersected cell"),
                ("los_t_mid", first["t_mid"], "dimensionless", "Midpoint parameter t for each intersected cell"),
                ("los_r_mid", first["r_mid"], "m", "R coordinate at each LoS cell midpoint"),
                ("los_z_mid", first["z_mid"], "m", "Z coordinate at each LoS cell midpoint"),
                ("los_length", first["los_length"], "m", "Length of the LoS intersection with each cell"),
            ]:
                dset = h5.create_dataset(name, data=data)
                dset.attrs["units"] = units
                dset.attrs["description"] = description

        for result in profile_results:
            table = np.column_stack([
                result["los_values"],
                result["r_mid"],
                result["z_mid"],
                result["los_length"],
            ])
            dset = h5.create_dataset(result["code"], data=table)
            dset.attrs["units"] = f"{result['units']}, m, m, m"
            dset.attrs["columns"] = "value, r_mid, z_mid, los_length"
            if result.get("time") is not None:
                dset.attrs["time"] = result["time"]
            dset.attrs["description"] = (
                f"Line of sight profile for run {run_name}; "
                f"TRAN variable {result['code']}: {result['description']}. "
                "Columns are value, R coordinate, Z coordinate, intersection length."
            )
            dset.attrs["tran_variable"] = result["code"]

    print(f"HDF5 export written to: {output_file}")
    return True


def plot_los_profile(profile_result):
    """
    Open one figure for one physical variable, preserving the existing plotting logic.
    """
    profile_name = profile_result["description"]
    profile_units = profile_result["units"]
    profile_time = profile_result.get("time")
    values = profile_result["values"]
    segments = profile_result["segments"]

    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    # --- Left: 2D mesh with line of sight ---
    ax = axes[0]
    title = f"{profile_name} + Line of Sight"
    if profile_time is not None:
        title += f" (TRAN t = {profile_time:.3f} s)"
    ax.set_title(title, fontsize=12, fontweight='bold')
    ax.set_aspect("equal")

    norm = plt.Normalize(vmin=values.min(), vmax=values.max())
    cmap = plt.cm.rainbow

    patches = []
    colors = []
    for qi, quad in enumerate(quads):
        poly = Polygon(nodes[quad], closed=True)
        patches.append(poly)
        colors.append(values[qi])

    pc = PatchCollection(patches, cmap=cmap, norm=norm, edgecolors="black", linewidths=0.3, alpha=0.9)
    pc.set_array(np.array(colors))
    ax.add_collection(pc)
    plt.colorbar(pc, ax=ax, label=f"{profile_name} ({profile_units})")

    # Draw vessel contour
    ax.plot(vessel_r, vessel_z, color='0.35', linewidth=2.0, label='Vessel contour', zorder=8)

    # Draw separatrix from RSEPX/ZSEPX.
    ax.plot(sep_r, sep_z, color='red', marker='.', markersize=4, linewidth=1.6, label='Separatrix', zorder=9)

    # Draw line of sight
    t_all = [s[0] for s in segments] + [s[1] for s in segments]
    if t_all:
        t_lo = min(t_all) - 0.1
        t_hi = max(t_all) + 0.1
    else:
        t_lo, t_hi = -0.5, 1.5

    ax.plot([point1[0], point2[0]], [point1[1], point2[1]],
            "r-", linewidth=2.5, label="Line of Sight", zorder=10)

    # Mark line endpoints
    ax.plot(*point1, "go", markersize=12, label="Point 1", zorder=11, markeredgecolor='darkgreen', markeredgewidth=2)
    ax.plot(*point2, "mo", markersize=12, label="Point 2", zorder=11, markeredgecolor='purple', markeredgewidth=2)

    # Mark intersected quad centers
    for t_enter, t_exit, val in segments:
        t_mid = 0.5 * (t_enter + t_exit)
        pt = origin + t_mid * direction
        ax.plot(*pt, "r.", markersize=8, zorder=12)

    ax.set_xlabel('R (m)')
    ax.set_ylabel('Z (m)')
    ax.legend(loc='upper right', fontsize=10)
    ax.grid(True, alpha=0.3)

    # --- Right: 1D profile ---
    ax2 = axes[1]
    title_1d = f"1D {profile_name} Along Line of Sight"
    if profile_time is not None:
        title_1d += f" (TRAN t = {profile_time:.3f} s)"
    ax2.set_title(title_1d, fontsize=12, fontweight='bold')

    for t_enter, t_exit, val in segments:
        color = cmap(norm(val))
        ax2.fill_between([t_enter, t_exit], [val, val], step="post", alpha=0.7, color=color)
        ax2.hlines(val, t_enter, t_exit, colors=color, linewidths=2.5)
        ax2.vlines([t_enter, t_exit], 0, val, colors="gray", linewidths=0.6, linestyles="--")

    ax2.set_xlabel("Parameter t along ray")
    ax2.set_ylabel(f"{profile_name} ({profile_units})")
    ax2.set_xlim(t_lo, t_hi)
    ax2.set_ylim(values.min() - 0.05 * (values.max() - values.min()),
                 values.max() + 0.1 * (values.max() - values.min()))
    ax2.grid(True, alpha=0.3)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    plt.colorbar(sm, ax=ax2, label=f"{profile_name} ({profile_units})")

    plt.tight_layout()


# ============================================================================
# Apply LoS to extract and plot all requested profiles
# ============================================================================
# Convert to parametric form: origin + t*direction
origin = point1
direction = point2 - point1

print(f"\n{'='*70}")
print("LINE OF SIGHT PROFILES")
print(f"{'='*70}")
print(f"Origin: {origin}")
print(f"Direction: {direction}")
print(f"Point 1: {point1}")
print(f"Point 2: {point2}")

profile_results = []
for profile_code in PROFILES:
    profile_result = load_profile_data(profile_code)
    if profile_result is None:
        continue
    segments = extract_profile(nodes, quads, profile_result["values"], origin, direction)
    t_enter, t_exit, t_mid, r_mid, z_mid, los_length, los_values = segment_arrays(segments, origin, direction)

    profile_result.update({
        "segments": segments,
        "t_enter": t_enter,
        "t_exit": t_exit,
        "t_mid": t_mid,
        "r_mid": r_mid,
        "z_mid": z_mid,
        "los_length": los_length,
        "los_values": los_values,
    })
    profile_results.append(profile_result)

    print(f"\n{profile_code}: number of intersected cells: {len(segments)}")
    print(f"{'t_enter':>10}  {'t_exit':>10}  {'length_m':>10}  {'value':>10}")
    print("-" * 50)
    for t_in, t_out, length_m, val in zip(t_enter, t_exit, los_length, los_values):
        print(f"{t_in:10.4f}  {t_out:10.4f}  {length_m:10.4f}  {val:10.4f}")

    plot_los_profile(profile_result)

if not profile_results:
    raise ValueError("No requested profiles were found in the TRAN file.")

export_los_profiles_hdf5(HDF5_OUTPUT, run_name, TRAN_FILE, point1, point2, profile_results)

plt.show()

print(f"\nLine of Sight analysis complete")
print(f"{'='*70}\n")
