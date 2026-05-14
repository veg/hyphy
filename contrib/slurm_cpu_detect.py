#!/usr/bin/env python3
import subprocess
import sys
import argparse

def get_slurm_nodes():
    try:
        # Get all unique nodes in the cluster
        output = subprocess.check_output(["sinfo", "-h", "-o", "%n"], text=True)
        return list(set(output.split()))
    except (subprocess.CalledProcessError, FileNotFoundError):
        return []

def get_node_flags(nodes):
    if not nodes:
        return None

    # Run grep on all nodes simultaneously to find CPU flags
    # We use -w to specify nodes, and --ntasks-per-node=1 to ensure one task per node
    # We join nodes with commas for -w
    nodes_str = ",".join(nodes)
    try:
        # Use a small timeout to avoid hanging
        cmd = ["srun", "-w", nodes_str, "--ntasks-per-node=1", "grep", "-m", "1", "flags", "/proc/cpuinfo"]
        output = subprocess.check_output(cmd, text=True, stderr=subprocess.STDOUT)
        
        node_flags = {}
        for line in output.splitlines():
            # Expected output format from srun: "node_name: flags : ..."
            if ":" in line:
                parts = line.split(":", 2)
                if len(parts) >= 3:
                    node_name = parts[0].strip()
                    flags = set(parts[2].strip().split())
                    node_flags[node_name] = flags
        return node_flags
    except subprocess.CalledProcessError as e:
        print(f"Error running srun: {e.output}", file=sys.stderr)
        return None

def find_common_flags(node_flags):
    if not node_flags:
        return set()
    
    common = None
    for flags in node_flags.values():
        if common is None:
            common = flags
        else:
            common = common.intersection(flags)
    return common if common else set()

def main():
    parser = argparse.ArgumentParser(description="Detect common CPU flags across SLURM nodes.")
    parser.add_argument("--cmake", action="store_true", help="Output in CMake -D format")
    args = parser.parse_args()

    nodes = get_slurm_nodes()
    if not nodes:
        if not args.cmake:
            print("No SLURM nodes found or sinfo not available.")
        sys.exit(0)

    node_flags = get_node_flags(nodes)
    if not node_flags:
        if not args.cmake:
            print("Failed to retrieve CPU flags from nodes.")
        sys.exit(1)

    common = find_common_flags(node_flags)

    # HyPhy specific flags
    # We want to determine if we should DISABLE certain features
    has_avx2 = "avx2" in common and "fma" in common
    has_avx = "avx" in common
    has_sse41 = "sse4_1" in common
    # NEON on ARM is usually indicated by 'neon' or 'asimd'
    has_neon = "neon" in common or "asimd" in common

    if args.cmake:
        # Output CMake definitions
        # If we DON'T have a feature, we set NO<FEATURE>=ON
        print(f"-DNOAVX2={'OFF' if has_avx2 else 'ON'}")
        print(f"-DNOAVX={'OFF' if has_avx else 'ON'}")
        print(f"-DNOSSE4={'OFF' if has_sse41 else 'ON'}")
        print(f"-DNONEON={'OFF' if has_neon else 'ON'}")
    else:
        print(f"Nodes surveyed: {len(node_flags)}")
        print(f"Common flags: {' '.join(sorted(list(common)))}")
        print(f"AVX2 support: {has_avx2}")
        print(f"AVX support: {has_avx}")
        print(f"SSE4.1 support: {has_sse41}")
        print(f"NEON support: {has_neon}")

if __name__ == "__main__":
    main()
