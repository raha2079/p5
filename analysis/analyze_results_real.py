#!/usr/bin/env python3
"""
Protein Folding - Real Data Analysis
Generates figures from actual CSV exports
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
import sys

def load_metadata(metadata_file):
    """Load metadata"""
    metadata = {}
    with open(metadata_file, 'r') as f:
        section = None
        for line in f:
            line = line.strip()
            if line.startswith('['):
                section = line[1:-1]
            elif '=' in line and section:
                key, value = line.split('=', 1)
                metadata[f"{section}.{key.strip()}"] = value.strip()
    return metadata

def plot_landscape(results_dir):
    """Generate landscape.png from 15s seed 42 vertices"""
    vertices = pd.read_csv(results_dir / 'char_15s_seed42_data_vertices.csv')

    fig, ax = plt.subplots(figsize=(10, 8))

    # Clip energy for visualization
    e_min, e_max = np.percentile(vertices['energy'], [1, 99])

    sc = ax.scatter(vertices['R_ee'], vertices['R_g'],
                   c=vertices['energy'], s=15, alpha=0.6,
                   cmap='viridis', vmin=e_min, vmax=e_max)

    # Mark q0 (vertex 0) and lowest energy state
    q0_idx = 0
    lowest_idx = vertices['energy'].idxmin()

    ax.scatter(vertices.loc[q0_idx, 'R_ee'], vertices.loc[q0_idx, 'R_g'],
              marker='*', s=400, c='red', edgecolors='black', linewidths=2,
              label=f'q₀ (E={vertices.loc[q0_idx, "energy"]:.2f})', zorder=10)
    ax.scatter(vertices.loc[lowest_idx, 'R_ee'], vertices.loc[lowest_idx, 'R_g'],
              marker='*', s=400, c='cyan', edgecolors='black', linewidths=2,
              label=f'E_min={vertices.loc[lowest_idx, "energy"]:.2f}', zorder=10)

    ax.set_xlabel('R_ee (Å)', fontsize=12)
    ax.set_ylabel('R_g (Å)', fontsize=12)
    ax.set_title('15s Characterization: Conformational Landscape (seed 42)', fontsize=14)
    ax.legend()
    ax.grid(alpha=0.3)

    plt.colorbar(sc, ax=ax, label='Energy')
    plt.tight_layout()
    plt.savefig(results_dir / 'landscape.png', dpi=300)
    print(f"Saved: landscape.png")
    print(f"  q₀: R_ee={vertices.loc[q0_idx, 'R_ee']:.2f}, R_g={vertices.loc[q0_idx, 'R_g']:.2f}, E={vertices.loc[q0_idx, 'energy']:.2f}")
    print(f"  E_min: R_ee={vertices.loc[lowest_idx, 'R_ee']:.2f}, R_g={vertices.loc[lowest_idx, 'R_g']:.2f}, E={vertices.loc[lowest_idx, 'energy']:.2f}")
    print(f"  Energy range: [{vertices['energy'].min():.2f}, {vertices['energy'].max():.2f}]")
    plt.close()

def plot_serial_trace(results_dir):
    """Generate serial_trace.png from 15s seed 42 splits"""
    splits = pd.read_csv(results_dir / 'char_15s_seed42_data_splits.csv')

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), sharex=True)

    # Leaf count progression
    ax1.plot(splits['iteration'], splits['total_leaves'], 'b-o', linewidth=2, markersize=6)
    ax1.set_ylabel('Total Leaf Cells', fontsize=12)
    ax1.set_title('Quadtree Subdivision Progression (15s seed 42)', fontsize=14)
    ax1.grid(alpha=0.3)

    # Depth progression
    ax2.plot(splits['iteration'], splits['depth'], 'r-s', linewidth=2, markersize=6)
    ax2.set_xlabel('Iteration', fontsize=12)
    ax2.set_ylabel('Split Depth', fontsize=12)
    ax2.set_title('Split Depth Over Time', fontsize=12)
    ax2.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(results_dir / 'serial_trace.png', dpi=300)
    print(f"Saved: serial_trace.png")
    print(f"  Total splits: {len(splits)}")
    print(f"  Final leaves: {splits.iloc[-1]['total_leaves']}")
    print(f"  Max depth: {splits['depth'].max()}")
    plt.close()

def plot_comparison_coverage(results_dir):
    """Generate comparison_coverage.png from 10s runs"""
    planners = ['quadtree_energy', 'quadtree_no_energy', 'rrt']
    names = ['Quadtree (energy)', 'Quadtree (no energy)', 'RRT']

    # Fixed grid [0, 40.41] x [0, 11.67] - 20x20 bins
    L = 40.41
    bins = 20
    x_edges = np.linspace(0, L, bins+1)
    y_edges = np.linspace(0, 11.67, bins+1)

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    for idx, (planner, name) in enumerate(zip(planners, names)):
        ax = axes[idx]

        # Combine all 3 seeds
        all_x, all_y = [], []
        for seed in [42, 43, 44]:
            vfile = results_dir / f'comp_10s_{planner}_seed{seed}_data_vertices.csv'
            if vfile.exists():
                v = pd.read_csv(vfile)
                all_x.extend(v['R_ee'])
                all_y.extend(v['R_g'])

        # 2D histogram
        H, _, _ = np.histogram2d(all_x, all_y, bins=[x_edges, y_edges])
        occupied = np.sum(H > 0)

        # Heatmap
        im = ax.imshow(H.T, origin='lower', extent=[0, L, 0, 11.67],
                      cmap='YlOrRd', aspect='auto', interpolation='nearest')
        ax.set_xlabel('R_ee (Å)')
        ax.set_ylabel('R_g (Å)')
        ax.set_title(f'{name}\\n{occupied}/{bins*bins} bins occupied')
        plt.colorbar(im, ax=ax, label='Count')

    fig.suptitle('Coverage Comparison (10s runs, 3 seeds combined)', fontsize=14)
    plt.tight_layout()
    plt.savefig(results_dir / 'comparison_coverage.png', dpi=300)
    print(f"Saved: comparison_coverage.png")
    plt.close()

def compute_aggregate(results_dir):
    """Generate aggregate.csv from 10s runs"""
    planners = ['quadtree_energy', 'quadtree_no_energy', 'rrt']
    seeds = [42, 43, 44]

    L = 40.41
    bins = 20
    x_edges = np.linspace(0, L, bins+1)
    y_edges = np.linspace(0, 11.67, bins+1)

    rows = []
    for planner in planners:
        occupied_bins = []
        num_states = []
        total_attempts = []
        acceptance_rates = []
        min_energies = []
        median_energies = []

        for seed in seeds:
            # Load vertices
            vfile = results_dir / f'comp_10s_{planner}_seed{seed}_data_vertices.csv'
            mfile = results_dir / f'comp_10s_{planner}_seed{seed}_metadata.txt'

            if vfile.exists() and mfile.exists():
                v = pd.read_csv(vfile)
                meta = load_metadata(mfile)

                # Occupied bins
                H, _, _ = np.histogram2d(v['R_ee'], v['R_g'], bins=[x_edges, y_edges])
                occupied_bins.append(np.sum(H > 0))

                # Metadata
                num_states.append(int(meta['Results.num_states']))
                total_attempts.append(int(meta['Results.total_attempts']))
                acceptance_rates.append(float(meta['Results.acceptance_rate']))
                min_energies.append(float(meta['Results.min_energy']))
                median_energies.append(float(meta['Results.median_energy']))

        if len(occupied_bins) == 3:
            rows.append({
                'planner': planner,
                'occupied_bins_mean': np.mean(occupied_bins),
                'occupied_bins_std': np.std(occupied_bins),
                'num_states_mean': np.mean(num_states),
                'num_states_std': np.std(num_states),
                'total_attempts_mean': np.mean(total_attempts),
                'total_attempts_std': np.std(total_attempts),
                'acceptance_rate_mean': np.mean(acceptance_rates),
                'acceptance_rate_std': np.std(acceptance_rates),
                'min_energy_mean': np.mean(min_energies),
                'min_energy_std': np.std(min_energies),
                'median_energy_mean': np.mean(median_energies),
                'median_energy_std': np.std(median_energies)
            })
            # Also save per-seed bins
            for i, seed in enumerate(seeds):
                rows[-1][f'seed{seed}_bins'] = occupied_bins[i]

    df = pd.DataFrame(rows)
    df.to_csv(results_dir / 'aggregate.csv', index=False, float_format='%.4f')
    print(f"Saved: aggregate.csv")
    print(df[['planner', 'occupied_bins_mean', 'num_states_mean', 'acceptance_rate_mean']])

def main():
    project_root = Path(__file__).parent.parent
    results_dir = project_root / 'results'

    print("Generating figures from real CSV data...")
    print("=" * 80)

    plot_landscape(results_dir)
    print()
    plot_serial_trace(results_dir)
    print()
    plot_comparison_coverage(results_dir)
    print()
    compute_aggregate(results_dir)
    print()
    print("=" * 80)
    print("Analysis complete!")

if __name__ == '__main__':
    main()
