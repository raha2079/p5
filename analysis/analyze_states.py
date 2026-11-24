#!/usr/bin/env python3
"""
Protein Folding State Analysis

This script reads exported state data from the protein folding planner
and provides basic verification and analysis capabilities.

For Milestone 1, this demonstrates that:
1. Python can successfully read the CSV format
2. The expected columns are present
3. The data has the correct structure

Future milestones will add:
- Energy landscape visualization (2D projection with energy coloring)
- Pathway analysis between conformations
- Statistical analysis of the conformational space
"""

import sys
import csv
import argparse
from pathlib import Path
from typing import List, Dict, Any


def read_states_csv(filename: str) -> tuple[List[Dict[str, Any]], List[str]]:
    """
    Read states from a CSV file exported by the C++ planner.

    Args:
        filename: Path to the CSV file

    Returns:
        Tuple of (list of state dictionaries, list of column names)
    """
    states = []
    columns = []

    with open(filename, 'r') as f:
        reader = csv.DictReader(f)
        columns = reader.fieldnames or []

        for row in reader:
            # Convert numeric fields
            state = {}
            for key, value in row.items():
                try:
                    # Try to convert to float
                    state[key] = float(value)
                except ValueError:
                    # Keep as string if not numeric
                    state[key] = value
            states.append(state)

    return states, columns


def verify_format(columns: List[str], expected_residues: int = 2) -> bool:
    """
    Verify that the CSV has the expected format.

    Args:
        columns: List of column names from the CSV
        expected_residues: Expected number of residues (default: 2 for test data)

    Returns:
        True if format is valid
    """
    print(f"Verifying CSV format...")
    print(f"  Columns found: {columns}")

    # Check for state_id
    if 'state_id' not in columns:
        print("  ✗ Missing 'state_id' column")
        return False
    print("  ✓ Found 'state_id' column")

    # Check for phi and psi columns for each residue
    for i in range(expected_residues):
        phi_col = f'phi_{i}'
        psi_col = f'psi_{i}'

        if phi_col not in columns:
            print(f"  ✗ Missing '{phi_col}' column")
            return False
        if psi_col not in columns:
            print(f"  ✗ Missing '{psi_col}' column")
            return False

        print(f"  ✓ Found '{phi_col}' and '{psi_col}' columns")

    # Check for optional energy column
    if 'energy' in columns:
        print("  ✓ Found optional 'energy' column")

    return True


def print_summary(states: List[Dict[str, Any]], columns: List[str]) -> None:
    """
    Print a summary of the loaded states.

    Args:
        states: List of state dictionaries
        columns: List of column names
    """
    print(f"\nData summary:")
    print(f"  Number of states: {len(states)}")
    print(f"  Columns: {', '.join(columns)}")

    if len(states) > 0:
        print(f"\nFirst state:")
        for key, value in states[0].items():
            print(f"    {key}: {value}")

        if len(states) > 1:
            print(f"\nLast state:")
            for key, value in states[-1].items():
                print(f"    {key}: {value}")

        # Basic statistics for numeric columns
        print(f"\nNumeric column ranges:")
        for col in columns:
            if col != 'state_id' and all(isinstance(s.get(col), (int, float)) for s in states):
                values = [s[col] for s in states]
                print(f"  {col}: [{min(values):.4f}, {max(values):.4f}]")


def main():
    parser = argparse.ArgumentParser(
        description='Analyze protein conformational states from CSV export'
    )
    parser.add_argument('csv_file', help='CSV file containing state data')
    parser.add_argument('--residues', type=int, default=2,
                        help='Expected number of residues (default: 2)')
    parser.add_argument('--verify-only', action='store_true',
                        help='Only verify format without detailed analysis')

    args = parser.parse_args()

    # Check file exists
    csv_path = Path(args.csv_file)
    if not csv_path.exists():
        print(f"Error: File not found: {args.csv_file}")
        return 1

    print(f"Reading states from: {args.csv_file}")
    print("=" * 60)

    try:
        # Read the CSV
        states, columns = read_states_csv(args.csv_file)

        # Verify format
        if not verify_format(columns, args.residues):
            print("\n✗ CSV format verification FAILED")
            return 1

        print("\n✓ CSV format verification PASSED")

        # Print summary unless verify-only mode
        if not args.verify_only:
            print_summary(states, columns)

        print("\n" + "=" * 60)
        print("Analysis complete!")

        # Note about future functionality
        if 'energy' in columns:
            print("\nNote: Energy landscape visualization will be added in future milestones.")
        else:
            print("\nNote: To visualize energy landscape, export states with energy values.")

        return 0

    except Exception as e:
        print(f"\n✗ Error: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == '__main__':
    sys.exit(main())
