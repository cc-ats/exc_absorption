"""
EXTRACT: Q-Chem Output Parser for Excited States, Interactions, and Geometries

This script processes Q-Chem output files (`.out`) to extract **excited state data**, 
**Coulomb & dipole-dipole interactions**, and **geometries** from the "Standard Nuclear 
Orientation" section. The extracted data is then **saved in JSON format** for further analysis.

---

## FUNCTIONALITY:

1. **Data Extraction**:
    - Reads **Q-Chem output files** (`.out`) and extracts:
        - **Excited states** (excitation energy, transition dipole moments, oscillator strengths).
        - **Interaction terms** (Coulomb & dipole-dipole interactions).
        - **Atomic positions** from the "Standard Nuclear Orientation" sections.

2. **Data Processing**:
    - Identifies and **parses excitation energies** for the first two excited states.
    - Extracts **transition dipole moments** along X, Y, and Z axes.
    - Retrieves **interaction matrices** for Coulomb and dipole-dipole interactions.
    - **Filters atomic geometries**, ensuring they belong to the proper `$rem` section.

3. **Data Output**:
    - Saves extracted **excited states, interactions, and geometries** into structured **JSON files**.
    - Each `.out` file generates a corresponding **`<filename>_output.json`**.

---

## USAGE:

### Command-line execution:
python3 extract.py <molecule_name> <functional>
#@author: Sanghita Sengupta, cc-@Shao Group
"""

import os
import glob
import argparse
import json
import re
import numpy as np

def parse_qchem_output(filename):
    """
    Parses a Q-Chem output file to extract excited states, interaction terms, and molecular geometries.

    This function processes the contents of a Q-Chem `.out` file and extracts:
    
    1. **Excited States**:
       - Excitation energies (in eV).
       - Transition dipole moments (X, Y, Z components).
       - Oscillator strengths.
       
    2. **Interaction Terms**:
       - Coulomb interaction (J_coul in eV).
       - Exchange interaction (J_ex in eV).
       - Dipole-dipole interaction (J_dip in eV).
       
    3. **Molecular Geometries**:
       - Extracts atomic positions from the "Standard Nuclear Orientation" section.
       - Filters valid geometries by ensuring they belong to the `$rem` block.

    ---
    
    Parameters:
    -----------
    filename : str
        Path to the Q-Chem `.out` file to be parsed.

    Returns:
    --------
    excited_states : list of dict
        A list of dictionaries, each containing:
        - 'Excited State' (int): Excited state number (1 or 2).
        - 'Excitation Energy (eV)' (float): Energy of the excited state.
        - 'Transition Dipole Moment' (list of floats): [Tx, Ty, Tz] in atomic units.
        - 'Oscillator Strength' (float): Strength of the transition.

    interactions : list of dict
        A list of dictionaries containing interaction terms between states:
        - 'State Pair' (tuple): Excited states interacting (e.g., (1,2)).
        - 'Coulomb Interaction (eV)' (float): Coulomb interaction strength.
        - 'Exchange Interaction (eV)' (float): Exchange interaction strength.
        - 'Dipole-Dipole Interaction (eV)' (float): Dipole-dipole interaction strength.

    geometries : dict
        A dictionary where keys are 'geometry N' (N being the index of the geometry),
        and values are lists of dictionaries with:
        - 'Atomic Index' (int): Index of the atom in the molecule.
        - 'Atom Type' (str): Element symbol (e.g., 'C', 'H', 'O').
        - 'Position Vector (Angstrom)' (list of floats): [x, y, z] coordinates.

    ---
    
    Notes:
    ------
    - Uses regular expressions to match patterns corresponding to excitation energies,
      transition dipole moments, and interactions in the Q-Chem output.
    - Only the first two excited states (1,2) are processed.
    - Handles potential parsing errors gracefully by printing warnings.
    - Calls `extract_geometries_with_rem_block()` to retrieve atomic coordinates.

    """
    excited_states = []
    interactions = []
    geometries = {}

    # Regular expression patterns for the data we need to extract
    excited_state_pattern = re.compile(r'Excited state\s+(\d+): excitation energy \(eV\)\s+=\s+([\d.]+)')
    trans_mom_pattern = re.compile(r'Trans\. Mom\.\:\s+([\-\d\.]+)\s+X\s+([\-\d\.]+)\s+Y\s+([\-\d\.]+)\s+Z')
    strength_pattern = re.compile(r'Strength\s+:\s+([\d.]+)')
    interaction_pattern = re.compile(r'\(\s*(\d+)\s*,\s*(\d+)\s*\):\s*([\d.\-eE]+)\s*([\d.\-eE]+)\s*([\d.\-eE]+)')
    
    with open(filename, 'r') as f:
        data = f.read()

    # Extract all matches for excited states
    excited_state_matches = excited_state_pattern.findall(data)
    trans_mom_matches = trans_mom_pattern.findall(data)
    strength_matches = strength_pattern.findall(data)
    interaction_matches = interaction_pattern.findall(data)

    # Process excited states and transition dipole moments
    for i, (state, excitation_energy) in enumerate(excited_state_matches):
        state_num = int(state)
        if state_num in [1, 2]:
            try:
                trans_mom_str = trans_mom_matches[i]
                trans_mom = np.array([float(trans_mom_str[0]), float(trans_mom_str[1]), float(trans_mom_str[2])])
                strength = float(strength_matches[i])
                
                sign = np.sign(trans_mom[np.argmax(np.absolute(trans_mom))])
                print (sign)
                excited_states.append({
                    'Excited State': state_num,
                    'Excitation Energy (eV)': float(excitation_energy),
                    'Transition Dipole Moment unmodified': trans_mom.tolist(),
                    'Transition Dipole Moment': (trans_mom * sign).tolist(),
                    'Oscillator Strength': strength,
                    'sign': np.sign(trans_mom[np.argmax(np.absolute(trans_mom))])
                })
            except Exception as e:
                print(f"Error parsing state {i+1}: {e}")

    # Process interactions (Coulomb, Dipole-Dipole, etc.)
    for match in interaction_matches:
        state1, state2, coulomb, exchange, interaction = match
        state_pair = (int(state1), int(state2))
        if state_pair in [(1,1), (1,2), (2,1), (2,2)]:
            sign = excited_states[state_pair[0]-1]['sign'] * excited_states[state_pair[1]+1]['sign']
            interactions.append({
                'State Pair': state_pair,
                'Coulomb Interaction (eV)': float(coulomb) * sign,
                'Exchange Interaction (eV)': float(exchange) * sign,
                'Dipole-Dipole Interaction (eV)': float(interaction) * sign
            })
        if int(state1)==2 and int(state2)==2:
            
            diff = interactions[-2]['Coulomb Interaction (eV)'] - interactions[-3]['Coulomb Interaction (eV)']
            print('diff', diff,  interactions[-2]['Coulomb Interaction (eV)'] )
    # Extract geometries (atomic positions from 'Standard Nuclear Orientation' sections)
    geometries = extract_geometries_with_rem_block(data)

    return excited_states, interactions, geometries

def extract_geometries_with_rem_block(data):
    """
    Extracts atomic positions from the 'Standard Nuclear Orientation' sections 
    in a Q-Chem output file, but only if they are preceded by a valid '$rem' block.

    This function identifies and extracts molecular geometries (atomic positions)
    from the 'Standard Nuclear Orientation' sections of a Q-Chem output file. 
    It ensures that the extracted geometries are part of a computation that 
    includes a '$rem' block (which defines calculation parameters) and ends with 
    'ECP $end'. 

    ---
    
    Parameters:
    -----------
    data : str
        The entire text content of a Q-Chem `.out` file.

    Returns:
    --------
    geometries : dict
        A dictionary where:
        - Keys are 'geometry N' (e.g., 'geometry 1', 'geometry 2', ...).
        - Values are lists of dictionaries, each containing:
            - 'Atomic Index' (int): The index of the atom in the molecule.
            - 'Atom Type' (str): The element symbol (e.g., 'C', 'H', 'O').
            - 'Position Vector (Angstrom)' (list of floats): [x, y, z] coordinates.

    ---
    
    Notes:
    ------
    - The `$rem ... ECP $end` block is used to verify the computation settings.
    - Uses regular expressions to extract atomic indices, element types, and 
      Cartesian coordinates in Angstroms.
    - Returns a dictionary of geometries where each entry corresponds to a different
      optimization step or scan in the Q-Chem calculation.

    """
    geometries = {}
    geometry_counter = 0

    # Regular expression to match the $rem ... ECP $end block
    rem_block_pattern = re.compile(r'\$rem.*?ECP\s*\$end', re.DOTALL)

    # Regular expression to extract atomic positions
    position_pattern = re.compile(r'\s*(\d+)\s+([A-Za-z]+)\s+([\d\.\-]+)\s+([\d\.\-]+)\s+([\d\.\-]+)')

    # Pattern to find 'Standard Nuclear Orientation'
    orientation_pattern = re.compile(r'Standard Nuclear Orientation \(Angstroms\)')

    # Find all occurrences of 'Standard Nuclear Orientation'
    orientation_matches = list(orientation_pattern.finditer(data))

    # Iterate over orientation matches to find atomic positions
    for i, orientation_match in enumerate(orientation_matches):
        # Skip the first block of positions (if needed)
        if i == 0:
            continue

        orientation_start = orientation_match.end()
        orientation_end = -1 if i == len(orientation_matches) - 1 else orientation_matches[i + 1].start()

        if orientation_start != -1:
            geometry_counter += 1
            geometry_label = f'geometry {geometry_counter}'
            atomic_positions = []

            # Extract atomic positions in the 'Standard Nuclear Orientation' section
            atomic_matches = position_pattern.findall(data[orientation_start:orientation_end])

            for match in atomic_matches:
                atomic_index = int(match[0])
                atom_type = match[1]
                position_vector = np.array([float(match[2]), float(match[3]), float(match[4])])

                atomic_positions.append({
                    'Atomic Index': atomic_index,
                    'Atom Type': atom_type,
                    'Position Vector (Angstrom)': position_vector.tolist()
                })

            if atomic_positions:
                geometries[geometry_label] = atomic_positions

    return geometries

def save_output(output_data, output_file):
    """Save the extracted data to a JSON file."""
    with open(output_file, 'w') as f:
        json.dump(output_data, f, indent=4)

def generate_output_filename(input_file):
    """Generate an output filename for the processed data."""
    base_name = os.path.basename(input_file)
    base_name = os.path.splitext(base_name)[0]
    output_filename = base_name + "_output.json"
    return output_filename

def main():
    """Main function to parse Q-Chem files and save the extracted data."""
    parser = argparse.ArgumentParser(description="Process Q-Chem output files for a specified molecule.")
    parser.add_argument('molecule_name', type=str, help='Name of the molecule folder within Qchem_files.')
    parser.add_argument('functional', type=str, help='Name of the molecule folder within Qchem_files.')

    args = parser.parse_args()

    molecule_folder = os.path.join('..', 'systems', 'Qchem_files', args.molecule_name, f'{args.functional}/*.out')

    # Find all .out files within the molecule folder
    matched_files = glob.glob(molecule_folder)

    if not matched_files:
        print(f"No .out files matched the pattern in: {molecule_folder}")
        return

    # Ensure the output folder exists
    base_output_folder = os.path.join('..', 'systems', 'Qchem_files', args.molecule_name, args.functional)
    os.makedirs(base_output_folder, exist_ok=True)

    for qchem_file in matched_files:
        # Parse the Q-Chem file to extract data
        excited_states, interactions, geometries = parse_qchem_output(qchem_file)

        # Prepare the output data
        output_data = {
            'excited_states': excited_states,
            'interactions': interactions,
            'geometries': geometries
        }

        # Generate the output filename
        output_filename = generate_output_filename(qchem_file)

        # Save the extracted data to a JSON file
        output_file_path = os.path.join(base_output_folder, output_filename)
        save_output(output_data, output_file_path)

        print(f"Processed {qchem_file}, results saved to {output_file_path}")

if __name__ == "__main__":
    main()
