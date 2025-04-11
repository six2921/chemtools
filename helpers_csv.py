import pandas as pd
from rdkit import Chem

# ---------------------
# SMILES Checker Functions
# ---------------------

def is_valid_smiles(s: str) -> bool:
    """Return True if SMILES is valid, else False."""
    mol = Chem.MolFromSmiles(s)
    return mol is not None

def standardize_smiles(s: str) -> str:
    """Return canonical SMILES as a standardized version (if valid)."""
    mol = Chem.MolFromSmiles(s)
    return Chem.MolToSmiles(mol, canonical=True) if mol else ""

def neutralize_smiles(s: str) -> str:
    """
    Return a neutralized SMILES string.
    (This is a placeholder implementation using canonical SMILES.)
    """
    mol = Chem.MolFromSmiles(s)
    return Chem.MolToSmiles(mol, canonical=True) if mol else ""

def strip_smiles(s: str) -> str:
    """
    Return a SMILES string with salts removed.
    Keeps only the largest fragment.
    """
    mol = Chem.MolFromSmiles(s)
    if not mol:
        return ""
    fragments = Chem.GetMolFrags(mol, asMols=True)
    if not fragments:
        return ""
    largest = max(fragments, key=lambda m: m.GetNumHeavyAtoms())
    return Chem.MolToSmiles(largest, canonical=True)

# ---------------------
# SMILES Sanitizer (Interactive)
# ---------------------

def sanitize_smiles(smiles_series: pd.Series, column_name: str) -> pd.Series:
    """
    Interactively sanitize SMILES values in a Pandas Series.
    
    The user is prompted with the following transformation options:
      (t) - Standardize (convert to canonical SMILES)
      (n) - Neutralize (apply neutralization; placeholder implementation)
      (p) - Strip salts (retain only the largest fragment)
      (s) - Save and finish
    
    Once an option is applied, that option will no longer be offered.
    
    Returns:
      Updated SMILES series.
    """
    current_series = smiles_series.copy()
    applied = {"standardized": False, "neutralized": False, "stripped": False}
    
    while True:
        print("\nCurrent sample SMILES: {}".format(current_series.iloc[0]))
        print("Which transformation would you like to perform next?")
        # Build options string based on transformations not yet applied.
        options = []
        if not applied["standardized"]:
            options.append("(t) Standardize")
        if not applied["neutralized"]:
            options.append("(n) Neutralize")
        if not applied["stripped"]:
            options.append("(p) Strip salts")
        options.append("(s) Save and finish")
        print("  " + "\n  ".join(options))
        
        choice = input("Your choice (t/n/p/s): ").strip().lower()
        if choice == 't' and not applied["standardized"]:
            new_series = current_series.apply(standardize_smiles)
            print("Standardization applied.")
            current_series = new_series
            applied["standardized"] = True
        elif choice == 'n' and not applied["neutralized"]:
            new_series = current_series.apply(neutralize_smiles)
            print("Neutralization applied.")
            current_series = new_series
            applied["neutralized"] = True
        elif choice == 'p' and not applied["stripped"]:
            new_series = current_series.apply(strip_smiles)
            print("Salt stripping applied.")
            current_series = new_series
            applied["stripped"] = True
        elif choice == 's':
            print("Finalizing changes.")
            break
        else:
            print("Invalid option or transformation already applied. Please choose a valid, not-yet-applied option.")
    return current_series

# ---------------------
# Column Listing and Selection Functions
# ---------------------

def list_columns(df: pd.DataFrame) -> None:
    """
    Print column information in a formatted table with headers.
    Format:
      idx (width=5), column (width=15), type (width=10), value(first) (width=15)
    """
    header = f"{'idx':<5} {'column':<15} {'type':<10} {'value(first)':<15}"
    print(header)
    print("-" * len(header))
    for idx, col in enumerate(df.columns):
        col_type = str(df[col].dtype)
        value_first = str(df[col].iloc[0]) if not df.empty else ""
        # Truncate strings to the desired width if they are too long.
        col_disp = col[:15] if len(col) > 15 else col
        type_disp = col_type[:10] if len(col_type) > 10 else col_type
        value_disp = value_first[:15] if len(value_first) > 15 else value_first
        print(f"{idx:<5d} {col_disp:<15} {type_disp:<10} {value_disp:<15}")

def select_columns(df: pd.DataFrame, *col_types) -> tuple:
    """
    Display the column information (using list_columns) once, then prompt the user to select one or more columns.
    
    Usage examples:
      select_columns(df, "name")
      select_columns(df, "smiles")
      select_columns(df, "name", "smiles")
      select_columns(df, "text", "numer")
    
    For each column type provided:
      - A prompt "Enter the number for the <col_type> column: " is displayed.
      - If col_type is "smiles", the first value of the column is validated to be a SMILES.
      - If col_type is "numer" or "number", the first value is validated to be numeric.
    
    Returns:
      A tuple containing the selected column names corresponding to each col_type.
    """
    list_columns(df)  # Print the column listing once.
    selected = []
    for col_type in col_types:
        prompt = f"{col_type} column: "
        while True:
            try:
                index = int(input(prompt).strip())
                if index < 0 or index >= len(df.columns):
                    raise ValueError
                col_name = df.columns[index]
                # If the column type is 'smiles', validate the first value.
                if col_type.lower() == "smiles":
                    test_val = str(df[col_name].iloc[0]).strip()
                    if not is_valid_smiles(test_val):
                        print("The selected column does not appear to contain valid SMILES. Please choose a different column.")
                        continue
                # If the column type is 'numer' or 'number', try to cast the first value to float.
                if col_type.lower() in ["numer", "number"]:
                    try:
                        float(df[col_name].iloc[0])
                    except ValueError:
                        print("The selected column does not appear to contain numeric values. Please choose a different column.")
                        continue
                selected.append(col_name)
                break
            except ValueError:
                print("Invalid input. Please enter a valid column number.")
    return tuple(selected)
