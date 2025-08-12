import pickle
import pandas as pd
import pyreadr
import os
import sys
from tqdm import tqdm

def convert_pickle_to_rds(pickle_path, rds_path=None):
    """
    Converts a pickle file to an RDS file.

    Args:
        pickle_path (str): Path to the pickle file.
        rds_path (str, optional): Path to save the RDS file. 
            If None, it will be generated from the pickle file name
            by replacing the extension with .rds.
            Defaults to None.

    Returns:
        pandas.DataFrame: The converted data as a pandas DataFrame.
                        Returns None if there is an error
    """
    print(f"Processing {pickle_path}...")

    # 1. Input Validation and Path Handling
    if not os.path.exists(pickle_path):
        print(f"Error: Pickle file not found at {pickle_path}")
        return None

    if not pickle_path.lower().endswith(".pkl") and not pickle_path.lower().endswith(".pickle"):
        print(f"Error: Invalid file extension.  Expected .pkl or .pickle, got {pickle_path}")
        return None
    
    # Determine output rds_path
    if rds_path is None:
        base_name = os.path.splitext(os.path.basename(pickle_path))[0]
        rds_path = os.path.join(os.path.dirname(pickle_path), f"{base_name}.rds")
    # Ensure the directory exists
    os.makedirs(os.path.dirname(rds_path), exist_ok=True)

    # 2. Load Pickle File with Progress
    file_size = os.path.getsize(pickle_path)
    try:
        with tqdm(total=file_size, unit='B', unit_scale=True, desc="Loading pickle") as pbar:
            with open(pickle_path, 'rb') as f:
                data = b''
                chunk_size = 1024 * 1024  # 1MB chunks
                while True:
                    chunk = f.read(chunk_size)
                    if not chunk:
                        break
                    data += chunk
                    pbar.update(len(chunk))
        data_loaded = pickle.loads(data)
    except Exception as e:
        print(f"Error: Failed to load pickle file: {e}")
        return None

    # 3. Convert to DataFrame if necessary
    try:
        with tqdm(desc="Converting to DataFrame", total=1) as pbar:
            if not isinstance(data_loaded, pd.DataFrame):
                data_loaded = pd.DataFrame(data_loaded)
            pbar.update(1)
    except Exception as e:
        print(f"Error: Failed to convert data to DataFrame: {e}")
        return None

    # 4. Save as RDS with Progress
    try:
        with tqdm(desc=f"Saving to {rds_path}", total=1) as pbar:
            pyreadr.write_rds(rds_path, data_loaded)
            pbar.update(1)
        print(f"✓ Successfully converted {pickle_path} to {rds_path}")
        return data_loaded  # Return the DataFrame
    except Exception as e:
        print(f"Error: Failed to save RDS file: {e}")
        return None

if __name__ == "__main__":
    # Check if the user provided a pickle file path
    if len(sys.argv) != 2:
        print("Usage: python pickle_to_rds.py <pickle_file_path>")
        sys.exit(1)

    pickle_file_path = sys.argv[1]
    
    # Convert the pickle file to RDS
    df = convert_pickle_to_rds(pickle_file_path) # No need to provide output path, will create default

    if df is not None:
        print("\nConversion Summary:")
        print(f"Data shape: {df.shape}")
        print("Column information:")
        print(df.head())  # Print the first few rows of the DataFrame
        print("\nConversion completed.")
    else:
        print("\nConversion failed.")