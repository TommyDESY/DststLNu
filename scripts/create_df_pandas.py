import pandas as pd
from importlib_resources import files
import uproot
import argparse
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.getcwd(), '..')))
import root, samples


def root_to_pandas(file_name: str, tree_name: str) -> pd.DataFrame:
    """
    Loads root file via uproot and converts to a dataframe.

    Args:
        file_name (str): path to file
        tree_name (str): name of tree in file

    Returns:
        pd.DataFrame: dataframe of the ntuple stored in the root file.
    """
    tree = uproot.open(file_name)[tree_name]
    try:
        df = tree.arrays(library='pd')[0]
    except:
        df = tree.arrays(library='pd')
    finally:
        assert isinstance(df, pd.DataFrame)
    # defragment with the copy.
    return df.copy()


def parse_cmdline_args() -> argparse.Namespace:
    """
    Creates the argument parser and returns command line arguments.
    """
    parser = argparse.ArgumentParser(description='Inclusive B->Xlv reconstruction steering file.')
    parser.add_argument(
        '-input',
        '--input-path',
        dest='input_path',
        help='Path to the file that has the input files stored',
    )

    parser.add_argument(
        '-output',
        '--output-file',
        dest='output_file',
        help='Full path of the output file',
    )
    return parser.parse_args()


def main():
    args = parse_cmdline_args()
    # Input files
    input_files = args.input_path
    # Output file
    output_filename = args.output_file
    # Save the data frame to pandas parquet format
    df = root_to_pandas(files(root).joinpath(input_files), 'InclusiveSLGen')
    df['__eventType__'] = df['__eventType__'].astype(str)
    df.to_parquet(files(samples).joinpath(output_filename), engine='pyarrow')
    print(f'Pandas data frame succesfully saved in DststLNu/{files(samples).name}!')


if __name__ == '__main__':
    main()
