
import argparse
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import glob
import os

def extract_participant_id(df, id_column='id'):
    """Extracts the participant ID from a column."""
    if id_column in df.columns:
        df['participant_id'] = df[id_column].str.split('_').str[0]
    return df

def consolidate_parquet_files(input_dir, output_file, id_column='id'):
    """Consolidates all parquet files in a directory into a single file using incremental writing."""
    parquet_files = glob.glob(os.path.join(input_dir, '*.parquet'))
    if not parquet_files:
        print("No parquet files found.")
        return

    writer = None
    for i, f in enumerate(parquet_files):
        print(f"Processing parquet file {i+1}/{len(parquet_files)}: {os.path.basename(f)}")
        df = pd.read_parquet(f)
        df = extract_participant_id(df, id_column)

        if writer is None:
            # Initialize ParquetWriter with the schema of the first DataFrame
            table = pa.Table.from_pandas(df)
            writer = pq.ParquetWriter(output_file, table.schema, compression='zstd')
        
        # Write the DataFrame to the Parquet file
        writer.write_table(pa.Table.from_pandas(df))

    if writer:
        writer.close()
        print(f"Successfully consolidated {len(parquet_files)} parquet files into {output_file}")
    else:
        print("No data was written to the parquet file.")

def consolidate_csv_files(input_dir, output_file, id_column='id'):
    """Consolidates all csv files in a directory into a single file."""
    csv_files = glob.glob(os.path.join(input_dir, '*.csv'))
    if not csv_files:
        print("No CSV files found.")
        return

    all_dfs = []
    for i, f in enumerate(csv_files):
        print(f"Processing CSV file {i+1}/{len(csv_files)}: {os.path.basename(f)}")
        df = pd.read_csv(f)
        df = extract_participant_id(df, id_column)
        all_dfs.append(df)

    if all_dfs:
        combined_df = pd.concat(all_dfs, ignore_index=True)
        combined_df.to_csv(output_file, index=False)
        print(f"Successfully consolidated {len(csv_files)} CSV files into {output_file}")


def main():
    parser = argparse.ArgumentParser(description='Consolidate parquet and csv/json files.')
    parser.add_argument('--input_dir', required=True, help='Directory containing the files to consolidate.')
    parser.add_argument('--output_dir', required=True, help='Directory to save the consolidated files.')
    parser.add_argument('--id_column', default='id', help='Name of the column containing the participant ID.')
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # Consolidate parquet files
    consolidate_parquet_files(args.input_dir, os.path.join(args.output_dir, 'combined_temperature_data.parquet'), args.id_column)

    # Consolidate csv files
    consolidate_csv_files(args.input_dir, os.path.join(args.output_dir, 'combined_covariates.csv'), args.id_column)


if __name__ == '__main__':
    main()
