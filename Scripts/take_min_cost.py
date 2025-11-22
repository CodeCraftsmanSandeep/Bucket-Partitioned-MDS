import os
import csv
import glob

# --- Configuration ---
outputs_dir = "Outputs"
sub_directories = ["XS", "S", "M", "L", "XL", "XXL", "XXXL"]
header_fields = ["input", "alpha", "rho", "time(sec)", "max_memory_diff(MB)", "cost"]

# ----------------------------------------------------------------------------------
## 1. Aggregation within each Size Directory (Level 1)
# ----------------------------------------------------------------------------------

def aggregate_size_directories(base_dir, sizes):
    """
    Processes each size directory, finds the best minimal 'cost' from all
    CSV files within its sub-subdirectories, and writes the results to a
    summary CSV (e.g., 'S.csv') inside the size directory.
    """
    print("--- Starting Level 1 Aggregation (Best 'cost' per input file) ---")
    
    for size_dir_name in sizes:
        size_path = os.path.join(base_dir, size_dir_name)
        if not os.path.isdir(size_path):
            print(f"Directory not found: {size_path}. Skipping.")
            continue
            
        # The target summary CSV path (e.g., Outputs/L/L.csv)
        summary_csv_path = os.path.join(size_path, f"{size_dir_name}.csv")
        all_best_rows = []
        
        # Walk through the immediate sub-directories (e.g., Outputs/L/sub_dir_1)
        # and search for CSVs within them.
        for root, dirs, files in os.walk(size_path):
            # Only process if we are NOT in the size_path itself (i.e., we are in a sub-sub-directory)
            # and if there are CSV files.
            if root != size_path and files:
                csv_files = [f for f in files if f.endswith(".csv")]
                
                for csv_file in csv_files:
                    input_name = os.path.splitext(csv_file)[0]
                    path = os.path.join(root, csv_file)

                    min_cost = float("inf")
                    best_row = None
                    
                    try:
                        with open(path, "r", newline="") as fin:
                            reader = csv.reader(fin)
                            # Skip header (assuming header is always present and in order)
                            try:
                                next(reader)
                            except StopIteration:
                                continue # Empty file

                            for row in reader:
                                # Ensure row has enough columns before accessing
                                if len(row) < 5:
                                    continue
                                
                                try:
                                    # Extract columns: alpha, rho, time, mem, cost
                                    cost = float(row[4])
                                    
                                    if cost < min_cost:
                                        min_cost = cost
                                        # Construct the full best_row tuple: (input, alpha, rho, time, mem, cost)
                                        # The 'input_name' is derived from the filename.
                                        best_row = (input_name, row[0], row[1], float(row[2]), float(row[3]), cost)

                                except ValueError:
                                    # Handle cases where non-float data is in a numeric column
                                    continue

                    except FileNotFoundError:
                        print(f"File not found: {path}")
                        continue
                        
                    if best_row:
                        all_best_rows.append(best_row)

        
        # Sort the results by 'input' name before writing the size summary CSV
        all_best_rows.sort(key=lambda x: x[0])
        
        if all_best_rows:
            print(f"Creating summary CSV: {summary_csv_path}")
            with open(summary_csv_path, "w", newline="") as fout:
                writer = csv.writer(fout)
                writer.writerow(header_fields)
                writer.writerows(all_best_rows)
            print(f"✔ Created: {summary_csv_path} with {len(all_best_rows)} entries.")
        else:
            print(f"No valid data found to create summary CSV for {size_dir_name}.")

# ----------------------------------------------------------------------------------
## 2. Final Aggregation (Level 2)
# ----------------------------------------------------------------------------------

def combine_final_output(base_dir, sizes):
    """
    Combines all size summary CSVs (e.g., XS.csv, S.csv, etc.) into the
    final Outputs.csv, sorted by 'input' name.
    """
    print("\n--- Starting Level 2 Aggregation (Final 'Outputs.csv') ---")
    final_output_path = os.path.join(base_dir, "Outputs.csv")
    all_final_rows = []
    
    # Iterate through the expected summary files
    for size_dir_name in sizes:
        summary_csv_path = os.path.join(base_dir, size_dir_name, f"{size_dir_name}.csv")
        
        if not os.path.exists(summary_csv_path):
            print(f"Summary CSV not found: {summary_csv_path}. Skipping.")
            continue
            
        print(f"Reading data from: {summary_csv_path}")
        try:
            with open(summary_csv_path, "r", newline="") as fin:
                reader = csv.reader(fin)
                try:
                    next(reader) # Skip header
                except StopIteration:
                    continue # Empty file
                    
                for row in reader:
                    # The structure of rows here is already (input, alpha, rho, time, mem, cost)
                    all_final_rows.append(row)
        except Exception as e:
            print(f"Error reading {summary_csv_path}: {e}")

    # Sort the combined results by 'input' name (the first element in the row)
    all_final_rows.sort(key=lambda x: x[0])

    if all_final_rows:
        print(f"Creating final combined CSV: {final_output_path}")
        with open(final_output_path, "w", newline="") as fout:
            writer = csv.writer(fout)
            # Add a 'size' column for clarity in the final output
            final_header = ["size"] + header_fields
            writer.writerow(final_header)
            
            # Reconstruct rows to include the 'size' prefix
            # This requires iterating through the size summary files again to associate size.
            # Rerun the collection loop for clean size association:
            final_data_with_size = []
            for size_dir_name in sizes:
                summary_csv_path = os.path.join(base_dir, size_dir_name, f"{size_dir_name}.csv")
                
                if not os.path.exists(summary_csv_path):
                    continue
                    
                with open(summary_csv_path, "r", newline="") as fin:
                    reader = csv.reader(fin)
                    try:
                        next(reader) # Skip header
                    except StopIteration:
                        continue
                        
                    for row in reader:
                        # Prepend the size directory name to the row
                        final_data_with_size.append([size_dir_name] + row)

            # Sort the fully-structured data by 'input' name (which is the second element now)
            final_data_with_size.sort(key=lambda x: x[1])
            
            # Write the final sorted rows
            writer.writerows(final_data_with_size)
        print(f"✨ Successfully created final output: {final_output_path} with {len(final_data_with_size)} entries.")
    else:
        print("No data found to create the final Outputs.csv.")

# ----------------------------------------------------------------------------------
## 3. Execution
# ----------------------------------------------------------------------------------

if __name__ == "__main__":
    aggregate_size_directories(outputs_dir, sub_directories)
    combine_final_output(outputs_dir, sub_directories)