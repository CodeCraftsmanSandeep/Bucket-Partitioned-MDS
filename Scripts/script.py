# -*- coding: utf-8 -*-

import os
import subprocess
import numpy as np

num_executions = 5
rho_values = [5000]

inputs_dir      = "Inputs"
output_base_dir = "Outputs"
binary_path     = "./Bin/bucket-partitioned-MDS"

# Size categories
SMALL_SIZES = ["XS", "S", "M"]
LARGE_SIZES = ["L", "XL", "XXL", "XXXL"]

# Ensure output base dir
os.makedirs(output_base_dir, exist_ok=True)

for root, _, files in os.walk(inputs_dir):
    # Skip root folder "Inputs"
    if root == inputs_dir:
        continue

    # Detect first folder after Inputs/
    relative = os.path.relpath(root, inputs_dir)
    first_folder = relative.split(os.sep)[0]

    # Determine alpha range
    if first_folder == "XS":
        alpha_values = list(range(1, 360 + 1))
    elif first_folder == "S":
        alpha_values = list(range(1, 360 + 1))
    elif first_folder == "M":
        alpha_values = list(range(1, 180 + 1))
    elif first_folder == "L":
        alpha_values = list(range(1, 60 + 1))
    elif first_folder == "XL":
        alpha_values = list(range(1, 60 + 1))
    elif first_folder == "XXL":
        alpha_values = list(np.arange(0.1, 2.0 + 1e-6, 0.1)) + list(range(3, 16))
        alpha_values = [round(a, 2) for a in alpha_values]
    elif first_folder == "XXXL":
        alpha_values = list(np.arange(0.1, 2.0 + 1e-6, 0.1)) + list(range(3, 16))
        alpha_values = [round(a, 2) for a in alpha_values]
    else:
        raise ValueError(f"Unknown size category: {first_folder}")

    for file in files:
        input_path    = os.path.join(root, file)
        relative_path = os.path.relpath(root, inputs_dir)
        output_subdir = os.path.join(output_base_dir, relative_path)

        os.makedirs(output_subdir, exist_ok=True)

        file_base = os.path.splitext(file)[0]
        csv_path  = os.path.join(output_subdir, file_base + ".csv")

        with open(csv_path, "w") as csv_file:
            csv_file.write("alpha,rho,time(sec),max_memory_diff(MB),cost\n")
            csv_file.flush()

            for rho in rho_values:
                for alpha in alpha_values:

                    total_time = 0.0
                    total_max_memory_diff = 0.0
                    total_cost = 0.0

                    for _ in range(num_executions):

                        temp_output = os.path.join(output_subdir, "temp_output.txt")

                        cmd = [
                            binary_path,
                            f"--alpha={alpha}",
                            f"--rho={rho}",
                            f"--input={input_path}",
                            f"--output={temp_output}"
                        ]

                        print("Running:", " ".join(cmd))
                        subprocess.call(cmd)

                        with open(temp_output, "r") as f:
                            line1 = f.readline().strip()
                            line2 = f.readline().strip()
                            line3 = f.readline().strip()

                        os.remove(temp_output)

                        time_val = float(line1.split(":")[1].strip())
                        max_memory_diff = float(line2.split(":")[1].strip())
                        cost_val = float(line3.split(":")[1].strip())

                        total_time += time_val
                        total_max_memory_diff += max_memory_diff
                        total_cost += cost_val

                    avg_time = total_time / num_executions
                    avg_max_memory_diff = total_max_memory_diff / num_executions;
                    avg_cost = total_cost / num_executions

                    csv_file.write(f"{alpha},{rho},{avg_time:.6f},{avg_max_memory_diff:.6f},{avg_cost:.6f}\n")
                    csv_file.flush()


        print(f"✔ CSV created: {csv_path}")
