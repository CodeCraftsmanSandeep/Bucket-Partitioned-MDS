# -*- coding: utf-8 -*-

import os
import subprocess

lambda_values = [1]
rho_values = [5000]
alpha_values = [15]

inputs_dir = "Inputs/XXL"
output_base_dir = "Scripts/Outputs"
binary_path = "./Bin/bucket-partitioned-MDS"

# Create output base directory if it doesn’t exist
if not os.path.exists(output_base_dir):
    os.makedirs(output_base_dir)

for root, _, files in os.walk(inputs_dir):
    for file in files:
        input_path = os.path.join(root, file)
        relative_path = os.path.relpath(root, inputs_dir)  # subfolder inside Inputs
        output_subdir = os.path.join(output_base_dir, relative_path)

        # Create corresponding output subdirectory
        if not os.path.exists(output_subdir):
            os.makedirs(output_subdir)

        file_base = os.path.splitext(file)[0]

        for lam in lambda_values:
            for rho in rho_values:
                for alpha in alpha_values:
                    output_file = "%s_alpha=%d_lambda=%d_rho=%d.txt" % (file_base, alpha, lam, rho)
                    output_path = os.path.join(output_subdir, output_file)

                    cmd = [
                        binary_path,
                        "--alpha=%d" % alpha,
                        "--lambda=%d" % lam,
                        "--rho=%d" % rho,
                        "--input=%s" % input_path,
                        "--output=%s" % output_path
                    ]

                    print("Running:", " ".join(cmd))
                    subprocess.call(cmd)
