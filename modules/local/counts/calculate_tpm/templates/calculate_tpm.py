#!/usr/bin/env python3

import pandas as pd
import platform

def format_yaml_like(data: dict, indent: int = 0) -> str:
    """Formats a dictionary to a YAML-like string.

    Args:
        data (dict): The dictionary to format.
        indent (int): The current indentation level.

    Returns:
        str: A string formatted as YAML.
    """
    yaml_str = ""
    for key, value in data.items():
        spaces = "  " * indent
        if isinstance(value, dict):
            yaml_str += f"{spaces}{key}:\\n{format_yaml_like(value, indent + 1)}"
        else:
            yaml_str += f"{spaces}{key}: {value}\\n"
    return yaml_str

df_counts = pd.read_csv("$counts", index_col=0, header=0, sep="\\t")
df_lengths = pd.read_csv("$lengths", index_col=0, header=0, sep="\\t", usecols=["gene", "merged"])
df_lengths.columns = ["length"]
df_genes = pd.read_csv("$gene_map", sep="\\t", index_col=0)

conversion_dict = df_genes["gene_name"].to_dict()
df_lengths.index = df_lengths.index.map(conversion_dict).str.upper()

# Mean length for each gene (in kb)
df_lengths = df_lengths / 1e3
df_lengths = df_lengths.groupby(df_lengths.index).mean()

# Subset gene lengths and counts to common index
shared_index = df_lengths.index.intersection(df_counts.index)
df_lengths = df_lengths.loc[shared_index]
df_counts = df_counts.loc[shared_index]

# Calculate TPM
df_rpk = df_counts.div(df_lengths["length"], axis=0)
df_scale = df_rpk.sum() / 1e6
df_tpm = df_rpk.div(df_scale, axis=1)

# Save to file
df_tpm.to_csv("${meta.id}.tpm.tsv", sep="\\t")

# Create version file
versions = {
    "${task.process}" : {
        "python": platform.python_version(),
        "pandas": pd.__version__,
    }
}

with open("versions.yml", "w") as f:
    f.write(format_yaml_like(versions))
