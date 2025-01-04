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

def remove_version(gene_id):
    return gene_id.split(".")[0]

df_differential = pd.read_csv("$differential".replace("\\\\", ""), sep='\\t', index_col=0)
df_affinities = pd.read_csv("$affinities".replace("\\\\", ""), sep='\\t', index_col=0)
df_coefficients = pd.read_csv("$regression_coefficients".replace("\\\\", ""), sep='\\t', index_col=0)

# Remove version from gene ids
df_differential.index = df_differential.index.map(remove_version)
df_affinities.index = df_affinities.index.map(remove_version)

# Make sure genes are in common between the differential expression and affinities files
gene_intersection = df_differential.index.intersection(df_affinities.index)
assert len(gene_intersection) > 0, "No genes found in common between the differential expression and affinities files"

df_affinities = df_affinities.loc[gene_intersection]
df_differential = df_differential.loc[gene_intersection]

# Aggregate duplicated genes from version clipping
df_affinities = df_affinities.groupby(df_affinities.index).mean()
df_differential = df_differential.groupby(df_differential.index).mean()

# Make sure TFs are in common between the affinities and coefficients files
tf_intersection = df_affinities.columns.intersection(df_coefficients.index)
assert len(tf_intersection) > 0, "No TFs found in common between the affinities and coefficients files"

df_affinities = df_affinities[tf_intersection]
df_coefficients = df_coefficients.loc[tf_intersection]


# Calculate the TF-TG scores

## Multiply the log2FC by the affinities
result = (df_affinities
            .mul(abs(df_differential["log2FoldChange"]), axis=0)
            .mul(abs(df_coefficients["value"]), axis=1))

## Make sure results are not empty
assert not result.empty, "No TF-TG scores were calculated"

# Save the result
result.to_csv("${meta.id}.score.tsv", sep='\\t')

# Create version file
versions = {
    "${task.process}" : {
        "python": platform.python_version(),
        "pandas": pd.__version__
    }
}

with open("versions.yml", "w") as f:
    f.write(format_yaml_like(versions))
