#!/usr/bin/env python3

# Based on https://github.com/SchulzLab/TEPIC/blob/master/MachineLearningPipelines/DYNAMITE/Scripts/integrateData.py

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

df_affinities = pd.read_csv("$affinity_ratio".replace("\\\\", ""), sep="\\t", index_col=0)
df_expression = pd.read_csv("$differential_expression".replace("\\\\", ""), sep="\\t", index_col=0)

def remove_version(gene_id):
    return gene_id.split(".")[0]

df_affinities.index = df_affinities.index.map(remove_version)
df_expression.index = df_expression.index.map(remove_version)

gene_intersection = df_affinities.index.intersection(df_expression.index)

df_affinities = df_affinities.loc[gene_intersection]
df_expression = df_expression.loc[gene_intersection]

# Aggregate duplicated genes from version clipping
df_affinities = df_affinities.groupby(df_affinities.index).mean()
df_expression = df_expression.groupby(df_expression.index).mean()

df_affinities["Expression"] = 0
df_affinities.loc[df_expression["log2FoldChange"] > 0, "Expression"] = 1

df_affinities.to_csv("${meta.id}.preprocessed.tsv", sep="\\t")

# Create version file
versions = {
    "${task.process}" : {
        "python": platform.python_version(),
        "pandas": pd.__version__
    }
}

with open("versions.yml", "w") as f:
    f.write(format_yaml_like(versions))
