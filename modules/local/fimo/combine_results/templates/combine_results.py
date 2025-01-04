#!/usr/bin/env python3

import os
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
        spaces = "    " * indent
        if isinstance(value, dict):
            yaml_str += f"{spaces}{key}:\\n{format_yaml_like(value, indent + 1)}"
        else:
            yaml_str += f"{spaces}{key}: {value}\\n"
    return yaml_str


output_dirs = [os.path.join('fimo', d) for d in os.listdir('fimo') if os.path.isdir(os.path.join('fimo', d))]

output_tsv = "${meta.id}.tsv"
output_gff = "${meta.id}.gff"

with open(output_tsv, 'w') as tsv_out, open(output_gff, 'w') as gff_out:
    tsv_out.write('motif_id\\tmotif_alt_id\\tsequence_name\\tstart\\tstop\\tstrand\\tscore\\tp-value\\tq-value\\tmatched_sequence\\n')

    for output in output_dirs:
        with open(f"{output}/fimo.tsv", "r") as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#') and not line.startswith('motif_id'):
                    tsv_out.write(line + "\\n")

        with open(f"{output}/fimo.gff", "r") as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    gff_out.write(line + "\\n")

# Create version file
versions = {
    "${task.process}" : {
        "python": platform.python_version()
    }
}

with open("versions.yml", "w") as f:
    f.write(format_yaml_like(versions))
