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

tsvs = []
gffs = []
for output in output_dirs:
    with open(f'{output}/fimo.tsv', 'r') as f:
        tsv = f.read().split('\\n')
    with open(f'{output}/fimo.gff', 'r') as f:
        gff = f.read().split('\\n')

    tsvs.extend(tsv)
    gffs.extend(gff)

tsvs = [line for line in tsvs if not line.startswith('#') and not line.startswith('motif_id') and not line == '']
gffs = [line for line in gffs if not line.startswith('#') and not line == '']

tsvs = ['motif_id\\tmotif_alt_id\\tsequence_name\\tstart\\tstop\\tstrand\\tscore\\tp-value\\tq-value\\tmatched_sequence'] + tsvs

with open('${meta.id}.tsv', 'w') as f:
    f.write('\\n'.join(tsvs))

with open('${meta.id}.gff', 'w') as f:
    f.write('\\n'.join(gffs))


# Create version file
versions = {
    "${task.process}" : {
        "python": platform.python_version()
    }
}

with open("versions.yml", "w") as f:
    f.write(format_yaml_like(versions))
