#!/usr/bin/env python3

from jinja2 import Environment, PackageLoader, select_autoescape
import jinja2
import os
import shutil
import json
import pandas as pd
import sys
from collections import defaultdict

module_app = os.path.abspath("$moduleDir/app")
app_dir = "app"
out_dir = "report"

# Copy app_dir to current directory
shutil.copytree(module_app, os.path.join(os.getcwd(), app_dir), dirs_exist_ok=True)

# Copy dependencies to report
shutil.copytree(os.path.join(app_dir, "dependencies"), os.path.join(out_dir, "dependencies"), dirs_exist_ok=True)

params = json.loads(r'$params_string')
schema_path = "$schema"
with open(schema_path) as f:
    schema = json.load(f)["\$defs"]

env = Environment(
    loader=PackageLoader(app_dir),
    autoescape=select_autoescape()
)

rankings = {
    key: pd.read_csv(path, sep="\\t", index_col=0, usecols=[0,1], names=["TF", key], header=0)
    for key, path in {
        path[:-len(".tf_ranking.tsv")]: path
        for path in r"$tf_ranking".split(" ")
    }.items()
}

df_ranking = pd.concat(rankings.values(), axis=1)

# Remove all NaN values
tf_ranking = {
    tf: {assay: rank for assay, rank in ranks.items() if not pd.isna(rank)}
    for tf, ranks in df_ranking.to_dict(orient="index").items()
}

assays = df_ranking.columns.tolist()
sorted(assays, reverse=True)

raw_tf_tg_ranking = {
    assay: pd.read_csv(path, sep="\\t", index_col=0, header=0)
    for assay, path in {
        path[:-len(".tg_ranking.tsv")]: path
        for path in r"$tg_ranking".split(" ")
    }.items()
}

tf_tg_ranking = defaultdict(lambda: defaultdict(dict))
tg_tf_ranking = defaultdict(lambda: defaultdict(dict))
tg_ranking = defaultdict(lambda: defaultdict(int))
for assay, ranking in raw_tf_tg_ranking.items():
    for tf, genes in ranking.to_dict().items():
        for gene, dcg in genes.items():
            tf_tg_ranking[tf][gene][assay] = dcg
            tg_tf_ranking[gene][tf][assay] = dcg
            tg_ranking[gene][assay] += dcg

df_tg_ranking = pd.DataFrame(tg_ranking).T.rank(ascending=False)
df_tg_ranking = 1 - df_tg_ranking.apply(lambda x: x / x.count())

tg_ranking = {
    gene: {assay: rank for assay, rank in ranks.items() if not pd.isna(rank)}
    for gene, ranks in df_tg_ranking.to_dict(orient="index").items()
}

raw_differential = {
    pairing: pd.read_csv(path, sep="\\t", index_col=0, header=0)["log2FoldChange"].to_dict()
    for pairing, path in {
        path[:-len(".deseq2.results.tsv")]: path
        for path in "$differential".replace("\\\\", "").split(" ")
    }.items()
}

differential = defaultdict(dict)
for pairing, values in raw_differential.items():
    for gene, value in values.items():
        differential[gene][pairing] = value

pairings = list(raw_differential.keys())
sorted(pairings)

tf = env.get_template("tf.html")
tg = env.get_template("tg.html")
network = env.get_template("network.html")
snp = env.get_template("snp.html")
configuration = env.get_template("configuration.html")
styles = env.get_template("styles.css")
ranking_js = env.get_template("ranking.js")

os.makedirs(out_dir, exist_ok=True)
with open(os.path.join(out_dir, "index.html"), "w") as f:
    f.write(tf.render(tf_ranking=tf_ranking,
                      assays=assays,
                      tf_tg_ranking=tf_tg_ranking,
                      differential=differential,
                      pairings=pairings))

with open(os.path.join(out_dir, "target_genes.html"), "w") as f:
    f.write(tg.render(tg_ranking=tg_ranking,
                      assays=assays,
                      tg_tf_ranking=tg_tf_ranking,
                      differential=differential,
                      pairings=pairings))

with open(os.path.join(out_dir, "network.html"), "w") as f:
    f.write(network.render(tf_tg_ranking=tf_tg_ranking,
                           assays=assays,
                           ))

with open(os.path.join(out_dir, "snps.html"), "w") as f:
    f.write(snp.render())

with open(os.path.join(out_dir, "styles.css"), "w") as f:
    f.write(styles.render())

with open(os.path.join(out_dir, "params.json"), "w") as f:
    json.dump(params, f, indent=4)

with open(os.path.join(out_dir, "configuration.html"), "w") as f:
    f.write(configuration.render(params=params, schema=schema))

with open(os.path.join(out_dir, "ranking.js"), "w") as f:
    f.write(ranking_js.render())

with open("versions.yml", "w") as f:
    f.write('"${task.process}":\\n')
    f.write(f'  python: "{sys.executable}"\\n')
    f.write(f'  pandas: "{pd.__version__}"\\n')
    f.write(f'  jinja2: "{jinja2.__version__}"\\n')
