import argparse
import os
import sys
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--mtist-root", required=True)
parser.add_argument("--output", required=True)
args = parser.parse_args()
sys.path.insert(0, args.mtist_root)
from mtist import mtist_utils as mu
from mtist import lvsimulator

base = os.path.join(args.mtist_root, "mtist1.0")
dataset = pd.read_csv(os.path.join(base, "mtist_datasets", "dataset_37.csv"))
A = np.loadtxt(os.path.join(base, "ground_truths", "interaction_coefficients", "10_sp_aij_1.csv"), delimiter=",")
growth = np.loadtxt(os.path.join(base, "ground_truths", "growth_rates", "10_sp_gr_1.csv"), delimiter=",")
eco, gro = mu.create_lv_dicts(A, growth)
species = [f"species_{i}" for i in range(10)]
rows = []
dense_rows = []
for seed in sorted(dataset.timeseries_id.unique()):
    rng = np.random.default_rng(int(seed))
    yinit = dict(zip(species, rng.integers(1, 10, 10) / 100))
    lv = lvsimulator.LV(ecosystem=eco.copy(), growth_rates=gro.copy())
    t, y, t_all, y_all, first = lv.run_lv(
        random_seed=int(seed), tend=30, dt=0.1, yinit_specific=yinit,
        noise=0, sample_freq=100)
    frame = pd.DataFrame(y, columns=species)
    frame.insert(0, "time", t[:, 0])
    frame.insert(0, "subject", str(int(seed)))
    rows.append(frame)
    dense = pd.DataFrame(y_all, columns=species)
    dense.insert(0, "time", t_all[:, 0])
    dense.insert(0, "subject", str(int(seed)))
    dense_rows.append(dense)
out = pd.concat(rows, ignore_index=True)
dense = pd.concat(dense_rows, ignore_index=True)
os.makedirs(args.output, exist_ok=True)
out.to_csv(os.path.join(args.output, "official_noiseless_100.csv"), index=False)
dense.to_csv(os.path.join(args.output, "official_noiseless_dense.csv"), index=False)
