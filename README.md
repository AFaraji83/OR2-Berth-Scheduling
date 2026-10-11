# Berth Allocation & Scheduling (MILP)

Course project for **Operations Research 2** at the Department of Industrial Engineering, Sharif University of Technology (Spring 2025, instructor Dr. Maryam Radman).

A port has 4 berths and 8 incoming ships. The goal is to decide **which berth each ship goes to and in what order ships are unloaded**, so that the total cost of unloading time, ship waiting time and extra berth capacity is as low as possible. The problem is formulated as a mixed-integer linear program (MILP), then solved in Python (Pyomo and Gurobipy) and in GAMS, with sensitivity analysis on key parameters.

## Project phases

| Phase | What it covers | Files |
|---|---|---|
| 1. Modeling | Mathematical formulation: decision variables, parameters, objective, constraints | [`phase1-modeling/`](phase1-modeling) |
| 2. Python | Pyomo model solved with GLPK, Gurobipy version (bonus), results export, visualisations, relaxing the berth 3 rule, sensitivity analysis | [`phase2-python/`](phase2-python) |
| 3. GAMS | Same model in GAMS, plus sensitivity analysis loops (bonus) | [`phase3-gams/`](phase3-gams) |

The reports for each phase (PDF, written in Persian) are inside each phase folder.

## Problem description

- **Berths:** berths 1, 2 and 4 have 2 anchorage spots each, so they can serve at most 2 ships. Port policy says berth 3 must work at full capacity and serve **exactly 3 ships**.
- **Ships** each have a number of containers, a length, an arrival time and an unloading priority between 0 and 1 (closer to 1 means more urgent).
- A ship only fits a berth whose maximum length is at least the ship's length.
- Every berth has a container capacity. It can be increased up to a limit, at a cost per extra unit.
- Unloading time depends on the ship and the berth. Between two ships at the same berth there is a setup (lag) time.
- At a shared berth the ship with higher priority is unloaded first. A ship's waiting time is the gap between arrival and start of unloading, and it is capped at `W_max * (1 - priority)` with `W_max = 80`.
- Waiting cost is piecewise linear:

  | Waiting time W | Cost C(W) |
  |---|---|
  | 0 ≤ W ≤ 10 | 25W |
  | 10 ≤ W ≤ 40 | 10W + 150 |
  | 40 ≤ W ≤ 70 | 5W + 350 |

- Unloading costs 5 units per unit of time.

All input data lives in [`data/Parameters.xlsx`](data/Parameters.xlsx) (one sheet per parameter).

## Model summary

**Decision variables**

| Variable | Meaning |
|---|---|
| `X[i,j]` ∈ {0,1} | ship *i* is assigned to berth *j* |
| `y[i,q,j]` ∈ {0,1} | ship *i* is unloaded before ship *q* at berth *j* |
| `s[i]` ≥ 0 | start time of unloading for ship *i* |
| `W[i]` ≥ 0 | waiting time of ship *i* (`W[i] = s[i] - a[i]`) |
| `u[i,k]` ≥ 0, `f[i,b]` ∈ {0,1} | split of `W[i]` over the three cost segments (piecewise linearisation) |
| `IC[j]` ∈ ℤ⁺ | extra capacity added to berth *j* |

**Objective**

$$
\min \; 25\sum_i u_{i1} + 10\sum_i u_{i2} + 5\sum_i u_{i3} \;+\; 5\sum_{i,j} Un_{ij}\,X_{ij} \;+\; \sum_j CIC_j\,IC_j
$$

**Main constraints**

- every ship is assigned to exactly one berth; berth 3 gets exactly 3 ships and the others at most 2
- containers at a berth ≤ capacity + extra capacity, and extra capacity ≤ its allowed maximum
- ship length ≤ berth length
- `W[i] ≤ 80 (1 - P[i])`
- sequencing with big-M: if *i* goes before *q* at berth *j*, then `s[i] + Un[i,j] + O[j] ≤ s[q]`
- linking constraints so that `y` defines exactly one order for any two ships that share a berth
- piecewise segments are filled in order (`u[i,2] > 0` only if `u[i,1] = 10`, and so on)

After Phase 1 the model was revised for coding: waiting time became `W[i]` instead of `W[i,j]`, the priority variable `z[i,q]` was dropped because it duplicated `y[i,q,j]`, and a start-time variable `s[i]` was added.

## Results

Base model (solved with GLPK, confirmed with HiGHS and Gurobi): **objective = 1803.27**

| Berth | Ships (in unloading order) |
|---|---|
| 1 | 4 |
| 2 | 1 → 6 |
| 3 | 7 → 3 → 5 |
| 4 | 2 → 8 (1 unit of extra capacity) |

Ships 3, 5, 6 and 8 wait (19.3, 39.4, 20.7 and 8.2 time units). All high-priority ships start unloading as soon as they arrive.

**Relaxing the berth 3 rule** (no longer forced to take exactly 3 ships) drops the objective to **1605.56**. Ship 5 moves from berth 3 to berth 1 and its waiting time falls from 39.4 to 20.5.

**Sensitivity analysis** was run in both Pyomo and GAMS on:

1. lag time between two ships at berth 2 (5 to 10, step 0.5)
2. unloading time of ship 7 at berth 3 (5 to 15, step 0.5)
3. capacity of berth 2 (5 to 40)
4. capacity of berth 3 (5 to 40)
5. arrival time of ship 6 (1 to 5, step 0.1)

Plots and discussion are in the Phase 2 notebook and reports. GAMS results for each analysis are in [`phase3-gams/sensitivity/`](phase3-gams/sensitivity).

## Repository structure

```
.
├── data/
│   └── Parameters.xlsx              # all model inputs
├── phase1-modeling/
│   └── Phase1_Report.pdf            # mathematical model
├── phase2-python/
│   ├── berth_scheduling.py          # standalone Pyomo script
│   ├── pyomo_model.ipynb            # full Phase 2 notebook (model, plots, relaxation, sensitivity)
│   ├── gurobipy_model.ipynb         # Gurobipy version
│   ├── results/                     # exported solutions (base and berth 3 relaxed)
│   └── Phase2_Report.pdf
└── phase3-gams/
    ├── berth_scheduling.gms         # GAMS model + sensitivity loops
    ├── sensitivity/                 # results of the 5 sensitivity analyses
    └── Phase3_Report.pdf
```

## How to run

### Python

```bash
pip install -r requirements.txt
```

You also need a MILP solver. The notebooks use **GLPK** (`conda install -c conda-forge glpk`, or `apt install glpk-utils` on Linux). **HiGHS** is installed by `requirements.txt` and works without any extra setup.

```bash
cd phase2-python
python berth_scheduling.py --solver appsi_highs           # base model
python berth_scheduling.py --solver appsi_highs --relax-berth-3
python berth_scheduling.py --solver glpk --out results/schedule.xlsx
```

Or open `pyomo_model.ipynb` / `gurobipy_model.ipynb` in Jupyter. Run them from inside `phase2-python/` so the relative path to `../data/Parameters.xlsx` works.

### GAMS

Open `phase3-gams/berth_scheduling.gms` in GAMS Studio and run it. The data is written directly in the file. The script solves the base model, then loops over each sensitivity scenario and writes the results to `Sensitivity_sa1.gdx` … `Sensitivity_sa5.gdx`. The Excel files in `sensitivity/` were exported from those GDX files.

## Tools

Python · Pyomo · GLPK · Gurobi (gurobipy) · pandas · matplotlib · seaborn · GAMS

## Team

- Alireza Faraji
- Mohammad Parsa Keshmiri
- Mohammad Erfan Kalmatian
- Amyar Siosemardeh
