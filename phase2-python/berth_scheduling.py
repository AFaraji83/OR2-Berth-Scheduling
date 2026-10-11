"""
Berth allocation and scheduling MILP (OR2 course project, Phase 2).

Assigns 8 ships to 4 berths and sequences them so that the total cost of
unloading time, piecewise-linear waiting cost and extra berth capacity is
minimised. This is a standalone version of the model in `pyomo_model.ipynb`.

Usage:
    python berth_scheduling.py                    # base model, GLPK
    python berth_scheduling.py --solver appsi_highs
    python berth_scheduling.py --relax-berth-3    # berth 3 no longer needs exactly 3 ships
"""

import argparse
import ast
from pathlib import Path

import pandas as pd
import pyomo.environ as pyo

DATA_FILE = Path(__file__).resolve().parent.parent / "data" / "Parameters.xlsx"

N_SHIPS = 8
N_BERTHS = 4
BIG_M = 80          # also the global upper bound on waiting time (W_max)


def read_vector(sheet, column):
    """Read a two-column sheet (index, value) into a dict."""
    df = pd.read_excel(DATA_FILE, sheet_name=sheet, index_col=0)
    return df[column].to_dict()


def build_model(relax_berth_3=False):
    m = pyo.ConcreteModel(name="berth_scheduling")

    # ---------------- Sets ----------------
    m.i = pyo.RangeSet(1, N_SHIPS)                        # ships
    m.q = pyo.RangeSet(1, N_SHIPS)                        # alias of ships
    m.iq = pyo.Set(dimen=2, initialize=[(i, q) for i in m.i for q in m.i if i != q])
    m.j = pyo.RangeSet(1, N_BERTHS)                       # berths
    m.k = pyo.RangeSet(1, 3)                              # segments of the waiting-cost function
    m.bin = pyo.RangeSet(1, 2)                            # segment activation binaries

    # ---------------- Parameters ----------------
    unload = pd.read_excel(DATA_FILE, sheet_name="Unload Time")
    unload["Ship_Berth"] = unload["Ship_Berth"].apply(lambda x: tuple(map(int, ast.literal_eval(x))))
    unload_dict = dict(zip(unload["Ship_Berth"], unload["Unload_Time"]))

    m.C = pyo.Param(m.j, initialize=read_vector("Capacity", "Capacity"), mutable=True)          # berth capacity
    m.CA = pyo.Param(m.i, initialize=read_vector("Containers", "Containers"), mutable=True)     # containers per ship
    m.L = pyo.Param(m.j, initialize=read_vector("Max Length", "Max_Length"), mutable=True)      # max ship length per berth
    m.LS = pyo.Param(m.i, initialize=read_vector("Length", "Length"), mutable=True)             # ship length
    m.P = pyo.Param(m.i, initialize=read_vector("Priority", "Priority"), mutable=True)          # unloading priority (0-1)
    m.O = pyo.Param(m.j, initialize=read_vector("Lag", "Lag"), mutable=True)                    # setup time between two ships
    m.Un = pyo.Param(m.i, m.j, initialize=unload_dict, mutable=True)                            # unloading time
    m.a = pyo.Param(m.i, initialize=read_vector("Arrival Time", "Arrival_Time"), mutable=True)  # arrival time
    m.CIC = pyo.Param(m.j, initialize=read_vector("Cost for Extra Capacity", "Cost for Extra Capacity"), mutable=True)
    m.ext = pyo.Param(m.j, initialize=read_vector("Max Extra Allowable Capacity", "Max Extra Allowable Capacity"), mutable=True)

    # ---------------- Variables ----------------
    m.X = pyo.Var(m.i, m.j, domain=pyo.Binary)            # ship i assigned to berth j
    m.W = pyo.Var(m.i, domain=pyo.NonNegativeReals)       # waiting time of ship i
    m.s = pyo.Var(m.i, domain=pyo.NonNegativeReals)       # start of unloading for ship i
    m.IC = pyo.Var(m.j, domain=pyo.NonNegativeIntegers)   # extra capacity added to berth j
    m.y = pyo.Var(m.iq, m.j, domain=pyo.Binary)           # ship i is unloaded before ship q at berth j
    m.u = pyo.Var(m.i, m.k, domain=pyo.NonNegativeReals)  # waiting time falling in cost segment k
    m.f = pyo.Var(m.i, m.bin, domain=pyo.Binary)          # segment k+1 of the cost function is active

    # ---------------- Objective ----------------
    # Waiting cost  C(W) = 25W (0-10), 10W+150 (10-40), 5W+350 (40-70)  -> slopes 25 / 10 / 5
    # + 5 cost units per unit of unloading time + cost of extra capacity
    m.obj = pyo.Objective(
        expr=25 * pyo.quicksum(m.u[i, 1] for i in m.i)
        + 10 * pyo.quicksum(m.u[i, 2] for i in m.i)
        + 5 * pyo.quicksum(m.u[i, 3] for i in m.i)
        + 5 * pyo.quicksum(m.Un[i, j] * m.X[i, j] for i in m.i for j in m.j)
        + pyo.quicksum(m.CIC[j] * m.IC[j] for j in m.j),
        sense=pyo.minimize,
    )

    # ---------------- Constraints ----------------
    # Each ship goes to exactly one berth
    m.c1 = pyo.Constraint(m.i, rule=lambda m, i: pyo.quicksum(m.X[i, j] for j in m.j) == 1)

    # Berth 3 must serve exactly 3 ships (unless relaxed); the others at most 2
    def ships_per_berth(m, j):
        n_assigned = pyo.quicksum(m.X[i, j] for i in m.i)
        if j == 3:
            return n_assigned <= 3 if relax_berth_3 else n_assigned == 3
        return n_assigned <= 2
    m.c2 = pyo.Constraint(m.j, rule=ships_per_berth)

    # Container capacity, with optional (bounded) extra capacity
    m.c3 = pyo.Constraint(m.j, rule=lambda m, j: pyo.quicksum(m.X[i, j] * m.CA[i] for i in m.i) <= m.C[j] + m.IC[j])
    m.c4 = pyo.Constraint(m.j, rule=lambda m, j: m.IC[j] <= m.ext[j])

    # Ship must fit the berth length
    m.c5 = pyo.Constraint(m.i, m.j, rule=lambda m, i, j: m.LS[i] * m.X[i, j] <= m.L[j])

    # Priority-dependent upper bound on waiting time: W_i <= W_max * (1 - priority_i)
    m.c6 = pyo.Constraint(m.i, rule=lambda m, i: m.W[i] <= BIG_M * (1 - m.P[i]))

    # Waiting time = start time - arrival time
    m.c7 = pyo.Constraint(m.i, rule=lambda m, i: m.W[i] == m.s[i] - m.a[i])

    # If i precedes q at berth j, q starts after i is unloaded plus the setup lag
    m.c8 = pyo.Constraint(
        m.iq, m.j,
        rule=lambda m, i, q, j: m.s[i] + m.Un[i, j] + m.O[j] <= m.s[q] + BIG_M * (1 - m.y[i, q, j]),
    )

    # Ordering logic: exactly one order for two ships sharing a berth, none otherwise
    m.c9 = pyo.Constraint(m.iq, m.j, rule=lambda m, i, q, j: m.y[i, q, j] + m.y[q, i, j] <= 1)
    m.c10 = pyo.Constraint(m.iq, m.j, rule=lambda m, i, q, j: m.y[i, q, j] + m.y[q, i, j] >= m.X[i, j] + m.X[q, j] - 1)
    m.c11 = pyo.Constraint(m.iq, m.j, rule=lambda m, i, q, j: m.y[i, q, j] + m.y[q, i, j] <= m.X[i, j])

    # A ship that is first in line does not wait (berth 3 holds up to 3 ships)
    def wait_zero_if_first(m, i, q, j):
        if j == 3:
            n_after = pyo.quicksum(m.y[i, qq, j] for qq in m.q if qq != i)
            return m.W[i] <= BIG_M * (2 - n_after)
        return m.W[i] <= BIG_M * (1 - m.y[i, q, j])
    m.c13 = pyo.Constraint(m.iq, m.j, rule=wait_zero_if_first)

    # Piecewise-linear waiting cost: W_i split into segments of length 10, 30, 30
    m.c14 = pyo.Constraint(m.i, rule=lambda m, i: m.W[i] == pyo.quicksum(m.u[i, k] for k in m.k))
    seg_len = {1: 10, 2: 30}
    m.c15 = pyo.Constraint(m.i, m.bin, rule=lambda m, i, b: seg_len[b] * m.f[i, b] <= m.u[i, b])

    def segment_upper(m, i, k):
        if k == 1:
            return m.u[i, k] <= 10
        return m.u[i, k] <= 30 * m.f[i, k - 1]
    m.c16 = pyo.Constraint(m.i, m.k, rule=segment_upper)
    m.c17 = pyo.Constraint(m.i, rule=lambda m, i: m.f[i, 1] >= m.f[i, 2])

    return m


def results_table(m):
    rows = []
    for i in m.i:
        berth = next(j for j in m.j if pyo.value(m.X[i, j]) > 0.5)
        start = pyo.value(m.s[i])
        rows.append({
            "Ship": i,
            "Berth": berth,
            "Priority": round(pyo.value(m.P[i]), 3),
            "Arrival": round(pyo.value(m.a[i]), 2),
            "Start": round(start, 2),
            "Finish": round(start + pyo.value(m.Un[i, berth]), 2),
            "Wait": round(pyo.value(m.W[i]), 2),
        })
    return pd.DataFrame(rows).sort_values(["Berth", "Start"])


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--solver", default="glpk", help="Pyomo solver name (glpk, appsi_highs, gurobi, cbc, ...)")
    parser.add_argument("--relax-berth-3", action="store_true", help="drop the 'exactly 3 ships at berth 3' rule")
    parser.add_argument("--out", help="optional path to save the schedule as .xlsx")
    args = parser.parse_args()

    model = build_model(relax_berth_3=args.relax_berth_3)
    result = pyo.SolverFactory(args.solver).solve(model)
    if result.solver.termination_condition != pyo.TerminationCondition.optimal:
        raise SystemExit(f"Solver finished with status: {result.solver.termination_condition}")

    schedule = results_table(model)
    extra = {j: int(round(pyo.value(model.IC[j]))) for j in model.j}
    print(f"Objective value: {pyo.value(model.obj):.4f}\n")
    print(schedule.to_string(index=False))
    print(f"\nExtra capacity per berth: {extra}")

    if args.out:
        schedule.to_excel(args.out, index=False)
        print(f"Schedule saved to {args.out}")


if __name__ == "__main__":
    main()
