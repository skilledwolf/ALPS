# SPDX-License-Identifier: MIT
"""Check pyalps/seeds.py against the task rows of seed_vectors.txt.

The module is loaded from its file, so the check needs neither an installed
pyalps nor its compiled extensions. Usage: check_task_seeds.py SEEDS_PY TABLE
"""
import importlib.util
import sys


def main(seeds_py, table):
    spec = importlib.util.spec_from_file_location("pyalps_seeds", seeds_py)
    seeds = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(seeds)

    checked = 0
    failures = []
    with open(table) as rows:
        for line in rows:
            fields = line.split()
            if not fields or fields[0] != "task":
                continue
            baseseed, index, ntasks, expected = (int(x) for x in fields[1:])
            seed = seeds.task_seed(baseseed, index, ntasks)
            if seed != expected:
                failures.append(f"task_seed({baseseed}, {index}, {ntasks}) = {seed}, "
                                f"expected {expected}")
            checked += 1
    if checked == 0:
        failures.append("no task rows in " + table)
    for failure in failures:
        print(failure, file=sys.stderr)
    print(f"checked {checked} task seeds")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:]))
