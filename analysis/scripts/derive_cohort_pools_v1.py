"""Derive durable all-cohort AF/HF pooled estimates from the frozen submission table.

Input:
    .cache/bone-hf-submission-v1/results/tables/14-cardiac-ladder/big_cohorts_all.csv
Output:
    results/tables/14-cardiac-ladder/big_cohorts_pooled_all_v1.csv

This script does not download data or refit individual cohort estimates. It only
persists the random-effects pooling already used by the submission-package figure
script, so the headline pooled values have a durable table source.
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / ".cache" / "bone-hf-submission-v1" / "results" / "tables" / "14-cardiac-ladder" / "big_cohorts_all.csv"
OUTPUT = ROOT / "results" / "tables" / "14-cardiac-ladder" / "big_cohorts_pooled_all_v1.csv"


def pool(beta: np.ndarray, se: np.ndarray) -> dict[str, float]:
    beta = np.asarray(beta, dtype=float)
    se = np.asarray(se, dtype=float)
    weights = 1.0 / se**2
    fixed = float(np.sum(weights * beta) / np.sum(weights))
    q = float(np.sum(weights * (beta - fixed) ** 2))
    df = len(beta) - 1
    c = float(np.sum(weights) - np.sum(weights**2) / np.sum(weights))
    tau2 = max(0.0, (q - df) / c) if df > 0 and c > 0 else 0.0
    random_weights = 1.0 / (se**2 + tau2)
    pooled = float(np.sum(random_weights * beta) / np.sum(random_weights))
    pooled_se = float(np.sqrt(1.0 / np.sum(random_weights)))
    i2 = max(0.0, (q - df) / q) if q > 0 and df > 0 else 0.0
    p = float(2.0 * stats.norm.sf(abs(pooled / pooled_se)))
    return {
        "beta": pooled,
        "se": pooled_se,
        "lo": pooled - 1.96 * pooled_se,
        "hi": pooled + 1.96 * pooled_se,
        "p": p,
        "q": q,
        "q_df": df,
        "q_p": float(stats.chi2.sf(q, df)) if df > 0 else np.nan,
        "tau2": tau2,
        "i2": i2,
    }


def main() -> None:
    if not INPUT.exists():
        raise FileNotFoundError(f"Frozen input not found: {INPUT}")
    data = pd.read_csv(INPUT)
    rows = []
    for endpoint in ("AF", "HF"):
        subset = data.loc[data["endpoint"].eq(endpoint)].copy()
        if subset.empty:
            raise ValueError(f"No rows for endpoint {endpoint}")
        result = pool(subset["b"].to_numpy(), subset["se"].to_numpy())
        result.update(
            {
                "endpoint": endpoint,
                "k": int(len(subset)),
                "cases_total_reported": int(subset["label"].str.extract(r"\(([^,]+)", expand=False).str.replace(",", "", regex=False).astype(int).sum()),
                "input_table": INPUT.relative_to(ROOT).as_posix(),
                "pooling_method": "random-effects inverse-variance pooling; DerSimonian-Laird tau2",
            }
        )
        rows.append(result)
    output = pd.DataFrame(rows)[
        [
            "endpoint",
            "k",
            "cases_total_reported",
            "beta",
            "se",
            "lo",
            "hi",
            "p",
            "i2",
            "q",
            "q_df",
            "q_p",
            "tau2",
            "input_table",
            "pooling_method",
        ]
    ]
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(OUTPUT, index=False)
    print(output.to_string(index=False))
    print(f"saved {OUTPUT.relative_to(ROOT).as_posix()}")


if __name__ == "__main__":
    main()

