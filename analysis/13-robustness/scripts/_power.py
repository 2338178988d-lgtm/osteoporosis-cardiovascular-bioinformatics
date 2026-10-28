"""Power helpers for the robustness module: how large an effect each design could detect.

Every estimate in this module is reported with the power it had to detect the reference
effect of |b| = 0.064 (the FinnGen eBMD estimate) and the smallest effect it could detect
at 80 % power. Without those two numbers a null cannot be read.
"""
import numpy as np
from scipy import stats

REFERENCE = 0.064


def power(se, true_b, alpha=0.05):
    """Two-sided power at the given standard error."""
    if not np.isfinite(se) or se <= 0:
        return np.nan
    z = stats.norm.ppf(1 - alpha / 2)
    return float(stats.norm.sf(z - abs(true_b) / se) + stats.norm.cdf(-z - abs(true_b) / se))


def mde(se, target_power=0.80, alpha=0.05):
    """Smallest |b| detectable at the given power."""
    if not np.isfinite(se) or se <= 0:
        return np.nan
    z_a = stats.norm.ppf(1 - alpha / 2)
    z_b = stats.norm.ppf(target_power)
    return float((z_a + z_b) * se)

