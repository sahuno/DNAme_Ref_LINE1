##compare geometric and arithmeric mean; simulations

from scipy.stats import norm, gmean
import numpy as np
import matplotlib.pyplot as plt

var = [10, 5, 2, 1]
nDist = 4

# Create a list of 4 normal distributions
dists = list()
for i in range(nDist):
    sample = norm.rvs(5, var[i], size=10)
    dists.append(sample)
    # Compute arithmetic mean
    arithmetic_mean = np.mean(sample)
    # Try to compute geometric mean; if any sample value is non-positive, gmean will not work.
    try:
        geometric_mean = gmean(sample)
    except ValueError:
        geometric_mean = np.nan  # or you can handle it as you wish
    print(f"Distribution {i+1}:")
    print("  Sample:", sample)
    print("  Arithmetic Mean:", arithmetic_mean)
    print("  Geometric Mean:", geometric_mean)
    print()
# Plot the distributions