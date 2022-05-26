import tskit
import numpy as np

tsChr = tskit.load("Chr16.trees")
window_breaks = [10, 100]
focal=None

def pairwise_coalescence_counts(ts, window_breaks, focal=None):
    # Work in progress
    # Calculate counts of pairwise coalescence events within
    # time windows for a set of pairs of tips.
    # Also returns the number of uncoalesced "pairs" at the beginning of
    # each window (e.g. number of trees for which both members of pair are
    # not missing and not coalesced).
    #@@ -1597,86 +1598,83 @@ def pairwise_coalescence_counts(ts, window_breaks, focal=None):
    if focal is None:
        focal = []
        for i in range(ts.num_samples - 1):
            for j in range(i + 1, ts.num_samples):
                focal += [(i, j)]
    else:
        if not np.all([len(x) == 2 for x in focal]):
            raise Exception("'focal' must contain tuples of pairs of sample indices")

    window_breaks = np.sort(np.unique(np.append(window_breaks, [0.0])))
    if np.any(window_breaks < 0.0):
        raise Exception("'window_breaks' must be non-negative")
    num_windows = len(window_breaks) - 1
    windows = [(window_breaks[i], window_breaks[i + 1]) for i in range(num_windows)]

    counts = {}
    for pair in focal:
        # TMRCAs for each pair
        i, j = pair
        pairwise_tmrca = np.sort([t.tmrca(i, j) for t in ts.trees()])
        not_missing = np.invert(np.isnan(pairwise_tmrca))
        nonmissing_trees = np.sum(not_missing)
        pairwise_tmrca = pairwise_tmrca[not_missing]
        # TODO: How is missing data handled by tree.tmrca()? Does it return np.nan?

        # tally coalesced trees during window, and the number of uncoalesced
        # trees at the start of the window per pair
        coalesced, _ = np.histogram(pairwise_tmrca, window_breaks)
        total = nonmissing_trees - np.append([0], np.cumsum(coalesced))[:num_windows]
        counts[pair] = {"coalesced": coalesced, "total": total}

    return {"counts": counts, "windows": windows}


def pairwise_coalescence_rates(raw_counts):
    # Work in progress.
    # Calculate effective population size from counts
    # of coalescing lineage pairs in time windows

    def ne_from_counts(x, k, d):
        # estimate haploid Ne from counts
        p = np.where(k == 0.0, np.nan, x / k)
        q = 1.0 - p
        n = np.where(q == 0.0, np.nan, -d / np.log(q))
        # standard error assuming binomial(x; p, k)
        se = 1.0 / np.sqrt(
            np.where(
                q == 0.0,
                np.nan,
                -2.0 * k * np.log(q) / n ** 2
                + x * (1.0 + q / p) * np.log(q) / n ** 2 * (2.0 + q / p * np.log(q)),
                )
        )
        return {"mle": n, "std_err": se}

    windows = raw_counts["windows"]
    counts = raw_counts["counts"]
    window_duration = np.array([x[1] - x[0] for x in windows])
    num_windows = len(windows)

    # calculate Ne estimates from pairwise counts,
    # and from counts summed over pairs
    pairwise_estimates = {}
    coalesced = np.zeros(num_windows)
    total = np.zeros(num_windows)
    for pair in counts.keys():
        pairwise_estimates[pair] = ne_from_counts(
            counts[pair]["coalesced"], counts[pair]["total"], window_duration
        )
        coalesced += counts[pair]["coalesced"]
        total += counts[pair]["total"]
    global_estimate = ne_from_counts(coalesced, total, window_duration)

    return {
        "global_estimate": global_estimate,
        "pairwise_estimates": pairwise_estimates,
    }

first = pairwise_coalescence_counts(ts = tsChr, window_breaks = [0, 10, 100])
second = pairwise_coalescence_rates(first)