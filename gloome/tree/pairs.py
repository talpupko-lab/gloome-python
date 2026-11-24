"""
Memory-bounded, vectorized engine for all-pairs site correlations and their simulated null.

The number of site pairs grows quadratically: 22,820 candidate sites already give 260 M pairs, and
100 simulated datasets give ~15 billion null correlations. Nothing here materializes a full
correlation matrix or a per-pair Python object. Correlations are computed block by block with BLAS,
and the null is kept either exactly (small inputs) or as a fixed-size per-rate-bin histogram (large
inputs, constant memory regardless of input size).

Ties: identical site patterns give correlations that are equal in exact arithmetic but differ by
~1e-15 in floating point (e.g. r = 1 against a null spike at 1 - 2e-16). "At least as extreme" is
therefore |null| >= |r| - TIE_TOLERANCE; a plain >= would drop tied null values at random and make
p-values anti-conservative exactly where significant pairs are called.

Rate-bin codes: a pair of sites with rate categories (a, b), a <= b, has code ``a * K + b`` where K
is the number of rate categories; its label is the tuple ``(a, b)``.
"""
from typing import Dict, Iterator, List, Optional, Tuple

import numpy as np
from scipy.stats import beta as sp_beta, kstwobign

EPS2 = 1e-10
TIE_TOLERANCE = 1e-12

# Histogram null: total cells over all rate-bins (int64 -> 256 MB). The |r| resolution 1 / H is the
# finest power of two that fits: 2^-20 (~1e-6) for 4 rate categories, 2^-16 for 16.
HIST_CELLS_BUDGET = 1 << 25
MAX_HIST_BINS = 1 << 20
MIN_HIST_BINS = 1 << 12
# Above this many null values the null switches from exact arrays to histograms (~8 bytes / value).
EXACT_NULL_LIMIT = 20_000_000
# Above this many observed pairs the pairs are streamed in two passes instead of held in memory.
PAIRS_IN_MEMORY_LIMIT = 2_000_000
# Maximum rows written to Coevolution.tsv when pairs are streamed: the most significant pairs are kept
# (q-values are still computed over all pairs). 2 M rows ~ 130 MB of TSV; each extra million rows
# costs ~0.15 GB of peak memory.
MAX_TABLE_ROWS = 2_000_000
# Correlation values computed per BLAS block (bounds the working memory to a few hundred MB).
BLOCK_ELEMENTS = 1 << 22
# Points used for plots and for the descriptive Beta fit of the observed correlations.
PLOT_SAMPLE = 200_000
STRIPPLOT_SAMPLE = 20_000


def centered_rows(matrix: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Row-centered matrix and its row norms: corr(i, j) = c_i . c_j / (|c_i| |c_j|)."""
    centered = np.asarray(matrix, dtype=np.float64)
    centered = centered - centered.mean(axis=1, keepdims=True)
    norms = np.sqrt((centered ** 2).sum(axis=1))

    return centered, norms


def _block_rows(n_columns: int, block_elements: int) -> int:
    return max(1, block_elements // max(1, n_columns))


def iter_upper_pairs(centered: np.ndarray, norms: np.ndarray, block_elements: int = BLOCK_ELEMENTS
                     ) -> Iterator[Tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """
    Yield (i, j, r) for every pair i < j of rows, in the order of ``np.triu_indices(n, k=1)``.
    i and j are int64 row indices into ``centered``; r is the Pearson correlation.
    """
    n = centered.shape[0]
    start = 0
    while start < n - 1:
        stop = min(n - 1, start + _block_rows(n - start, block_elements))
        block = centered[start:stop] @ centered[start:].T
        block /= norms[start:stop, None] * norms[None, start:]
        ii, jj = np.nonzero(np.arange(n - start)[None, :] > np.arange(stop - start)[:, None])
        yield ii + start, jj + start, block[ii, jj]
        start = stop


def iter_bin_blocks(centered: np.ndarray, norms: np.ndarray, categories: np.ndarray, categories_quantity: int,
                    block_elements: int = BLOCK_ELEMENTS) -> Iterator[Tuple[int, np.ndarray]]:
    """
    Yield (code, r_values) for every pair of rows, grouped by rate-bin. Rows are grouped by rate
    category so each BLAS block belongs to a single rate-bin: no per-pair index or code arrays.
    """
    groups = [np.flatnonzero(categories == c) for c in range(categories_quantity)]
    for a in range(categories_quantity):
        rows_a = groups[a]
        if not rows_a.size:
            continue
        za, na = centered[rows_a], norms[rows_a]
        for b in range(a, categories_quantity):
            rows_b = groups[b]
            if not rows_b.size or (a == b and rows_b.size < 2):
                continue
            zb, nb = (za, na) if a == b else (centered[rows_b], norms[rows_b])
            code = a * categories_quantity + b
            n_b = zb.shape[0]
            start = 0
            last_row = n_b - 1 if a == b else za.shape[0]
            while start < last_row:
                col0 = start if a == b else 0
                stop = min(last_row, start + _block_rows(n_b - col0, block_elements))
                block = za[start:stop] @ zb[col0:].T
                block /= na[start:stop, None] * nb[None, col0:]
                if a == b:
                    block = block[np.arange(n_b - col0)[None, :] > np.arange(stop - start)[:, None]]
                yield code, block.ravel()
                start = stop


def pair_codes(categories_i: np.ndarray, categories_j: np.ndarray, categories_quantity: int) -> np.ndarray:
    low = np.minimum(categories_i, categories_j).astype(np.int32)
    high = np.maximum(categories_i, categories_j).astype(np.int32)

    return low * categories_quantity + high


def code_label(code: int, categories_quantity: int) -> Tuple[int, int]:
    return int(code // categories_quantity), int(code % categories_quantity)


def histogram_bins(categories_quantity: int) -> int:
    """Cells per side of zero for the histogram null of ``categories_quantity`` rate categories."""
    used_codes = categories_quantity * (categories_quantity + 1) // 2
    bins = 1 << int(np.floor(np.log2(max(1, HIST_CELLS_BUDGET // (2 * used_codes)))))

    return int(min(MAX_HIST_BINS, max(MIN_HIST_BINS, bins)))


def abs_index(r: np.ndarray, hist_bins: int) -> np.ndarray:
    """Histogram cell of |r| - TIE_TOLERANCE: cell h holds |r| in [h / H, (h + 1) / H)."""
    h = np.floor((np.abs(r) - TIE_TOLERANCE) * hist_bins).astype(np.int64)
    np.clip(h, 0, hist_bins - 1, out=h)

    return h


class NullDistribution:
    """
    Simulated null correlations per rate-bin code.

    Exact mode keeps every value (identical p-values to an unbounded list) until ``exact_limit``
    values; then everything is folded into per-code histograms of signed r: 2 * H cells of width
    1 / H over [-1, 1]. Memory is then fixed at n_codes * 2H * 8 bytes whatever the input size.
    """

    def __init__(self, categories_quantity: int, exact_limit: int = EXACT_NULL_LIMIT,
                 hist_bins: Optional[int] = None):
        self.categories_quantity = categories_quantity
        self.n_codes = categories_quantity * categories_quantity
        self.exact_limit = exact_limit
        self.hist_bins = hist_bins or histogram_bins(categories_quantity)
        # Histogram rows exist only for valid codes a * K + b with a <= b.
        valid = [a * categories_quantity + b for a in range(categories_quantity)
                 for b in range(a, categories_quantity)]
        self.valid_codes = np.asarray(valid, dtype=np.int64)
        self.row_of_code = np.full(self.n_codes, -1, dtype=np.int64)
        self.row_of_code[self.valid_codes] = np.arange(len(valid))
        self.exact = True
        self.size = 0
        self._values: Dict[int, List[np.ndarray]] = {}
        self.counts: Optional[np.ndarray] = None
        if exact_limit <= 0:
            self._to_histogram()

    def add(self, code: int, r: np.ndarray) -> None:
        if not r.size:
            return
        if self.exact and self.size + r.size > self.exact_limit:
            self._to_histogram()
        if self.exact:
            self._values.setdefault(code, []).append(np.array(r, dtype=np.float64))
        else:
            self._add_to_histogram(code, r)
        self.size += r.size

    def _signed_index(self, r: np.ndarray) -> np.ndarray:
        k = np.floor((np.asarray(r, dtype=np.float64) + 1.0) * self.hist_bins).astype(np.int64)
        np.clip(k, 0, 2 * self.hist_bins - 1, out=k)

        return k

    def _add_to_histogram(self, code: int, r: np.ndarray) -> None:
        self.counts[self.row_of_code[code]] += np.bincount(self._signed_index(r), minlength=2 * self.hist_bins)

    def _to_histogram(self) -> None:
        self.exact = False
        self.counts = np.zeros((len(self.valid_codes), 2 * self.hist_bins), dtype=np.int64)
        for code, chunks in self._values.items():
            for chunk in chunks:
                self._add_to_histogram(code, chunk)
        self._values = {}

    def samples(self) -> Dict[int, np.ndarray]:
        """Exact mode only: signed null values per code."""
        return {code: np.concatenate(chunks) for code, chunks in self._values.items()}

    # ----- histogram-mode helpers ---------------------------------------------------------------
    def tail_counts(self) -> np.ndarray:
        """(valid codes, H): number of null values with |r| >= h / H, for every cell h."""
        h = self.hist_bins
        positive = self.counts[:, h:]
        negative = self.counts[:, :h][:, ::-1]

        return (np.cumsum(positive[:, ::-1], axis=1)[:, ::-1] +
                np.cumsum(negative[:, ::-1], axis=1)[:, ::-1])

    def empirical_p_table(self) -> np.ndarray:
        """
        (n_codes, H), indexed by code: (1 + #null with |r| >= h / H) / (1 + n). Conservative: null
        values inside the observed value's own cell (width 1 / H) count as at least as extreme.
        """
        totals = self.counts.sum(axis=1).astype(np.float64)
        rows = (1.0 + self.tail_counts()) / (1.0 + totals[:, None])
        rows[totals == 0] = 1.0
        table = np.ones((self.n_codes, self.hist_bins), dtype=np.float64)
        table[self.valid_codes] = rows

        return table

    def fit_beta(self, code: Optional[int] = None, min_fit: int = 100, gof_alpha: float = 0.01
                 ) -> Optional[Tuple[float, float]]:
        """Histogram analogue of ``Tree.fit_beta_null`` (``code=None`` pools every code)."""
        counts = self.counts.sum(axis=0) if code is None else self.counts[self.row_of_code[code]]

        return fit_beta_histogram(counts, self.hist_bins, min_fit=min_fit, gof_alpha=gof_alpha)


def fit_beta_histogram(counts: np.ndarray, hist_bins: int, min_fit: int = 100, gof_alpha: float = 0.01
                       ) -> Optional[Tuple[float, float]]:
    """
    Method-of-moments Beta(a, b, loc=-1, scale=2) fit from a signed-r histogram, with the same
    rejection rules as ``Tree.fit_beta_null``. The Kolmogorov-Smirnov distance is measured at
    the cell edges (asymptotic distribution, n is large in histogram mode).
    """
    counts = np.asarray(counts, dtype=np.float64)
    n = counts.sum()
    if n < max(min_fit, 2):
        return None
    edges = -1.0 + np.arange(2 * hist_bins + 1) / hist_bins
    mids = np.clip(((edges[:-1] + edges[1:]) / 2 + 1.0) / 2.0, EPS2, 1 - EPS2)
    mean = (counts * mids).sum() / n
    var = (counts * (mids - mean) ** 2).sum() / (n - 1)
    max_var = mean * (1 - mean)
    if not (0.0 < var < max_var):
        return None
    common = max_var / var - 1.0
    a, b = mean * common, (1 - mean) * common
    if not (np.isfinite(a) and np.isfinite(b) and a > 0 and b > 0):
        return None
    ecdf = np.concatenate(([0.0], np.cumsum(counts) / n))
    distance = float(np.max(np.abs(ecdf - sp_beta.cdf(edges, a, b, loc=-1, scale=2))))
    p_gof = float(kstwobign.sf(distance * np.sqrt(n)))

    return None if p_gof < gof_alpha else (float(a), float(b))


def beta_tail(abs_r: np.ndarray, params: Tuple[float, float]) -> np.ndarray:
    """Vectorized ``Tree.beta_p``: P(|R| >= |r|) under Beta(a, b, loc=-1, scale=2)."""
    a, b = params
    x = np.minimum(np.abs(abs_r), 1 - EPS2)

    return sp_beta.cdf(-x, a, b, loc=-1, scale=2) + sp_beta.sf(x, a, b, loc=-1, scale=2)


def empirical_p_exact(abs_r: np.ndarray, null_values: Optional[np.ndarray]) -> np.ndarray:
    """(1 + #|null| >= |r|) / (1 + n) for every r, by binary search (no per-pair loop)."""
    if null_values is None or not null_values.size:
        return np.ones(abs_r.shape, dtype=np.float64)
    sorted_abs = np.sort(np.abs(null_values))
    at_least = sorted_abs.size - np.searchsorted(sorted_abs, np.asarray(abs_r) - TIE_TOLERANCE, side='left')

    return (1.0 + at_least) / (1.0 + sorted_abs.size)


def benjamini_hochberg_weighted(p_values: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """
    BH q-values when p_values[k] is shared by weights[k] tests. Identical to running BH on the
    expanded list of sum(weights) p-values (tied p-values get the same q either way).
    """
    order = np.argsort(p_values, kind='stable')
    ranked = p_values[order]
    rank = np.cumsum(weights[order]).astype(np.float64)
    q_sorted = np.minimum.accumulate((ranked * rank[-1] / rank)[::-1])[::-1]
    q = np.empty_like(q_sorted)
    q[order] = np.clip(q_sorted, 0, 1)

    return q
