from __future__ import annotations

import warnings
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
import pandas as pd
import scipy.optimize
import scipy.stats as st
from scipy.ndimage import convolve1d
from scipy.special import logsumexp, ndtr, ndtri_exp

from genome_tools.plotting.modular_plot.api import DataBundle, PlotDataLoader
from footprint_tools.stats import differential
from footprint_tools.stats.distributions import invchi2
from collections.abc import Mapping
from pathlib import Path
from typing import Any


__all__ = [
    "DifferentialModelConfig",
    "DifferentialFit",
    "DifferentialLikelihood",
    "DifferentialDiagnostics",
    "WindowedScore",
    "DifferentialAnalysis",
    "DifferentialFootprintModel",
    "WindowedStoufferScorer",
    "DifferentialFitLoader",
    "WindowedScoreLoader",
]


@dataclass(frozen=True, slots=True)
class DifferentialModelConfig:
    theta_x: np.ndarray | None = None
    mu_x: np.ndarray | None = None
    sig2_x: np.ndarray | None = None
    theta_step: float = 0.1
    theta_tail_sd: float = 5.0
    min_theta_mass: float = 0.999
    base_chunk_size: int = 32
    max_boundary_fraction: float = 0.01


@dataclass(frozen=True, slots=True)
class DifferentialFit:
    """Pointwise fit arrays; group-dependent arrays are group x position."""

    statistic: np.ndarray
    df: int
    common_mu: np.ndarray
    group_mu: np.ndarray
    group_conditional_sig2: np.ndarray

    @property
    def llr(self) -> np.ndarray:
        return self.statistic / 2.0

    @property
    def log_pvalue(self) -> np.ndarray:
        return st.chi2.logsf(self.statistic, self.df)

    @property
    def pvalue(self) -> np.ndarray:
        floor = np.log(np.finfo(np.float64).tiny)
        return np.exp(np.maximum(self.log_pvalue, floor))

    @property
    def neglog10_pvalue(self) -> np.ndarray:
        return -self.log_pvalue / np.log(10.0)

    @property
    def group_conditional_sigma(self) -> np.ndarray:
        return np.sqrt(self.group_conditional_sig2)


@dataclass(frozen=True, slots=True)
class DifferentialLikelihood:
    """Likelihood grids retained from the fit.

    Attributes
    ----------
    mu_x
        Mean grid, shape ``mu``.
    sig2_x
        Variance grid, shape ``sig2``.
    group_loglik_mu
        Sigma-squared-integrated log likelihood,
        shape ``group x position x mu``.
    group_log_sig2_prior_mass
        Normalized finite-grid inverse-chi-squared log prior masses,
        shape ``group x sig2``. These include sigma-squared quadrature weights.
    group_loglik_mu_sig2
        Optional raw log likelihood before the variance prior is applied,
        shape ``group x position x mu x sig2``. Theta is already integrated
        and sample log likelihoods are summed within each group.
    """

    mu_x: np.ndarray
    sig2_x: np.ndarray
    group_loglik_mu: np.ndarray
    group_log_sig2_prior_mass: np.ndarray
    group_loglik_mu_sig2: np.ndarray | None = None

    @property
    def sigma_x(self) -> np.ndarray:
        return np.sqrt(self.sig2_x)

    @property
    def group_sig2_prior_mass(self) -> np.ndarray:
        return np.exp(self.group_log_sig2_prior_mass)

    @property
    def common_loglik_mu(self) -> np.ndarray:
        """Common-mean log likelihood, shape position x mu."""
        return self.group_loglik_mu.sum(axis=0)

    @property
    def raw_surface_nbytes(self) -> int:
        value = self.group_loglik_mu_sig2
        return 0 if value is None else value.nbytes


@dataclass(frozen=True, slots=True)
class DifferentialDiagnostics:
    theta_x: np.ndarray
    group_invchi2_params: np.ndarray  # group x (nu, prior_sig2)
    theta_coverage: float
    mu_boundary_fraction: float


@dataclass(frozen=True, slots=True)
class WindowedScore:
    """Local Stouffer score with nominal p-values under independence."""

    radius: int
    score: np.ndarray

    @property
    def nominal_log_pvalue(self) -> np.ndarray:
        return st.norm.logsf(self.score)

    @property
    def nominal_pvalue(self) -> np.ndarray:
        floor = np.log(np.finfo(np.float64).tiny)
        return np.exp(np.maximum(self.nominal_log_pvalue, floor))

    @property
    def nominal_neglog10_pvalue(self) -> np.ndarray:
        return -self.nominal_log_pvalue / np.log(10.0)


@dataclass(slots=True)
class DifferentialAnalysis:
    group_names: tuple[str, ...]
    fit: DifferentialFit
    likelihood: DifferentialLikelihood
    diagnostics: DifferentialDiagnostics | None = None
    windowed: WindowedScore | None = None

    def group_index(self, group: str) -> int:
        try:
            return self.group_names.index(group)
        except ValueError:
            raise KeyError(f"Unknown group {group!r}") from None

    def mu(self, group: str) -> np.ndarray:
        return self.fit.group_mu[self.group_index(group)]

    def conditional_sig2(self, group: str) -> np.ndarray:
        return self.fit.group_conditional_sig2[self.group_index(group)]

    def conditional_sigma(self, group: str) -> np.ndarray:
        return np.sqrt(self.conditional_sig2(group))

    def loglik_mu(self, group: str) -> np.ndarray:
        """Sigma-squared-integrated log likelihood, position x mu."""
        return self.likelihood.group_loglik_mu[self.group_index(group)]

    def loglik_mu_sig2(self, group: str) -> np.ndarray:
        """Raw log likelihood before the variance prior, position x mu x sig2."""
        value = self.likelihood.group_loglik_mu_sig2
        if value is None:
            raise ValueError(
                "The raw mu/sig2 likelihood was not retained; rerun with "
                "keep_sig2_loglik=True"
            )
        return value[self.group_index(group)]

    def log_sig2_prior_mass(self, group: str) -> np.ndarray:
        return self.likelihood.group_log_sig2_prior_mass[
            self.group_index(group)
        ]

    def sig2_prior_mass(self, group: str) -> np.ndarray:
        return np.exp(self.log_sig2_prior_mass(group))

    def to_npz(
        self,
        path,
        *,
        compressed: bool = True,
        require_sig2_loglik: bool = True,
        extra: Mapping[str, Any] = None,
    ) -> Path:
        """Save results needed for plotting and downstream segmentation."""
        path = Path(path)

        if path.suffix != ".npz":
            path = Path(f"{path}.npz")

        path.parent.mkdir(parents=True, exist_ok=True)

        fit = self.fit
        likelihood = self.likelihood
        raw_loglik = likelihood.group_loglik_mu_sig2

        if require_sig2_loglik and raw_loglik is None:
            raise ValueError(
                "The position × mu × sig2 likelihood was not retained. "
                "Rerun DifferentialFitLoader with keep_sig2_loglik=True."
            )

        arrays: dict[str, np.ndarray] = {
            # File schema
            "format_version": np.array(1, dtype=np.int16),
            "group_names": np.asarray(
                self.group_names,
                dtype=np.str_,
            ),
            "n_position": np.array(
                fit.statistic.size,
                dtype=np.int64,
            ),

            # Parameter grids
            "mu_x": np.asarray(
                likelihood.mu_x,
                dtype=np.float64,
            ),
            "sig2_x": np.asarray(
                likelihood.sig2_x,
                dtype=np.float64,
            ),

            # Required input for segmentation:
            # group × position × mu
            "group_loglik_mu": np.asarray(
                likelihood.group_loglik_mu,
                dtype=np.float64,
            ),

            # Normalized finite-grid inverse-chi-squared prior masses:
            # group × sig2
            "group_log_sig2_prior_mass": np.asarray(
                likelihood.group_log_sig2_prior_mass,
                dtype=np.float64,
            ),

            # Pointwise differential results
            "lrt_statistic": np.asarray(
                fit.statistic,
                dtype=np.float64,
            ),
            "lrt_df": np.array(
                fit.df,
                dtype=np.int16,
            ),
            "lrt_log_pvalue": np.asarray(
                fit.log_pvalue,
                dtype=np.float64,
            ),
            "common_mu": np.asarray(
                fit.common_mu,
                dtype=np.float64,
            ),
            "group_mu": np.asarray(
                fit.group_mu,
                dtype=np.float64,
            ),
            "group_conditional_sig2": np.asarray(
                fit.group_conditional_sig2,
                dtype=np.float64,
            ),
        }

        # Raw likelihood before applying the sigma² prior:
        # group × position × mu × sig2
        if raw_loglik is not None:
            arrays["group_loglik_mu_sig2"] = np.asarray(
                raw_loglik,
                dtype=np.float64,
            )

        if self.windowed is not None:
            arrays.update({
                "window_radius": np.array(
                    self.windowed.radius,
                    dtype=np.int32,
                ),
                "window_score": np.asarray(
                    self.windowed.score,
                    dtype=np.float64,
                ),
                "window_nominal_log_pvalue": np.asarray(
                    self.windowed.nominal_log_pvalue,
                    dtype=np.float64,
                ),
            })

        if self.diagnostics is not None:
            arrays.update({
                "theta_x": np.asarray(
                    self.diagnostics.theta_x,
                    dtype=np.float64,
                ),
                "group_invchi2_params": np.asarray(
                    self.diagnostics.group_invchi2_params,
                    dtype=np.float64,
                ),
                "theta_coverage": np.array(
                    self.diagnostics.theta_coverage,
                    dtype=np.float64,
                ),
                "mu_boundary_fraction": np.array(
                    self.diagnostics.mu_boundary_fraction,
                    dtype=np.float64,
                ),
            })

        if extra is not None:
            for key, value in extra.items():
                if not isinstance(key, str) or not key:
                    raise ValueError(
                        "Extra metadata keys must be nonempty strings"
                    )

                if key in arrays:
                    raise KeyError(
                        f"Extra metadata key {key!r} conflicts with "
                        "a reserved result key"
                    )

                array = np.asarray(value)

                if array.dtype == object:
                    flat = array.ravel()

                    if all(
                        isinstance(item, (str, np.str_))
                        for item in flat
                    ):
                        array = array.astype(np.str_)
                    else:
                        raise TypeError(
                            f"Extra value {key!r} has object dtype. "
                            "Pass a numeric, boolean, or string array instead."
                        )

                arrays[key] = array

        save = np.savez_compressed if compressed else np.savez
        save(path, **arrays)

        return path


@dataclass(frozen=True, slots=True)
class _AlignedData:
    group_names: tuple[str, ...]
    group_indices: tuple[np.ndarray, ...]
    obs: np.ndarray
    exp: np.ndarray
    disp_models: np.ndarray


@dataclass(frozen=True, slots=True)
class _Grids:
    theta_x: np.ndarray
    mu_x: np.ndarray
    sig2_x: np.ndarray
    sig2_weights: np.ndarray
    theta_kernel: np.ndarray
    theta_coverage: float


class DifferentialFootprintModel:
    """Fit the integrated group-mean model and retain requested surfaces."""

    def __init__(self, config: DifferentialModelConfig | None = None):
        self.config = config or DifferentialModelConfig()

    def fit(
        self,
        groups_data: pd.Series,
        obs: pd.DataFrame,
        exp: pd.DataFrame,
        disp_models: pd.Series,
        selected_groups: Sequence[str] | None = None,
        nb: np.ndarray | None = None,
        keep_diagnostics: bool = False,
        keep_sig2_loglik: bool = False,
    ) -> DifferentialAnalysis:
        aligned = self._align(
            groups_data,
            obs,
            exp,
            disp_models,
            selected_groups,
        )
        grids = self._make_grids()
        nb = self._get_nb(aligned, grids, nb)
        fit, likelihood, prior, boundary_fraction = self._fit_lrt(
            aligned,
            grids,
            nb,
            keep_sig2_loglik,
        )
        diagnostics = (
            DifferentialDiagnostics(
                theta_x=grids.theta_x,
                group_invchi2_params=prior,
                theta_coverage=grids.theta_coverage,
                mu_boundary_fraction=boundary_fraction,
            )
            if keep_diagnostics
            else None
        )
        return DifferentialAnalysis(
            group_names=aligned.group_names,
            fit=fit,
            likelihood=likelihood,
            diagnostics=diagnostics,
        )

    @staticmethod
    def _align(
        groups: pd.Series,
        obs: pd.DataFrame,
        exp: pd.DataFrame,
        disp_models: pd.Series,
        selected: Sequence[str] | None,
    ) -> _AlignedData:
        if not isinstance(groups, pd.Series) or not groups.index.is_unique:
            raise ValueError("groups_data must be a Series with a unique index")
        if groups.isna().any():
            raise ValueError("Group labels contain missing values")
        if not obs.columns.equals(exp.columns):
            raise ValueError("obs and exp must have identical columns")

        if selected is None:
            names = tuple(pd.unique(groups))
        else:
            names = tuple(dict.fromkeys(selected))
            available = set(groups)
            missing = [name for name in names if name not in available]
            if missing:
                raise ValueError(f"Selected groups are absent: {missing}")
            groups = groups[groups.isin(names)]

        if len(names) < 2:
            raise ValueError("At least two groups are required")

        index = groups.index
        try:
            obs_array = np.ascontiguousarray(obs.loc[index].to_numpy(float))
            exp_array = np.ascontiguousarray(exp.loc[index].to_numpy(float))
            disp_array = disp_models.loc[index].to_numpy()
        except KeyError as error:
            raise ValueError("A grouped sample is missing from the input data") from error

        if obs_array.shape != exp_array.shape:
            raise ValueError("obs and exp must have the same shape")
        if np.any(~np.isfinite(obs_array)) or np.any(obs_array < 0):
            raise ValueError("obs must be finite and nonnegative")
        if np.any(~np.isfinite(exp_array)) or np.any(exp_array < 0):
            raise ValueError("exp must be finite and nonnegative")

        labels = groups.to_numpy()
        group_indices = tuple(
            np.flatnonzero(labels == name) for name in names
        )
        small = [
            name
            for name, index_ in zip(names, group_indices)
            if index_.size < 2
        ]
        if small:
            raise ValueError(f"Groups need at least two samples: {small}")

        return _AlignedData(
            group_names=names,
            group_indices=group_indices,
            obs=obs_array,
            exp=exp_array,
            disp_models=disp_array,
        )

    def _make_grids(self) -> _Grids:
        config = self.config
        mu_x = self._grid(
            config.mu_x,
            np.linspace(-6.0, 6.0, 481),
            "mu_x",
        )
        sig2_x = self._grid(
            config.sig2_x,
            np.geomspace(1e-3, 10.0, 181),
            "sig2_x",
        )

        if np.any(sig2_x <= 0):
            raise ValueError("sig2_x must be positive")
        if (
            config.theta_step <= 0
            or config.theta_tail_sd <= 0
            or config.base_chunk_size < 1
        ):
            raise ValueError(
                "theta_step, theta_tail_sd, and base_chunk_size must be positive"
            )
        if not 0 < config.min_theta_mass <= 1:
            raise ValueError("min_theta_mass must be in (0, 1]")
        if not 0 <= config.max_boundary_fraction <= 1:
            raise ValueError("max_boundary_fraction must be in [0, 1]")

        if config.theta_x is None:
            margin = config.theta_tail_sd * np.sqrt(sig2_x[-1])
            lo = mu_x[0] - margin
            hi = mu_x[-1] + margin
            n_theta = int(np.ceil((hi - lo) / config.theta_step)) + 1
            theta_x = np.linspace(lo, hi, n_theta)
        else:
            theta_x = self._grid(config.theta_x, None, "theta_x")

        theta_coverage = self._theta_coverage(mu_x, sig2_x, theta_x)
        if theta_coverage < config.min_theta_mass:
            warnings.warn(
                f"The theta grid covers as little as {theta_coverage:.2%} "
                "before residual tail mass is assigned to its edge points; "
                "widen theta_x.",
                RuntimeWarning,
            )

        return _Grids(
            theta_x=theta_x,
            mu_x=mu_x,
            sig2_x=sig2_x,
            sig2_weights=self._quadrature_weights(sig2_x),
            theta_kernel=self._theta_kernel(mu_x, sig2_x, theta_x),
            theta_coverage=theta_coverage,
        )

    @staticmethod
    def _get_nb(
        data: _AlignedData,
        grids: _Grids,
        nb: np.ndarray | None,
    ) -> np.ndarray:
        if nb is None:
            step = np.diff(grids.theta_x)
            if not np.allclose(step, step[0]):
                raise ValueError(
                    "compute_logpmf_values requires equally spaced theta_x"
                )
            nb = differential.compute_logpmf_values(
                data.disp_models,
                data.obs,
                data.exp,
                float(grids.theta_x[0]),
                float(grids.theta_x[-1]),
                grids.theta_x.size,
            )

        nb = np.asarray(nb, dtype=np.float64)
        expected = (
            grids.theta_x.size,
            data.obs.shape[1],
            data.obs.shape[0],
        )
        if nb.shape != expected:
            raise ValueError(f"Expected nb shape {expected}, got {nb.shape}")
        if np.any(np.isnan(nb)) or np.any(np.isposinf(nb)):
            raise ValueError("nb contains NaN or +inf")
        if np.any(np.all(np.isneginf(nb), axis=0)):
            raise ValueError(
                "A position/sample has -inf likelihood at every theta"
            )
        return nb

    def _fit_lrt(
        self,
        data: _AlignedData,
        grids: _Grids,
        nb: np.ndarray,
        keep_sig2_loglik: bool,
    ) -> tuple[
        DifferentialFit,
        DifferentialLikelihood,
        np.ndarray,
        float,
    ]:
        n_group = len(data.group_names)
        n_position = data.obs.shape[1]
        n_mu = grids.mu_x.size
        n_sig2 = grids.sig2_x.size
        rows = np.arange(n_position)

        group_loglik_mu = np.empty((n_group, n_position, n_mu))
        group_loglik_mu_sig2 = (
            np.empty((n_group, n_position, n_mu, n_sig2))
            if keep_sig2_loglik
            else None
        )
        group_log_sig2_prior_mass = np.empty((n_group, n_sig2))
        group_mu = np.empty((n_group, n_position))
        group_conditional_sig2 = np.empty((n_group, n_position))
        group_invchi2_params = np.empty((n_group, 2))

        common_loglik_mu = np.zeros((n_position, n_mu))
        alternative_loglik = np.zeros(n_position)
        boundary_count = 0
        log_ratio = np.log1p(data.obs) - np.log1p(data.exp)

        for group_index, sample_index in enumerate(data.group_indices):
            nu, prior_sig2 = self._fit_prior(log_ratio[sample_index])
            log_prior_mass = self._sig2_log_prior_mass(
                grids.sig2_x,
                grids.sig2_weights,
                nu,
                prior_sig2,
            )
            loglik_mu, conditional_sig2_index, raw_loglik = (
                self._integrate_group(
                    theta_kernel=grids.theta_kernel,
                    log_sig2_prior_mass=log_prior_mass,
                    nb=nb,
                    sample_index=sample_index,
                    keep_sig2_loglik=keep_sig2_loglik,
                )
            )

            mu_index = np.argmax(loglik_mu, axis=1)
            group_loglik_mu[group_index] = loglik_mu
            group_log_sig2_prior_mass[group_index] = log_prior_mass
            if group_loglik_mu_sig2 is not None:
                group_loglik_mu_sig2[group_index] = raw_loglik

            common_loglik_mu += loglik_mu
            alternative_loglik += loglik_mu[rows, mu_index]
            group_mu[group_index] = grids.mu_x[mu_index]
            group_conditional_sig2[group_index] = grids.sig2_x[
                conditional_sig2_index
            ]
            group_invchi2_params[group_index] = nu, prior_sig2
            boundary_count += np.count_nonzero(
                (mu_index == 0) | (mu_index == n_mu - 1)
            )

        common_mu_index = np.argmax(common_loglik_mu, axis=1)
        null_loglik = common_loglik_mu[rows, common_mu_index]
        gain = alternative_loglik - null_loglik
        minimum_gain = float(gain.min())
        if minimum_gain < -1e-8:
            raise RuntimeError(
                "Nested likelihood violation: minimum alternative-null "
                f"difference was {minimum_gain:.3e}"
            )

        boundary_fraction = boundary_count / (n_group * n_position)
        if boundary_fraction > self.config.max_boundary_fraction:
            warnings.warn(
                f"{boundary_fraction:.2%} of group/position means hit a "
                "mu-grid boundary.",
                RuntimeWarning,
            )

        fit = DifferentialFit(
            statistic=2.0 * np.maximum(gain, 0.0),
            df=n_group - 1,
            common_mu=grids.mu_x[common_mu_index],
            group_mu=group_mu,
            group_conditional_sig2=group_conditional_sig2,
        )
        likelihood = DifferentialLikelihood(
            mu_x=grids.mu_x,
            sig2_x=grids.sig2_x,
            group_loglik_mu=group_loglik_mu,
            group_log_sig2_prior_mass=group_log_sig2_prior_mass,
            group_loglik_mu_sig2=group_loglik_mu_sig2,
        )
        return fit, likelihood, group_invchi2_params, boundary_fraction

    def _integrate_group(
        self,
        theta_kernel: np.ndarray,
        log_sig2_prior_mass: np.ndarray,
        nb: np.ndarray,
        sample_index: np.ndarray,
        keep_sig2_loglik: bool,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray | None]:
        """Calculate one group's raw and variance-integrated likelihoods."""
        n_sig2, n_mu, n_theta = theta_kernel.shape
        kernel = np.ascontiguousarray(
            theta_kernel.reshape(n_sig2 * n_mu, n_theta)
        )
        n_position = nb.shape[1]
        n_sample = sample_index.size

        loglik_mu = np.empty((n_position, n_mu))
        conditional_sig2_index = np.empty(n_position, dtype=np.int32)
        loglik_mu_sig2 = (
            np.empty((n_position, n_mu, n_sig2))
            if keep_sig2_loglik
            else None
        )

        for lo in range(0, n_position, self.config.base_chunk_size):
            hi = min(lo + self.config.base_chunk_size, n_position)
            n_chunk = hi - lo

            x = np.take(nb[:, lo:hi], sample_index, axis=2)
            shift = np.max(x, axis=0)
            if not np.all(np.isfinite(shift)):
                raise ValueError(
                    "A position/sample has no finite theta likelihood"
                )

            scaled = np.exp(x - shift[None]).reshape(n_theta, -1)
            with np.errstate(divide="ignore"):
                raw = np.log(kernel @ scaled).reshape(
                    n_sig2,
                    n_mu,
                    n_chunk,
                    n_sample,
                )
            raw = raw.sum(axis=-1) + shift.sum(axis=-1)[None, None]

            integrated = logsumexp(
                raw + log_sig2_prior_mass[:, None, None],
                axis=0,
            ).T
            mu_index = np.argmax(integrated, axis=1)
            selected_raw = raw[:, mu_index, np.arange(n_chunk)]

            loglik_mu[lo:hi] = integrated
            conditional_sig2_index[lo:hi] = np.argmax(
                selected_raw,
                axis=0,
            )
            if loglik_mu_sig2 is not None:
                loglik_mu_sig2[lo:hi] = raw.transpose(2, 1, 0)

        if np.any(~np.isfinite(np.max(loglik_mu, axis=1))):
            raise RuntimeError(
                "No finite integrated likelihood for at least one position"
            )
        return loglik_mu, conditional_sig2_index, loglik_mu_sig2

    @staticmethod
    def _fit_prior(values: np.ndarray) -> tuple[float, float]:
        variance = np.nanvar(values, axis=0, ddof=1)
        variance = variance[np.isfinite(variance) & (variance > 0)]
        if variance.size == 0:
            raise ValueError("No finite positive variances available")

        def objective(log_parameters: np.ndarray) -> float:
            nu, sig2 = np.exp(log_parameters)
            return -float(invchi2.log_likelihood(variance, nu, sig2))

        result = scipy.optimize.minimize(
            objective,
            np.log([5.0, np.median(variance)]),
            method="L-BFGS-B",
            bounds=[(-8.0, 12.0), (-20.0, 12.0)],
        )
        if not result.success:
            raise RuntimeError(
                f"Inverse-chi-squared prior fit failed: {result.message}"
            )
        nu, sig2 = np.exp(result.x)
        return float(nu), float(sig2)

    @staticmethod
    def _sig2_log_prior_mass(
        sig2_x: np.ndarray,
        quadrature_weights: np.ndarray,
        nu: float,
        sig2: float,
    ) -> np.ndarray:
        log_mass = np.fromiter(
            (invchi2.logpdf(x, nu, sig2) for x in sig2_x),
            dtype=np.float64,
            count=sig2_x.size,
        ) + np.log(quadrature_weights)
        if not np.any(np.isfinite(log_mass)):
            raise ValueError("Variance prior has no finite mass on sig2_x")
        return log_mass - logsumexp(log_mass)

    @staticmethod
    def _theta_kernel(
        mu_x: np.ndarray,
        sig2_x: np.ndarray,
        theta_x: np.ndarray,
    ) -> np.ndarray:
        boundaries = np.r_[
            -np.inf,
            (theta_x[:-1] + theta_x[1:]) / 2.0,
            np.inf,
        ]
        z = (
            boundaries[None, None, :] - mu_x[None, :, None]
        ) / np.sqrt(sig2_x)[:, None, None]
        mass = np.maximum(np.diff(ndtr(z), axis=-1), 0.0)
        mass /= mass.sum(axis=-1, keepdims=True)
        return np.ascontiguousarray(mass)

    @staticmethod
    def _theta_coverage(
        mu_x: np.ndarray,
        sig2_x: np.ndarray,
        theta_x: np.ndarray,
    ) -> float:
        lo = theta_x[0] - (theta_x[1] - theta_x[0]) / 2.0
        hi = theta_x[-1] + (theta_x[-1] - theta_x[-2]) / 2.0
        scale = np.sqrt(sig2_x)[:, None]
        coverage = (
            ndtr((hi - mu_x[None]) / scale)
            - ndtr((lo - mu_x[None]) / scale)
        )
        return float(np.min(coverage))

    @staticmethod
    def _quadrature_weights(x: np.ndarray) -> np.ndarray:
        weights = np.empty_like(x)
        weights[0] = 0.5 * (x[1] - x[0])
        weights[-1] = 0.5 * (x[-1] - x[-2])
        weights[1:-1] = 0.5 * (x[2:] - x[:-2])
        return weights

    @staticmethod
    def _grid(
        value: np.ndarray | None,
        default: np.ndarray | None,
        name: str,
    ) -> np.ndarray:
        x = np.asarray(default if value is None else value, dtype=np.float64)
        if (
            x.ndim != 1
            or x.size < 2
            or not np.all(np.isfinite(x))
            or np.any(np.diff(x) <= 0)
        ):
            raise ValueError(
                f"{name} must be finite, one-dimensional, and increasing"
            )
        return x


class WindowedStoufferScorer:
    """Compute a local score and nominal normal p-value under independence."""

    def __init__(self, radius: int = 3):
        if radius < 0:
            raise ValueError("radius must be nonnegative")
        self.radius = radius

    def compute(self, fit: DifferentialFit) -> WindowedScore:
        statistic = np.asarray(fit.statistic, dtype=np.float64)
        valid = np.isfinite(statistic)
        if statistic.ndim != 1 or fit.df <= 0:
            raise ValueError("Invalid LRT statistic or degrees of freedom")
        if np.any(statistic[valid] < 0):
            raise ValueError("LRT statistics must be nonnegative")

        point_score = np.zeros_like(statistic)
        point_score[valid] = self._point_score(statistic[valid], fit.df)
        kernel = np.ones(2 * self.radius + 1)
        total = convolve1d(
            point_score,
            kernel,
            mode="constant",
            cval=0.0,
        )
        count = convolve1d(
            valid.astype(np.float64),
            kernel,
            mode="constant",
            cval=0.0,
        )
        score = np.divide(
            total,
            np.sqrt(count),
            out=np.full_like(total, np.nan),
            where=count > 0,
        )
        return WindowedScore(radius=self.radius, score=score)

    @staticmethod
    def _point_score(statistic: np.ndarray, df: int) -> np.ndarray:
        score = -ndtri_exp(st.chi2.logsf(statistic, df))
        score[statistic == 0] = 0.0
        return score


class DifferentialFitLoader(PlotDataLoader):
    """Stage 1: fit the model and attach ``data.differential``."""

    def apply(self, data: DataBundle, **kwargs) -> DataBundle:
        return self._load(data, **kwargs)

    def _load(
        self,
        data: DataBundle,
        selected_groups: Sequence[str] | None = None,
        nb: np.ndarray | None = None,
        config: DifferentialModelConfig | None = None,
        keep_diagnostics: bool = False,
        keep_sig2_loglik: bool = False,
    ) -> DataBundle:
        data.differential = DifferentialFootprintModel(config).fit(
            groups_data=data.groups_data,
            obs=data.obs,
            exp=data.exp,
            disp_models=data.disp_models,
            selected_groups=selected_groups,
            nb=nb,
            keep_diagnostics=keep_diagnostics,
            keep_sig2_loglik=keep_sig2_loglik,
        )
        return data


class WindowedScoreLoader(PlotDataLoader):
    """Stage 2: attach ``data.differential.windowed``."""

    def apply(self, data: DataBundle, **kwargs) -> DataBundle:
        return self._load(data, **kwargs)

    def _load(
        self,
        data: DataBundle,
        radius: int = 3,
    ) -> DataBundle:
        analysis = getattr(data, "differential", None)
        if not isinstance(analysis, DifferentialAnalysis):
            raise ValueError("Run DifferentialFitLoader first")
        analysis.windowed = WindowedStoufferScorer(radius).compute(analysis.fit)
        return data
