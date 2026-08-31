"""Helpers for the polarization-CHSH analysis.

Only relative phase separations enter the Bell correlations, so an unknown
common shift of the sampled phase origin does not affect the result.  All
errors returned by this module are one-standard-deviation (1 sigma) errors.
"""

from __future__ import annotations

from pathlib import Path
from typing import Sequence

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit, minimize_scalar


N_CURVES = 8
N_PARAMETERS = 4


def plot_coinc(phases: np.ndarray, counts: np.ndarray, label: str) -> None:
    """Plot mean coincidence counts with 1-sigma error bars."""
    plt.errorbar(phases, counts[:, 0], yerr=counts[:, 1], fmt="none",
                 color="k", capsize=2, zorder=1)
    plt.scatter(phases, counts[:, 0], s=12, zorder=2, label=label)


def load_coinc(filenames: Sequence[str | Path], number_of_points: int):
    """Load repeats and return coincidence/Alice/Bob means and standard errors.

    Every row is ``Alice singles, Bob singles, coincidences``.  The error on a
    mean is the sample standard deviation (ddof=1) divided by sqrt(repeats).
    This is preferable to ddof=0 for the three repeats in this experiment.
    """
    if len(filenames) < 2:
        raise ValueError("At least two repeated files are required.")

    repeated = np.empty((len(filenames), number_of_points, 3), dtype=float)
    for repeat_index, filename in enumerate(filenames):
        values = np.loadtxt(filename, delimiter="\t", ndmin=2)
        if values.shape != (number_of_points, 3):
            raise ValueError(f"{filename} has shape {values.shape}; expected "
                             f"({number_of_points}, 3).")
        repeated[repeat_index] = values

    means = np.mean(repeated, axis=0)
    errors = np.std(repeated, axis=0, ddof=1) / np.sqrt(len(filenames))
    # Background subtraction can yield a negative fluctuation.  Clip only its
    # mean; keep the measured error so it retains a finite statistical weight.
    means[:, 2] = np.clip(means[:, 2], 0.0, None)

    alice = np.column_stack((means[:, 0], errors[:, 0]))
    bob = np.column_stack((means[:, 1], errors[:, 1]))
    coincidences = np.column_stack((means[:, 2], errors[:, 2]))
    return coincidences, alice, bob


def cos_f(x, offset, amplitude, angular_frequency, phase):
    """Cosine fringe: offset + amplitude*cos(frequency*x + phase)."""
    return offset + amplitude * np.cos(angular_frequency * x + phase)


def cos_f2(x, offset, amplitude, phase):
    """Unit-frequency cosine fringe retained for backwards compatibility."""
    return offset + amplitude * np.cos(x + phase)


def cos_jacobian(x: float, parameters: np.ndarray) -> np.ndarray:
    """Derivative of the fringe value with respect to all four parameters."""
    _, amplitude, angular_frequency, phase = parameters
    argument = angular_frequency * x + phase
    return np.array([1.0, np.cos(argument),
                     -amplitude * x * np.sin(argument),
                     -amplitude * np.sin(argument)])


def interf_curve(phases, curve, initial_guess):
    """Fit one fringe, treating the measured errors as absolute 1-sigma values."""
    sigma = np.asarray(curve[:, 1], dtype=float)
    if np.any(~np.isfinite(sigma)) or np.any(sigma <= 0):
        raise ValueError("Every coincidence point needs a finite positive error.")
    return curve_fit(cos_f, phases, curve[:, 0], sigma=sigma,
                     absolute_sigma=True, p0=initial_guess, maxfev=50_000)


def visibility_with_error(parameters, covariance):
    """Return |amplitude/offset| with covariance-propagated uncertainty."""
    offset, amplitude = parameters[:2]
    if offset == 0 or amplitude == 0:
        raise ZeroDivisionError("Visibility needs non-zero offset and amplitude.")
    visibility = abs(amplitude / offset)
    gradient = np.array([-visibility / offset, visibility / amplitude, 0.0, 0.0])
    variance = float(gradient @ covariance @ gradient)
    return visibility, np.sqrt(max(variance, 0.0))


def prediction_with_error(phase_a, parameters, covariance):
    """Evaluate a fitted fringe with its propagated fit uncertainty."""
    value = float(cos_f(phase_a, *parameters))
    gradient = cos_jacobian(phase_a, parameters)
    variance = float(gradient @ covariance @ gradient)
    return value, np.sqrt(max(variance, 0.0))


def block_covariance(covariances):
    """Assemble the covariance of eight independent four-parameter fits."""
    if covariances.shape != (N_CURVES, N_PARAMETERS, N_PARAMETERS):
        raise ValueError("covariances must have shape (8, 4, 4).")
    result = np.zeros((N_CURVES * N_PARAMETERS,) * 2)
    for index, covariance in enumerate(covariances):
        block = slice(N_PARAMETERS * index, N_PARAMETERS * (index + 1))
        result[block, block] = covariance
    return result


def correlation_with_error(phase_a, bob_setting, parameters, covariances):
    """Calculate E and propagate the covariance of its four fitted fringes.

    E = (N11 + N22 - N12 - N21) / sum(Nij).  A 32-element gradient is also
    returned so correlations sharing fit parameters can be combined correctly.
    """
    if bob_setting not in (0, 1):
        raise ValueError("bob_setting must be 0 or 1.")
    first = 4 * bob_setting
    indices = np.arange(first, first + 4)
    signs = np.array([1.0, -1.0, -1.0, 1.0])
    predicted = np.array([cos_f(phase_a, *parameters[i]) for i in indices])
    denominator = float(np.sum(predicted))
    if np.isclose(denominator, 0.0):
        raise ZeroDivisionError("The four predicted coincidence rates sum to zero.")

    correlation = float(signs @ predicted / denominator)
    d_e_d_n = (signs - correlation) / denominator
    gradient = np.zeros(N_CURVES * N_PARAMETERS)
    for local_index, curve_index in enumerate(indices):
        block = slice(N_PARAMETERS * curve_index,
                      N_PARAMETERS * (curve_index + 1))
        gradient[block] = d_e_d_n[local_index] * cos_jacobian(
            phase_a, parameters[curve_index])

    covariance = block_covariance(covariances)
    variance = float(gradient @ covariance @ gradient)
    return correlation, np.sqrt(max(variance, 0.0)), gradient


def chsh_with_error(reference_phase, parameters, covariances):
    """Return S, its full propagated error, and the four E components.

    S = E_xx + E_xy + E_yx - E_yy.  The absolute phase zero is absorbed in
    ``reference_phase``; only the pi/2 separation of the analyser settings and
    the pi separation between the two output ports matter.
    """
    total, total_gradient, components = _chsh_value_gradient(
        reference_phase, parameters, covariances)
    covariance = block_covariance(covariances)
    variance = float(total_gradient @ covariance @ total_gradient)
    return total, np.sqrt(max(variance, 0.0)), components


def _chsh_value_gradient(reference_phase, parameters, covariances):
    """Internal CHSH value/gradient calculation used for error propagation."""
    definitions = {
        "E_xx": (reference_phase + np.pi / 2, 0, +1.0),
        "E_xy": (reference_phase + np.pi / 2, 1, +1.0),
        "E_yx": (reference_phase, 0, +1.0),
        "E_yy": (reference_phase, 1, -1.0),
    }
    total = 0.0
    total_gradient = np.zeros(N_CURVES * N_PARAMETERS)
    components = {}
    for name, (phase_a, bob_setting, sign) in definitions.items():
        value, error, gradient = correlation_with_error(
            phase_a, bob_setting, parameters, covariances)
        components[name] = (value, error)
        total += sign * value
        total_gradient += sign * gradient

    return total, total_gradient, components


def chsh_extremum_with_error(parameters, covariances, kind="maximum",
                             phase_bounds=(-np.pi, np.pi)):
    """Locate a CHSH extremum and propagate errors to both S and its phase.

    The uncertainty of the optimal phase follows from implicit differentiation
    of dS/dphase = 0.  Mixed derivatives and curvature are evaluated by stable
    central differences.  The S error uses the envelope theorem, so it already
    includes the fact that the optimum shifts when fit parameters vary.
    """
    if kind not in {"maximum", "minimum"}:
        raise ValueError("kind must be 'maximum' or 'minimum'.")
    sign = -1.0 if kind == "maximum" else 1.0

    def objective(phase):
        value = _chsh_value_gradient(phase, parameters, covariances)[0]
        return sign * value

    optimum = minimize_scalar(objective, bounds=phase_bounds, method="bounded",
                              options={"xatol": 1e-12})
    phase = float(optimum.x)
    value, gradient, _ = _chsh_value_gradient(phase, parameters, covariances)
    covariance = block_covariance(covariances)
    value_error = np.sqrt(max(float(gradient @ covariance @ gradient), 0.0))

    step = 1e-4
    value_minus, gradient_minus, _ = _chsh_value_gradient(
        phase - step, parameters, covariances)
    value_plus, gradient_plus, _ = _chsh_value_gradient(
        phase + step, parameters, covariances)
    curvature = (value_plus - 2.0 * value + value_minus) / step**2
    mixed_derivative = (gradient_plus - gradient_minus) / (2.0 * step)
    phase_gradient = -mixed_derivative / curvature
    phase_error = np.sqrt(max(float(phase_gradient @ covariance @ phase_gradient), 0.0))
    return phase, phase_error, value, value_error


def compute_corr(phase_a, bob_setting, parameters):
    """Backwards-compatible correlation value without an error."""
    zero_covariance = np.zeros((N_CURVES, N_PARAMETERS, N_PARAMETERS))
    return correlation_with_error(
        phase_a, bob_setting, parameters, zero_covariance)[0]
