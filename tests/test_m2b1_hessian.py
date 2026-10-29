import numpy as np

from iontrap_rng.hessian import fit_local_hessian


def test_unconstrained_normalized_fit_recovers_harmonic_hessian():
    rng = np.random.default_rng(271828)
    H = np.array([[1.0, .2, -.1], [.2, -.4, .15], [-.1, .15, -.6]])
    points = rng.uniform(-.4, .4, (500, 3))
    values = .5 * np.einsum('ni,ij,nj->n', points, H, points)
    fit = fit_local_hessian(points, values, radius=.5, order=2)
    assert fit.rank == 10
    assert np.isfinite(fit.condition_number) and fit.condition_number < 20
    assert np.linalg.norm(fit.hessian - H) < 1e-12
    assert abs(np.trace(fit.hessian)) < 1e-12


def test_harmonic_fit_trace_constraint_is_not_independent_validation():
    rng = np.random.default_rng(314159)
    points = rng.uniform(-.4, .4, (500, 3))
    values = points[:, 0] ** 2 + points[:, 1] ** 2 + points[:, 2] ** 2
    unconstrained = fit_local_hessian(points, values, radius=.5, order=2)
    constrained = fit_local_hessian(points, values, radius=.5, order=2,
                                    harmonic_quadratic=True)
    assert np.trace(unconstrained.hessian) > 5.9
    assert np.trace(constrained.hessian) == 0.0
    assert constrained.relative_residual > unconstrained.relative_residual


def test_rank_deficiency_is_reported_for_insufficient_samples():
    points = np.array([[0., 0., 0.], [.1, 0., 0.], [0., .1, 0.]])
    fit = fit_local_hessian(points, np.array([0., .1, .1]), radius=.2, order=2)
    assert fit.rank < 10
    assert np.isfinite(fit.condition_number)
