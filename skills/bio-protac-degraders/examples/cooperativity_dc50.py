#!/usr/bin/env python3
'''Cooperativity alpha and DC50/Dmax/hook-effect fitting for PROTAC dose-response data.

Implements the two prose-only workflows from SKILL.md's "Cooperativity (Alpha)" and
"DC50 / Dmax Characterization" sections as runnable code. Reference: numpy, scipy
(no version pin required beyond a `curve_fit` with `bounds=` support, scipy 1.0+).

Alpha and DC50/Dmax are always computed from binding/degradation MEASUREMENTS, not
from a structural prediction -- this script does not touch ternary-complex modeling.
The synthetic dataset below exists only to demonstrate the fitting code; do not treat
its numbers as real assay data.
'''

import numpy as np
from scipy.optimize import curve_fit


def cooperativity_alpha(kd_binary, kd_ternary):
    '''alpha = Kd_binary,target / Kd_ternary,target (SKILL.md, "Cooperativity (Alpha)").

    alpha > 1: positive cooperativity (ternary stronger than binary)
    alpha = 1: no cooperativity
    alpha < 1: negative cooperativity
    '''
    if kd_binary <= 0 or kd_ternary <= 0:
        raise ValueError('Kd values must be positive')
    alpha = kd_binary / kd_ternary
    if alpha > 1.0:
        label = 'positive cooperativity'
    elif alpha < 1.0:
        label = 'negative cooperativity'
    else:
        label = 'no cooperativity'
    return {'alpha': alpha, 'label': label}


def _hill(conc, dmax, dc50, hill):
    '''Standard ascending Hill/logistic degradation curve.'''
    return dmax / (1.0 + (dc50 / conc) ** hill)


def _hook_curve(conc, dmax, dc50, hill, hook_k, hook_hill):
    '''Bell-shaped dose-response: ascending Hill divided by a descending Hill (hook).

    Used only to generate the synthetic demonstration dataset below -- fitting uses
    `_hill` on the ascending arm alone, per SKILL.md's hook-effect guidance.
    '''
    return _hill(conc, dmax, dc50, hill) / (1.0 + (conc / hook_k) ** hook_hill)


def synthetic_dose_response(seed=42, n_points=14):
    '''SYNTHETIC data for demonstration only -- not a real cellular degradation assay.'''
    rng = np.random.default_rng(seed)
    true_params = dict(dmax=80.0, dc50=40.0, hill=1.3, hook_k=5000.0, hook_hill=1.6)
    conc = np.logspace(0, 4.7, n_points)  # ~1 nM to ~50 uM
    clean = _hook_curve(conc, **true_params)
    noise = rng.normal(0, 2.0, size=conc.shape)
    observed = np.clip(clean + noise, 0, 100)
    return conc, observed, true_params


def detect_hook(conc, degradation, downturn_threshold_pp=15.0):
    '''Flag a hook effect: peak degradation followed by a >threshold-point downturn
    at higher concentration (SKILL.md's "DC50 / Dmax Characterization" definition).
    '''
    peak_idx = int(np.argmax(degradation))
    peak = degradation[peak_idx]
    final = degradation[-1]
    hook = bool((peak - final) > downturn_threshold_pp) and peak_idx < len(degradation) - 1
    return {
        'peak_idx': peak_idx,
        'peak_conc': conc[peak_idx],
        'peak_pct': peak,
        'final_pct': final,
        'hook_effect': hook,
    }


PLATEAU_MIN_RATIO = 10.0  # peak conc / fitted DC50 needed for the ascending arm to reach its plateau


def fit_dc50(conc, degradation, peak_idx=None):
    '''Fit DC50/Dmax with a 3-parameter Hill curve.

    If `peak_idx` is given (from `detect_hook`), the fit is restricted to the
    ascending arm (through the point just past the peak) so a hook-affected
    high-dose tail does not distort the fit to a monotonic sigmoid -- this is the
    exclusion SKILL.md's hook-effect section requires.

    With a hook, the result also carries `plateau_reached` and `caveat`. The peak is
    where the hook already suppresses the curve, so when the peak concentration is
    less than PLATEAU_MIN_RATIO x the fitted DC50 (hook onset close to DC50) the arm
    never reaches its plateau and Dmax is underestimated (audit: 15 pp on a noise-free
    curve). The caveat then says to treat Dmax as a lower bound.
    '''
    peak_conc = None
    if peak_idx is not None:
        peak_conc = conc[peak_idx]
        sl = slice(0, peak_idx + 2)
        conc, degradation = conc[sl], degradation[sl]
    p0 = [max(degradation), float(np.median(conc)), 1.0]
    bounds = ([1.0, 1e-3, 0.1], [200.0, 1e6, 10.0])
    popt, _ = curve_fit(_hill, conc, degradation, p0=p0, bounds=bounds, maxfev=20000)
    out = {'dmax_fit': popt[0], 'dc50_fit': popt[1], 'hill_fit': popt[2],
           'plateau_reached': None, 'caveat': None}
    if peak_conc is not None:
        ratio = peak_conc / popt[1]
        out['plateau_reached'] = bool(ratio >= PLATEAU_MIN_RATIO)
        if not out['plateau_reached']:
            out['caveat'] = (f'peak at {peak_conc:.0f} nM is only {ratio:.1f}x the fitted DC50 '
                             f'(< {PLATEAU_MIN_RATIO:.0f}x): ascending-arm data may not reach the true '
                             'plateau; treat Dmax as a lower bound')
    return out


if __name__ == '__main__':
    print('--- Cooperativity alpha (synthetic Kd pairs, nM) ---')
    for kd_binary, kd_ternary in [(50.0, 3.6), (20.0, 35.0), (10.0, 10.0)]:
        result = cooperativity_alpha(kd_binary, kd_ternary)
        print(f'  Kd_binary={kd_binary}  Kd_ternary={kd_ternary}  '
              f'alpha={result["alpha"]:.2f} ({result["label"]})')

    print('\n--- DC50 / Dmax / hook effect (SYNTHETIC dose-response) ---')
    conc, degradation, truth = synthetic_dose_response()
    hook_info = detect_hook(conc, degradation)
    print(f'  Peak degradation: {hook_info["peak_pct"]:.1f}% at {hook_info["peak_conc"]:.0f} nM')
    print(f'  Final-point degradation: {hook_info["final_pct"]:.1f}% at {conc[-1]:.0f} nM')
    print(f'  Hook effect flag (>15pp downturn): {hook_info["hook_effect"]}')

    fit = fit_dc50(conc, degradation, peak_idx=hook_info['peak_idx'])
    print(f'  Fitted (ascending-arm only) DC50 = {fit["dc50_fit"]:.1f} nM, '
          f'Dmax = {fit["dmax_fit"]:.1f}%, Hill = {fit["hill_fit"]:.2f}')
    print(f'  Plateau reached: {fit["plateau_reached"]}; caveat: {fit["caveat"]}')

    # This dataset's true generating parameters are known (they built the synthetic
    # curve above); recovering them from the noisy points is the check that the fit
    # code works, not evidence about any real PROTAC.
    dc50_err = abs(fit['dc50_fit'] - truth['dc50']) / truth['dc50']
    dmax_err = abs(fit['dmax_fit'] - truth['dmax'])
    print(f'  [verification] planted DC50={truth["dc50"]} nM, relative error={dc50_err:.1%}')
    print(f'  [verification] planted Dmax={truth["dmax"]}%, absolute error={dmax_err:.1f}pp')
    assert hook_info['hook_effect'] is True, 'expected the synthetic curve to show a hook effect'
    assert dc50_err < 0.30, 'fitted DC50 too far from the planted synthetic truth'
    assert dmax_err < 12.0, 'fitted Dmax too far from the planted synthetic truth'
    assert fit['plateau_reached'] is True and fit['caveat'] is None, 'well-separated hook should not be caveated'
    print('  Recovery-against-planted-truth check: PASS')

    print('\n--- Narrow-separation hook (SYNTHETIC, hook_k/DC50 = 20, noise-free) ---')
    narrow = dict(dmax=65.0, dc50=120.0, hill=1.8, hook_k=2500.0, hook_hill=1.1)
    conc_n = np.logspace(0, 4.3, 18)
    deg_n = _hook_curve(conc_n, **narrow)
    hook_n = detect_hook(conc_n, deg_n)
    fit_n = fit_dc50(conc_n, deg_n, peak_idx=hook_n['peak_idx'])
    print(f'  Fitted Dmax = {fit_n["dmax_fit"]:.1f}% (planted {narrow["dmax"]}%), '
          f'DC50 = {fit_n["dc50_fit"]:.1f} nM (planted {narrow["dc50"]} nM)')
    print(f'  Caveat: {fit_n["caveat"]}')
    assert hook_n['hook_effect'] is True and fit_n['plateau_reached'] is False and fit_n['caveat'],         'narrow-separation hook must carry the lower-bound caveat'
    assert fit_n['dmax_fit'] < narrow['dmax'], 'caveated Dmax should sit below the planted asymptote'
    print('  Narrow-hook caveat check: PASS')
