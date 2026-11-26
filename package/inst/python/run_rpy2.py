"""Call coherentSurv's R engine through rpy2; no statistical reimplementation.

Install coherentSurv in R first, then install rpy2, pandas and numpy in Python.
This example is checked with rpy2 3.6 on Linux; use run_forecast.py for the
process-based route on Windows. See https://rpy2.github.io/doc/latest/html/overview.html
for rpy2 platform prerequisites. No thread-safety or speed advantage is promised;
use separate processes for independent concurrent R sessions. The optional --input CSV must follow
the demonstrated exact-history schema. Adapt formula, marks, units and horizons
to your study. The default fits once and retains the R fit for prediction.
--verify-parity performs a second, direct-R fit as an explicit transport check.
--timings reports measured stages, excluding imports and R initialization.
This example is not held-out predictive validation.
"""
import argparse
import json
from pathlib import Path
from time import perf_counter


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, help='Optional exact-history CSV; default: bundled synthetic data.')
    parser.add_argument('--output', type=Path, default=Path('python_forecasts.csv'))
    parser.add_argument('--verify-parity', action='store_true', help='Refit directly in R and compare coefficients, covariance and forecasts.')
    parser.add_argument('--timings', action='store_true', help='Print conversion, fit, prediction and output timings in seconds.')
    args = parser.parse_args()
    if args.input is not None and not args.input.is_file():
        parser.error('The input CSV does not exist: ' + str(args.input))
    # Parse --help before importing optional dependencies or initializing R.
    try:
        import numpy as np
        import pandas as pd
        import rpy2.robjects as ro
        from rpy2.robjects import pandas2ri
        from rpy2.robjects.packages import importr, PackageNotInstalledError
    except (ImportError, OSError, RuntimeError, ValueError) as exc:
        parser.error('Cannot initialize the Python/R bridge. Install compatible R, '
                     'rpy2, pandas and numpy; see the platform prerequisites in --help. '
                     'Original error: ' + str(exc))
    try:
        engine = importr('coherentSurv')
    except PackageNotInstalledError as exc:
        parser.error('Install coherentSurv in the R library used by rpy2. ' + str(exc))
    stats = importr('stats')
    base = importr('base')
    timings = {}
    started = perf_counter()
    if args.input:
        # Keep label text, including leading zeros and literal "NA". An empty
        # event label means no event; empty numeric fields remain missing.
        try:
            frame = pd.read_csv(args.input, dtype={'id': 'string', 'event_type': 'string'},
                                keep_default_na=False,
                                na_values={name: [''] for name in ('start', 'stop', 'terminal', 'x')})
        except (OSError, ValueError) as exc:
            parser.error('Cannot read the input CSV: ' + str(exc))
    else:
        ro.r('data("sim_perinatal", package = "coherentSurv")')
        reference_data = ro.globalenv['sim_perinatal']
        with (ro.default_converter + pandas2ri.converter).context():
            frame = ro.conversion.get_conversion().rpy2py(reference_data)
    timings['input_read_seconds'] = perf_counter() - started
    required = {'id', 'start', 'stop', 'event_type', 'terminal', 'x'}
    missing = required.difference(frame.columns)
    if missing:
        parser.error('Missing required interval columns: ' + ', '.join(sorted(missing)))
    if frame.empty:
        parser.error('The input CSV contains no interval rows.')
    started = perf_counter()
    with (ro.default_converter + pandas2ri.converter).context():
        r_data = ro.conversion.get_conversion().py2rpy(frame)
        r_history = ro.conversion.get_conversion().py2rpy(frame.loc[frame.id == frame.id.iloc[0]])
    timings['input_conversion_seconds'] = perf_counter() - started
    started = perf_counter()
    fit = engine.cs_fit(ro.Formula('~ x'), data=r_data,
                        event_types=ro.StrVector(['1', '2']), time_unit='months')
    timings['fit_seconds'] = perf_counter() - started
    horizons = ro.FloatVector([0.25, 0.5, 1.0])
    started = perf_counter()
    prediction = stats.predict(fit, newdata=r_history, landmark=0,
                               horizon=horizons, type='all')
    timings['prediction_seconds'] = perf_counter() - started
    # Keep model objects in R; convert only the requested result table.
    started = perf_counter()
    with (ro.default_converter + pandas2ri.converter).context():
        result = base.as_data_frame(prediction)
    timings['output_conversion_seconds'] = perf_counter() - started
    if result.empty or not np.isfinite(result['estimate']).all():
        raise RuntimeError('R returned empty or nonfinite forecasts; inspect the fit diagnostics.')
    wide = result.pivot(index=['unit', 'horizon'], columns='target', values='estimate')
    if not np.allclose(wide['risk'] + wide['survival'], 1, atol=1e-12, rtol=0):
        raise AssertionError('Forecast risk and survival do not sum to one.')
    started = perf_counter()
    audit_text = str(engine.cs_audit_report(fit, predictions=prediction, format='json')[0])
    timings['audit_seconds'] = perf_counter() - started
    errors = None
    if args.verify_parity:
        started = perf_counter()
        if args.input:
            # Load independently of pandas, with the same declared label schema.
            reference_data = ro.r('''function(path) utils::read.csv(path,
                colClasses=c(id="character", event_type="character"),
                na.strings=character(), stringsAsFactors=FALSE, check.names=FALSE)''')(str(args.input.resolve()))
        reference = ro.r('''function(d) {
        f <- coherentSurv::cs_fit(~ x, d, event_types=c("1","2"), time_unit="months")
        p <- predict(f, d[d$id == d$id[1], ], landmark=0,
                     horizon=c(.25,.5,1), type="all")
        list(coefficients=coef(f), covariance=vcov(f), estimates=p$estimate)
        }''')(reference_data)
        errors = {
            'coefficients': float(np.max(np.abs(np.asarray(stats.coef(fit)) - np.asarray(reference.rx2('coefficients'))))),
            'covariance': float(np.max(np.abs(np.asarray(stats.vcov(fit)) - np.asarray(reference.rx2('covariance'))))),
            'forecasts': float(np.max(np.abs(result['estimate'].to_numpy() - np.asarray(reference.rx2('estimates'))))),
        }
        if not all(value < 1e-9 for value in errors.values()):
            raise AssertionError(errors)
        timings['parity_verification_seconds'] = perf_counter() - started
    report = {'package': 'coherentSurv', 'version': str(ro.r('as.character(packageVersion("coherentSurv"))')[0]),
              'verification_status': 'passed' if args.verify_parity else 'not_requested',
              'rows': len(result), 'maximum_absolute_difference_from_direct_r': errors,
              'r_convergence': json.loads(audit_text)['convergence'],
              'scope': ('Transport and numerical parity on the supplied dataset, not predictive validation.'
                        if args.verify_parity else 'Forecast export; direct-R parity was not requested. Not predictive validation.')}
    # Validate before publishing results. Always replace the parity sidecar so a
    # previous verified run is not mistaken for this run when verification is off.
    started = perf_counter()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, index=False)
    Path(str(args.output) + '.audit.json').write_text(audit_text, encoding='utf-8')
    Path(str(args.output) + '.parity.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    timings['output_write_seconds'] = perf_counter() - started
    print(result[['unit','horizon','target','estimate']].to_string(index=False))
    print(json.dumps(report, indent=2))
    if args.timings:
        print(json.dumps({'timings_seconds': timings,
                          'scope': 'One run; imports and R initialization excluded. Not a comparative benchmark.'}, indent=2))


if __name__ == '__main__':
    main()
