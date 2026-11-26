"""Call coherentSurv's R engine through rpy2; no statistical reimplementation.

Install coherentSurv in R first, then install rpy2, pandas and numpy in Python.
Run this file on a supported rpy2 platform. The optional --input CSV must follow
the demonstrated exact-history schema. Adapt formula, marks, units and horizons
to your study. This example is not held-out predictive validation.
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path)
    parser.add_argument('--output', type=Path, default=Path('python_forecasts.csv'))
    args = parser.parse_args()
    engine = importr('coherentSurv')
    stats = importr('stats')
    base = importr('base')
    utils = importr('utils')
    if args.input:
        frame = pd.read_csv(args.input)
        # The independently loaded R frame is the conversion-parity reference.
        reference_data = utils.read_csv(str(args.input.resolve()))
    else:
        ro.r('data("sim_perinatal", package = "coherentSurv")')
        reference_data = ro.globalenv['sim_perinatal']
        with (ro.default_converter + pandas2ri.converter).context():
            frame = ro.conversion.get_conversion().rpy2py(reference_data)
    required = {'id', 'start', 'stop', 'event_type', 'terminal', 'x'}
    if not required.issubset(frame.columns):
        raise ValueError('Required interval columns: ' + ', '.join(sorted(required)))
    with (ro.default_converter + pandas2ri.converter).context():
        r_data = ro.conversion.get_conversion().py2rpy(frame)
        r_history = ro.conversion.get_conversion().py2rpy(frame.loc[frame.id == frame.id.iloc[0]])
    fit = engine.cs_fit(ro.Formula('~ x'), data=r_data,
                        event_types=ro.StrVector(['1', '2']), time_unit='months')
    horizons = ro.FloatVector([0.25, 0.5, 1.0])
    prediction = stats.predict(fit, newdata=r_history, landmark=0,
                               horizon=horizons, type='all')
    # Keep model objects in R; convert only the requested result table.
    with (ro.default_converter + pandas2ri.converter).context():
        result = ro.conversion.get_conversion().rpy2py(base.as_data_frame(prediction))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, index=False)
    audit_text = str(engine.cs_audit_report(fit, predictions=prediction, format='json')[0])
    Path(str(args.output) + '.audit.json').write_text(audit_text, encoding='utf-8')
    # The reference fit is evaluated entirely in R, independently of pandas.
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
    wide = result.pivot(index=['unit','horizon'], columns='target', values='estimate')
    assert np.allclose(wide['risk'] + wide['survival'], 1, atol=1e-12, rtol=0)
    report = {'package': 'coherentSurv', 'version': str(ro.r('as.character(packageVersion("coherentSurv"))')[0]),
              'rows': len(result), 'maximum_absolute_difference_from_direct_r': errors,
              'r_convergence': json.loads(audit_text)['convergence'],
              'scope': 'Transport and numerical parity on the supplied dataset, not predictive validation.'}
    Path(str(args.output) + '.parity.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(result[['unit','horizon','target','estimate']].to_string(index=False))
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
