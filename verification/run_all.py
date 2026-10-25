#!/usr/bin/env python3
"""Run the three documented scientific checks and compare committed references."""
from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys


HERE = Path(__file__).resolve().parent


def execute(script, report):
    path = HERE/'actual'/report
    path.parent.mkdir(exist_ok=True)
    run = subprocess.run([sys.executable, str(HERE/script), '--output', str(path)],
                         capture_output=True, text=True, timeout=180)
    if run.returncode:
        raise RuntimeError(f'{script} failed:\n{run.stdout}\n{run.stderr}')
    return json.loads(path.read_text()), json.loads((HERE/'expected'/report).read_text())


def close(a, b, tolerance):
    assert abs(float(a)-float(b)) <= tolerance, (a, b, tolerance)


def main():
    actual, expected = execute('bures_verify.py', 'bures_verification.json')
    assert actual['status'] == expected['status'] == 'passed'
    for key in ('peeling_case_count', 'amplification_case_count', 'maximum_peeling_rounds'):
        assert actual['summary'][key] == expected['summary'][key]
    assert actual['summary']['minimum_captured_fraction'] >= .5-1e-10
    assert actual['summary']['maximum_normalization_identity_residual'] < 1e-10
    assert actual['summary']['maximum_rank_tail'] < 1e-10
    print('PASS: Bures peeling, normalization, factor rounding, and amplification rank checks')

    actual, expected = execute('verify_certified_optimization.py', 'certified_optimization_checks.json')
    assert actual == expected
    assert actual['all_checks_passed']
    print('PASS: exact rational physical-fit and primal/dual certificates (4 and 20/9)')

    actual, expected = execute('verify_holevo_circuit.py', 'holevo_circuit.json')
    assert actual['model'] == expected['model']
    certificate, reference = actual['certificate'], expected['certificate']
    assert certificate['constraint_residual_exact'] == '0'
    assert certificate['fisher_lower_bound_exact'] == '4/5'
    assert float(certificate['certified_width_upper']) < 1e-6
    for key in ('holevo_lower', 'holevo_upper', 'sld_trace_lower', 'sld_trace_upper'):
        close(certificate[key], reference[key], 2e-7)
    assert actual['compression_check']['compressed_dimension'] == 6
    assert actual['compression_check']['tangent_projection_residual_float'] < 1e-12
    print('PASS: four-qubit exact interval, incompatibility penalty, and dimension reduction')
    print('All three scientific checks passed. Reports are in verification/actual/.')


if __name__ == '__main__':
    main()
