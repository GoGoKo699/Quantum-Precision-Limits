#!/usr/bin/env python3
"""Exact rational certificates for a four-qubit, rank-two sensing example.

Floating point optimization only proposes witnesses. All affine constraints,
positive quadratic dual forms, and final lower/upper bounds are checked with
exact rational arithmetic. This is a classical check, not quantum tomography.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import scipy.linalg
from scipy.optimize import minimize_scalar
import sympy as sp


I2 = sp.eye(2)
X = sp.Matrix([[0, 1], [1, 0]])
Y = sp.Matrix([[0, -sp.I], [sp.I, 0]])
Z = sp.diag(1, -1)


def pauli(n, entries):
    factors = [entries.get(i, I2) for i in range(n)]
    return sp.kronecker_product(*factors)


def rational_circuit(n=4):
    """Hardcoded gates with rational cosines/sines and exact even parity."""
    psi = sp.zeros(2**n, 1)
    psi[0] = 1
    gates = [
        ({0: X, 1: X}, 3, 4, 5),
        ({1: Z}, 5, 12, 13),
        ({1: X, 2: X}, 8, 15, 17),
        ({0: Z, 1: Z}, 7, 24, 25),
        ({2: X, 3: X}, 20, 21, 29),
        ({2: Z}, 12, 35, 37),
        ({1: Z, 2: Z}, 9, 40, 41),
        ({0: X, 3: X}, 11, 60, 61),
        ({0: Z}, 28, 45, 53),
    ]
    for locations, c, s, den in gates:
        assert c*c + s*s == den*den
        psi = (sp.Rational(c, den)*psi
               - sp.I*sp.Rational(s, den)*pauli(n, locations)*psi)
    psi = psi.applyfunc(sp.expand)
    assert sp.expand((psi.conjugate().T*psi)[0]) == 1
    parity = pauli(n, {i: Z for i in range(n)})
    assert parity*psi == psi
    partner = pauli(n, {0: X})*psi
    assert sp.expand((psi.conjugate().T*partner)[0]) == 0
    return psi, partner, gates


def exact_completion(psi, partner):
    """Complete rational complex orthonormal vectors without square roots."""
    d = len(psi)
    eye = sp.eye(d)

    def send_basis_to(vector, index):
        e = eye[:, index]
        if vector == e:
            return eye
        v = e-vector
        return (eye-v*v.conjugate().T/(1-sp.conjugate(vector[index]))).applyfunc(sp.cancel)

    q1 = send_basis_to(psi, 0)
    z = (q1.conjugate().T*partner).applyfunc(sp.cancel)
    assert z[0] == 0
    q2 = send_basis_to(z, 1)
    q = (q1*q2).applyfunc(sp.cancel)
    assert q[:, 0] == psi and q[:, 1] == partner
    return q


def exact_model():
    n = 4
    psi, partner, gates = rational_circuit(n)
    q = exact_completion(psi, partner)
    d = 2**n
    weights = [sp.Rational(3, 4), sp.Rational(1, 4)] + [sp.S.Zero]*(d-2)
    generators = [pauli(n, {1: X})/2, pauli(n, {1: Y})/2]
    ge = [(q.conjugate().T*g*q).applyfunc(sp.cancel) for g in generators]
    pairs = [(i, j) for i in range(2) for j in range(i+1, d)]
    # E_R=|i><j|+|j><i| and E_I=-i|i><j|+i|j><i|.
    a = sp.zeros(2, 2*len(pairs))
    diag_r, pair_i = [], []
    for k, (i, j) in enumerate(pairs):
        w, z = weights[i]+weights[j], weights[i]-weights[j]
        diag_r += [w, w]
        pair_i.append(z)
        for ell in range(2):
            dij = sp.expand(-sp.I*ge[ell][i, j]*(weights[j]-weights[i]))
            a[ell, 2*k] = 2*sp.re(dij)
            a[ell, 2*k+1] = -2*sp.im(dij)
    a = a.applyfunc(sp.cancel)
    # Diagonal estimator entries have zero derivative constraints and only add
    # nonnegative real covariance. Kernel-kernel entries are immaterial.
    return a, diag_r, pair_i, {
        'qubits': n, 'physical_dimension': d, 'state_rank': 2,
        'occupied_weights': ['3/4', '1/4'], 'parameters': 2,
        'gates': [{ 'paulis': {str(k): ('X' if v == X else 'Z') for k,v in loc.items()},
                    'cos': f'{c}/{den}', 'sin': f'{s}/{den}'} for loc,c,s,den in gates],
        'generators': ['X_1/2', 'Y_1/2'],
        'real_coordinates_per_estimator': len(diag_r),
    }


def hinv_apply(v, beta, r, pair_i):
    """Apply H_beta^-1 by exact 2x2 blocks; works for NumPy or SymPy."""
    size = len(r)
    out = [0]*len(v)
    for k, z in enumerate(pair_i):
        a, b = 2*k, 2*k+1
        w, t = r[a], beta*z
        den = w*w-t*t
        assert den > 0
        out[a] = (w*v[a]-t*v[size+b])/den
        out[size+b] = (w*v[size+b]-t*v[a])/den
        out[b] = (w*v[b]+t*v[size+a])/den
        out[size+a] = (w*v[size+a]+t*v[b])/den
    return out


def propose_witnesses(a, r, pair_i):
    af = np.array(a, dtype=float)
    rf = np.array(r, dtype=float)
    iz = np.array(pair_i, dtype=float)
    n = len(r)
    aa = scipy.linalg.block_diag(af, af)
    target = np.array([1., 0., 0., 1.])

    def at(beta):
        inv_at = np.column_stack([hinv_apply(v, beta, rf, iz) for v in aa])
        gram = aa@inv_at
        mult = np.linalg.solve(gram, target)
        coord = inv_at@mult
        return float(target@mult), coord, mult

    opt = minimize_scalar(lambda beta: -at(beta)[0], bounds=(-.999999, .999999),
                          method='bounded', options={'xatol': 1e-12})
    assert opt.success
    low, coord, mult = at(float(opt.x))
    assert np.max(np.abs(aa@coord-target)) < 1e-8
    _, _, pivot = scipy.linalg.qr(af, mode='economic', pivoting=True)
    return float(opt.x), coord, mult, list(map(int, pivot[:2]))


def rationalize(x, digits=11):
    scale = 10**digits
    return sp.Rational(round(float(x)*scale), scale)


def certify(a, r, pair_i, beta_float, coord, mult, pivot):
    n = len(r)
    beta = rationalize(beta_float)
    assert -1 < beta < 1
    y = sp.Matrix([rationalize(x) for x in mult])
    xs = [sp.Matrix([rationalize(x) for x in coord[k*n:(k+1)*n]]) for k in range(2)]
    ap = a[:, pivot]
    assert ap.det() != 0
    for k, x in enumerate(xs):
        rhs = sp.eye(2)[:, k]-a*x
        correction = ap.inv()*rhs
        for j, index in enumerate(pivot):
            x[index] += correction[j]
        assert a*x == sp.eye(2)[:, k]

    variance = sum(r[i]*(xs[0][i]**2+xs[1][i]**2) for i in range(n))
    imag = sum(z*(xs[0][2*k]*xs[1][2*k+1]-xs[0][2*k+1]*xs[1][2*k])
               for k,z in enumerate(pair_i))
    upper = sp.cancel(variance+2*abs(imag))
    aty = list(a.T*y[:2, :])+list(a.T*y[2:, :])
    inv_aty = hinv_apply(aty, beta, r, pair_i)
    lower = sp.cancel(2*(y[0]+y[3])-sum(x*z for x,z in zip(aty, inv_aty)))
    assert lower <= upper
    assert upper-lower < sp.Rational(1, 10**6)

    fisher = a*sp.diag(*[1/x for x in r])*a.T
    sld = sp.cancel(sp.trace(fisher.inv()))
    assert lower > sld
    fisher_margin = fisher-sp.Rational(4, 5)*sp.eye(2)
    assert fisher_margin[0, 0] > 0 and fisher_margin.det() > 0

    def enclosing_decimal(value, upper=False, digits=9):
        scale = 10**digits
        numerator, denominator = sp.fraction(value)
        scaled_num = int(numerator)*scale
        scaled_den = int(denominator)
        integer = (-((-scaled_num)//scaled_den)) if upper else scaled_num//scaled_den
        return f'{integer//scale}.{integer%scale:0{digits}d}'

    witness_digest = hashlib.sha256(str((beta, y, xs)).encode()).hexdigest()
    return {
        'holevo_lower': enclosing_decimal(lower),
        'holevo_upper': enclosing_decimal(upper, upper=True),
        'certified_width_upper': enclosing_decimal(upper-lower, upper=True),
        'sld_trace_lower': enclosing_decimal(sld),
        'sld_trace_upper': enclosing_decimal(sld, upper=True),
        'incompatibility_gap_lower': enclosing_decimal(lower-sld),
        'dual_beta': str(beta),
        'fisher_matrix_float': [[float(fisher[i,j]) for j in range(2)] for i in range(2)],
        'fisher_min_eigenvalue_float': float(np.linalg.eigvalsh(np.array(fisher,dtype=float))[0]),
        'fisher_lower_bound_exact': '4/5',
        'constraint_residual_exact': '0',
        'arithmetic': 'exact rational primal/dual checks; floating point only proposes witnesses',
        'witness_sha256': witness_digest,
    }


def check_compression(reference_certificate):
    """Numerical crosscheck of the exact tangent-support reduction."""
    psi, partner, _ = rational_circuit()
    psi = np.array(psi, dtype=complex).reshape(-1)
    partner = np.array(partner, dtype=complex).reshape(-1)
    factor = np.column_stack([np.sqrt(.75)*psi, .5*partner])
    gs = [np.array(pauli(4, {1: g})/2, dtype=complex) for g in (X, Y)]
    tangent = np.column_stack([factor]+[-1j*g@factor for g in gs])
    u, singular, _ = np.linalg.svd(tangent, full_matrices=False)
    rank = int(np.count_nonzero(singular > 1e-12))
    support = u[:, :rank]
    rho = support.conj().T@factor@factor.conj().T@support
    values, basis = np.linalg.eigh(rho)
    order = np.argsort(values)[::-1]
    values, basis = values[order], basis[:, order]
    assert np.max(np.abs(values[2:])) < 1e-12
    values[2:] = 0.
    pairs = [(i,j) for i in range(2) for j in range(i+1, rank)]
    a = np.zeros((2, 2*len(pairs)))
    r, iz = [], []
    ge = [basis.conj().T@support.conj().T@g@support@basis for g in gs]
    for k,(i,j) in enumerate(pairs):
        r += [values[i]+values[j]]*2
        iz.append(values[i]-values[j])
        for ell in range(2):
            dij = -1j*ge[ell][i,j]*(values[j]-values[i])
            a[ell,2*k], a[ell,2*k+1] = 2*dij.real, -2*dij.imag
    beta, coord, mult, _ = propose_witnesses(a, r, iz)
    size = len(r)
    lower = mult[0]+mult[3]
    variance = np.dot(r, coord[:size]**2+coord[size:]**2)
    imag = sum(z*(coord[2*k]*coord[size+2*k+1]-coord[2*k+1]*coord[size+2*k])
               for k,z in enumerate(iz))
    upper = variance+2*abs(imag)
    exact_lower = float(reference_certificate['holevo_lower'])
    exact_upper = float(reference_certificate['holevo_upper'])
    assert lower <= upper+1e-10
    assert abs(lower-exact_lower) < 2e-8
    assert abs(upper-exact_upper) < 2e-7
    residual = float(np.linalg.norm(tangent-support@support.conj().T@tangent))
    assert residual < 1e-12
    return {'physical_dimension': 16, 'compressed_dimension': rank,
            'compressed_lower_float': float(lower), 'compressed_upper_float': float(upper),
            'tangent_projection_residual_float': residual,
            'scope': 'Numerical dimension-reduction check; the full-model interval is certified separately by exact rational arithmetic.'}


def run():
    a, r, pair_i, model = exact_model()
    beta, coord, mult, pivot = propose_witnesses(a, r, pair_i)
    certificate = certify(a, r, pair_i, beta, coord, mult, pivot)
    return {'model': model, 'certificate': certificate, 'compression_check': check_compression(certificate),
            'scope': 'One fixed four-qubit example; no hardware, query-complexity, or generic-SDP benchmark.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = run()
    payload = json.dumps(result, indent=2, sort_keys=True)+'\n'
    if args.output:
        args.output.write_text(payload)
    print(payload)
