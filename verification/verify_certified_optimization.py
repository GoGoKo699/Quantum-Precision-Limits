#!/usr/bin/env python3
"""Exact rational checks for the finite Holevo certification pipeline.

This is a certificate verifier and two analytic witness generators, not a
generic SDP optimizer. All accepted inequalities and objective values use
SymPy rational arithmetic. No floating-point tolerance decides acceptance.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import sympy as sp


I = sp.I
Q = sp.Rational


def clean(matrix):
    return sp.Matrix(matrix).applyfunc(sp.simplify)


def is_zero(matrix):
    return all(sp.simplify(x) == 0 for x in matrix)


def must_reject(check):
    try:
        check()
    except AssertionError:
        return True
    raise AssertionError("An invalid certificate was accepted")


def psd_factor(matrix):
    """Verify PSD exactly and return rational K=R.H*diag(d)*R.

    The recursive Schur-complement test handles singular matrices exactly.
    Its pivots are positive rational numbers; a zero diagonal with a nonzero
    remaining matrix is rejected when no positive pivot remains.
    """
    matrix = clean(matrix)
    n = matrix.rows
    assert matrix.cols == n
    assert all(sp.re(z).is_Rational and sp.im(z).is_Rational for z in matrix), (
        "Certificate entries must be exact complex rationals"
    )
    assert is_zero(matrix - matrix.H), "Matrix is not Hermitian"
    if n == 0 or is_zero(matrix):
        return [], sp.zeros(0, n)
    diagonal = [sp.simplify(matrix[k, k]) for k in range(n)]
    assert all(x.is_real is True and x.is_nonnegative is True for x in diagonal), (
        "Negative or non-real diagonal"
    )
    pivot = next((k for k, x in enumerate(diagonal) if x.is_positive is True), None)
    assert pivot is not None, "Zero diagonal but nonzero Hermitian matrix"
    order = [pivot] + [k for k in range(n) if k != pivot]
    permuted = matrix.extract(order, order)
    a = permuted[0, 0]
    row = sp.Matrix([[1] + [sp.simplify(permuted[0, j] / a) for j in range(1, n)]])
    schur = clean(permuted[1:, 1:] - permuted[1:, :1] * permuted[:1, 1:] / a)
    tail_d, tail_r = psd_factor(schur)
    rows = row.col_join(sp.zeros(len(tail_d), 1).row_join(tail_r))
    restored = sp.zeros(rows.rows, n)
    for j, original in enumerate(order):
        restored[:, original] = rows[:, j]
    pivots = [a] + tail_d
    assert is_zero(matrix - restored.H * sp.diag(*pivots) * restored)
    return pivots, restored


def real_vector(matrix):
    flat = [matrix[i, j] for j in range(matrix.cols) for i in range(matrix.rows)]
    return sp.Matrix([sp.re(z) for z in flat] + [sp.im(z) for z in flat])


def action_constraints(R, D, rank, parameters):
    """Return exact independent A,b with A*realvec(Y)=b.

    Y columns are vec(B_i), where the physical action is x_i=sqrt(D)*B_i.
    Constraints are Hermitian realizability, centering, and local unbiasedness.
    """
    dimension = D.rows
    n = dimension * rank
    variables = sp.symbols(f"z0:{2*n*parameters}", real=True)
    Y = sp.zeros(n, parameters)
    for j in range(parameters):
        for i in range(n):
            ell = j * n + i
            Y[i, j] = variables[ell] + I * variables[n * parameters + ell]
    R0 = R[:, :rank]
    equations = []
    for i in range(parameters):
        B = sp.Matrix(dimension, rank, lambda a, b: Y[a + dimension*b, i])
        H = R0.H * D * B
        for a in range(rank):
            equations.append(sp.im(H[a, a]))
            for b in range(a + 1, rank):
                entry = H[a, b] - sp.conjugate(H[b, a])
                equations.extend([sp.re(entry), sp.im(entry)])
        equations.append(sp.re(sp.trace(H)))
        for j in range(parameters):
            Rj = R[:, (j+1)*rank:(j+2)*rank]
            equations.append(2*sp.re(sp.trace(Rj.H * D * B)) - int(i == j))
    A, b = sp.linear_eq_to_matrix(equations, variables)
    indices = list(A.T.rref()[1])
    A_reduced = A.extract(indices, range(A.cols))
    b_reduced = b.extract(indices, [0])
    assert A_reduced.rank() == A_reduced.rows
    return A_reduced, b_reduced, A, b


def verify_fit(K, Z, S, rank, parameters, nu, L_squared, gamma):
    psd_factor(K)
    psd_factor(K[:rank, :rank] - nu * sp.eye(rank) / 2)
    assert sp.trace(K[:rank, :rank]) == 1
    for i in range(parameters):
        block = K[:rank, (i+1)*rank:(i+2)*rank]
        assert sp.re(sp.trace(block)) == 0
    derivative_trace = sum(K[a, a] for a in range(rank, K.rows))
    assert derivative_trace <= 2 * L_squared
    psd_factor(K.row_join(Z).col_join(Z.H.row_join(S)))
    objective = sp.simplify(sp.trace(K) + sp.trace(S) - 2*sp.re(sp.trace(Z)))
    assert objective <= gamma**2 / 4
    return objective


def verify_holevo_certificate(D, R, rank, parameters, Y, V, dual_top):
    n = D.rows * rank
    inverse_metric = sp.kronecker_product(sp.eye(rank), D.inv())
    metric = inverse_metric.inv()
    A, b, all_A, all_b = action_constraints(R, D, rank, parameters)
    y = real_vector(Y)
    assert is_zero(all_A*y - all_b)
    assert V == V.T and all(sp.im(x) == 0 for x in V)
    primal = V.row_join(Y.H).col_join(Y.row_join(inverse_metric))
    psd_factor(primal)

    # This analytic construction produces a dual proposal. Every dual equation
    # and the full PSD condition are checked separately below.
    dual_bottom_left = -metric * Y * dual_top
    dual_bottom = metric * Y * dual_top * Y.H * metric
    dual = dual_top.row_join(dual_bottom_left.H).col_join(
        dual_bottom_left.row_join(dual_bottom)
    )
    psd_factor(dual)
    assert clean(sp.re(dual_top)) == sp.eye(parameters)
    rhs = -2 * real_vector(dual_bottom_left)
    multipliers, free = A.T.gauss_jordan_solve(rhs)
    assert free.rows == 0
    assert is_zero(A.T*multipliers + 2*real_vector(dual_bottom_left))
    upper = sp.trace(V)
    lower = sp.simplify((multipliers.T*b)[0] - sp.trace(dual_bottom*inverse_metric))
    assert lower <= upper
    return {
        "upper": str(upper),
        "lower": str(lower),
        "gap": str(sp.simplify(upper-lower)),
        "effective_hilbert_dimension": D.rows,
        "schur_order": parameters + n,
        "independent_affine_constraints": A.rows,
        "primal_psd_exact": True,
        "dual_psd_exact": True,
        "affine_constraints_exact": True,
    }


def example(c, s):
    assert c*c + s*s == 1
    rank, parameters = 1, 2
    if s == 0:
        V_columns = sp.Matrix([[1, 0, 0], [0, Q(1,2), I/2]])
        actions = sp.Matrix([[0, 0], [1, I]])
    else:
        V_columns = sp.Matrix([[1,0,0], [0,Q(1,2),I*c/2], [0,0,s/2]])
        t = c/(1+s)
        actions = sp.Matrix([[0,0], [1,I*t], [I*t,1]])
    K = clean(V_columns.H*V_columns)
    pivots, R = psd_factor(K)
    D = sp.diag(*pivots)
    square_root_D = sp.diag(*[sp.sqrt(d) for d in pivots])
    assert is_zero(square_root_D*R - V_columns)
    Y = clean(square_root_D.inv()*actions)
    t = c/(1+s)
    dual_top = sp.Matrix([[1,I*t], [-I*t,1]])
    cost = 4/(1+s)
    V_cov = (cost/2)*sp.eye(2)

    gamma = Q(1,100)
    m = K.rows
    regularizer = gamma**2/(64*m)
    S = K + regularizer*sp.eye(m)
    fit_objective = verify_fit(
        K, K, S, rank, parameters, Q(1), Q(1,2), gamma
    )
    certificate = verify_holevo_certificate(D, R, rank, parameters, Y, V_cov, dual_top)
    assert certificate["upper"] == str(cost)
    assert certificate["gap"] == "0"
    must_reject(lambda: verify_holevo_certificate(
        D, R, rank, parameters, Y, V_cov-Q(1,100)*sp.eye(2), dual_top
    ))
    must_reject(lambda: verify_fit(
        K, sp.zeros(m), S, rank, parameters, Q(1), Q(1,2), gamma
    ))
    certificate.update({
        "c": str(c), "s": str(s),
        "tangent_gram_order": m,
        "tangent_gram_rank": len(pivots),
        "tangent_gram_singular": len(pivots) < m,
        "fit_objective": str(fit_objective),
        "fit_objective_threshold": str(gamma**2/4),
        "fit_feasible_exact": True,
        "expected_holevo": str(cost),
        "sld_cost": "2",
        "false_upper_bound_rejected": True,
        "excess_fit_objective_rejected": True,
    })
    return certificate


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    must_reject(lambda: psd_factor(sp.Matrix([[0,1],[1,0]])))
    report = {
        "arithmetic": "exact complex rationals (SymPy)",
        "implements_generic_optimizer": False,
        "examples": [example(Q(1),Q(0)), example(Q(3,5),Q(4,5))],
        "all_checks_passed": True,
        "indefinite_zero_diagonal_matrix_rejected": True,
    }
    text = json.dumps(report, indent=2)
    if args.output:
        args.output.write_text(text + "\n", encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
