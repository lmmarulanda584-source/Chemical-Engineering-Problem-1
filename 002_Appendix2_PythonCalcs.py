import numpy as np

def solve_matrix():
    # ---------------- Given ----------------
    V4_L_per_h = 112_075.0
    Vm_STP_L_per_mol = 22.414   # adjust if your class uses 22.4, etc.
    purge_frac = 0.05
    ratio_B_to_E = 8.0
    selectivity = 20.0

    n7_total = 100.0
    y7_EB = 0.99
    y7_T  = 0.01

    y1_T = 0.01
    y1_B = 0.99

    # Known from specs
    n4_E = (V4_L_per_h / Vm_STP_L_per_mol) / 1000.0  # kmol/h
    n7_EB = y7_EB * n7_total
    n7_T  = y7_T  * n7_total

    # ---------------- Unknown vector x ----------------
    # 0  n1     total fresh feed
    # 1  n1B
    # 2  n1T
    # 3  n2E
    # 4  n10B
    # 5  n10T
    # 6  n3B
    # 7  n3E
    # 8  n3EB
    # 9  n3DEB
    # 10 n5B
    # 11 n5T
    # 12 n9B
    # 13 n9T
    # 14 n6EB
    # 15 xi1
    # 16 xi2
    #

    nvars = 17
    A = np.zeros((17, 17), dtype=float)
    b = np.zeros(17, dtype=float)

    def eq(row, coeffs, rhs=0.0):
        for j, v in coeffs:
            A[row, j] = v
        b[row] = rhs

    r = 0

    # (1) Stream 1 composition: n1B = 0.99 n1
    eq(r, [(1, 1.0), (0, -y1_B)], 0.0); r += 1
    # (2) Stream 1 composition: n1T = 0.01 n1
    eq(r, [(2, 1.0), (0, -y1_T)], 0.0); r += 1

    # (3) Flash+ D-101 ethylene: n3E = n4E (known)
    eq(r, [(7, 1.0)], n4_E); r += 1

    # (4) D-102 / product: EB in stream 6 equals EB product
    eq(r, [(14, 1.0)], n7_EB); r += 1

    # (5) Flash+ D-101 EB routing: n3EB = n6EB
    eq(r, [(8, 1.0), (14, -1.0)], 0.0); r += 1

    # (6) Reactor DEB balance (recycle DEB = 0): n3DEB = xi2
    eq(r, [(9, 1.0), (16, -1.0)], 0.0); r += 1

    # (7) Selectivity: xi1 - 20 xi2 = 0
    eq(r, [(15, 1.0), (16, -selectivity)], 0.0); r += 1

    # (8) Reactor EB balance (recycle EB = 0): -n3EB + xi1 - xi2 = 0
    eq(r, [(8, -1.0), (15, 1.0), (16, -1.0)], 0.0); r += 1

    # (9) Reactor ethylene balance: n2E - n3E - xi1 - xi2 = 0
    eq(r, [(3, 1.0), (7, -1.0), (15, -1.0), (16, -1.0)], 0.0); r += 1

    # (10) Reactor benzene balance: n1B + n10B - n3B - xi1 = 0
    eq(r, [(1, 1.0), (4, 1.0), (6, -1.0), (15, -1.0)], 0.0); r += 1

    # (11) Reactor feed ratio: (n1B + n10B) - 8 n2E = 0
    eq(r, [(1, 1.0), (4, 1.0), (3, -ratio_B_to_E)], 0.0); r += 1

    # (12) Flash+ D-101 benzene routing: n5B = n3B
    eq(r, [(10, 1.0), (6, -1.0)], 0.0); r += 1

    # (13) Split benzene: n9B - 0.05 n5B = 0
    eq(r, [(12, 1.0), (10, -purge_frac)], 0.0); r += 1

    # (14) Split benzene: n10B - 0.95 n5B = 0
    eq(r, [(4, 1.0), (10, -(1.0 - purge_frac))], 0.0); r += 1

    # (15) Split toluene: n9T - 0.05 n5T = 0
    eq(r, [(13, 1.0), (11, -purge_frac)], 0.0); r += 1

    # (16) Split toluene: n10T - 0.95 n5T = 0
    eq(r, [(5, 1.0), (11, -(1.0 - purge_frac))], 0.0); r += 1

    # (17) Overall T balance: n1T - n9T = n7T
    eq(r, [(2, 1.0), (13, -1.0)], n7_T); r += 1

    assert r == 17

    # ---------- Solve ----------
    # Optional debugging:
    # print("Rank(A) =", np.linalg.matrix_rank(A), " of ", nvars)

    x = np.linalg.solve(A, b)

    names = [
        "n1", "n1B", "n1T", "n2E", "n10B", "n10T",
        "n3B", "n3E", "n3EB", "n3DEB", "n5B", "n5T",
        "n9B", "n9T", "n6EB", "xi1", "xi2"
    ]
    sol = dict(zip(names, x))

    # totals and mole fractions for purge & recycle
    n9 = sol["n9B"] + sol["n9T"]
    n10 = sol["n10B"] + sol["n10T"]

    y10_T = sol["n10T"] / n10
    y10_B = sol["n10B"] / n10
    y9_T  = sol["n9T"] / n9
    y9_B  = sol["n9B"] / n9

    return sol, n9, n10, y10_T, y10_B, y9_T, y9_B

if __name__ == "__main__":
    sol, n9, n10, y10_T, y10_B, y9_T, y9_B = solve_matrix()

    print("=== Matrix solution (A x = b) ===")
    print(f"(i)  n1  = {sol['n1']:.3f} kmol/h")
    print(f"(ii) n9  = {n9:.3f} kmol/h")
    print(f"(iii)n10 = {n10:.3f} kmol/h")
    print("Compositions (mole fractions):")
    print(f"y10,T = {y10_T:.5f} ; y10,B = {y10_B:.5f}")
    print(f"y9,T  = {y9_T:.5f} ; y9,B  = {y9_B:.5f}")
    print("\n(Extra)")
    print(f"xi1 = {sol['xi1']:.5f}, xi2 = {sol['xi2']:.5f}")
    print(f"n3E (should equal off-gas) = {sol['n3E']:.6f} kmol/h")

