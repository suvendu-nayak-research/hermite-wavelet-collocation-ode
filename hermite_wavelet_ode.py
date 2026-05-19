import math
import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt


# ============================================================
# 0. USER SETTINGS: CHANGE ONLY THIS BLOCK
# ============================================================

# Hermite wavelet parameters
M_x = 8
k = 1
N_x = 2 ** (k - 1)

# Domain
a = 0.0
b = 1.0

# Initial conditions at x = a
y0 = 1.0          # y(a)
dy0 = 0.0         # y'(a)
d2y0 = -2.0       # y''(a)

# Printing and plotting control
PRINT_MATRICES = False
MAKE_PLOTS = True


# ============================================================
# 1. Hermite polynomial H_j(z)
# ============================================================

def hermite_polynomial(j, z):
    """
    Compute the physicists' Hermite polynomial H_j(z).

    H_0(z) = 1
    H_1(z) = 2z

    For j >= 2:

        H_j(z) = 2z H_{j-1}(z) - 2(j-1) H_{j-2}(z)
    """

    if j < 0:
        raise ValueError("j must be non-negative.")

    if j == 0:
        return 1.0

    if j == 1:
        return 2.0 * z

    H_j_minus_2 = 1.0
    H_j_minus_1 = 2.0 * z

    for n in range(2, j + 1):
        H_j_current = (
            2.0 * z * H_j_minus_1
            - 2.0 * (n - 1) * H_j_minus_2
        )

        H_j_minus_2 = H_j_minus_1
        H_j_minus_1 = H_j_current

    return H_j_minus_1


# ============================================================
# 2. Normalization constant A_{i,j}
# ============================================================

def A_ij(j):
    """
    Compute A_{i,j}.

    General Hermite wavelet normalization:

        A_{i,j} = sqrt( 2 N_x / (2^j j! sqrt(pi)) )

    For N_x = 1:

        A_{i,j} = A_{1,j}
                = sqrt( 2 / (2^j j! sqrt(pi)) )
    """

    numerator = 2.0 * N_x
    denominator = (2.0 ** j) * math.factorial(j) * math.sqrt(math.pi)

    return math.sqrt(numerator / denominator)


# ============================================================
# 3. Support interval I_i
# ============================================================

def interval_support(i):
    """
    Compute the support interval I_i of psi_{i,j}(x).

    For i = 1,2,...,N_x:

        I_i = [a + (i-1)(b-a)/N_x, a + i(b-a)/N_x]
    """

    left = a + (i - 1) * (b - a) / N_x
    right = a + i * (b - a) / N_x

    return left, right


# ============================================================
# 4. Local coordinate z_i(x)
# ============================================================

def local_coordinate(i, x):
    """
    Compute the local coordinate z_i(x).

    The general mapping from I_i = [left,right] to [-1,1] is:

        z_i = 2(x-left)/(right-left) - 1

    For [0,1] with N_x subintervals:

        z_i = 2 N_x x - 2i + 1

    For N_x = 1 and i = 1:

        z_1 = 2x - 1
    """

    left, right = interval_support(i)

    z = 2.0 * (x - left) / (right - left) - 1.0

    return z


# ============================================================
# 5. Hermite wavelet basis psi_{i,j}(x)
# ============================================================

def psi_ij(i, j, x):
    """
    Compute psi_{i,j}(x).

    Definition:

        psi_{i,j}(x) = A_{i,j} H_j(z_i(x))

    on the interval I_i.

    Outside I_i:

        psi_{i,j}(x) = 0
    """

    left, right = interval_support(i)

    if x < left or x > right:
        return 0.0

    z = local_coordinate(i, x)

    return A_ij(j) * hermite_polynomial(j, z)


# ============================================================
# 6. Basis index list
# ============================================================

def build_basis_indices():
    """
    Build all basis index pairs (i,j).

    For N_x = 1, M_x = 8:

        [(1,0), (1,1), ..., (1,7)]

    For N_x = 2, M_x = 8:

        [(1,0), ..., (1,7), (2,0), ..., (2,7)]

    Each pair (i,j) corresponds to one coefficient c_{i,j}.
    """

    basis_indices = []

    for i in range(1, N_x + 1):
        for j in range(M_x):
            basis_indices.append((i, j))

    return basis_indices


# ============================================================
# 7. Collocation points x_r
# ============================================================

def collocation_points():
    """
    Compute collocation points.

    Total unknown coefficients:

        total_basis = N_x M_x

    Therefore we use total_basis collocation points:

        x_r = a + (2r - 1)(b-a)/(2 total_basis),
        r = 1,2,...,total_basis.

    For the benchmark:

        N_x = 1, M_x = 8

    so:

        x_r = (2r - 1)/16.
    """

    total_basis = N_x * M_x

    r = np.arange(1, total_basis + 1)

    points = a + (2.0 * r - 1.0) * (b - a) / (2.0 * total_basis)

    return points


# ============================================================
# 8. Build collocation matrix Psi
# ============================================================

def build_Psi(points, basis_indices):
    """
    Build the collocation matrix Psi.

    General definition:

        Psi[row, col] = psi_{i,j}(x_r)

    where column col corresponds to the basis pair (i,j).
    """

    num_points = len(points)
    num_basis = len(basis_indices)

    Psi = np.zeros((num_points, num_basis))

    for row, x_r in enumerate(points):
        for col, (i, j) in enumerate(basis_indices):
            Psi[row, col] = psi_ij(i, j, x_r)

    return Psi


# ============================================================
# 9. Right-hand side of the ODE
# ============================================================

def rhs_function(x):
    """
    Right-hand side of the ODE:

        y'''(x) = 3 sin(x)
    """

    return 3.0 * np.sin(x)


def build_F(points):
    """
    Build F using:

        F_r = 3 sin(x_r)
    """

    F = rhs_function(points)

    return F


# ============================================================
# 10. Exact solution
# ============================================================

def y_exact(x):
    """
    Exact solution:

        y(x) = 3 cos(x) + x^2/2 - 2
    """

    return 3.0 * np.cos(x) + 0.5 * x ** 2 - 2.0


# ============================================================
# 11. Initial-condition contribution
# ============================================================

def initial_part(x):
    """
    Contribution from initial conditions after integrating y''' three times:

        y(a) + (x-a)y'(a) + ((x-a)^2/2)y''(a)

    For this benchmark, a = 0:

        1 + x(0) + x^2/2 (-2) = 1 - x^2
    """

    return y0 + (x - a) * dy0 + 0.5 * (x - a) ** 2 * d2y0


# ============================================================
# 12. One entry of P_3
# ============================================================

def P_3_entry(x_r, i, j):
    """
    Compute one entry of P_3.

    Mathematically:

        (P_3)_{r,(i,j)}
        =
        integral_a^{x_r}
        [ ((x_r - s)^2)/2 * psi_{i,j}(s) ] ds

    But psi_{i,j}(s) is nonzero only on I_i = [left, right].

    Therefore we integrate only over the overlap:

        [a, x_r] intersection [left, right]

    Since in this code left >= a, the useful interval is:

        [left, min(x_r, right)]
    """

    left, right = interval_support(i)

    integration_left = left
    integration_right = min(x_r, right)

    if integration_right <= integration_left:
        return 0.0

    def integrand(s):
        kernel = 0.5 * (x_r - s) ** 2
        basis_value = psi_ij(i, j, s)
        return kernel * basis_value

    value, estimated_error = quad(
        integrand,
        integration_left,
        integration_right,
        epsabs=1e-12,
        epsrel=1e-12
    )

    return value


# ============================================================
# 13. Build third-integration matrix P_3
# ============================================================

def build_P_3(points, basis_indices):
    """
    Build the third-integration matrix P_3.

    Each entry is:

        (P_3)_{r,(i,j)}
        =
        integral_a^{x_r}
        [ ((x_r - s)^2)/2 * psi_{i,j}(s) ] ds

    For N_x = 1, this reduces to:

        (P_3)_{r,j}
        =
        integral_0^{x_r}
        [ ((x_r - s)^2)/2 * psi_{1,j}(s) ] ds
    """

    num_points = len(points)
    num_basis = len(basis_indices)

    P_3 = np.zeros((num_points, num_basis))

    for row, x_r in enumerate(points):
        for col, (i, j) in enumerate(basis_indices):
            P_3[row, col] = P_3_entry(x_r, i, j)

    return P_3


# ============================================================
# 14. Approximate solution y_h(x)
# ============================================================

def y_h_value(x, C, basis_indices):
    """
    Compute y_h(x).

    General formula:

        y_h(x)
        =
        initial_part(x)
        +
        sum_{i,j} c_{i,j}
        integral_a^x [ ((x-s)^2)/2 * psi_{i,j}(s) ] ds

    For N_x = 1:

        y_h(x)
        =
        1 - x^2
        +
        sum_{j=0}^{M_x-1} c_{1,j}
        integral_0^x [ ((x-s)^2)/2 * psi_{1,j}(s) ] ds
    """

    value = initial_part(x)

    for coefficient, (i, j) in zip(C, basis_indices):
        value += coefficient * P_3_entry(x, i, j)

    return value


# ============================================================
# 15. Main program
# ============================================================

def main():
    np.set_printoptions(precision=10, suppress=True)

    print("\nHermite wavelet collocation method")
    print("----------------------------------")
    print(f"M_x = {M_x}")
    print(f"k   = {k}")
    print(f"N_x = 2^(k-1) = {N_x}")
    print(f"Domain = [{a}, {b}]")

    total_basis = N_x * M_x
    print(f"Total unknown coefficients = N_x * M_x = {total_basis}")

    # --------------------------------------------------------
    # Step 1: Build basis index list
    # --------------------------------------------------------

    basis_indices = build_basis_indices()

    print("\nBasis index mapping:")
    for col, (i, j) in enumerate(basis_indices):
        print(f"column {col:2d} -> c_{{{i},{j}}}, psi_{{{i},{j}}}(x)")

    # --------------------------------------------------------
    # Step 2: Build collocation points
    # --------------------------------------------------------

    points = collocation_points()

    print("\nCollocation points x_r:")
    print(points)

    # --------------------------------------------------------
    # Step 3: Build Psi
    # --------------------------------------------------------

    Psi = build_Psi(points, basis_indices)

    if PRINT_MATRICES:
        print("\nCollocation matrix Psi:")
        print(Psi)

    print("\nShape of Psi:")
    print(Psi.shape)

    condition_number = np.linalg.cond(Psi)
    print("\nCondition number of Psi:")
    print(condition_number)

    # --------------------------------------------------------
    # Step 4: Build F
    # --------------------------------------------------------

    F = build_F(points)

    print("\nRight-hand side vector F:")
    print(F)

    # --------------------------------------------------------
    # Step 5: Solve Psi C = F
    # --------------------------------------------------------

    C = np.linalg.solve(Psi, F)

    print("\nCoefficient vector C:")
    print(C)

    print("\nCoefficient mapping:")
    for coefficient, (i, j) in zip(C, basis_indices):
        print(f"c_{{{i},{j}}} = {coefficient:.16e}")

    # --------------------------------------------------------
    # Step 6: Build P_3
    # --------------------------------------------------------

    P_3 = build_P_3(points, basis_indices)

    if PRINT_MATRICES:
        print("\nThird-integration matrix P_3:")
        print(P_3)

    print("\nShape of P_3:")
    print(P_3.shape)

    # --------------------------------------------------------
    # Step 7: Compute numerical solution vector Y_h
    # --------------------------------------------------------

    initial_vector = initial_part(points)

    Y_h = initial_vector + P_3 @ C

    Y_exact = y_exact(points)

    absolute_error = np.abs(Y_exact - Y_h)

    print("\nApproximate solution vector Y_h:")
    print(Y_h)

    print("\nNumerical solution y_h(x) at collocation points:")
    print("-" * 60)
    print(f"{'x':>15} {'y_h(x)':>25}")
    print("-" * 60)

    for x_value, numerical_value in zip(points, Y_h):
        print(f"{x_value:15.8f} {numerical_value:25.14f}")

    print("-" * 60)

    print("\nExact solution values:")
    print(Y_exact)

    print("\nAbsolute errors:")
    print(absolute_error)

    # --------------------------------------------------------
    # Step 8: Clean result table
    # --------------------------------------------------------

    print("\nTable: Hermite wavelet collocation result")
    print("-" * 82)
    print(
        f"{'x':>12} "
        f"{'exact y(x)':>22} "
        f"{'Hermite wavelet y_h(x)':>28} "
        f"{'absolute error':>18}"
    )
    print("-" * 82)

    for x_value, exact_value, approx_value, error_value in zip(
        points,
        Y_exact,
        Y_h,
        absolute_error
    ):
        print(
            f"{x_value:12.8f} "
            f"{exact_value:22.14f} "
            f"{approx_value:28.14f} "
            f"{error_value:18.6e}"
        )

    print("-" * 82)

    # --------------------------------------------------------
    # Step 9: Plots
    # --------------------------------------------------------

    if MAKE_PLOTS:
        x_plot = np.linspace(a, b, 400)

        y_exact_plot = y_exact(x_plot)

        y_h_plot = np.array([
            y_h_value(x, C, basis_indices)
            for x in x_plot
        ])

        plt.figure(figsize=(8, 5))
        plt.plot(x_plot, y_exact_plot, label="Exact solution")
        plt.plot(x_plot, y_h_plot, "--", label="Hermite wavelet solution")
        plt.scatter(points, Y_h, label="Collocation values")
        plt.xlabel("x")
        plt.ylabel("y(x)")
        plt.title("Exact solution vs Hermite wavelet collocation solution")
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.savefig("solution_comparison.png", dpi=300)
        plt.show()

        error_plot = np.abs(y_exact_plot - y_h_plot)

        plt.figure(figsize=(8, 5))
        plt.plot(x_plot, error_plot)
        plt.scatter(points, absolute_error)
        plt.xlabel("x")
        plt.ylabel("absolute error")
        plt.title("Absolute error of Hermite wavelet collocation solution")
        plt.grid(True)
        plt.tight_layout()
        plt.savefig("absolute_error.png", dpi=300)
        plt.show()


# ============================================================
# 16. Run the program
# ============================================================

if __name__ == "__main__":
    main()
