# Method Derivation

This document describes the Hermite wavelet collocation method used in this repository.

## Model Problem

We solve the third-order ordinary differential equation

$$
y'''(x)=3\sin(x), \qquad 0 \le x \le 1,
$$

with initial conditions

$$
y(0)=1, \qquad y'(0)=0, \qquad y''(0)=-2.
$$

The exact solution is

$$
y(x)=3\cos(x)+\frac{x^2}{2}-2.
$$

## Hermite Wavelet Approximation

The third derivative is approximated using Hermite wavelet basis functions:

$$
y_h'''(x)
=
\sum_{i=1}^{N_x}
\sum_{j=0}^{M_x-1}
c_{i,j}\psi_{i,j}(x).
$$

Here \(c_{i,j}\) are the unknown coefficients.

## Collocation System

At the collocation points \(x_r\), we impose

$$
y_h'''(x_r)=3\sin(x_r).
$$

Therefore,

$$
\sum_{i=1}^{N_x}
\sum_{j=0}^{M_x-1}
c_{i,j}\psi_{i,j}(x_r)
=
3\sin(x_r).
$$

This gives the algebraic system

$$
\Psi C = F,
$$

where

$$
\Psi_{r,(i,j)}=\psi_{i,j}(x_r),
$$

and

$$
F_r=3\sin(x_r).
$$

Solving this system gives the coefficient vector \(C\).

## Reconstruction of the Numerical Solution

The coefficients \(c_{i,j}\) approximate the third derivative \(y'''(x)\), not the solution \(y(x)\) directly.

To recover \(y_h(x)\), we integrate three times:

$$
y_h(x)
=$$
$$
y(0)+xy'(0)+\frac{x^2}{2}y''(0)
+
\int_0^x
\frac{(x-s)^2}{2}y_h'''(s)\,ds.
$$

Using the initial conditions,

$$
y(0)+xy'(0)+\frac{x^2}{2}y''(0)
=$$
$$
1-x^2.
$$

Thus,

$$
y_h(x)$$
$$
=
1-x^2
+
\sum_{i=1}^{N_x}
\sum_{j=0}^{M_x-1}
c_{i,j}
\int_0^x
\frac{(x-s)^2}{2}\psi_{i,j}(s)\,ds.
$$

## Third-Integration Matrix

At the collocation points, define the third-integration matrix \(P_3\) by

$$
(P_3)_{r,(i,j)}$$

$$ =
\int_0^{x_r}
\frac{(x_r-s)^2}{2}\psi_{i,j}(s)\,ds.
$$

Then the numerical solution vector is reconstructed as

$$
Y_h
=$$
$$
Y_{\mathrm{initial}} + P_3 C.
$$

This is the main reconstruction formula implemented in the code.
