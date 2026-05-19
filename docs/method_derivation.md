# Method Derivation

This document describes the Hermite wavelet collocation method used in this repository.

## Model Problem

We solve the third-order ordinary differential equation

$$
y'''(x)=3\sin(x), \qquad 0 \le x \le 1.
$$

The initial conditions are

$$
y(0)=1.
$$

$$
y'(0)=0.
$$

$$
y''(0)=-2.
$$

The exact solution is

$$
y(x)=3\cos(x)+\frac{x^2}{2}-2.
$$

## Hermite Wavelet Approximation

The third derivative is approximated using Hermite wavelet basis functions:

$$
y_h'''(x)$$

$$=
\sum_{i=1}^{N_x}
\sum_{j=0}^{M_x-1}
c_{i,j}\psi_{i,j}(x).
$$

Here,

$$
i
$$

is the interval index,

$$
j
$$

is the Hermite polynomial degree,

$$
c_{i,j}
$$

is the unknown wavelet coefficient, and

$$
\psi_{i,j}(x)
$$

is the Hermite wavelet basis function.

The goal of the collocation procedure is to determine the unknown coefficients

$$
c_{i,j}.
$$

## Collocation System

At the collocation points

$$
x_r,
$$

we impose the differential equation directly:

$$
y_h'''(x_r)=3\sin(x_r).
$$

Using the Hermite wavelet approximation, this becomes

$$
\sum_{i=1}^{N_x}
\sum_{j=0}^{M_x-1}
c_{i,j}\psi_{i,j}(x_r)$$

$$
= 3\sin(x_r).
$$

This gives the linear algebraic system

$$
\Psi C = F.
$$

The collocation matrix is defined by

$$
\Psi_{r,(i,j)}$$


$$=
\psi_{i,j}(x_r).
$$

The coefficient vector is

$$
C
=$$
$$
\begin{bmatrix}
c_{1,0} \\
c_{1,1} \\
\vdots \\
c_{i,j}
\end{bmatrix}.
$$

The right-hand side vector is

$$
F_r
=$$
$$
3\sin(x_r).
$$

Solving the system

$$
\Psi C = F
$$

gives the coefficient vector

$$
C.
$$

These coefficients approximate the third derivative

$$
y_h'''(x),
$$

not the solution

$$
y_h(x)
$$

directly.

## Reconstruction of the Numerical Solution

Since the coefficients approximate the third derivative, the numerical solution must be recovered by integrating three times.

From calculus, we have

$$
y_h(x)
=$$
$$
y(0)
+
xy'(0)
+
\frac{x^2}{2}y''(0)
+
\int_0^x
\frac{(x-s)^2}{2}
y_h'''(s)\,ds.
$$

Using the initial conditions,

$$
y(0)=1,
$$

$$
y'(0)=0,
$$

and

$$
y''(0)=-2,
$$

the initial-condition part becomes

$$
y(0)
+
xy'(0)
+
\frac{x^2}{2}y''(0)
=$$
$$
1-x^2.
$$

Therefore,

$$
y_h(x)
=$$
$$
1-x^2
+
\int_0^x
\frac{(x-s)^2}{2}
y_h'''(s)\,ds.
$$

Substituting the Hermite wavelet expansion of the third derivative gives

$$
y_h(x)
=$$
$$
1-x^2
+
\sum_{i=1}^{N_x}
\sum_{j=0}^{M_x-1}
c_{i,j}
\int_0^x
\frac{(x-s)^2}{2}
\psi_{i,j}(s)\,ds.
$$

This is the formula used to reconstruct the approximate solution.

## Third-Integration Matrix

At the collocation points, define the third-integration matrix

$$
P_3
$$

by

$$
(P_3)_{r,(i,j)}
=$$

$$
\int_0^{x_r}
\frac{(x_r-s)^2}{2}
\psi_{i,j}(s)\,ds.
$$

The initial-condition vector is

$$
Y_{\mathrm{initial}}
=$$
$$
\begin{bmatrix}
1-x_1^2 \\
1-x_2^2 \\
\vdots \\
1-x_r^2
\end{bmatrix}.
$$

The numerical solution vector is

$$
Y_h
=$$
$$
\begin{bmatrix}
y_h(x_1) \\
y_h(x_2) \\
\vdots \\
y_h(x_r)
\end{bmatrix}.
$$

The reconstruction formula is

$$
Y_h
=$$
$$
Y_{\mathrm{initial}}
+
P_3C.
$$

This is the main formula implemented in the code.

## Summary of the Computational Flow

The complete computational flow is

$$
y'''(x)=3\sin(x),
$$

$$
y_h'''(x)
=$$
$$
\sum_{i=1}^{N_x}
\sum_{j=0}^{M_x-1}
c_{i,j}\psi_{i,j}(x),
$$

$$
\Psi C = F,
$$

$$
C = \Psi^{-1}F,
$$

and finally

$$
Y_h
=$$
$$
Y_{\mathrm{initial}}
+
P_3C.
$$

Thus, the method first solves for the wavelet coefficients of the third derivative and then reconstructs the numerical solution by triple integration.
# Method Derivation

This document describes the Hermite wavelet collocation method used in this repository.

## Model Problem

We solve the third-order ordinary differential equation

$$y'''(x)=3\sin(x), \qquad 0 \le x \le 1.$$

The initial conditions are

$$y(0)=1, \qquad y'(0)=0, \qquad y''(0)=-2.$$

The exact solution is

$$y(x)=3\cos(x)+\frac{x^2}{2}-2.$$

## Hermite Wavelet Approximation

The third derivative is approximated using Hermite wavelet basis functions:

$$y_h'''(x)=\sum_{i=1}^{N_x}\sum_{j=0}^{M_x-1}c_{i,j}\psi_{i,j}(x).$$

Here, the interval index is

$$i.$$

The Hermite polynomial degree is

$$j.$$

The unknown wavelet coefficient is

$$c_{i,j}.$$

The Hermite wavelet basis function is

$$\psi_{i,j}(x).$$

The goal of the collocation procedure is to determine the unknown coefficients

$$c_{i,j}.$$

## Collocation System

At the collocation points

$$x_r,$$

we impose the differential equation directly:

$$y_h'''(x_r)=3\sin(x_r).$$

Using the Hermite wavelet approximation, this becomes

$$\sum_{i=1}^{N_x}\sum_{j=0}^{M_x-1}c_{i,j}\psi_{i,j}(x_r)=3\sin(x_r).$$

This gives the linear algebraic system

$$\Psi C=F.$$

The collocation matrix is defined by

$$\Psi_{r,(i,j)}=\psi_{i,j}(x_r).$$

The coefficient vector is

$$C=[c_{1,0},c_{1,1},\ldots,c_{i,j}]^T.$$

The right-hand side vector is

$$F_r=3\sin(x_r).$$

Solving the system

$$\Psi C=F$$

gives the coefficient vector

$$C.$$

These coefficients approximate the third derivative

$$y_h'''(x),$$

not the solution

$$y_h(x)$$

directly.

## Reconstruction of the Numerical Solution

Since the coefficients approximate the third derivative, the numerical solution must be recovered by integrating three times.

From calculus,

$$y_h(x)=y(0)+xy'(0)+\frac{x^2}{2}y''(0)+\int_0^x\frac{(x-s)^2}{2}y_h'''(s)\,ds.$$

Using the initial conditions,

$$y(0)=1, \qquad y'(0)=0, \qquad y''(0)=-2,$$

the initial-condition part becomes

$$y(0)+xy'(0)+\frac{x^2}{2}y''(0)=1-x^2.$$

Therefore,

$$y_h(x)=1-x^2+\int_0^x\frac{(x-s)^2}{2}y_h'''(s)\,ds.$$

Substituting the Hermite wavelet expansion of the third derivative gives

$$y_h(x)=1-x^2+\sum_{i=1}^{N_x}\sum_{j=0}^{M_x-1}c_{i,j}\int_0^x\frac{(x-s)^2}{2}\psi_{i,j}(s)\,ds.$$

This is the formula used to reconstruct the approximate solution.

## Third-Integration Matrix

At the collocation points, define the third-integration matrix

$$P_3.$$

Its entries are

$$(P_3)_{r,(i,j)}=\int_0^{x_r}\frac{(x_r-s)^2}{2}\psi_{i,j}(s)\,ds.$$

The initial-condition vector is

$$Y_{\mathrm{initial}}=[1-x_1^2,1-x_2^2,\ldots,1-x_r^2]^T.$$

The numerical solution vector is

$$Y_h=[y_h(x_1),y_h(x_2),\ldots,y_h(x_r)]^T.$$

The reconstruction formula is

$$Y_h=Y_{\mathrm{initial}}+P_3C.$$

This is the main formula implemented in the code.

## Summary of the Computational Flow

The complete computational flow is

$$y'''(x)=3\sin(x),$$

$$y_h'''(x)=\sum_{i=1}^{N_x}\sum_{j=0}^{M_x-1}c_{i,j}\psi_{i,j}(x),$$

$$\Psi C=F,$$

$$C=\Psi^{-1}F,$$

and finally

$$Y_h=Y_{\mathrm{initial}}+P_3C.$$

Thus, the method first solves for the wavelet coefficients of the third derivative and then reconstructs the numerical solution by triple integration.

