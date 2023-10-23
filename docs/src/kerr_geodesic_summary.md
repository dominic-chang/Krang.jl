# The Geodesics of the Kerr Metric

### Geodesic Equations

Geodesics are extremal solutions to the local optimization problem of the action:

```math 
\begin{equation}
S = \int_{\tilde\tau_i}^{\tilde\tau_f} d\tilde\tau\sqrt{-g_{\mu\nu}\dot x^\mu \dot x^\nu},
\end{equation}
```

defined on a manifold with metric $g_{\mu\nu}$, which in the case of the Kerr, can be written in Boyer-Lindquist coordinates[^BL] as:

```math
\begin{equation}
g_{\mu\nu}=\left(
\begin{array}{cccc}
 \frac{r_s r}{\Sigma }-1 & 0 & 0 & -\frac{a^2   r_s r\sin ^2\theta}{\Sigma } \\
 0 & \frac{\Sigma }{\Delta } & 0 & 0 \\
 0 & 0 & \Sigma  & 0 \\
 -\frac{a^2 r  r_s\sin ^2\theta }{\Sigma } & 0 & 0 & \sin ^2\theta  \left(\frac{a^2
    r_s r\sin ^2\theta }{\Sigma }+a^2+r^2\right) \\
\end{array}
\right),
\end{equation}
```
where $\Sigma = r^2 +a^2\cos(\theta)^2$ and $\Delta=r^2- r_s r+a^2$. 
In general, geodesics are constrained to satisfy the geodesic equations,

```math
\begin{equation}
\frac{d^2x^\mu}{d\tilde\tau^2} =-\Gamma^\mu{}_{\alpha\beta}\dot x^\alpha \dot x^\beta,
\end{equation}
```

where,
```math
\begin{equation}
\Gamma^\mu{}_{\alpha\beta} = \frac{1}{2}g^{\mu\nu}(g_{\alpha\nu, \beta} + g_{\beta\nu, \alpha} - g_{\alpha\beta, \nu})
\end{equation}
```
are the Christoffel symbols.

### Hamilton-Jacobi Equations

The Kerr family of metrics belongs to the class of type-D (type-{2,2}) spacetimes.
These spacetimes admit a Killing-Yano tensor [^PW] from which it is possible to define a Carter constant and Penrose-Walker constant. 
The Carter constant, in particular, allows solutions of the geodesic equations to be expressed in terms of quadratures through the Hamilton-Jacobi approach[^Carter], whose equations are given by:

```math
\begin{align}
H
    &=\frac{1}{2}g^{\mu\nu}p_{\mu}p_{\nu},\\
K
    &=\frac{\partial S}{\partial \tilde\tau} + H\left(x, \frac{\partial S}{\partial x}, \tilde\tau\right)=0,\\
p_\mu
    &=\frac{\partial S}{\partial x^\mu}.
\end{align}
```

The Carter constant is useful since it lets the Hamilton-Jacobi equations be seperable, resulting in a Hamiltonian principle function that takes the form,

```math
\begin{equation}
S =  \int \frac{\sqrt{\mathcal R}}{\Delta}dr + \int\sqrt\Theta d\theta-Et +L\phi+ \frac{1}{2}\kappa\tilde\tau
\end{equation}
```

where $\kappa$ is the conserved four momentum magnitude, $E$ is the conserved energy, $L$ is the conserved azimuthal angular momentum, and

```math
\begin{align}
\mathcal R
    &=\Delta[-C+\kappa r - (L-aE)^2]+[E(r^2+a^2)-La]^2,\\
\Theta
    &=C+\cos^2\theta\left((\kappa+E^2)a^2-\frac{L^2}{\sin^2\theta}\right)
\end{align}
```
are the radial and inclination potentials that depend on the previously mentioned conserved quantities along with the Carter constant $C$.

### Integral solutions

A noteworthy fact is that null geodesics, $\kappa=0$, only depend on the reduced Carter constant and angular momentum

```math
\begin{align}
    \lambda
        &=\frac{L}{E}\\
    \eta
        &=\frac{C}{E^2}.
\end{align}
```
We use this freedom to fix $E=1$ and write the 
radial and inclination potentials for null geodesics as,
```math
\begin{align}
\mathcal R
    &=-\Delta[\eta + (\lambda-a)^2]+[(r^2+a^2)-\lambda a]^2,\\
\Theta
    &=\eta+\cos^2\theta\left(a^2-\frac{\lambda^2}{\sin^2\theta}\right).
\end{align}
```
The Hamilton-Jacobi Equations then reduce to four first order differential equations: 

```math
\begin{align}
\frac{\Sigma}{E}\dot r
    &= \sqrt{\mathcal R},\\
\frac{\Sigma}{E}\dot \theta
    &= \sqrt{\Theta},\\
\frac{\Sigma}{E}\dot \phi
    &= \frac{a}{\Delta}(r^2+a^2-a\lambda)+\frac{\lambda}{\sin^2\theta}-a,\\
\frac{\Sigma}{E}\dot t
    &= \frac{r^2+a^2}{\Delta}(r^2+a^2-a\lambda)+a(\lambda-a\sin^2\theta),\\

\end{align}
```
a fact which follows from the relationship, $\dot x^\mu=p^\mu=g^{\mu\nu}\partial_\nu S$.
Equations (15)-(18) define a 'Kepler-like' problem where a new parameter, known as the Mino time $\tau$, is defined in terms of an integral in $r$ or $\theta$, and acts similar to the Keplerian coordinate time in that problem.
The Mino time is defined by its relationship to the affine parameter by $d\tau=d\tilde\tau/\Sigma$.
Furthermore, Equations (15) and (16) imply that the Mino time can be written as a function of $r$ or $\theta$ through the relationships,
```math
\begin{align}
\tau=\tau_r(r_f,r_i)
    &=\oint\sqrt{R}\;dr,\text{ and}\\
\tau=\tau_\theta(\theta_f,\theta_i)
    &=\oint\sqrt{\Theta}\;d\theta.
\end{align}
```

These equations are useful since $\tau_\theta$ and $\tau_r$ can be inverted to write,
```math
\begin{align}
\theta_i
    &=\tau_\theta^{-1}(\tau,\theta_f)\\
r_i
    &=\tau_r^{-1}(\tau,r_f)\\

\end{align}
```

Thus for a given $t_f, r_f, \theta_f, \phi_f, \lambda$ and $\eta$, one can choose to parameterize the solutions of equations (15)-(18) by either $\Delta\tau,r_i$ or $\theta_i$.

More details on this formalism and the solutions can be found in this review artivle by [Gralla & Lupsasca](https://doi.org/10.1103/PhysRevD.101.044032)[^GL]


## References
[^PW]: Walker, M., Penrose, R. On quadratic first integrals of the geodesic equations for type {22} spacetimes. Commun.Math. Phys. 18, 265–274 (1970).

[^Carter]: Brandon Carter Phys. Rev. 174, 1559 – Published 25 October 1968

[^BL]: Robert H. Boyer, Richard W. Lindquist; Maximal Analytic Extension of the Kerr Metric. J. Math. Phys. 1 February 1967; 8 (2): 265–281.

[^GL]: Samuel E. Gralla and Alexandru Lupsasca Phys. Rev. D 101, 044032