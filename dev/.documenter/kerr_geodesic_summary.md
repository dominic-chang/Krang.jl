
# The Geodesics of the Kerr Metric {#The-Geodesics-of-the-Kerr-Metric}

### Geodesic Equations {#Geodesic-Equations}

Geodesics are extremal solutions to the local optimization problem of the action:

$$\begin{equation}
S = \int_{\tilde\tau_i}^{\tilde\tau_f} d\tilde\tau\sqrt{-g_{\mu\nu}\dot x^\mu \dot x^\nu},
\end{equation}$$

defined on a manifold with metric $g_{\mu\nu}$, which in the case of the Kerr, can be written in Boyer-Lindquist coordinates[^BL] as:

$$\begin{equation}
g_{\mu\nu}=\left(
\begin{array}{cccc}
 \frac{r_s r}{\Sigma }-1 & 0 & 0 & -\frac{a^2   r_s r\sin ^2\theta}{\Sigma } \\
 0 & \frac{\Sigma }{\Delta } & 0 & 0 \\
 0 & 0 & \Sigma  & 0 \\
 -\frac{a^2 r  r_s\sin ^2\theta }{\Sigma } & 0 & 0 & \sin ^2\theta  \left(\frac{a^2
    r_s r\sin ^2\theta }{\Sigma }+a^2+r^2\right) \\
\end{array}
\right),
\end{equation}$$

where $\Sigma = r^2 +a^2\cos(\theta)^2$ and $\Delta=r^2- r_s r+a^2$.  In general, geodesics are constrained to satisfy the geodesic equations,

$$\begin{equation}
\frac{d^2x^\mu}{d\tilde\tau^2} =-\Gamma^\mu{}_{\alpha\beta}\dot x^\alpha \dot x^\beta,
\end{equation}$$

where,

$$\begin{equation}
\Gamma^\mu{}_{\alpha\beta} = \frac{1}{2}g^{\mu\nu}(g_{\alpha\nu, \beta} + g_{\beta\nu, \alpha} - g_{\alpha\beta, \nu})
\end{equation}$$

are the Christoffel symbols.

### Hamilton-Jacobi Equations {#Hamilton-Jacobi-Equations}

The Kerr family of metrics belongs to the class of type-D (type-{2,2}) spacetimes. These spacetimes admit a Killing tensor from which it is possible to define a Carter constant and Penrose-Walker constant.[^PW]  The Carter constant, in particular, allows solutions of the geodesic equations to be expressed in terms of quadratures through the Hamilton-Jacobi approach,[^Carter] with equations:

$$\begin{align}
H
    &=\frac{1}{2}g^{\mu\nu}p_{\mu}p_{\nu},\\
K
    &=\frac{\partial S}{\partial \tilde\tau} + H\left(x, \frac{\partial S}{\partial x}, \tilde\tau\right)=0,\\
p_\mu
    &=\frac{\partial S}{\partial x^\mu}.
\end{align}$$

The Carter constant is useful since it allows for seperable solutions to the Hamilton-Jacobi equations, resulting in a Hamiltonian principle function that takes the form,

$$\begin{equation}
S =  \int \frac{\sqrt{\mathcal R}}{\Delta}dr + \int\sqrt\Theta d\theta-Et +L\phi+ \frac{1}{2}\kappa\tilde\tau
\end{equation}$$

where $\kappa$ is the conserved four momentum magnitude, $E$ is the conserved energy, $L$ is the conserved azimuthal angular momentum, and

$$\begin{align}
\mathcal R
    &=\Delta[-C+\kappa r - (L-aE)^2]+[E(r^2+a^2)-La]^2,\\
\Theta
    &=C+\cos^2\theta\left((\kappa+E^2)a^2-\frac{L^2}{\sin^2\theta}\right)
\end{align}$$

are the radial and inclination potentials that depend on the previously mentioned conserved quantities along with the Carter constant $C$.

### Solutions by Quadrature {#Solutions-by-Quadrature}

Null geodesics have $\kappa=0$; thus they are only dependent on the reduced Carter constant and angular momentum:

$$\begin{align}
    \eta
        &=\frac{C}{E^2}\\
    \lambda
        &=\frac{L}{E}.
\end{align}$$

We use this fact to fix $E=1$ and write the  radial and inclination potentials for null geodesics as,

$$\begin{align}
\mathcal R
    &=-\Delta[\eta + (\lambda-a)^2]+[(r^2+a^2)-\lambda a]^2,\\
\Theta
    &=\eta+\cos^2\theta\left(a^2-\frac{\lambda^2}{\sin^2\theta}\right).
\end{align}$$

The Hamilton-Jacobi Equations then reduce to four first order differential equations: 

$$\begin{align}
\Sigma\dot r =\frac{d r}{d\tau}
    &= \sqrt{\mathcal R},\\
\Sigma\dot \theta=\frac{d \theta}{d\tau}
    &= \sqrt{\Theta},\\
\Sigma\dot \phi=\frac{d \phi}{d\tau}
    &= \frac{a}{\Delta}(r^2+a^2-a\lambda)+\frac{\lambda}{\sin^2\theta}-a,\\
\Sigma\dot t=\frac{d t}{d\tau}
    &= \frac{r^2+a^2}{\Delta}(r^2+a^2-a\lambda)+a(\lambda-a\sin^2\theta),
\end{align}$$

where we have defined the Mino time by its relationship to the affine parameter;  $d\tau=d\tilde\tau/\Sigma.$

Equations (15)-(18) define a 'Kepler-like' problem for photons in the Kerr spacetime.   The Mino time here acts similar to the Keplerian coordinate time in that problem, and can be expressed either in terms of $r$ or $\theta$ using equations (15) and (16);

$$\begin{align}
\tau=\tau_r(r_f,r_i)
    &=\oint\frac{dr}{\sqrt{R}},\text{ and}\\
\tau=\tau_\theta(\theta_f,\theta_i)
    &=\oint\frac{d\theta}{\sqrt{\Theta}}.
\end{align}$$

These equations are useful since $\tau_\theta$ and $\tau_r$ can be inverted to write,

$$\begin{align}
\theta_i
    &=\tau_\theta^{-1}(\tau,\theta_f)\\
r_i
    &=\tau_r^{-1}(\tau,r_f)\\

\end{align}$$

Thus for a given observer location, $(t_f, r_f, \theta_f, \phi_f)$, and observer screen coordinate $(\alpha(\lambda, \theta_f), \;\beta(\lambda,\eta,\theta_f))$, one can choose to parameterize the solutions of equations (15)-(18) by either $\Delta\tau,r_i$ or $\theta_i$.

Details on the solutions to equations (15)-(18) can be found in this review article by [Gralla & Lupsasca](https://doi.org/10.1103/PhysRevD.101.044032)[^GL].

## References {#References}

[^PW]: Walker, Martin, and Roger Penrose. ‘On Quadratic First Integrals of the Geodesic Equations for Type 22 Spacetimes’. Communications in Mathematical Physics, vol. 18, no. 4, Springer Science and Business Media LLC, Dec. 1970, pp. 265–274, https://doi.org10.1007/bf01649445.


[^Carter]: Carter, Brandon. ‘Global Structure of the Kerr Family of Gravitational Fields’. The Physical Review, vol. 174, no. 5, American Physical Society (APS), Oct. 1968, pp. 1559–1571, https://doi.org10.1103/physrev.174.1559.


[^BL]: Robert H. Boyer, Richard W. Lindquist; Maximal Analytic Extension of the Kerr Metric. J. Math. Phys. 1 February 1967; 8 (2): 265–281. https://doi.org/10.1063/1.1705193


[^GL]: Gralla, Samuel E., and Alexandru Lupsasca. ‘Null Geodesics of the Kerr Exterior’. Physical Review. D. (2016), vol. 101, no. 4, American Physical Society (APS), Feb. 2020, https://doi.org10.1103/physrevd.101.044032.

