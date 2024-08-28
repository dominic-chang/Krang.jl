
# Linear-Polarization

## Linear Polarization of null geodesics in the Newmann-Penrose Formalism {#Linear-Polarization-of-null-geodesics-in-the-Newmann-Penrose-Formalism}

The Newmann-Penrose formalism allows the linear polarization of a null vector to be modelled with the definition of a flag plane. Here, we provide a summarized construction.[^PR]

Let $k^\mu$ be a null vector. Then there exist a family of spin vectors $\kappa^A$ that can be embedded into $k^\mu$ as,

$$k^\mu=\sigma^\mu{}_{AA'}\kappa^A\bar\kappa^{A'}.$$

There exists a set of equivalent spin vectors for each vector. These spin vectors are related to each other through a complex phase $\kappa\rightarrow e^{i\theta}\kappa$. The goal of this construction is to show that the orientation of the linear polarization can be encoded in the phase.

The relationship between linear polarization and the phase is first acheived by constructing a polarization tensor from $\kappa^A$ and $\epsilon^{AB}$ as,

$$P^{\mu\nu}=\sigma^{\mu}{}_{AA'}\sigma^{\nu}{}_{BB'}\left(
    \kappa^{A}\kappa^B\epsilon^{A'B'}+
\epsilon^{AB}\bar\kappa^{A'}\bar\kappa^{B'}\right).$$

Both $\kappa^{A}\kappa^B$ and $\bar\kappa^{A}\bar\kappa^B$ are included in the defition of $P^{\mu\nu}$ to ensure that it is real. Then in general, $P^{\mu\nu}$ can be written as the antisymmetric outer product of $k^\mu=\sigma^\mu_{AA'}\kappa^A\bar\kappa^{A'}$ with another vector spacelike vector $f^\mu = \sigma^{\mu}{}_{AA'}(\kappa^A\bar\tau^{A'}+\tau^A\bar\kappa^{A'})$,

$$P^{\mu\nu}
=k^\mu f^\nu-f^\mu k^\nu,$$

a fact which follows from expressing $\epsilon^{AB}$ as,

$$\epsilon^{AB}=\kappa^A\tau^{B}-\tau^A\kappa^{B}.$$

From this consideration, we can define another spacelike vector, $p^\mu= \sigma^{\mu}{}_{AA'}\left(i\,\kappa^A\bar\tau^{A'}-i\,\tau^A\bar\kappa^{A'}\right)$, which satisfies $p_\mu f^\mu=p_\mu k^\mu=0$. The vectors $f^\mu$ and $p^\mu$ span a flag plane that lies orthogonal to null geodesics with $k^\mu$ as their tangent vector.

Reproducing the above construction with a phase rotated spinor, $\kappa\rightarrow e^{i\theta}\kappa$, induces a phase rotation into the definition of $\tau$, and physical rotation of the flag plane, $(f,p)$,

$$f^\mu\rightarrow \cos(2\theta) f^\mu+\sin(2\theta) p^\mu.$$

Thus a polarization orientation can be encoded as a phase rotation in $\kappa^A=e^{i\theta}\kappa_0$ relative to some reference spinor $\kappa_0$.

## Walker-Penrose Constant and Killing Spinors of type D spacetimes {#Walker-Penrose-Constant-and-Killing-Spinors-of-type-D-spacetimes}

It can be shown that type-$D$ spacetimes admits a Killing spinor,[^NP]

$$\chi_{AB}=\psi_2^{-1/3}\sigma_{(A}\iota_{B)},$$

with

$$\psi_2 = -\frac{M}{(r-ia\cos\theta)^3},$$

which satisfies the spinoral Killing equation,

$$\nabla^{A'}_{(A}\chi_{BC)} = 0.$$

The Walker-Penrose constant is then defined by first constructing the complex tensor,

$$J_{\mu\nu}
=\sigma_\mu^{AA'}\sigma_\nu^{BB'}\chi_{AB}\epsilon_{A'B'},$$

which in Boyer-Lindquist coordinates takes the form[^AN][^GKL]

$$J_{\mu\nu}=2(r-ia\cos\theta)\left(l_{[\mu}n_{\nu]}-m_{[\mu}\bar m_{\nu]}\right).$$

Contracting $J_{\mu\nu}$ with a null vector $k^\mu$ and its polarization vector $f^\nu$ results in a conserved quantity known as the Walker-Penrose constant.[^CS]

$$\kappa=(A+iB)(r-ia\cos\theta)$$

where

$$\begin{align}
    A
        &=k^tf^r-k^rf^t+a\sin^2\theta(k^r f^\phi-k^\phi f^r)\\
    B
        &=\left[(r^2+a^2)(k^\phi f^\theta-k^\theta f^\phi)-a*(k^t f^\theta - k^\theta f^t)\right]\sin\theta.
\end{align}$$

## References

[^PR]: Penrose, R., &amp; Rindler, W. (1984). Spinors and Space-Time (Cambridge Monographs on Mathematical Physics). Cambridge: Cambridge University Press. doi:10.1017/CBO9780511564048


[^NP]: Ezra Newman, Roger Penrose; An Approach to Gravitational Radiation by a Method of Spin Coefficients. J. Math. Phys. 3, 566–578 (1962)


[^CS]: CONNORS, P. A., &amp; STARK, R. F. (1977). Observable gravitational effects on polarised radiation coming from near a black hole. In Nature (Vol. 269, Issue 5624, pp. 128–129). Springer Science and Business Media LLC. https://doi.org/10.1038/269128a0


[^AN]: Adamo, T., &amp; Newman, E. T. (2014). The Kerr-Newman metric: A Review. arXiv. https://doi.org/10.48550/ARXIV.1410.6626


[^GKL]: Gates, D., Kapec, D., Lupsasca, A., Shi, Y., &amp; Strominger, A. (2018). Polarization Whorls from M87* at the Event Horizon Telescope. arXiv. https://doi.org/10.48550/ARXIV.1809.09092

