
# Newmann-Penrose Formalism {#Newmann-Penrose-Formalism}

The Newmann-Penrose (NP) Formalism is a set of notation that is ideal for expressing spacetime quantities in terms of a null-tetrad/null-spinoral basis.

## Notation {#Notation}

We will use Greek indices, $\{\alpha,\beta,\dots\}$, to index through space time indices in some coordinate system. The first half of the lower-case Latin alphabet, $\{a,\dots,l\}$, will be used to index the components of the tetrad associated Dyad basis, while the first half of the upper-case Latin alphabet, $\{A,\dots,L\}$, will be used to index through any other Dyad basis. Lastly, the second of of the lower-case Latin alphabet,$\{m,\dots,z\}$ will be used to index through the Tetrad basis components.

## Tetrad Calculus {#Tetrad-Calculus}

We follow [^NP] and define the tetrads as a basis of null vectors, $z_m{}^\mu=\langle l^\mu,n^\mu, m^\mu,\bar m^\mu\rangle_m$, which satisfy;

$$g^{\mu\nu}=\tilde\eta^{mn}z^{\mu}{}_mz^{\nu}{}_n,$$

where,

$$\tilde\eta^{mn}
    =\begin{pmatrix}
        0 & -1 & 0 & 0\\
        -1 & 0 & 0 & 0\\
        0 & 0 & 0 & 1\\
        0 & 0 & 1 & 0
    \end{pmatrix}^{mn}.$$

The $l^\mu$ tetrad component is usually taken to be the principal null direction of the Weyl tensor.

## Spinoral Calculus {#Spinoral-Calculus}

Spinors are defined interms of unitary projective representations of the Lorentz group. Spinors act on a 2D Complex vector space of spin vectors whose basis elements are known as dyads.  The spinors are often defined by their embeddings in vector representations of various spacetime tensors. In particular, the spinors, $\sigma^\mu{}_{AB}$, define an embedding into spacetime tensors through its action on the metric in the following way;[^PR]

$$g_{\mu\nu}=\epsilon_{AB}\epsilon_{A'B'}\sigma_{\mu}{}^{AA'}\sigma_{\nu}{}^{BB'},$$

or

$$g_{\mu\nu}\sigma^{\mu}{}_{AA'}\sigma^{\nu}{}_{BB'}=\epsilon_{AB}\epsilon_{A'B'}.$$

where,

$$\epsilon_{AB}=
\begin{pmatrix}
0 & 1\\
-1 & 0
\end{pmatrix}
\qquad\text{ and }\qquad
\epsilon_{A'B'}=
\begin{pmatrix}
0 & 1\\
-1 & 0
\end{pmatrix}
$$

are Levi-Civita symbols.  The $\epsilon_{AB}$ and $\epsilon_{A'B'}$ serve the purpose of the metric on spinor space, and act to raise and lower spinoral indices.

The notation is chosen to make a distinction between primed and unprimed indices.  This choice is conventient since it allows spinoral quantities that have pairs of primed and unprimed indices to be easily embedded into their tensoral equivalents with the use of the $\sigma^\mu{}_{AA'}$. In particular, spin vectors are embedded into vectors as;

$$\zeta^\mu=\sigma^{\mu}_{AA'}z^A\bar z^{A'},$$

and spinors into tensors as

$$\Chi^{\mu_1,\cdots,\mu_{N}}_{\nu_1,\cdots,\nu_{N}}=\left(\sigma^{\mu_1}_{A_1A_1'}\cdots\sigma^{\mu_n}_{A_n A_n'}\right)\left(\sigma_{\nu_1}^{B_1 B_1'}\cdots\sigma_{\nu_n}^{B_nB_n'}\right)\chi^{A_1\cdots A_n\,A_1'\cdots A_n'}_{B_1\cdots B_n\,B_1'\cdots B_n'}.$$

A canonical dyad basis, $\mathcal\zeta_{a}{}^A=(o^A,\iota^A)_a$, is normally chosen such that

$$\begin{align}
    l^\mu
        &=\sigma^\mu_{AA'}o^A\bar o^{A'}\\
    n^\mu
        &=\sigma^\mu_{AA'}\iota^A\bar \iota^{A'}\\
    m^\mu
        &=\sigma^\mu_{AA'}o^A\bar \iota^{A'}\\
    \bar m^\mu
        &=\sigma^\mu_{AA'}\iota^A\bar o^{A'}
\end{align}$$

where $\bar o^{A'}$, $\bar \iota^{A'}$ are the complex conjugates of $o^A$ and $\iota^A$.

## Type D SpaceTimes {#Type-D-SpaceTimes}

It is possible to perform a spacetime classification based on the eigen bi-vector space of the Weyl tensor, [^Petrov]

$$\frac{1}{2}C^{\alpha\beta}{}_{\gamma\delta}X^{\gamma\delta}=\lambda \, X^{\alpha\beta}.$$

The number eigen bi-vectors multiplicities define the spacetimes classes. type D spacetimes are of particular interst since this is the class that the Kerr family of metrics belong.

Certain properties of spacetimes are better ellucidated through the Newmann-Penrose formalism. In particular, one can show that that the Weyl spinor has only 1 scalar degree of freedom for the class of vacuum type D spacetimes.

The Weyl scalar $\Psi_{ABCD}$ is related to the Weyl tensor $C_{\alpha\beta\gamma\delta}$ through the relationship,

$$C_{\alpha\beta\gamma\delta}=\sigma_{\alpha}{}^{AA'}\sigma_{\beta}{}^{BB'}\sigma_{\gamma}{}^{CC'}\sigma_{\delta}{}^{DD'}(\Psi_{ABCD}\epsilon_{A'B'}\epsilon_{C'D'}+\bar\Psi_{A'B'C'D'}\epsilon_{AB}\epsilon_{CD}),$$

which for vacuum type D spacetimes, can be written in terms of a scalar $\psi_2$ through,

$$\Psi_{ABCD}=\psi_2o_{(A}o_{B}\iota_{C}\iota_{D)}.$$

The Kerr space time has [^AN]

$$\psi_2 = -\frac{M}{(r-ia\cos\theta)^3}.$$

## Ricci Rotation Coefficients {#Ricci-Rotation-Coefficients}

It is important to be able to define a connection and derivatives on tetrads and dyads. The covariant derivative on tetrads is adopted from its traditional definition on vectors. The covariant derivative on the tetrad basis is used to define the Ricci Rotation coefficients:

$$\gamma_m{}^{np}=z_{m\mu;\nu}z^{n\mu}z^{p\nu}$$

The covariant derivative on dyads however requires specification. The covariant derivative on dyads is traditionally defined to annihilate $\epsilon_{AC}$, $\epsilon_{B'D'}$ and $\sigma^{\mu}{}_{AB'}$. The resulting covariant derivative on the dyad basis is then,

$$\zeta_{aB;\mu}=\Gamma_{aB\mu}=\zeta_{a}{}^{A}\Gamma_{AB\mu}=\zeta_{a}{}^{C}\epsilon_{CA}\epsilon^{AD}\Gamma_{DB\mu}=-\zeta_{aA}\Gamma^{A}{}_{B\mu}.$$

The dyadic equivalent of the Ricci rotation coefficients are;

$$\Gamma_{abcd'}=\zeta_{aB;\mu}\zeta_{b}{}^B\sigma^{\mu}{}_{cd'}.$$

## References {#References}

[^NP]: Newman, Ezra, and Roger Penrose. “An Approach to Gravitational Radiation by a Method of Spin Coefficients.” Journal of Mathematical Physics, vol. 3, no. 3, AIP Publishing, 1 May 1962, pp. 566–578. Crossref, doi:10.1063/1.1724257.


[^PR]: Penrose, R., & Rindler, W. (1984). Spinors and Space-Time (Cambridge Monographs on Mathematical Physics). Cambridge: Cambridge University Press. doi:10.1017/CBO9780511564048


[^Petrov]: Petrov, A. Z. “The Classification of Spaces Defining Gravitational Fields.” General Relativity and Gravitation, vol. 32, no. 8, Springer Science and Business Media LLC, Aug. 2000, pp. 1665–1685. Crossref, doi:10.1023/a:1001910908054.


[^AN]: Adamo, T., & Newman, E. T. (2014). The Kerr-Newman metric: A Review. arXiv. https://doi.org/10.48550/ARXIV.1410.6626

