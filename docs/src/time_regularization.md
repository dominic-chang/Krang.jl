# Regularized time

The Kerr spacetime admits a unique Killing vector field, $\xi^t$, that is timelike at $r=\infty$, and whose integral curves serve as a natural time coordinate for observers sitting far away from the black hole.
This is the typical time coordinate used for the Kerr when expressed in Boyer-Lindquist form.

It is then possible to calculate the total coordinate time elapsed time over the trajectory of any null geodesic $\gamma(r,\theta)$ as,

```math
\begin{align}
    \Delta t
        &= \Delta I_t+a^2\Delta G_t
\end{align}
```
where,
```math
\begin{align}
    \Delta I_t(r_o,r_s)
        &= I_t(r_o) - I_t(r_s)
        =\int
    _{r_s}^{r_o}\frac{r^2\Delta(r)+2Mr(r^2+a^2-a\lambda)}{\Delta(r)\sqrt{\mathcal R(r)}}\,dr,\\
    \Delta G_t(\theta_o, \theta_s)
        &=G_t(\theta_o)-G_t(\theta_s)
        =\int_{\theta_s}^{\theta_o}\frac{\cos^2\theta}{\sqrt{\Theta(\theta)}}\,d\theta,\\
    \Delta(r)
        &=r^2-2Mr+a^2,\\
    \mathcal R(r)
        &= (r^2+a^2-a\lambda)^2-\Delta(r)[\eta+(\lambda-a)^2], \text{ and}\\
    \Theta(\theta)
        &=\eta+a^2\cos^2\theta-\lambda^2\cot^s\theta.
\end{align}
```

One idiosyncracy of this coordinate is that it introduces ambiguities when evaluating the elapsed time of any geodesic terminating at the asymptotic observer.
The ambiguity is due to logarithmic and linear divergences in the integral $I_t$ at $r_o=\infty$ which take the form,[^CAZ]
```math
\begin{align}
    \lim_{r\rightarrow\infty} \frac{dt}{dr}\approx \frac{2M}{r}+1,
\end{align}
```
These divergences are indepedent of the nature of the geodesicas.

One way around these ambiguities is for the observer to ask questions that have a regulating effect on the divergences.
An example of a valid question is; "what is the relative time delay of the arival of a photon, $i$, with respect to a reference photon, $j$";
```math
\begin{align}
    t_i-t_j
        &=\Delta t_i - \Delta t_j + t_{si} - t_{sj},
\end{align}
```
where $t_{si}$ and $t_{sj}$ are the emission times of photons $i$ and $j$ respectively. 

Another approch is to avoid the divergent quality of $I_t$ altogether by introducing a regularized version of the integral, $\tilde I_t$, defined as;
```math
\begin{align}
    \Delta \tilde I_t(r_s,\infty)
        &=\lim_{r_o\rightarrow\infty} \int_{r_s}^{r_o}\frac{dt}{dr}\,dr-2M\log\left(r_o e^{\frac{r_o}{2M}}\right),
\end{align}
```

The definition of this integral is possible in practice since $I_t$ takes the forms of the sums of Elliptic Integrals of the first $(F)$, second $(E)$ and third $(\Pi)$ kind, as well as logarithms and linear functions of $r$.

```math
\begin{align}
    I_t(r\rightarrow\infty)
        &\simeq \#_1 F(\phi(r)\mid k) + \#_2E(\phi(r)\mid k)+ \#_3\Pi(n;\phi(r)\mid k) + \#_4\log(\#_5 r^2) + r + \#_6
\end{align}
```
Some of the logarithmic divergence resides in the contribution from the Elliptic Integral of the third kind; a fact which follows from the logarithmic divergence of $\Pi(n;\phi\,|\,k)$ as $n\rightarrow \csc^2\phi$.
This divergence can be extracted with the connection formula,[^DLMF]
```math
\Pi\left(\alpha^{2};\phi\mid k\right)+\Pi\left(\omega^{2}; \phi\mid k\right)=F\left(%
\phi\mid k\right)+\sqrt{c}\,R_{C}\left((c-1)(c-k),(c-\alpha^{2})(c-\omega^{2})%
\right),
```
where,
```math
\alpha^2\omega^2=k\text{ and }c=\csc^2\phi
```
and
```math
    R_c(x,y)=\frac{1}{\sqrt{x-y}}\log\left\lvert\frac{\sqrt{x} + \sqrt{x-y}}{\sqrt{y}}\right\rvert.
```

>**Our convention here differs from [^DLMF] to be consistent with [^GL]. The two are related by,**
```math
\begin{align}
 F(\phi , \sqrt{k})
    &=F\left(\phi \mid k\right)\text{,}\\ 
E(\phi ,\sqrt{k})
    &=E\left(\phi \mid k\right)\text{, and}\\ 
\Pi(\phi, n ,\sqrt{k})
    &=\Pi\left(n;\phi \mid k\right).
\end{align}
```

Thus, our regularized time takes the form 
```math
\begin{align}
    \Delta \tilde I_t(r_s,\infty)
        &=\tilde I_t(\infty) - \tilde I_t(r_s),
\end{align}
```
with $\tilde I_t(r_s)=I_t(r_s)$ for finite $r_s$ and,
```math
\begin{align}
    \Delta \tilde I_t(\infty)
        &= \#_1 F(\phi\mid k) + \#_2E(\phi\mid k)+ \#_3\tilde\Pi(n;\phi\mid k) + \tilde\#_4\log(\tilde\#_5) + \#_6,
\end{align}
```
where we have used the connection formula to define a regularized Elliptic integral of the third kind by,
```math
\tilde \Pi(\alpha^2; \phi\mid k)
    =F(\phi\mid k) - \Pi(\omega^2, \phi\mid k).
```

## References
[^CAZ]: Alejandro Cárdenas-Avendaño, Alexandru Lupsasca, and Hengrui Zhu Phys. Rev. D 107, 043030 – Published 22 February 2023

[^DLMF]: NIST Digital Library of Mathematical Functions. https://dlmf.nist.gov/, Release 1.1.11 of 2023-09-15. F. W. J. Olver, A. B. Olde Daalhuis, D. W. Lozier, B. I. Schneider, R. F. Boisvert, C. W. Clark, B. R. Miller, B. V. Saunders, H. S. Cohl, and M. A. McClain, eds.

[^GL]: Samuel E. Gralla and Alexandru Lupsasca Phys. Rev. D 101, 044032