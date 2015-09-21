---
layout: page
title: Algorithms
---


#### General Relativistic Extended Magnetohydrodynamics (EMHD)
---
The model, derived by _Chandra et al._ (2015)[^EMHD_model_paper] is a one-fluid
model of a plasma consisting of electrons and ions. It considers the following
number current vector $$N^\mu$$ for the ions (set to be the same for electrons)
and total (electrons+ions) stress-tensor $$T^{\mu \nu}$$,<br>
$$\begin{align}
N^\mu & = n u^\mu \\
T^{\mu\nu} & =  (\rho + u + \frac{1}{2} b^2)u^\mu u^\nu 
              + (P + \frac{1}{2} b^2) h^{\mu\nu} 
              - b^\mu b^\nu + q^\mu u^\nu + q^\nu u^\mu 
              + \Pi^{\mu\nu}
\end{align}$$,<br>
where<br>
* $$n$$ is the number density of ions, which is equal to the number density of electrons.
* $$\rho = (m_i+m_e)n \approx m_i n$$ is the total rest mass density.<br>
* $$u$$ is the total internal energy.<br>
* $$P = (\gamma-1)u$$ is the total pressure approximated by a Gamma-law equation of state.<br>
* $$u^\mu$$ is the fluid four-velocity.<br>
* $$B^i$$ is the magnetic field three-vector.<br>
* $$q^\mu$$ is the heat flux four-vector.<br>
* $$\Pi^{\mu \nu}$$ is the shear stress.

The forms of $$q^\mu$$ and $$\Pi^{\mu \nu}$$ are obtained by assuming that the Larmor radius
of the particles is much lesser than the system scale,
![larmor_radius](../larmor_radius.png){: style="max-width: 400px; height: auto;"}
thus leading to a gyrotropic
distribution function $$f(t, x^i, p^i) = f(t, x^i, p_\parallel, p_\perp)$$.
Taking moments of this distribution functions, one gets
$$\begin{align}
q^\mu & = q \hat{b}^\mu \\
\Pi^{\mu \nu} & = - \Delta P \left(\hat{b}^\mu \hat{b}^\nu - \frac{1}{3}h^{\mu \nu} \right)
\end{align}$$,<br>
where $$\hat{b}^\mu$$ is a unit vector along the magnetic field lines and
$$h^{\mu \nu} = u^\mu u^\nu + g^{\mu \nu}$$ is a projection tensor. Thus heat
flows only along field lines and shear stress leads to a pressure anisotropy
$$\Delta P$$

The heat flux along $$q$$ and $$\Delta P$$ are visualized in momentum space in the figure below
![ideal_MHD_vs_GRMHD](../ideal_MHD_vs_GRMHD.png){: style="max-width: 400px; height: auto;"}

The evolution equations of the full system are
$$\begin{align}
\partial_\mu\left(\sqrt{-g}\rho u^\mu\right)  & = 0, \\
\partial_\mu\left(\sqrt{-g} T^\mu_\nu\right) & =
                                               \sqrt{-g}T^\kappa_\lambda
                                               \Gamma^\lambda_{\nu\kappa}, \\
\nabla_\mu (\tilde{q} u^\mu) &= - \frac{\tilde{q}-\tilde{q}_0}{\tau_R} + \frac{\tilde{q}}{2} \nabla_\mu u^\mu,\\
\nabla_\mu (\Delta \tilde{P} u^\mu) &= -\frac{\Delta \tilde{P} - \Delta \tilde{P}_0}{\tau_R} + 
\frac{\Delta\tilde{P}}{2} \nabla_\mu u^\mu.
\end{align}$$<br>
where<br>
$$\begin{align}
\tilde{q} &= q \left(\frac{\tau_R}{\chi \rho \Theta^2}\right)^{1/2},\\
\Delta \tilde{P} &= \Delta P \left(\frac{\tau_R}{\nu \rho\Theta}\right)^{1/2},\\
q_0 &= -\rho \chi \hat b^\mu (\nabla_\mu \Theta + \Theta u^\nu \nabla_\nu u_\mu),\\
\Delta P_0 &= 3\rho\nu (\hat{b}^\mu \hat{b}^\nu\nabla_\mu u_\nu - \frac{1}{3} \nabla_\mu u^\mu),
\end{align}$$.<br>
The parameters $$\tau_R$$ is a relaxation time scale over which $$\tilde{q},
\Delta\tilde{P}$$ relax to $$\tilde{q}_0$$ and $$\Delta\tilde{P}_0$$
respectively. The parameters $$\chi$$ and $$\nu$$ are the transport
coefficients, which along with $$\tau_R$$ are inputs to the model.

In $$\mathtt{grim}$$, the equations are evolved as
$$\begin{align}
\partial_t {\bf U} + \partial_i {\bf F}^i = {\bf S}
\end{align}$$,
where $${\bf U}=\sqrt{-g}(\rho
u^t,T^t_\mu,B^i,\tilde{q}u^t,\Delta\tilde{P}u^t)$$
is the vector of 10 conservative variables, $$g$$ is the determinant of the
4-metric $$g_{\mu\nu}$$, $${\bf F}^i$$ are the numerical fluxes and $${\bf S}$$
contains the source terms.

[^EMHD_model_paper]: Chandra, Gammie, Foucart, Quataert (2015), [arXiv:1508.00878](http://arxiv.org/abs/1508.00878)
