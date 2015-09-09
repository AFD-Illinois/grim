︠fbddea82-3103-42dd-bbf9-04c4a571eee5as︠
%auto
typeset_mode(True)
︡9100ba9e-5e05-492d-8113-c5431822b7b4︡{"auto":true}︡
︠cfc00840-33c2-4987-9665-f39a6125d569i︠
%html
<p style="text-align: center;"><img style="vertical-align: middle;" src="http://art.ngfiles.com/images/284000/284257_emutoons_derpemon-bulbasaur.png" alt="" width="320" height="180" /></p>
<h2 style="text-align: center;">Balbusaur</h2>

︡c46734ed-1981-4b41-8b91-ca14e9dd1ff4︡{"html":"<p style=\"text-align: center;\"><img style=\"vertical-align: middle;\" src=\"http://art.ngfiles.com/images/284000/284257_emutoons_derpemon-bulbasaur.png\" alt=\"\" width=\"320\" height=\"180\" /></p>\n<h2 style=\"text-align: center;\">Balbusaur</h2>\n\n"}︡
︠82740d77-8853-4cb1-a194-df8f1673c123i︠
%md

#### Variables involved:
1. Density $\rho$
2. Internal energy $u$
3. Velocity in $x^1$ direction $u^1$
4. Velocity in $x^2$ direction $u^2$
5. Magnetic field in $x^1$ direction $B^1$
6. Magnetic field in $x^2$ direction $B^2$
7. Heat flux magnitude $\phi$
8. Shear stress magnitude $\psi$

#### Problem descriptions:

1. CONDUCTION_1D:
   * Mean variables      : $\rho_0$, $u_0$, $B^1_0$
   * Perturbed variables : $\rho$, $u$, $u^1$, $\phi$
   * Wavenumbers         : $k_1$

2. VISCOSITY_1D:
   * Mean variables      : $\rho_0$, $u_0$, $B^1_0$
   * Perturbed variables : $\rho$, $u$, $u^1$, $\psi$
   * Wavenumbers         : $k_1$

3. CONDUCTION_2D:
   * Mean variables      : $\rho_0$, $u_0$, $B^1_0$, $B^2_0$
   * Perturbed variables : $\rho$, $u$, $u^1$, $u^2$, $\phi$
   * Wavenumbers         : $k_1$, $k_2$

4. VISCOSITY_2D:
   * Mean variables      : $\rho_0$, $u_0$, $B^1_0$, $B^2_0$
   * Perturbed variables : $\rho$, $u$, $u^1$, $u^2$, $B^1$, $B^2$, $\psi$
   * Wavenumbers         : $k_1$, $k_2$

5. CONDUCTION_AND_VISCOSITY_1D:
   * Mean variables      : $\rho_0$, $u_0$, $B^1_0$
   * Perturbed variables : $\rho$, $u$, $u^1$, $\phi$, $\psi$
   * Wavenumbers         : $k_1$

3. CONDUCTION_AND_VISCOSITY_2D:
   * Mean variables      : $\rho_0$, $u_0$, $B^1_0$, $B^2_0$
   * Perturbed variables : $\rho$, $u$, $u^1$, $u^2$, $\phi$, $\psi$
   * Wavenumbers         : $k_1$, $k_2$
︡7d01e076-a4fd-40cf-acbb-a223b84f62b9︡{"md":"\n#### Variables involved:\n1. Density $\\rho$\n2. Internal energy $u$\n3. Velocity in $x^1$ direction $u^1$\n4. Velocity in $x^2$ direction $u^2$\n5. Magnetic field in $x^1$ direction $B^1$\n6. Magnetic field in $x^2$ direction $B^2$\n7. Heat flux magnitude $\\phi$\n8. Shear stress magnitude $\\psi$\n\n#### Problem descriptions:\n\n1. CONDUCTION_1D: \n   * Mean variables      : $\\rho_0$, $u_0$, $B^1_0$\n   * Perturbed variables : $\\rho$, $u$, $u^1$, $\\phi$\n   * Wavenumbers         : $k_1$\n   \n2. VISCOSITY_1D: \n   * Mean variables      : $\\rho_0$, $u_0$, $B^1_0$\n   * Perturbed variables : $\\rho$, $u$, $u^1$, $\\psi$\n   * Wavenumbers         : $k_1$\n\n3. CONDUCTION_2D: \n   * Mean variables      : $\\rho_0$, $u_0$, $B^1_0$, $B^2_0$\n   * Perturbed variables : $\\rho$, $u$, $u^1$, $u^2$, $\\phi$\n   * Wavenumbers         : $k_1$, $k_2$\n   \n4. VISCOSITY_2D: \n   * Mean variables      : $\\rho_0$, $u_0$, $B^1_0$, $B^2_0$\n   * Perturbed variables : $\\rho$, $u$, $u^1$, $u^2$, $B^1$, $B^2$, $\\psi$\n   * Wavenumbers         : $k_1$, $k_2$\n\n5. CONDUCTION_AND_VISCOSITY_1D: \n   * Mean variables      : $\\rho_0$, $u_0$, $B^1_0$\n   * Perturbed variables : $\\rho$, $u$, $u^1$, $\\phi$, $\\psi$\n   * Wavenumbers         : $k_1$\n\n3. CONDUCTION_AND_VISCOSITY_2D: \n   * Mean variables      : $\\rho_0$, $u_0$, $B^1_0$, $B^2_0$\n   * Perturbed variables : $\\rho$, $u$, $u^1$, $u^2$, $\\phi$, $\\psi$\n   * Wavenumbers         : $k_1$, $k_2$\n"}︡
︠c953ed04-b6ab-4a98-a57f-2711f7d0d62b︠
# Inputs:

# Choose problem here:
problem = "IDEALMHD"

# Inputs for numerical diagonalization for finite k modes
rho0_num  = 1.
u0_num    = 2.
u1_num    = 0.
B10_num   = 1.
B20_num   = 0.0

Gamma_num = 4./3
P0_num    = (Gamma_num - 1.)*u0_num
T0_num    = P0_num/rho0_num
psi0_num  = 0.
k1_num    = 2.*pi
k2_num    = 0
kappa_num = 1.
eta_num   = 0.1
tau_num   = 1.01818181818182

  
︡77754425-39d4-4f31-bccd-1ef5bfd1a786︡
︠a7637a36-f643-4769-aed3-98f0e4e300b9︠
# Spatiotemporal variables
t, omega, k1, k2 = var('t, omega, k1, k2')

# Constants:
# Gamma : Adiabatic index
# kappa : Heat conductivity
# eta   : shear viscosity
# tau   : relaxation time scale
Gamma, kappa, eta, tau, fE = var('Gamma, kappa, eta, tau, fE')

# Background mean values
rho0, u0, u10, B10, B20, psi0, phi0 = var('rho0, u0, u10, B10, B20, psi0, phi0')

# Perturbations in space
delta_rho, delta_u, delta_u1, delta_u2, delta_u3, delta_B1, delta_B2, delta_B3, delta_phi, delta_psi = \
    var('delta_rho, delta_u, delta_u1, delta_u2, delta_u3, delta_B1, delta_B2, delta_B3, delta_phi, delta_psi')

# Perturbations in time
delta_rho_dt, delta_u_dt, delta_u1_dt, delta_u2_dt, delta_u3_dt, delta_B1_dt, delta_B2_dt, delta_B3_dt, delta_phi_dt, delta_psi_dt = \
    var('delta_rho_dt, delta_u_dt, delta_u1_dt, delta_u2_dt, delta_u3_dt, delta_B1_dt, delta_B2_dt, delta_B3_dt, delta_phi_dt, delta_psi_dt')

#eta = 0
if (problem=="IDEALMHD"):
    rho = rho0 + delta_rho
    u   = u0 + delta_u
    u1  = delta_u1
    u2  = delta_u2
    u3  = 0
    B1  = B10+delta_B1
    B2  = delta_B2
    B3  = 0
    phi = 0
    psi = 0
    
if (problem=="CONDUCTION_1D"):
    rho = rho0 + delta_rho
    u   = u0 + delta_u
    u1  = delta_u1
    u2  = 0
    u3  = 0
    B1  = B10
    B2  = 0
    B3  = 0
    phi = phi0 + delta_phi
    psi = 0

elif (problem=="CONDUCTION_2D"):
    rho = rho0 + delta_rho
    u   = u0 + delta_u
    u1  = delta_u1
    u2  = delta_u2
    u3  = 0
    B1  = B10
    B2  = B20
    B3  = 0
    phi = delta_phi
    psi = 0
elif (problem=="CONDUCTION_2D_TRANSVERSE"):
    rho = rho0 + delta_rho
    u   = u0 + delta_u
    u1  = delta_u1
    u2  = delta_u2
    u3  = 0
    B1  = B10+delta_B1
    B2  = 0+delta_B2
    B3  = 0
    phi = delta_phi
    psi = 0
elif (problem=="VISCOSITY_1D"):
    rho = rho0 + delta_rho
    u   = u0 + delta_u
    u1  = delta_u1
    u2  = 0
    u3  = 0
    B1  = B10
    B2  = 0
    B3  = 0
    phi = 0
    psi = delta_psi
elif (problem=="VISCOSITY_1D_MOVING"):
    rho = rho0 + delta_rho
    u   = u0
    u1  = u10 + delta_u1
    u2  = 0
    u3  = 0
    B1  = B10
    B2  = 0
    B3  = 0
    phi = 0
    psi = delta_psi
elif (problem=="CONDUCTION_1D_MOVING"):
    rho = rho0 + delta_rho
    u   = u0 + delta_u
    u1  = u10 + delta_u1
    u2  = 0
    u3  = 0
    B1  = B10
    B2  = 0
    B3  = 0
    phi = delta_phi
    psi = 0
elif (problem=="VISCOSITY_TRANS"):
    rho = rho0 
    u   = u0 
    u1  = u10
    u2  = delta_u2
    u3  = 0
    B1  = B10
    B2  = B20+delta_B2
    B3  = 0
    phi = 0
    psi = delta_psi
elif (problem=="VISCOSITY_2D"):
    rho = rho0 + delta_rho
    u   = u0 + delta_u
    u1  = delta_u1
    u2  = delta_u2
    u3  = 0
    B1  = B10  + delta_B1
    B2  = B20 + delta_B2
    B3  = 0
    phi = 0
    psi = delta_psi
elif (problem=="VISCOSITY_2D_TRANSVERSE"):
    rho = rho0 + delta_rho
    u   = u0 + delta_u
    u1  = delta_u1
    u2  = delta_u2
    u3  = 0
    B1  = B10  + delta_B1
    B2  = delta_B2
    B3  = 0
    phi = 0
    psi = delta_psi
elif (problem=="CONDUCTION_AND_VISCOSITY_1D"):
    rho = rho0 + delta_rho
    u   = u0 + delta_u
    u1  = delta_u1
    u2  = 0
    u3  = 0
    B1  = B10
    B2  = 0
    B3  = 0
    phi = delta_phi
    psi = delta_psi

elif (problem=="CONDUCTION_AND_VISCOSITY_2D"):
    rho = rho0 + delta_rho
    u   = u0 + delta_u
    u1  = delta_u1
    u2  = delta_u2
    u3  = 0
    B1  = B10
    B2  = B20
    B3  = 0
    phi = delta_phi
    psi = delta_psi
elif (problem=="PRINCIPAL_COMPONENTS"):
    rho = rho0 + delta_rho
    u   = u0 + delta_u
    u1  = delta_u1
    u2  = delta_u2
    u3  = 0
    B1  = B10+delta_B1
    B2  = delta_B2
    B3  = 0
    phi = delta_phi
    psi = delta_psi
elif (problem=="FIREHOSE"):
    rho = rho0 + delta_rho
    u   = u0 + delta_u
    u1  = 0
    u2  = delta_u2
    u3  = 0
    B1  = B10
    B2  = delta_B2
    B3  = 0
    phi = 0
    psi = psi0 + delta_psi


gcon = Matrix([ [-1, 0, 0, 0],
                [0, 1, 0, 0],
                [0, 0, 1, 0],
                [0, 0, 0, 1]
              ]
             )

gcov = gcon.inverse()

gamma = sqrt(1 +    gcov[1][1]*u1*u1 + gcov[2][2]*u2*u2 + gcov[3][3]*u3*u3
               + 2*(gcov[1][2]*u1*u2 + gcov[1][3]*u1*u3 + gcov[2][3]*u2*u3)
            )

ucon = [gamma, u1, u2, u3]
ucov = [-gamma, u1, u2, u3]

bcon0 = B1*ucov[1] + B2*ucov[2] + B3*ucov[3]
bcon1 = (B1 + bcon0*ucon[1])/ucon[0]
bcon2 = (B2 + bcon0*ucon[2])/ucon[0]
bcon3 = (B3 + bcon0*ucon[3])/ucon[0]

bcon = [bcon0, bcon1, bcon2, bcon3]
bcov = [-bcon0, bcon1, bcon2, bcon3]

bsqr = bcon[0]*bcov[0] + bcon[1]*bcov[1] + bcon[2]*bcov[2] + bcon[3]*bcov[3]

P = (Gamma - 1)*u
T = P/rho


def linearize(term):
    return taylor(term, (delta_rho, 0), \
                        (delta_u, 0),   \
                        (delta_u1, 0),  \
                        (delta_u2, 0),  \
                        (delta_u3, 0),  \
                        (delta_B1, 0),  \
                        (delta_B2, 0),  \
                        (delta_B3, 0),  \
                        (delta_phi, 0), \
                        (delta_psi, 0), \
                        (delta_rho_dt, 0), \
                        (delta_u_dt, 0),   \
                        (delta_u1_dt, 0),  \
                        (delta_u2_dt, 0),  \
                        (delta_u3_dt, 0),  \
                        (delta_B1_dt, 0),  \
                        (delta_B2_dt, 0),  \
                        (delta_B3_dt, 0),  \
                        (delta_phi_dt, 0), \
                        (delta_psi_dt, 0), 1 \
                 ).simplify_full()

def d_dX1(term):
    term  = Expression(SR, linearize(term))

    expr =   term.coefficient(delta_rho) * I * k1 * delta_rho \
           + term.coefficient(delta_u)   * I * k1 * delta_u   \
           + term.coefficient(delta_u1)  * I * k1 * delta_u1  \
           + term.coefficient(delta_u2)  * I * k1 * delta_u2  \
           + term.coefficient(delta_u3)  * I * k1 * delta_u3  \
           + term.coefficient(delta_B1)  * I * k1 * delta_B1  \
           + term.coefficient(delta_B2)  * I * k1 * delta_B2  \
           + term.coefficient(delta_B3)  * I * k1 * delta_B3  \
           + term.coefficient(delta_phi) * I * k1 * delta_phi \
           + term.coefficient(delta_psi) * I * k1 * delta_psi

    return expr

def d_dX2(term):
    term  = Expression(SR, linearize(term))

    expr =   term.coefficient(delta_rho) * I * k2 * delta_rho \
           + term.coefficient(delta_u)   * I * k2 * delta_u   \
           + term.coefficient(delta_u1)  * I * k2 * delta_u1  \
           + term.coefficient(delta_u2)  * I * k2 * delta_u2  \
           + term.coefficient(delta_u3)  * I * k2 * delta_u3  \
           + term.coefficient(delta_B1)  * I * k2 * delta_B1  \
           + term.coefficient(delta_B2)  * I * k2 * delta_B2  \
           + term.coefficient(delta_B3)  * I * k2 * delta_B3  \
           + term.coefficient(delta_phi) * I * k2 * delta_phi \
           + term.coefficient(delta_psi) * I * k2 * delta_psi

    return expr


def d_dt(term):
    term  = Expression(SR, linearize(term))

    expr =   term.coefficient(delta_rho) * delta_rho_dt \
           + term.coefficient(delta_u)   * delta_u_dt   \
           + term.coefficient(delta_u1)  * delta_u1_dt  \
           + term.coefficient(delta_u2)  * delta_u2_dt  \
           + term.coefficient(delta_u3)  * delta_u3_dt  \
           + term.coefficient(delta_B1)  * delta_B1_dt  \
           + term.coefficient(delta_B2)  * delta_B2_dt  \
           + term.coefficient(delta_B3)  * delta_B3_dt  \
           + term.coefficient(delta_phi) * delta_phi_dt \
           + term.coefficient(delta_psi) * delta_psi_dt

    return expr


def delta(mu, nu):
    if (mu==nu):
        return 1
    else:
        return 0

FakeEMHD = 0

def TUpDown(mu, nu):

    return  (rho + u + P + bsqr + FakeEMHD*bsqr/6)*ucon[mu]*ucov[nu] + (P + bsqr/2 + FakeEMHD*bsqr/6)*delta(mu, nu) - (1+0.5*FakeEMHD)*bcon[mu]*bcov[nu] \
          + phi/sqrt(bsqr)*(bcon[mu]*ucov[nu] + ucon[mu]*bcov[nu]) + psi/bsqr*(bcon[mu]*bcov[nu]) - psi/3*(ucon[mu]*ucov[nu] + delta(mu, nu))

def acon(mu):
    return linearize(ucon[0]*d_dt(ucon[mu]) + ucon[1]*d_dX1(ucon[mu]) + ucon[2]*d_dX2(ucon[mu]))

def qconEckart(mu):
    acov = [-acon(0), acon(1), acon(2), acon(3)]

    ans = -kappa*(ucon[mu]*ucon[0] + gcon[mu, 0])*(d_dt(T) +  T*acov[0]) \
          -kappa*(ucon[mu]*ucon[1] + gcon[mu, 1])*(d_dX1(T) + T*acov[1]) \
          -kappa*(ucon[mu]*ucon[2] + gcon[mu, 2])*(d_dX2(T) + T*acov[2])

    return linearize(ans)

Eqn_rho = linearize(d_dt(rho*ucon[0])   + d_dX1(rho*ucon[1])   + d_dX2(rho*ucon[2]))
Eqn_u   = linearize(d_dt(TUpDown(0, 0)) + d_dX1(TUpDown(1, 0)) + d_dX2(TUpDown(2, 0)))
Eqn_u1  = linearize(d_dt(TUpDown(0, 1)) + d_dX1(TUpDown(1, 1)) + d_dX2(TUpDown(2, 1)))
Eqn_u2  = linearize(d_dt(TUpDown(0, 2)) + d_dX1(TUpDown(1, 2)) + d_dX2(TUpDown(2, 2)))

Eqn_B1  = linearize(d_dt(B1) + d_dX2(bcon[1]*ucon[2] - bcon[2]*ucon[1]) )
Eqn_B2  = linearize(d_dt(B2) + d_dX1(bcon[2]*ucon[1] - bcon[1]*ucon[2]) )

if (problem=='VISCOSITY_1D_MOVING'):
    beta1 = 1
else:
    beta1 = tau/(kappa*T)

beta2 = tau/(2*eta)
phi_relaxed = (bcov[0]*qconEckart(0) + bcov[1]*qconEckart(1) + bcov[2]*qconEckart(2) + bcov[3]*qconEckart(3) )/sqrt(bsqr)
psi_relaxed = 0
for nu in xrange(4):
    psi_relaxed = psi_relaxed - 3*eta/bsqr * (bcon[nu]* (bcon[0]*d_dt(ucov[nu]) + bcon[1]*d_dX1(ucov[nu]) + bcon[2]*d_dX2(ucov[nu])) )

psi_relaxed = psi_relaxed + eta*(d_dt(ucon[0]) + d_dX1(ucon[1]) + d_dX2(ucon[2]) )

princcoeff = 1
if (problem=='PRINCIPAL_COMPONENTS'):
    princcoeff = 0

Eqn_phi = linearize(   ucon[0]*d_dt(phi) + ucon[1]*d_dX1(phi) + ucon[2]*d_dX2(phi) + (princcoeff*phi - phi_relaxed)/tau)
                    #+ (phi*T/(2*beta1))*(d_dt(beta1*ucon[0]/T) + d_dX1(beta1*ucon[1]/T) + d_dX2(beta1*ucon[2]/T))
                   
Eqn_psi = linearize(   ucon[0]*d_dt(psi) + ucon[1]*d_dX1(psi) + ucon[2]*d_dX2(psi) + (princcoeff*psi - psi_relaxed)/tau)
                    #+ (psi*T/(2*beta2))*(d_dt(beta2*ucon[0]/T) + d_dX1(beta2*ucon[1]/T) + d_dX2(beta2*ucon[2]/T))

if (problem=='IDEALMHD'):

    Eqns          = [Eqn_rho==0, Eqn_u==0, Eqn_u1==0, Eqn_u2==0, Eqn_B1==0, Eqn_B2==0]
    delta_vars    = [delta_rho, delta_u, delta_u1, delta_u2, delta_B1, delta_B2]
    delta_vars_dt = [delta_rho_dt, delta_u_dt, delta_u1_dt, delta_u2_dt, delta_B1_dt, delta_B2_dt]

if (problem=='CONDUCTION_1D'):

    Eqns          = [Eqn_rho==0, Eqn_u==0, Eqn_u1==0, Eqn_phi==0]
    delta_vars    = [delta_rho, delta_u, delta_u1, delta_phi]
    delta_vars_dt = [delta_rho_dt, delta_u_dt, delta_u1_dt, delta_phi_dt]

elif (problem=='CONDUCTION_2D'):

    Eqns          = [Eqn_rho==0, Eqn_u==0, Eqn_u1==0, Eqn_u2==0, Eqn_phi==0]
    delta_vars    = [delta_rho, delta_u, delta_u1, delta_u2, delta_phi]
    delta_vars_dt = [delta_rho_dt, delta_u_dt, delta_u1_dt, delta_u2_dt, delta_phi_dt]
elif (problem=='CONDUCTION_2D_TRANSVERSE'):

    Eqns          = [Eqn_rho==0, Eqn_u==0, Eqn_u1==0, Eqn_u2==0, Eqn_B1==0, Eqn_B2==0, Eqn_phi==0]
    delta_vars    = [delta_rho, delta_u, delta_u1, delta_u2, delta_B1, delta_B2, delta_phi]
    delta_vars_dt = [delta_rho_dt, delta_u_dt, delta_u1_dt, delta_u2_dt, delta_B1_dt, delta_B2_dt, delta_phi_dt]
elif (problem=='VISCOSITY_1D'):

    Eqns          = [Eqn_rho==0, Eqn_u==0, Eqn_u1==0, Eqn_psi==0]
    delta_vars    = [delta_rho, delta_u, delta_u1, delta_psi]
    delta_vars_dt = [delta_rho_dt, delta_u_dt, delta_u1_dt, delta_psi_dt]
elif (problem=='VISCOSITY_1D_MOVING'):
    
    Eqns          = [Eqn_rho==0, Eqn_u1==0, Eqn_psi==0]
    delta_vars    = [delta_rho, delta_u1, delta_psi]
    delta_vars_dt = [delta_rho_dt, delta_u1_dt, delta_psi_dt]
elif (problem=='CONDUCTION_1D_MOVING'):

    Eqns          = [Eqn_rho==0,Eqn_u==0, Eqn_u1==0, Eqn_phi==0]
    delta_vars    = [delta_rho,delta_u, delta_u1, delta_phi]
    delta_vars_dt = [delta_rho_dt,delta_u_dt, delta_u1_dt, delta_phi_dt]
elif (problem=='VISCOSITY_TRANS'):

    Eqns          = [Eqn_u2==0, Eqn_B2==0, Eqn_psi==0]
    delta_vars    = [delta_u2, delta_B2, delta_psi]
    delta_vars_dt = [delta_u2_dt, delta_B2_dt, delta_psi_dt]
elif (problem=='VISCOSITY_2D'):

    Eqns          = [Eqn_rho==0, Eqn_u==0, Eqn_u1==0, Eqn_u2==0, Eqn_B1==0, Eqn_B2==0, Eqn_psi==0]
    delta_vars    = [delta_rho, delta_u, delta_u1, delta_u2, delta_B1, delta_B2, delta_psi]
    delta_vars_dt = [delta_rho_dt, delta_u_dt, delta_u1_dt, delta_u2_dt, delta_B1_dt, delta_B2_dt, delta_psi_dt]
elif (problem=='VISCOSITY_2D_TRANSVERSE'):

    Eqns          = [Eqn_rho==0, Eqn_u==0, Eqn_u1==0, Eqn_u2==0, Eqn_B1==0, Eqn_B2==0, Eqn_psi==0]
    delta_vars    = [delta_rho, delta_u, delta_u1, delta_u2, delta_B1, delta_B2, delta_psi]
    delta_vars_dt = [delta_rho_dt, delta_u_dt, delta_u1_dt, delta_u2_dt, delta_B1_dt, delta_B2_dt, delta_psi_dt]
elif (problem=='CONDUCTION_AND_VISCOSITY_1D'):

    Eqns          = [Eqn_rho==0, Eqn_u==0, Eqn_u1==0, Eqn_phi==0, Eqn_psi==0]
    delta_vars    = [delta_rho, delta_u, delta_u1, delta_phi, delta_psi]
    delta_vars_dt = [delta_rho_dt, delta_u_dt, delta_u1_dt, delta_phi_dt, delta_psi_dt]

elif (problem=='CONDUCTION_AND_VISCOSITY_2D'):

    Eqns          = [Eqn_rho==0, Eqn_u==0, Eqn_u1==0, Eqn_u2==0, Eqn_phi==0, Eqn_psi==0]
    delta_vars    = [delta_rho, delta_u, delta_u1, delta_u2, delta_phi, delta_psi]
    delta_vars_dt = [delta_rho_dt, delta_u_dt, delta_u1_dt, delta_u2_dt, delta_phi_dt, delta_psi_dt]

elif (problem=='PRINCIPAL_COMPONENTS'):

    Eqns          = [Eqn_rho==0, Eqn_u==0, Eqn_u1==0, Eqn_u2==0, Eqn_B1==0, Eqn_B2==0, Eqn_phi==0, Eqn_psi==0]
    delta_vars    = [delta_rho, delta_u, delta_u1, delta_u2, delta_B1, delta_B2, delta_phi, delta_psi]
    delta_vars_dt = [delta_rho_dt, delta_u_dt, delta_u1_dt, delta_u2_dt, delta_B1_dt, delta_B2_dt, delta_phi_dt, delta_psi_dt]
elif (problem=='FIREHOSE'):

    Eqns          = [Eqn_rho==0, Eqn_u==0, Eqn_u2==0, Eqn_B2==0, Eqn_psi==0]
    delta_vars    = [delta_rho, delta_u, delta_u2, delta_B2, delta_psi]
    delta_vars_dt = [delta_rho_dt, delta_u_dt, delta_u2_dt, delta_B2_dt, delta_psi_dt]

solutions = solve(Eqns, delta_vars_dt, solution_dict=True)

solns_delta_vars_dt = []
for dvar_dt in delta_vars_dt:
    solns_delta_vars_dt.append(solutions[0][dvar_dt])

M = jacobian(solns_delta_vars_dt, delta_vars)
M = M.apply_map(lambda x : x.simplify_full())

pretty_print("Linearized system : ", )
print("\n")
pretty_print(Matrix(delta_vars_dt).transpose(), " = ", M, Matrix(delta_vars).transpose())
print("\n\n")
pretty_print("Eigenvalues and eigenvectors in the $k \\rightarrow 0$ limit : ", )
M.subs(k1=0, k2=0).eigenvectors_right()


# Numerical diagonalization:

M_numerical = M.subs(rho0=rho0_num, u0=u0_num, u10=u1_num, B10=B10_num, B20=B20_num, psi0=psi0_num, Gamma=Gamma_num, kappa=kappa_num, eta=eta_num, tau=tau_num, k1=k1_num, k2=k2_num)
M_numerical = M_numerical.change_ring(CDF)
eigenvecs   = M_numerical.eigenvectors_right()

print "Numerical eigenvalues and eigenvectors for k > 0:\n"
print "--------------------------\n"

if (problem=='IDEALMHD'):

    print "--------------------------"
    for i in xrange(len(eigenvecs)):
        print "Eigenvalue = ",     eigenvecs[i][0]
        print "delta_rho  = ",     eigenvecs[i][1][0][0]
        print "delta_u    = ",     eigenvecs[i][1][0][1]
        print "delta_u1   = ",     eigenvecs[i][1][0][2]
        print "delta_u2  = ",     eigenvecs[i][1][0][3]
        print "delta_B1   = ",     eigenvecs[i][1][0][4]
        print "delta_B2  = ",     eigenvecs[i][1][0][5]
        print "--------------------------"

                    
if (problem=='CONDUCTION_1D'):

    print "kappa = ", kappa_num
    print "tau   = ", tau_num
    print "--------------------------"
    for i in xrange(len(eigenvecs)):
        print "Eigenvalue = ",     eigenvecs[i][0]
        print "delta_rho  = ",     eigenvecs[i][1][0][0]
        print "delta_u    = ",     eigenvecs[i][1][0][1]
        print "delta_u1   = ",     eigenvecs[i][1][0][2]
        print "delta_phi  = ",     eigenvecs[i][1][0][3]
        print "--------------------------"

if (problem=='CONDUCTION_2D'):

    print "kappa = ", kappa_num
    print "tau   = ", tau_num
    print "--------------------------"
    for i in xrange(len(eigenvecs)):
        print "Eigenvalue = ",     eigenvecs[i][0]
        print "delta_rho  = ",     eigenvecs[i][1][0][0]
        print "delta_u    = ",     eigenvecs[i][1][0][1]
        print "delta_u1   = ",     eigenvecs[i][1][0][2]
        print "delta_u2   = ",     eigenvecs[i][1][0][3]
        print "delta_phi  = ",     eigenvecs[i][1][0][4]
        print "--------------------------"
if (problem=='CONDUCTION_2D_TRANSVERSE'):

    print "eta   = ", eta_num
    print "tau   = ", tau_num
    print "--------------------------"
    for i in xrange(len(eigenvecs)):
        print "Eigenvalue = ",     eigenvecs[i][0]
        print "delta_rho  = ",     eigenvecs[i][1][0][0]
        print "delta_u    = ",     eigenvecs[i][1][0][1]
        print "delta_u1   = ",     eigenvecs[i][1][0][2]
        print "delta_u2   = ",     eigenvecs[i][1][0][3]
        print "delta_B1   = ",     eigenvecs[i][1][0][4]
        print "delta_B2   = ",     eigenvecs[i][1][0][5]
        print "delta_phi  = ",     eigenvecs[i][1][0][6]
        print "--------------------------"
if (problem=='VISCOSITY_1D'):

    print "eta   = ", eta_num
    print "tau   = ", tau_num
    print "--------------------------"
    for i in xrange(len(eigenvecs)):
        print "Eigenvalue = ",     eigenvecs[i][0]
        print "delta_rho  = ",     eigenvecs[i][1][0][0]
        print "delta_u    = ",     eigenvecs[i][1][0][1]
        print "delta_u1   = ",     eigenvecs[i][1][0][2]
        print "delta_psi  = ",     eigenvecs[i][1][0][3]
        print "--------------------------"
if (problem=='VISCOSITY_1D_MOVING'):

    print "eta   = ", eta_num
    print "tau   = ", tau_num
    print "--------------------------"
    for i in xrange(len(eigenvecs)):
        print "Eigenvalue = ",     eigenvecs[i][0]
        print "delta_rho  = ",     eigenvecs[i][1][0][0]
        print "delta_u1   = ",     eigenvecs[i][1][0][1]
        print "delta_psi  = ",     eigenvecs[i][1][0][2]
        print "--------------------------"
if (problem=='CONDUCTION_1D_MOVING'):

    print "kappa   = ", kappa_num
    print "tau   = ", tau_num
    print "--------------------------"
    for i in xrange(len(eigenvecs)):
        print "Eigenvalue = ",     eigenvecs[i][0]
        print "delta_rho  = ",     eigenvecs[i][1][0][0]
        print "delta_u    = ",     eigenvecs[i][1][0][1]
        print "delta_u1   = ",     eigenvecs[i][1][0][2]
        print "delta_phi  = ",     eigenvecs[i][1][0][3]
        print "--------------------------"
if (problem=='VISCOSITY_TRANS'):

    print "eta   = ", eta_num
    print "tau   = ", tau_num
    print "--------------------------"
    for i in xrange(len(eigenvecs)):
        print "Eigenvalue = ",     eigenvecs[i][0]
        print "delta_u2   = ",     eigenvecs[i][1][0][0]
        print "delta_B2    = ",    eigenvecs[i][1][0][1]
        print "delta_psi  = ",     eigenvecs[i][1][0][2]
        print "--------------------------"
if (problem=='VISCOSITY_2D'):

    print "eta   = ", eta_num
    print "tau   = ", tau_num
    print "--------------------------"
    for i in xrange(len(eigenvecs)):
        print "Eigenvalue = ",     eigenvecs[i][0]
        print "delta_rho  = ",     eigenvecs[i][1][0][0]
        print "delta_u    = ",     eigenvecs[i][1][0][1]
        print "delta_u1   = ",     eigenvecs[i][1][0][2]
        print "delta_u2   = ",     eigenvecs[i][1][0][3]
        print "delta_B1   = ",     eigenvecs[i][1][0][4]
        print "delta_B2   = ",     eigenvecs[i][1][0][5]
        print "delta_psi  = ",     eigenvecs[i][1][0][6]
        print "--------------------------"
if (problem=='VISCOSITY_2D_TRANSVERSE'):

    print "eta   = ", eta_num
    print "tau   = ", tau_num
    print "--------------------------"
    for i in xrange(len(eigenvecs)):
        print "Eigenvalue = ",     eigenvecs[i][0]
        print "delta_rho  = ",     eigenvecs[i][1][0][0]
        print "delta_u    = ",     eigenvecs[i][1][0][1]
        print "delta_u1   = ",     eigenvecs[i][1][0][2]
        print "delta_u2   = ",     eigenvecs[i][1][0][3]
        print "delta_B1   = ",     eigenvecs[i][1][0][4]
        print "delta_B2   = ",     eigenvecs[i][1][0][5]
        print "delta_psi  = ",     eigenvecs[i][1][0][6]
        print "--------------------------"
if (problem=='CONDUCTION_AND_VISCOSITY_1D'):

    print "kappa = ", kappa_num
    print "eta   = ", eta_num
    print "tau   = ", tau_num
    print "--------------------------"
    for i in xrange(len(eigenvecs)):
        print "Eigenvalue = ",     eigenvecs[i][0]
        print "delta_rho  = ",     eigenvecs[i][1][0][0]
        print "delta_u    = ",     eigenvecs[i][1][0][1]
        print "delta_u1   = ",     eigenvecs[i][1][0][2]
        print "delta_phi  = ",     eigenvecs[i][1][0][3]
        print "delta_psi  = ",     eigenvecs[i][1][0][4]
        print "--------------------------"

if (problem=='CONDUCTION_AND_VISCOSITY_2D'):

    print "kappa = ", kappa_num
    print "eta   = ", eta_num
    print "tau   = ", tau_num
    print "--------------------------"
    for i in xrange(len(eigenvecs)):
        print "Eigenvalue = ",     eigenvecs[i][0]
        print "delta_rho  = ",     eigenvecs[i][1][0][0]
        print "delta_u    = ",     eigenvecs[i][1][0][1]
        print "delta_u1   = ",     eigenvecs[i][1][0][2]
        print "delta_u2   = ",     eigenvecs[i][1][0][3]
        print "delta_phi  = ",     eigenvecs[i][1][0][4]
        print "delta_psi  = ",     eigenvecs[i][1][0][5]
        print "--------------------------"
if (problem=='PRINCIPAL_COMPONENTS'):

    print "kappa = ", kappa_num
    print "eta   = ", eta_num
    print "tau   = ", tau_num
    print "--------------------------"
    for i in xrange(len(eigenvecs)):
        print "Eigenvalue = ",     eigenvecs[i][0]
        print "delta_rho  = ",     eigenvecs[i][1][0][0]
        print "delta_u    = ",     eigenvecs[i][1][0][1]
        print "delta_u1   = ",     eigenvecs[i][1][0][2]
        print "delta_u2   = ",     eigenvecs[i][1][0][3]
        print "delta_b1   = ",     eigenvecs[i][1][0][4]
        print "delta_b2   = ",     eigenvecs[i][1][0][5]
        print "delta_phi  = ",     eigenvecs[i][1][0][6]
        print "delta_psi  = ",     eigenvecs[i][1][0][7]
        print "--------------------------"

if (problem=='FIREHOSE'):

    print "eta   = ", eta_num
    print "tau   = ", tau_num
    print "--------------------------"
    for i in xrange(len(eigenvecs)):
        print "Eigenvalue = ",     eigenvecs[i][0]
        print "delta_rho  = ",     eigenvecs[i][1][0][0]
        print "delta_u    = ",     eigenvecs[i][1][0][1]
        print "delta_u1   = ",     eigenvecs[i][1][0][2]
        print "delta_B1   = ",     eigenvecs[i][1][0][3]
        #print "delta_psi  = ",     eigenvecs[i][1][0][4]
        print "--------------------------"




︡dfbf1812-f200-4936-a88a-7c1117ae8553︡{"html":"<div align='center'>Linearized system : </div>"}︡{"stdout":"\n\n"}︡{"html":"<div align='center'>$\\displaystyle \\left(\\begin{array}{r}\n\\delta_{\\rho_{\\mathit{dt}}} \\\\\n\\delta_{u_{\\mathit{dt}}} \\\\\n\\delta_{\\mathit{u1}_{\\mathit{dt}}} \\\\\n\\delta_{\\mathit{u2}_{\\mathit{dt}}} \\\\\n\\delta_{\\mathit{B1}_{\\mathit{dt}}} \\\\\n\\delta_{\\mathit{B2}_{\\mathit{dt}}}\n\\end{array}\\right)$  =  $\\displaystyle \\left(\\begin{array}{rrrrrr}\n0 &amp; 0 &amp; -i \\, k_{1} \\rho_{0} &amp; -i \\, k_{2} \\rho_{0} &amp; 0 &amp; 0 \\\\\n0 &amp; 0 &amp; -i \\, \\Gamma k_{1} u_{0} &amp; -i \\, \\Gamma k_{2} u_{0} &amp; 0 &amp; 0 \\\\\n0 &amp; \\frac{{\\left(-i \\, \\Gamma + i\\right)} k_{1}}{\\Gamma u_{0} + \\rho_{0}} &amp; 0 &amp; 0 &amp; \\frac{i \\, B_{10} k_{1}}{\\Gamma u_{0} + \\rho_{0}} &amp; \\frac{i \\, B_{10} k_{2}}{\\Gamma u_{0} + \\rho_{0}} \\\\\n0 &amp; \\frac{{\\left(-i \\, \\Gamma + i\\right)} k_{2}}{B_{10}^{2} + \\Gamma u_{0} + \\rho_{0}} &amp; 0 &amp; 0 &amp; -\\frac{i \\, B_{10} k_{2}}{B_{10}^{2} + \\Gamma u_{0} + \\rho_{0}} &amp; \\frac{i \\, B_{10} k_{1}}{B_{10}^{2} + \\Gamma u_{0} + \\rho_{0}} \\\\\n0 &amp; 0 &amp; 0 &amp; -i \\, B_{10} k_{2} &amp; 0 &amp; 0 \\\\\n0 &amp; 0 &amp; 0 &amp; i \\, B_{10} k_{1} &amp; 0 &amp; 0\n\\end{array}\\right)$ $\\displaystyle \\left(\\begin{array}{r}\n\\delta_{\\rho} \\\\\n\\delta_{u} \\\\\n\\delta_{u_{1}} \\\\\n\\delta_{u_{2}} \\\\\n\\delta_{B_{1}} \\\\\n\\delta_{B_{2}}\n\\end{array}\\right)$</div>"}︡{"stdout":"\n\n\n"}︡{"html":"<div align='center'>Eigenvalues and eigenvectors in the $k \\rightarrow 0$ limit : </div>"}︡{"tex":{"tex":"\\left[\\left(0, \\left[\\left(1,\\,0,\\,0,\\,0,\\,0,\\,0\\right), \\left(0,\\,1,\\,0,\\,0,\\,0,\\,0\\right), \\left(0,\\,0,\\,1,\\,0,\\,0,\\,0\\right), \\left(0,\\,0,\\,0,\\,1,\\,0,\\,0\\right), \\left(0,\\,0,\\,0,\\,0,\\,1,\\,0\\right), \\left(0,\\,0,\\,0,\\,0,\\,0,\\,1\\right)\\right], 6\\right)\\right]","display":true}}︡{"stdout":"Numerical eigenvalues and eigenvectors for k > 0:\n\n"}︡{"stdout":"--------------------------\n\n"}︡{"stdout":"--------------------------\nEigenvalue =  0.0\ndelta_rho  =  1.0\ndelta_u    =  0.0\ndelta_u1   =  0.0\ndelta_u2  =  0.0\ndelta_B1   =  0.0\ndelta_B2  =  0.0\n--------------------------\nEigenvalue =  3.09362659024*I\ndelta_rho  =  0.345991032308\ndelta_u    =  0.922642752822\ndelta_u1   =  -0.170354208129\ndelta_u2  =  0.0\ndelta_B1   =  0.0\ndelta_B2  =  0.0\n--------------------------\nEigenvalue =  -3.09362659024*I\ndelta_rho  =  0.345991032308\ndelta_u    =  0.922642752822\ndelta_u1   =  0.170354208129 - 1.38777878078e-17*I\ndelta_u2  =  0.0\ndelta_B1   =  0.0\ndelta_B2  =  0.0\n--------------------------\nEigenvalue =  2.90854962399*I\ndelta_rho  =  0.0\ndelta_u    =  0.0\ndelta_u1   =  0.0\ndelta_u2  =  0.420084025208\ndelta_B1   =  0.0\ndelta_B2  =  0.907485212973\n--------------------------\nEigenvalue =  -1.11022302463e-16 - 2.90854962399*I\ndelta_rho  =  0.0\ndelta_u    =  0.0\ndelta_u1   =  0.0\ndelta_u2  =  -0.420084025208 + 2.77555756156e-17*I\ndelta_B1   =  0.0\ndelta_B2  =  0.907485212973\n--------------------------\nEigenvalue =  0.0\ndelta_rho  =  1.0\ndelta_u    =  -4.09547243262e-275 - 3.06162672181e-291*I\ndelta_u1   =  -7.57814631902e-292 - 2.33379417107e-292*I\ndelta_u2  =  0.0\ndelta_B1   =  -1.36515747754e-275\ndelta_B2  =  0.0\n--------------------------\n"}︡
︠d39e6e67-8c79-4c5a-8606-03aa2a4281dc︠
Mk1  = (M -I*omega*identity_matrix(6))
dispersion = Mk1.subs(k2=0).determinant().numerator().simplify_full()
soln = solve(dispersion==0, omega, solution_dict=True)
soln
︡10032e23-36d0-4a20-9c10-e5dd722b11c3︡{"tex":{"tex":"\\left[\\left\\{\\omega : -\\sqrt{\\frac{B_{10}^{2} \\Gamma \\mathit{fE}}{B_{10}^{2} \\mathit{fE} - 3 \\, \\Gamma u_{0} - 3 \\, \\rho_{0}} - \\frac{B_{10}^{2} \\mathit{fE}}{B_{10}^{2} \\mathit{fE} - 3 \\, \\Gamma u_{0} - 3 \\, \\rho_{0}} - \\frac{3 \\, \\Gamma^{2} u_{0}}{B_{10}^{2} \\mathit{fE} - 3 \\, \\Gamma u_{0} - 3 \\, \\rho_{0}} + \\frac{3 \\, \\Gamma u_{0}}{B_{10}^{2} \\mathit{fE} - 3 \\, \\Gamma u_{0} - 3 \\, \\rho_{0}}} k_{1}\\right\\}, \\left\\{\\omega : \\sqrt{\\frac{B_{10}^{2} \\Gamma \\mathit{fE}}{B_{10}^{2} \\mathit{fE} - 3 \\, \\Gamma u_{0} - 3 \\, \\rho_{0}} - \\frac{B_{10}^{2} \\mathit{fE}}{B_{10}^{2} \\mathit{fE} - 3 \\, \\Gamma u_{0} - 3 \\, \\rho_{0}} - \\frac{3 \\, \\Gamma^{2} u_{0}}{B_{10}^{2} \\mathit{fE} - 3 \\, \\Gamma u_{0} - 3 \\, \\rho_{0}} + \\frac{3 \\, \\Gamma u_{0}}{B_{10}^{2} \\mathit{fE} - 3 \\, \\Gamma u_{0} - 3 \\, \\rho_{0}}} k_{1}\\right\\}, \\left\\{\\omega : -B_{10} k_{1} \\sqrt{\\frac{3 \\, \\mathit{fE}}{B_{10}^{2} \\mathit{fE} + 6 \\, B_{10}^{2} + 6 \\, \\Gamma u_{0} + 6 \\, \\rho_{0}} + \\frac{6}{B_{10}^{2} \\mathit{fE} + 6 \\, B_{10}^{2} + 6 \\, \\Gamma u_{0} + 6 \\, \\rho_{0}}}\\right\\}, \\left\\{\\omega : B_{10} k_{1} \\sqrt{\\frac{3 \\, \\mathit{fE}}{B_{10}^{2} \\mathit{fE} + 6 \\, B_{10}^{2} + 6 \\, \\Gamma u_{0} + 6 \\, \\rho_{0}} + \\frac{6}{B_{10}^{2} \\mathit{fE} + 6 \\, B_{10}^{2} + 6 \\, \\Gamma u_{0} + 6 \\, \\rho_{0}}}\\right\\}, \\left\\{\\omega : 0\\right\\}\\right]","display":true}}︡
︠1cd9e8c5-212b-4073-a7ff-f63cfac7c73b︠
Mk1  = (M -I*omega*identity_matrix(5))
dispersion = Mk1.subs(B10=0, k2=0).determinant().numerator().simplify_full()
soln = solve(dispersion==0, omega, solution_dict=True)
soln
#soln[2][omega]
︡c2237abd-b5ff-4019-a8c4-ade5950b55cf︡{"tex":{"tex":"\\left[\\left\\{\\omega : -\\sqrt{\\frac{B_{20}^{2} \\Gamma}{B_{20}^{2} + \\Gamma u_{0} + \\rho_{0}} + \\frac{\\Gamma^{2} u_{0}}{B_{20}^{2} + \\Gamma u_{0} + \\rho_{0}} - \\frac{B_{20}^{2}}{B_{20}^{2} + \\Gamma u_{0} + \\rho_{0}} - \\frac{\\Gamma u_{0}}{B_{20}^{2} + \\Gamma u_{0} + \\rho_{0}}} k_{1}\\right\\}, \\left\\{\\omega : \\sqrt{\\frac{B_{20}^{2} \\Gamma}{B_{20}^{2} + \\Gamma u_{0} + \\rho_{0}} + \\frac{\\Gamma^{2} u_{0}}{B_{20}^{2} + \\Gamma u_{0} + \\rho_{0}} - \\frac{B_{20}^{2}}{B_{20}^{2} + \\Gamma u_{0} + \\rho_{0}} - \\frac{\\Gamma u_{0}}{B_{20}^{2} + \\Gamma u_{0} + \\rho_{0}}} k_{1}\\right\\}, \\left\\{\\omega : -\\frac{\\Gamma \\rho_{0} u_{0} + \\rho_{0}^{2}}{i \\, \\rho_{0}^{2} \\tau + {\\left(i \\, \\Gamma \\rho_{0} \\tau + {\\left(-i \\, \\Gamma + i\\right)} \\kappa\\right)} u_{0}}\\right\\}, \\left\\{\\omega : 0\\right\\}\\right]","display":true}}︡
︠9e508615-92b7-482f-a992-24af261ea9d5︠
#For firehose
Mk1        = (M -I*omega*identity_matrix(5))
dispersion = Mk1.determinant().numerator().simplify_full()
solve(dispersion.subs(k2=0, eta=0), omega)
︡5380d3a4-d73d-4cee-bbd6-4400c1da5199︡{"tex":{"tex":"\\left[\\omega = -\\sqrt{\\frac{3 \\, B_{10}^{2}}{3 \\, B_{10}^{2} + 3 \\, \\Gamma u_{0} - \\psi_{0} + 3 \\, \\rho_{0}} - \\frac{3 \\, \\psi_{0}}{3 \\, B_{10}^{2} + 3 \\, \\Gamma u_{0} - \\psi_{0} + 3 \\, \\rho_{0}}} k_{1}, \\omega = \\sqrt{\\frac{3 \\, B_{10}^{2}}{3 \\, B_{10}^{2} + 3 \\, \\Gamma u_{0} - \\psi_{0} + 3 \\, \\rho_{0}} - \\frac{3 \\, \\psi_{0}}{3 \\, B_{10}^{2} + 3 \\, \\Gamma u_{0} - \\psi_{0} + 3 \\, \\rho_{0}}} k_{1}, \\omega = \\frac{i}{\\tau}, \\omega = 0\\right]","display":true}}︡
︠e4ce921b-aa0b-42fa-9b88-56a263235393︠
#For VISCOSITY_2D_TRANSVERSE
Mk1=(M -I*omega*identity_matrix(7)).subs(k2=0)
dispersion1=factor( (Mk1.determinant().numerator()/(-omega**2+B10**2*k1**2/(B10**2+Gamma*u0+rho0))/omega**2/(B10**2+Gamma*u0+rho0)).simplify_full())
Mk2=(M -I*omega*identity_matrix(7)).subs(k1=0)
dispersion2=factor((Mk2.determinant().numerator()/omega**4).subs(tau=0).simplify_full())
dispersion1
dispersion2
︡0a43489e-9967-4eea-b5ca-f2c22894afb6︡{"tex":{"tex":"3 i \\, \\Gamma^{2} k_{1}^{2} \\omega \\tau u_{0} - 3 i \\, \\Gamma k_{1}^{2} \\omega \\tau u_{0} - 3 i \\, \\Gamma \\omega^{3} \\tau u_{0} - 3 i \\, \\omega^{3} \\rho_{0} \\tau + 3 \\, \\Gamma^{2} k_{1}^{2} u_{0} + 4 i \\, \\eta k_{1}^{2} \\omega - 3 \\, \\Gamma k_{1}^{2} u_{0} - 3 \\, \\Gamma \\omega^{2} u_{0} - 3 \\, \\omega^{2} \\rho_{0}","display":true}}︡{"tex":{"tex":"-3 \\, \\Gamma^{2} k_{2}^{2} u_{0} - 3 \\, B_{10}^{2} k_{2}^{2} - i \\, \\eta k_{2}^{2} \\omega + 3 \\, B_{10}^{2} \\omega^{2} + 3 \\, \\Gamma k_{2}^{2} u_{0} + 3 \\, \\Gamma \\omega^{2} u_{0} + 3 \\, \\omega^{2} \\rho_{0}","display":true}}︡
︠97926afc-8799-41c0-89ed-ff806ce932c6︠
#For CONDUCTION_2D_TRANSVERSE
Mk1=(M -I*omega*identity_matrix(7)).subs(k2=0)
dispersion1=factor(Mk1.determinant().numerator()/omega/(B10**2*k1**2-B10**2*omega**2-Gamma*omega**2*u0-omega**2*rho0))
dispersion1
Mk2=(M -I*omega*identity_matrix(7)).subs(k1=0)
dispertion2=factor((Mk2.determinant().numerator()/omega**4).simplify_full())
dispertion2

︡b8e2058b-4c3a-4593-8180-2e1adf05e18c︡{"tex":{"tex":"i \\, \\Gamma^{2} k_{1}^{2} \\omega^{2} \\rho_{0} \\tau u_{0} - i \\, \\Gamma^{2} k_{1}^{4} \\kappa u_{0} - i \\, \\Gamma^{2} k_{1}^{2} \\kappa \\omega^{2} u_{0} - i \\, \\Gamma k_{1}^{2} \\omega^{2} \\rho_{0} \\tau u_{0} - i \\, \\Gamma \\omega^{4} \\rho_{0} \\tau u_{0} + i \\, \\Gamma k_{1}^{2} \\kappa \\omega^{2} \\rho_{0} - i \\, \\omega^{4} \\rho_{0}^{2} \\tau + 2 i \\, \\Gamma k_{1}^{4} \\kappa u_{0} + 3 i \\, \\Gamma k_{1}^{2} \\kappa \\omega^{2} u_{0} + i \\, \\Gamma \\kappa \\omega^{4} u_{0} + \\Gamma^{2} k_{1}^{2} \\omega \\rho_{0} u_{0} - i \\, k_{1}^{2} \\kappa \\omega^{2} \\rho_{0} - i \\, k_{1}^{4} \\kappa u_{0} - 2 i \\, k_{1}^{2} \\kappa \\omega^{2} u_{0} - i \\, \\kappa \\omega^{4} u_{0} - \\Gamma k_{1}^{2} \\omega \\rho_{0} u_{0} - \\Gamma \\omega^{3} \\rho_{0} u_{0} - \\omega^{3} \\rho_{0}^{2}","display":true}}︡{"tex":{"tex":"-{\\left(\\Gamma^{2} k_{2}^{2} u_{0} + B_{10}^{2} k_{2}^{2} - B_{10}^{2} \\omega^{2} - \\Gamma k_{2}^{2} u_{0} - \\Gamma \\omega^{2} u_{0} - \\omega^{2} \\rho_{0}\\right)} {\\left(i \\, \\Gamma \\omega \\rho_{0} \\tau u_{0} + i \\, \\omega \\rho_{0}^{2} \\tau - i \\, \\Gamma \\kappa \\omega u_{0} + i \\, \\kappa \\omega u_{0} + \\Gamma \\rho_{0} u_{0} + \\rho_{0}^{2}\\right)}","display":true}}︡
︠0e073419-5b35-48b1-a142-2fa8f582b47e︠
soln = solve(factor(dispersion1.subs(u0=0))==0, omega, solution_dict=True)
︡2708f88a-a1fb-4056-81ea-cc4547479f52︡
︠f70ddff2-e85a-409d-b616-9ed90c1b5748︠
var('epsilon')
soln[0][omega].subs(k1=1/epsilon).taylor(epsilon, 0, 1)
︡3e30a511-64ea-4469-8371-dd14c0edbf47︡{"tex":{"tex":"\\epsilon","display":true}}︡{"tex":{"tex":"\\frac{\\sqrt{\\Gamma - 1} \\epsilon \\sqrt{\\rho_{0}}}{8 \\, {\\left(\\Gamma \\sqrt{\\tau} - \\sqrt{\\tau}\\right)} \\sqrt{\\kappa} \\tau} + \\frac{i}{2 \\, \\tau} - \\frac{\\sqrt{\\Gamma - 1} \\sqrt{\\kappa}}{\\epsilon \\sqrt{\\rho_{0}} \\sqrt{\\tau}}","display":true}}︡
︠75a452bb-2c88-46b9-9b13-9978f7054988︠

dispersion_relation = (M -I*omega*identity_matrix(7)).determinant()#.simplify_full()
dispersion_relation.numerator()
#solns = solve(dispersion_relation.numerator()==0, omega, solution_dict=True)
︡73769226-231a-4772-9507-62c6ce482dc8︡{"tex":{"tex":"3 i \\, B_{10}^{2} \\Gamma^{2} k_{1}^{4} \\omega^{3} \\tau u_{0} - 3 i \\, B_{10}^{2} \\Gamma^{2} k_{1}^{2} \\omega^{5} \\tau u_{0} - 3 i \\, \\Gamma^{3} k_{1}^{2} \\omega^{5} \\tau u_{0}^{2} - 3 i \\, B_{10}^{2} \\Gamma k_{1}^{4} \\omega^{3} \\tau u_{0} + 3 i \\, B_{10}^{2} \\Gamma \\omega^{7} \\tau u_{0} - 3 i \\, \\Gamma^{2} k_{1}^{2} \\omega^{5} \\rho_{0} \\tau u_{0} + 3 i \\, \\Gamma^{2} k_{1}^{2} \\omega^{5} \\tau u_{0}^{2} + 3 i \\, \\Gamma^{2} \\omega^{7} \\tau u_{0}^{2} - 3 i \\, B_{10}^{2} k_{1}^{2} \\omega^{5} \\rho_{0} \\tau + 3 i \\, B_{10}^{2} \\omega^{7} \\rho_{0} \\tau + 3 \\, B_{10}^{2} \\Gamma^{2} k_{1}^{4} \\omega^{2} u_{0} - 3 \\, B_{10}^{2} \\Gamma^{2} k_{1}^{2} \\omega^{4} u_{0} + 3 i \\, \\Gamma k_{1}^{2} \\omega^{5} \\rho_{0} \\tau u_{0} + 6 i \\, \\Gamma \\omega^{7} \\rho_{0} \\tau u_{0} - 3 \\, \\Gamma^{3} k_{1}^{2} \\omega^{4} u_{0}^{2} + 4 i \\, B_{10}^{2} \\eta k_{1}^{4} \\omega^{3} - 4 i \\, B_{10}^{2} \\eta k_{1}^{2} \\omega^{5} + 3 i \\, \\omega^{7} \\rho_{0}^{2} \\tau - 3 \\, B_{10}^{2} \\Gamma k_{1}^{4} \\omega^{2} u_{0} - 4 i \\, \\Gamma \\eta k_{1}^{2} \\omega^{5} u_{0} + 3 \\, B_{10}^{2} \\Gamma \\omega^{6} u_{0} - 3 \\, \\Gamma^{2} k_{1}^{2} \\omega^{4} \\rho_{0} u_{0} + 3 \\, \\Gamma^{2} k_{1}^{2} \\omega^{4} u_{0}^{2} + 3 \\, \\Gamma^{2} \\omega^{6} u_{0}^{2} - 3 \\, B_{10}^{2} k_{1}^{2} \\omega^{4} \\rho_{0} - 4 i \\, \\eta k_{1}^{2} \\omega^{5} \\rho_{0} + 3 \\, B_{10}^{2} \\omega^{6} \\rho_{0} + 3 \\, \\Gamma k_{1}^{2} \\omega^{4} \\rho_{0} u_{0} + 6 \\, \\Gamma \\omega^{6} \\rho_{0} u_{0} + 3 \\, \\omega^{6} \\rho_{0}^{2}","display":true}}︡{"tex":{"tex":"3 i \\, B_{10}^{2} \\Gamma^{2} k_{1}^{4} \\omega^{3} \\tau u_{0} + 3 i \\, B_{10}^{2} \\Gamma^{2} k_{1}^{2} k_{2}^{2} \\omega^{3} \\tau u_{0} - 3 i \\, B_{10}^{2} \\Gamma^{2} k_{1}^{2} \\omega^{5} \\tau u_{0} - 3 i \\, \\Gamma^{3} k_{1}^{2} \\omega^{5} \\tau u_{0}^{2} - 3 i \\, \\Gamma^{3} k_{2}^{2} \\omega^{5} \\tau u_{0}^{2} - 3 i \\, B_{10}^{2} \\Gamma k_{1}^{4} \\omega^{3} \\tau u_{0} - 3 i \\, B_{10}^{2} \\Gamma k_{1}^{2} k_{2}^{2} \\omega^{3} \\tau u_{0} - 3 i \\, B_{10}^{2} \\Gamma k_{2}^{2} \\omega^{5} \\tau u_{0} + 3 i \\, B_{10}^{2} \\Gamma \\omega^{7} \\tau u_{0} - 3 i \\, \\Gamma^{2} k_{1}^{2} \\omega^{5} \\rho_{0} \\tau u_{0} - 3 i \\, \\Gamma^{2} k_{2}^{2} \\omega^{5} \\rho_{0} \\tau u_{0} + 3 i \\, \\Gamma^{2} k_{1}^{2} \\omega^{5} \\tau u_{0}^{2} + 3 i \\, \\Gamma^{2} k_{2}^{2} \\omega^{5} \\tau u_{0}^{2} + 3 i \\, \\Gamma^{2} \\omega^{7} \\tau u_{0}^{2} - 3 i \\, B_{10}^{2} k_{1}^{2} \\omega^{5} \\rho_{0} \\tau - 3 i \\, B_{10}^{2} k_{2}^{2} \\omega^{5} \\rho_{0} \\tau + 3 i \\, B_{10}^{2} \\omega^{7} \\rho_{0} \\tau + 3 \\, B_{10}^{2} \\Gamma^{2} k_{1}^{4} \\omega^{2} u_{0} + 3 \\, B_{10}^{2} \\Gamma^{2} k_{1}^{2} k_{2}^{2} \\omega^{2} u_{0} + 9 i \\, \\Gamma^{2} \\eta k_{1}^{2} k_{2}^{2} \\omega^{3} u_{0} - 3 \\, B_{10}^{2} \\Gamma^{2} k_{1}^{2} \\omega^{4} u_{0} + 3 i \\, \\Gamma k_{1}^{2} \\omega^{5} \\rho_{0} \\tau u_{0} + 3 i \\, \\Gamma k_{2}^{2} \\omega^{5} \\rho_{0} \\tau u_{0} + 6 i \\, \\Gamma \\omega^{7} \\rho_{0} \\tau u_{0} - 3 \\, \\Gamma^{3} k_{1}^{2} \\omega^{4} u_{0}^{2} - 3 \\, \\Gamma^{3} k_{2}^{2} \\omega^{4} u_{0}^{2} + 4 i \\, B_{10}^{2} \\eta k_{1}^{4} \\omega^{3} + 4 i \\, B_{10}^{2} \\eta k_{1}^{2} k_{2}^{2} \\omega^{3} - 4 i \\, B_{10}^{2} \\eta k_{1}^{2} \\omega^{5} + 3 i \\, \\omega^{7} \\rho_{0}^{2} \\tau - 3 \\, B_{10}^{2} \\Gamma k_{1}^{4} \\omega^{2} u_{0} - 3 \\, B_{10}^{2} \\Gamma k_{1}^{2} k_{2}^{2} \\omega^{2} u_{0} - 9 i \\, \\Gamma \\eta k_{1}^{2} k_{2}^{2} \\omega^{3} u_{0} - 3 \\, B_{10}^{2} \\Gamma k_{2}^{2} \\omega^{4} u_{0} - 4 i \\, \\Gamma \\eta k_{1}^{2} \\omega^{5} u_{0} - i \\, \\Gamma \\eta k_{2}^{2} \\omega^{5} u_{0} + 3 \\, B_{10}^{2}[...]","display":true}}︡
︠518c09f4-7ca8-4540-96db-a4f1558aaa52︠

︠ef4cef1c-1ec7-4c20-8608-cb624e26d0f4︠
factor((dispersion_relation.numerator().subs(rho0=0)).simplify_full() )
︡e8a7ffb1-b7a6-4214-b860-27c81ad9a741︡{"tex":{"tex":"{\\left(-3 i \\, \\Gamma k_{1}^{4} \\omega \\tau u_{0} - 3 i \\, \\Gamma k_{1}^{2} \\omega^{3} \\tau u_{0} + 3 i \\, k_{1}^{4} \\omega \\tau u_{0} + 6 i \\, k_{1}^{2} \\omega^{3} \\tau u_{0} + 3 i \\, \\omega^{5} \\tau u_{0} - 4 i \\, \\eta k_{1}^{4} \\omega - 3 \\, \\Gamma k_{1}^{4} u_{0} - 3 \\, \\Gamma k_{1}^{2} \\omega^{2} u_{0} + 3 \\, k_{1}^{4} u_{0} + 6 \\, k_{1}^{2} \\omega^{2} u_{0} + 3 \\, \\omega^{4} u_{0}\\right)} {\\left(\\Gamma - 1\\right)} \\kappa","display":true}}︡
︠aac2272f-808d-433f-b576-0520b5b4c610︠
M_k0 = M.subs(k1 = 0)
M_k0
︡a7636137-0dfe-4813-81cc-ef03dbadd25f︡{"tex":{"tex":"\\left(\\begin{array}{rrr}\n0 & 0 & -\\frac{2 \\, \\sqrt{u_{10}^{2} + 1} \\rho_{0} u_{10}^{2}}{{\\left(6 \\, \\Gamma \\tau u_{0} + 3 \\, \\rho_{0} \\tau - 4 \\, \\eta\\right)} u_{10}^{4} + 3 \\, \\Gamma \\tau u_{0} + {\\left(9 \\, \\Gamma \\tau u_{0} + 6 \\, \\rho_{0} \\tau - 4 \\, \\eta\\right)} u_{10}^{2} + 3 \\, \\rho_{0} \\tau} \\\\\n0 & 0 & \\frac{2 \\, \\sqrt{u_{10}^{2} + 1} u_{10}}{3 \\, \\Gamma \\tau u_{0} + {\\left(6 \\, \\Gamma \\tau u_{0} + 3 \\, \\rho_{0} \\tau - 4 \\, \\eta\\right)} u_{10}^{2} + 3 \\, \\rho_{0} \\tau} \\\\\n0 & 0 & -\\frac{3 \\, {\\left({\\left(2 \\, \\Gamma u_{0} + \\rho_{0}\\right)} u_{10}^{2} + \\Gamma u_{0} + \\rho_{0}\\right)} \\sqrt{u_{10}^{2} + 1}}{{\\left(6 \\, \\Gamma \\tau u_{0} + 3 \\, \\rho_{0} \\tau - 4 \\, \\eta\\right)} u_{10}^{4} + 3 \\, \\Gamma \\tau u_{0} + {\\left(9 \\, \\Gamma \\tau u_{0} + 6 \\, \\rho_{0} \\tau - 4 \\, \\eta\\right)} u_{10}^{2} + 3 \\, \\rho_{0} \\tau}\n\\end{array}\\right)","display":true}}︡
︠61034d72-0dad-4f8c-9ffb-d16ef40a5f73︠
Eqn_rho_omega = Eqn_rho.subs(delta_rho_dt=(I*omega*delta_rho), delta_u1_dt=(I*omega*delta_u1), delta_psi_dt=(I*omega*delta_psi))
Eqn_u1_omega  = Eqn_u1.subs(delta_rho_dt=(I*omega*delta_rho), delta_u1_dt=(I*omega*delta_u1), delta_psi_dt=(I*omega*delta_psi))
Eqn_psi_omega = Eqn_psi.subs(delta_rho_dt=(I*omega*delta_rho), delta_u1_dt=(I*omega*delta_u1), delta_psi_dt=(I*omega*delta_psi))

M = jacobian([Eqn_rho_omega, Eqn_u1_omega, Eqn_psi_omega], delta_vars)
M = M.apply_map(lambda x : x.simplify_full())
M_k0 = M.subs(k1 = 0)
M_k0
︡d923ecd4-f272-4902-a25b-f6d96aad8789︡{"tex":{"tex":"\\left(\\begin{array}{rrr}\n\\frac{i \\, \\omega u_{10}^{2} + i \\, \\omega}{\\sqrt{u_{10}^{2} + 1}} & \\frac{i \\, \\omega \\rho_{0} u_{10}}{\\sqrt{u_{10}^{2} + 1}} & 0 \\\\\n\\frac{i \\, \\omega u_{10}^{3} + i \\, \\omega u_{10}}{\\sqrt{u_{10}^{2} + 1}} & \\frac{2 i \\, \\omega \\rho_{0} u_{10}^{2} + i \\, \\omega \\rho_{0}}{\\sqrt{u_{10}^{2} + 1}} & \\frac{2 i \\, \\omega u_{10}^{3} + 2 i \\, \\omega u_{10}}{3 \\, \\sqrt{u_{10}^{2} + 1}} \\\\\n0 & \\frac{2 i \\, \\sqrt{u_{10}^{2} + 1} \\eta \\omega u_{10}}{\\tau u_{10}^{2} + \\tau} & \\frac{i \\, \\sqrt{u_{10}^{2} + 1} \\omega \\tau + 1}{\\tau}\n\\end{array}\\right)","display":true}}︡
︠dd503676-6a9a-459e-96d6-4f71ea8e2ffe︠
︠9a131943-5b6c-4546-a441-9d35661b548e︠

dis_princ_k2_full = (M.subs(k1=0,B20=0) - I*omega*identity_matrix(8)).determinant().numerator()
solve(dis_princ_k2_full/omega**6==0, omega)

︡449f7268-307b-4a64-bc4e-c52f9d99fc04︡{"tex":{"tex":"\\left[\\omega = -\\sqrt{\\frac{\\Gamma^{2} \\tau u_{0}}{B_{10}^{2} \\tau + \\Gamma \\tau u_{0} + \\rho_{0} \\tau} + \\frac{B_{10}^{2} \\tau}{B_{10}^{2} \\tau + \\Gamma \\tau u_{0} + \\rho_{0} \\tau} - \\frac{\\Gamma \\tau u_{0}}{B_{10}^{2} \\tau + \\Gamma \\tau u_{0} + \\rho_{0} \\tau} + \\frac{\\eta}{3 \\, {\\left(B_{10}^{2} \\tau + \\Gamma \\tau u_{0} + \\rho_{0} \\tau\\right)}}} k_{2}, \\omega = \\sqrt{\\frac{\\Gamma^{2} \\tau u_{0}}{B_{10}^{2} \\tau + \\Gamma \\tau u_{0} + \\rho_{0} \\tau} + \\frac{B_{10}^{2} \\tau}{B_{10}^{2} \\tau + \\Gamma \\tau u_{0} + \\rho_{0} \\tau} - \\frac{\\Gamma \\tau u_{0}}{B_{10}^{2} \\tau + \\Gamma \\tau u_{0} + \\rho_{0} \\tau} + \\frac{\\eta}{3 \\, {\\left(B_{10}^{2} \\tau + \\Gamma \\tau u_{0} + \\rho_{0} \\tau\\right)}}} k_{2}\\right]","display":true}}︡
︠907d8f00-e1e4-459c-954d-ccfa89b73424︠
Ma = M.subs(k2=0,B20=0,B10=0)
Mb = I*omega*identity_matrix(8)
dis_princ_k1_full = (Ma - Mb).determinant().numerator()
factor(dis_princ_k1_full/omega**4)
solve(factor(dis_princ_k1_full/omega**4).subs(omega=k1, eta=0)==0, tau)
solve(factor(dis_princ_k1_full/omega**4).subs(omega=k1, kappa=0)==0, tau)
solve(factor(dis_princ_k1_full/omega**4).subs(omega=k1)==0, tau)

︡8b4944da-e245-4584-b56e-1e53e936287f︡{"tex":{"tex":"-3 \\, \\Gamma^{2} k_{1}^{2} \\omega^{2} \\rho_{0} \\tau^{2} u_{0} + 3 \\, \\Gamma^{2} k_{1}^{4} \\kappa \\tau u_{0} + 3 \\, \\Gamma^{2} k_{1}^{2} \\kappa \\omega^{2} \\tau u_{0} + 3 \\, \\Gamma k_{1}^{2} \\omega^{2} \\rho_{0} \\tau^{2} u_{0} + 3 \\, \\Gamma \\omega^{4} \\rho_{0} \\tau^{2} u_{0} - 3 \\, \\Gamma k_{1}^{2} \\kappa \\omega^{2} \\rho_{0} \\tau + 3 \\, \\omega^{4} \\rho_{0}^{2} \\tau^{2} - 6 \\, \\Gamma k_{1}^{4} \\kappa \\tau u_{0} - 9 \\, \\Gamma k_{1}^{2} \\kappa \\omega^{2} \\tau u_{0} - 3 \\, \\Gamma \\kappa \\omega^{4} \\tau u_{0} + 4 \\, \\Gamma \\eta k_{1}^{4} \\kappa - 4 \\, \\eta k_{1}^{2} \\omega^{2} \\rho_{0} \\tau + 3 \\, k_{1}^{2} \\kappa \\omega^{2} \\rho_{0} \\tau + 3 \\, k_{1}^{4} \\kappa \\tau u_{0} + 6 \\, k_{1}^{2} \\kappa \\omega^{2} \\tau u_{0} + 3 \\, \\kappa \\omega^{4} \\tau u_{0} - 4 \\, \\eta k_{1}^{4} \\kappa","display":true}}︡{"tex":{"tex":"\\left[\\tau = -\\frac{{\\left(\\Gamma - 1\\right)} \\kappa \\rho_{0} - 2 \\, {\\left(\\Gamma^{2} - 3 \\, \\Gamma + 2\\right)} \\kappa u_{0}}{{\\left(\\Gamma^{2} - 2 \\, \\Gamma\\right)} \\rho_{0} u_{0} - \\rho_{0}^{2}}, \\tau = 0\\right]","display":true}}︡{"tex":{"tex":"\\left[\\tau = -\\frac{4 \\, \\eta}{3 \\, {\\left({\\left(\\Gamma^{2} - 2 \\, \\Gamma\\right)} u_{0} - \\rho_{0}\\right)}}, \\tau = 0\\right]","display":true}}︡{"tex":{"tex":"\\left[\\tau = \\frac{6 \\, {\\left(\\Gamma^{2} - 3 \\, \\Gamma + 2\\right)} \\kappa u_{0} - {\\left(3 \\, {\\left(\\Gamma - 1\\right)} \\kappa + 4 \\, \\eta\\right)} \\rho_{0} - \\sqrt{36 \\, {\\left(\\Gamma^{4} - 6 \\, \\Gamma^{3} + 13 \\, \\Gamma^{2} - 12 \\, \\Gamma + 4\\right)} \\kappa^{2} u_{0}^{2} - {\\left(24 \\, {\\left(\\Gamma - 1\\right)} \\eta \\kappa - 9 \\, {\\left(\\Gamma^{2} - 2 \\, \\Gamma + 1\\right)} \\kappa^{2} - 16 \\, \\eta^{2}\\right)} \\rho_{0}^{2} + 12 \\, {\\left(4 \\, {\\left(\\Gamma^{3} - 4 \\, \\Gamma^{2} + 5 \\, \\Gamma - 2\\right)} \\eta \\kappa - 3 \\, {\\left(\\Gamma^{3} - 4 \\, \\Gamma^{2} + 5 \\, \\Gamma - 2\\right)} \\kappa^{2}\\right)} \\rho_{0} u_{0}}}{6 \\, {\\left({\\left(\\Gamma^{2} - 2 \\, \\Gamma\\right)} \\rho_{0} u_{0} - \\rho_{0}^{2}\\right)}}, \\tau = \\frac{6 \\, {\\left(\\Gamma^{2} - 3 \\, \\Gamma + 2\\right)} \\kappa u_{0} - {\\left(3 \\, {\\left(\\Gamma - 1\\right)} \\kappa + 4 \\, \\eta\\right)} \\rho_{0} + \\sqrt{36 \\, {\\left(\\Gamma^{4} - 6 \\, \\Gamma^{3} + 13 \\, \\Gamma^{2} - 12 \\, \\Gamma + 4\\right)} \\kappa^{2} u_{0}^{2} - {\\left(24 \\, {\\left(\\Gamma - 1\\right)} \\eta \\kappa - 9 \\, {\\left(\\Gamma^{2} - 2 \\, \\Gamma + 1\\right)} \\kappa^{2} - 16 \\, \\eta^{2}\\right)} \\rho_{0}^{2} + 12 \\, {\\left(4 \\, {\\left(\\Gamma^{3} - 4 \\, \\Gamma^{2} + 5 \\, \\Gamma - 2\\right)} \\eta \\kappa - 3 \\, {\\left(\\Gamma^{3} - 4 \\, \\Gamma^{2} + 5 \\, \\Gamma - 2\\right)} \\kappa^{2}\\right)} \\rho_{0} u_{0}}}{6 \\, {\\left({\\left(\\Gamma^{2} - 2 \\, \\Gamma\\right)} \\rho_{0} u_{0} - \\rho_{0}^{2}\\right)}}\\right]","display":true}}︡
︠6b8cbd48-871f-4288-b4a2-2dc791995954︠
︡a805685a-da21-4b96-8a4b-166be22ecfe7︡
︠7a9800c8-4aa3-4115-868c-effddef35907︠
Ma = M.subs(k2=0,B20=0)
Mb = I*omega*identity_matrix(4)
dis_princ_k1_full = (Ma - Mb).determinant().numerator()
factor(dis_princ_k1_full)
solve(factor(dis_princ_k1_full).subs(omega=k1)==0, tau)
︡fd858af1-78bc-4633-ab31-dda4074a05c3︡{"tex":{"tex":"-\\Gamma^{2} k_{1}^{2} \\omega^{2} \\rho_{0} \\tau u_{0} + \\Gamma^{2} k_{1}^{4} \\kappa u_{0} + \\Gamma^{2} k_{1}^{2} \\kappa \\omega^{2} u_{0} + \\Gamma k_{1}^{2} \\omega^{2} \\rho_{0} \\tau u_{0} + \\Gamma \\omega^{4} \\rho_{0} \\tau u_{0} - \\Gamma k_{1}^{2} \\kappa \\omega^{2} \\rho_{0} + \\omega^{4} \\rho_{0}^{2} \\tau - 2 \\, \\Gamma k_{1}^{4} \\kappa u_{0} - 3 \\, \\Gamma k_{1}^{2} \\kappa \\omega^{2} u_{0} - \\Gamma \\kappa \\omega^{4} u_{0} + k_{1}^{2} \\kappa \\omega^{2} \\rho_{0} + k_{1}^{4} \\kappa u_{0} + 2 \\, k_{1}^{2} \\kappa \\omega^{2} u_{0} + \\kappa \\omega^{4} u_{0}","display":true}}︡{"tex":{"tex":"\\left[\\tau = -\\frac{{\\left(\\Gamma - 1\\right)} \\kappa \\rho_{0} - 2 \\, {\\left(\\Gamma^{2} - 3 \\, \\Gamma + 2\\right)} \\kappa u_{0}}{{\\left(\\Gamma^{2} - 2 \\, \\Gamma\\right)} \\rho_{0} u_{0} - \\rho_{0}^{2}}\\right]","display":true}}︡


︠49b43220-918b-4bbd-a412-638da5235ceb︠
︠c70cdcd3-9461-48fe-82f8-2c1451457837︠
Eqn_phi
︡669e5c39-e4e2-4736-be4b-beabc2de335f︡{"tex":{"tex":"\\frac{{\\left(i \\, \\Gamma - i\\right)} \\delta_{u} k_{1} \\kappa \\rho_{0} + \\delta_{\\phi_{\\mathit{dt}}} \\rho_{0}^{2} \\tau + \\delta_{\\phi} \\rho_{0}^{2} + {\\left({\\left(-i \\, \\Gamma + i\\right)} \\delta_{\\rho} k_{1} \\kappa + {\\left(\\Gamma - 1\\right)} \\delta_{\\mathit{u1}_{\\mathit{dt}}} \\kappa \\rho_{0}\\right)} u_{0}}{\\rho_{0}^{2} \\tau}","display":true}}︡
︠1f58881c-0700-4115-a8fa-931d73fca61a︠
M_HL=Matrix([[omega*rho0**2/u0/(Gamma-1), omega*(rho0+u0)*rho0/u0/(Gamma-1), -rho0, 0],[omega*(rho0+u0)*rho0/u0/(Gamma-1), omega/(Gamma-1)/u0*(rho0**2+2*rho0*u0+Gamma*u0**2), -(rho0+Gamma*u0), -1],[-rho0, -(rho0+Gamma*u0), omega*(rho0+Gamma*u0), omega],[0, -kappa/tau*(Gamma-1)*u0/rho0, kappa/tau*(Gamma-1)*u0/rho0*omega, omega]])
M_HL
solve(M_HL.determinant().subs(omega=1)==0, tau)

︡961ee614-2347-4e0b-91cd-279d8cfae1a4︡{"stdout":"[                           omega*rho0^2/((Gamma - 1)*u0)                  omega*(rho0 + u0)*rho0/((Gamma - 1)*u0)                                                    -rho0                                                        0]\n[                 omega*(rho0 + u0)*rho0/((Gamma - 1)*u0) (Gamma*u0^2 + rho0^2 + 2*rho0*u0)*omega/((Gamma - 1)*u0)                                         -Gamma*u0 - rho0                                                       -1]\n[                                                   -rho0                                         -Gamma*u0 - rho0                                  (Gamma*u0 + rho0)*omega                                                    omega]\n[                                                       0                         -(Gamma - 1)*kappa*u0/(rho0*tau)                    (Gamma - 1)*kappa*omega*u0/(rho0*tau)                                                    omega]\n"}︡{"stdout":"[tau == -((Gamma - 1)*kappa*rho0 - 2*(Gamma^2 - 3*Gamma + 2)*kappa*u0)/((Gamma^2 - 2*Gamma)*rho0*u0 - rho0^2)]\n"}︡
︠84abc116-dfcd-414e-8e63-08f717a86f66︠
Ma = M.subs(k2=0,B20=0)
Mb = I*omega*identity_matrix(5)
dis_princ_k1_full = (Ma - Mb).determinant().numerator()
factor(dis_princ_k1_full/omega)
solve(factor(dis_princ_k1_full).subs(omega=k1)==0, tau)
︡0651c3be-9fdf-43e5-bab9-e2e6e7017aee︡{"stdout":"3*I*Gamma^2*k1^2*omega^2*rho0*tau^2*u0 - 3*I*Gamma^2*k1^4*kappa*tau*u0 - 3*I*Gamma^2*k1^2*kappa*omega^2*tau*u0 - 3*I*Gamma*k1^2*omega^2*rho0*tau^2*u0 - 3*I*Gamma*omega^4*rho0*tau^2*u0 + 3*I*Gamma*k1^2*kappa*omega^2*rho0*tau - 3*I*omega^4*rho0^2*tau^2 + 6*I*Gamma*k1^4*kappa*tau*u0 + 9*I*Gamma*k1^2*kappa*omega^2*tau*u0 + 3*I*Gamma*kappa*omega^4*tau*u0 - 4*I*Gamma*eta*k1^4*kappa + 4*I*eta*k1^2*omega^2*rho0*tau - 3*I*k1^2*kappa*omega^2*rho0*tau - 3*I*k1^4*kappa*tau*u0 - 6*I*k1^2*kappa*omega^2*tau*u0 - 3*I*kappa*omega^4*tau*u0 + 4*I*eta*k1^4*kappa\n"}︡{"stdout":"[tau == 1/6*(6*(Gamma^2 - 3*Gamma + 2)*kappa*u0 - (3*(Gamma - 1)*kappa + 4*eta)*rho0 - sqrt(36*(Gamma^4 - 6*Gamma^3 + 13*Gamma^2 - 12*Gamma + 4)*kappa^2*u0^2 - (24*(Gamma - 1)*eta*kappa - 9*(Gamma^2 - 2*Gamma + 1)*kappa^2 - 16*eta^2)*rho0^2 + 12*(4*(Gamma^3 - 4*Gamma^2 + 5*Gamma - 2)*eta*kappa - 3*(Gamma^3 - 4*Gamma^2 + 5*Gamma - 2)*kappa^2)*rho0*u0))/((Gamma^2 - 2*Gamma)*rho0*u0 - rho0^2), tau == 1/6*(6*(Gamma^2 - 3*Gamma + 2)*kappa*u0 - (3*(Gamma - 1)*kappa + 4*eta)*rho0 + sqrt(36*(Gamma^4 - 6*Gamma^3 + 13*Gamma^2 - 12*Gamma + 4)*kappa^2*u0^2 - (24*(Gamma - 1)*eta*kappa - 9*(Gamma^2 - 2*Gamma + 1)*kappa^2 - 16*eta^2)*rho0^2 + 12*(4*(Gamma^3 - 4*Gamma^2 + 5*Gamma - 2)*eta*kappa - 3*(Gamma^3 - 4*Gamma^2 + 5*Gamma - 2)*kappa^2)*rho0*u0))/((Gamma^2 - 2*Gamma)*rho0*u0 - rho0^2)]\n"}︡
︠3a57fea8-9f33-4198-8b4e-2b940e1cac0e︠









