︠fbddea82-3103-42dd-bbf9-04c4a571eee5ao︠
%auto
typeset_mode(True)
︡9df0a4aa-9d61-4c0e-82d3-4d0c824af3f7︡{"auto":true}︡
︠0d1a5418-4348-4f75-9af5-685a94ee8334i︠
%md
## $\mathtt{balbusaur}$<br> ##
#### A framework for automated linear analysis

#### Initialize the operators $\mathtt{d\_dt()}$, $\mathtt{d\_dX1()}$ and $\mathtt{d\_dX2()}$ acting on variables:
1. Density $\rho$
2. Internal energy $u$
3. Velocity in $X^1$ direction $u^1$
4. Velocity in $X^2$ direction $u^2$
5. Velocity in $X^3$ direction $u^3$
6. Magnetic field in $X^1$ direction $B^1$
7. Magnetic field in $X^2$ direction $B^2$
8. Magnetic field in $X^3$ direction $B^3$
9. Heat flux magnitude $q$
10. Pressure anisotropy magnitude $dp$

* Mean variables      : $\rho_0$, $u_0$, $u^1_0$, $u^2_0$, $u^3_0$, $B^1_0$, $B^2_0$, $B^3_0$, $q_0$, $dp_0$
* Perturbed variables : $\delta_\rho$, $\delta_u$, $\delta_{u^1}$, $\delta_{u^2}$, $\delta_{u^3}$, $\delta_{B^1}$, $\delta_{B^2}$, $\delta_{B^3}$, $\delta_q$, $\delta_{dp}$


︡de554342-ffe2-4092-97f1-65f5f7e57804︡{"md":"## $\\mathtt{balbusaur}$<br> ##\n#### A framework for automated linear analysis\n\n#### Initialize the operators $\\mathtt{d\\_dt()}$, $\\mathtt{d\\_dX1()}$ and $\\mathtt{d\\_dX2()}$ acting on variables:\n1. Density $\\rho$\n2. Internal energy $u$\n3. Velocity in $X^1$ direction $u^1$\n4. Velocity in $X^2$ direction $u^2$\n5. Velocity in $X^3$ direction $u^3$\n6. Magnetic field in $X^1$ direction $B^1$\n7. Magnetic field in $X^2$ direction $B^2$\n8. Magnetic field in $X^3$ direction $B^3$\n9. Heat flux magnitude $q$\n10. Pressure anisotropy magnitude $dp$\n\n* Mean variables      : $\\rho_0$, $u_0$, $u^1_0$, $u^2_0$, $u^3_0$, $B^1_0$, $B^2_0$, $B^3_0$, $q_0$, $dp_0$\n* Perturbed variables : $\\delta_\\rho$, $\\delta_u$, $\\delta_{u^1}$, $\\delta_{u^2}$, $\\delta_{u^3}$, $\\delta_{B^1}$, $\\delta_{B^2}$, $\\delta_{B^3}$, $\\delta_q$, $\\delta_{dp}$\n\n\n"}︡
︠81d0b45a-c4a1-4a7c-8b52-5a3198afcfbe︠
def linearize(term):
    return taylor(term, (delta_rho, 0), \
                        (delta_u, 0),   \
                        (delta_u1, 0),  \
                        (delta_u2, 0),  \
                        (delta_u3, 0),  \
                        (delta_B1, 0),  \
                        (delta_B2, 0),  \
                        (delta_B3, 0),  \
                        (delta_q, 0), \
                        (delta_dp, 0), \
                        (delta_rho_dt, 0), \
                        (delta_u_dt, 0),   \
                        (delta_u1_dt, 0),  \
                        (delta_u2_dt, 0),  \
                        (delta_u3_dt, 0),  \
                        (delta_B1_dt, 0),  \
                        (delta_B2_dt, 0),  \
                        (delta_B3_dt, 0),  \
                        (delta_q_dt, 0), \
                        (delta_dp_dt, 0), 1 \
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
           + term.coefficient(delta_q)   * I * k1 * delta_q \
           + term.coefficient(delta_dp)  * I * k1 * delta_dp

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
           + term.coefficient(delta_q)   * I * k2 * delta_q \
           + term.coefficient(delta_dp)  * I * k2 * delta_dp

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
           + term.coefficient(delta_q)   * delta_q_dt \
           + term.coefficient(delta_dp)  * delta_dp_dt

    return expr

︡09b9a633-2cf1-4611-8519-e7f347063589︡
︠82740d77-8853-4cb1-a194-df8f1673c123i︠
%md
#### Options:

1. $\mathtt{EVOLVE\_B\_FIELDS}$       : 0 or 1
2. $\mathtt{CONDUCTION}$              : 0 or 1
3. $\mathtt{VISCOSITY}$               : 0 or 1
4. $\mathtt{FAKE\_EMHD}$              : 0 or 1
5. $\mathtt{TURN\_OFF\_MEAN\_B2}$      : 0 or 1
6. $\mathtt{TURN\_OFF\_K2\_PERTURBATIONS}$ : 0 or 1
7. $\mathtt{PRINCIPAL\_COEFFICIENTS}$ : 0 or 1
︡d8bdd1ba-2212-41fb-8002-d1cb2a07f866︡{"md":"#### Options:\n\n1. $\\mathtt{EVOLVE\\_B\\_FIELDS}$       : 0 or 1\n2. $\\mathtt{CONDUCTION}$              : 0 or 1\n3. $\\mathtt{VISCOSITY}$               : 0 or 1\n4. $\\mathtt{FAKE\\_EMHD}$              : 0 or 1\n5. $\\mathtt{TURN\\_OFF\\_MEAN\\_B2}$      : 0 or 1\n6. $\\mathtt{TURN\\_OFF\\_K2\\_PERTURBATIONS}$ : 0 or 1\n7. $\\mathtt{PRINCIPAL\\_COEFFICIENTS}$ : 0 or 1\n"}︡
︠c953ed04-b6ab-4a98-a57f-2711f7d0d62b︠
# Spatiotemporal variables
t, omega, k1, k2 = var('t, omega, k1, k2')

# Constants:
# Gamma : Adiabatic index
# kappa : Heat conductivity
# eta   : shear viscosity
# tau   : relaxation time scale
# phi   : dimensionless coefficient for conduction
# psi   : dimensionless coefficient for viscosity
Gamma, kappa, eta, tau, phi, psi = var('Gamma, kappa, eta, tau, phi, psi')

# Background mean values: Symbolic variables
rho0, u0, u10, u20, u30, B10, B20, B30, q0, dp0 = var('rho0, u0, u10, u20, u30, B10, B20, B30, psi0, phi0')

# Perturbations in space
delta_rho, delta_u, delta_u1, delta_u2, delta_u3, delta_B1, delta_B2, delta_B3, delta_q, delta_dp = \
    var('delta_rho, delta_u, delta_u1, delta_u2, delta_u3, delta_B1, delta_B2, delta_B3, delta_q, delta_dp')

# Perturbations in time
delta_rho_dt, delta_u_dt, delta_u1_dt, delta_u2_dt, delta_u3_dt, delta_B1_dt, delta_B2_dt, delta_B3_dt, delta_q_dt, delta_dp_dt = \
    var('delta_rho_dt, delta_u_dt, delta_u1_dt, delta_u2_dt, delta_u3_dt, delta_B1_dt, delta_B2_dt, delta_B3_dt, delta_q_dt, delta_dp_dt')

# Inputs:

FAKE_EMHD                 = 0
EVOLVE_B_FIELDS           = 1
CONDUCTION                = 1
VISCOSITY                 = 1
TURN_OFF_MEAN_B2          = 0
TURN_OFF_K2_PERTURBATIONS = 0
PRINCIPAL_COEFFICIENTS    = 1

# Equilibrium states:
# rho0 and u0 are NEVER set to zero
u10 = 0
u20 = 0
u30 = 0
B30 = 0
q0  = 0
dp0 = 0
if (TURN_OFF_MEAN_B2):
    B20 = 0
if (TURN_OFF_K2_PERTURBATIONS):
    k2 = 0

rho = rho0 + delta_rho
u   = u0   + delta_u
u1  = u10  + delta_u1
u2  = u20  + delta_u2
u3  = u30  + delta_u3
B1  = B10  + delta_B1
B2  = B20  + delta_B2
B3  = B30  + delta_B3
q   = q0   + delta_q
dp  = dp0  + delta_dp

# Introducing:
# chi : Thermal diffusivity
# nu  : kinematic viscosity
# cs  : sound speed
P   = (Gamma - 1)*u
T   = P/rho
cs  = sqrt(Gamma * P/ (rho + Gamma*u) )
chi   = phi * cs**2 * tau
nu    = psi * cs**2 * tau
kappa = rho * chi
eta   = rho * nu

# Inputs for numerical diagonalization for finite k modes
k1_num    = 2*pi
k2_num    = 4*pi

rho0_num = 1
u0_num   = 2
u10_num  = 0
u20_num  = 0
u30_num  = 0
B10_num  = 0.1
B20_num  = 0.3
B30_num  = 0
q0_num   = 0
dp0_num  = 0

Gamma_num = 4/3
tau_num   = 1
phi_num   = 1
psi_num   = 1
︡0d7a9f74-2651-46a2-ac61-fc8962ff650b︡
︠c3616227-8585-41a4-ba9d-fe429350abe3i︠
%md

#### All the physics is below
︡6af4ab04-815e-44d4-bc62-363b012381ec︡{"md":"\n#### All the physics is below\n"}︡
︠a7637a36-f643-4769-aed3-98f0e4e300b9︠
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

def delta(mu, nu):
    if (mu==nu):
        return 1
    else:
        return 0

def TUpDown(mu, nu):

    return   (rho + u + P + bsqr + FAKE_EMHD*bsqr/6)*ucon[mu]*ucov[nu] \
           + (P + bsqr/2 + FAKE_EMHD*bsqr/6)*delta(mu, nu) - (1 + FAKE_EMHD/2)*bcon[mu]*bcov[nu] \
           + q/sqrt(bsqr) * (bcon[mu]*ucov[nu] + ucon[mu]*bcov[nu]) \
           + dp/bsqr      * (bcon[mu]*bcov[nu]) - dp/3*(ucon[mu]*ucov[nu] + delta(mu, nu))

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
Eqn_u3  = linearize(d_dt(TUpDown(0, 3)) + d_dX1(TUpDown(1, 3)) + d_dX2(TUpDown(2, 3)))

Eqn_B1  = linearize(d_dt(B1) + d_dX2(bcon[1]*ucon[2] - bcon[2]*ucon[1]) )
Eqn_B2  = linearize(d_dt(B2) + d_dX1(bcon[2]*ucon[1] - bcon[1]*ucon[2]) )
Eqn_B3  = linearize(d_dt(B3) + d_dX1(bcon[3]*ucon[1] - bcon[1]*ucon[3]) + d_dX2(bcon[3]*ucon[2] - bcon[2]*ucon[3]) )

q_relaxed  = (bcov[0]*qconEckart(0) + bcov[1]*qconEckart(1) + bcov[2]*qconEckart(2) + bcov[3]*qconEckart(3) )/sqrt(bsqr)
dp_relaxed = 0
for nu in xrange(4):
    dp_relaxed = dp_relaxed - 3*eta/bsqr * (bcon[nu]* (bcon[0]*d_dt(ucov[nu]) + bcon[1]*d_dX1(ucov[nu]) + bcon[2]*d_dX2(ucov[nu])) )

dp_relaxed = dp_relaxed + eta*(d_dt(ucon[0]) + d_dX1(ucon[1]) + d_dX2(ucon[2]) )

Eqn_q  = linearize(ucon[0]*d_dt(q)  + ucon[1]*d_dX1(q)  + ucon[2]*d_dX2(q)  + (PRINCIPAL_COEFFICIENTS*q  - q_relaxed)/tau)
Eqn_dp = linearize(ucon[0]*d_dt(dp) + ucon[1]*d_dX1(dp) + ucon[2]*d_dX2(dp) + (PRINCIPAL_COEFFICIENTS*dp - dp_relaxed)/tau)

Eqns          = [Eqn_rho==0, Eqn_u==0, Eqn_u1==0, Eqn_u2==0, Eqn_u3==0]
delta_vars    = [delta_rho, delta_u, delta_u1, delta_u2, delta_u3]
delta_vars_dt = [delta_rho_dt, delta_u_dt, delta_u1_dt, delta_u2_dt, delta_u3_dt]

if (EVOLVE_B_FIELDS):
    Eqns          = Eqns          + [Eqn_B1==0,   Eqn_B2==0,   Eqn_B3==0]
    delta_vars    = delta_vars    + [delta_B1,    delta_B2,    delta_B3]
    delta_vars_dt = delta_vars_dt + [delta_B1_dt, delta_B2_dt, delta_B3_dt]

if (CONDUCTION):
    Eqns          = Eqns          + [Eqn_q==0]
    delta_vars    = delta_vars    + [delta_q]
    delta_vars_dt = delta_vars_dt + [delta_q_dt]

if (VISCOSITY):
    Eqns          = Eqns          + [Eqn_dp==0]
    delta_vars    = delta_vars    + [delta_dp]
    delta_vars_dt = delta_vars_dt + [delta_dp_dt]

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
pretty_print("Analytic eigenvalues and eigenvectors in the $k_1, k_2 \\rightarrow 0$ limit : ", )
M.subs(k1=0, k2=0).eigenvectors_right()

# Numerical diagonalization:

M_numerical = M.subs(rho0=rho0_num, u0=u0_num, u10=u10_num, u20=u20_num, u30=u30_num, \
                     B10=B10_num, B20=B20_num, B30=B30_num, q0=q0_num,   dp0=dp0_num, \
                     Gamma=Gamma_num, tau=tau_num, phi=phi_num, psi=psi_num, \
                     k1=k1_num, k2=k2_num \
                    )

M_numerical = M_numerical.change_ring(CDF)
eigenvecs   = M_numerical.eigenvectors_right()

pretty_print("Numerical eigenvalues and eigenvectors for $k_1 = $", k1_num, " , $k_2$ = ", k2_num, ":\n")

rho_index = 0
u_index   = 1
u1_index  = 2
u2_index  = 3
u3_index  = 4
if (CONDUCTION==1 and VISCOSITY==0):
    q_index = 5
if (CONDUCTION==0 and VISCOSITY==1):
    dp_index = 5
if (CONDUCTION==1 and VISCOSITY==1):
    q_index = 5
    dp_index = 6

if (EVOLVE_B_FIELDS):
    b1_index  = 5
    b2_index  = 6
    b3_index  = 7

    if (CONDUCTION==1 and VISCOSITY==0):
        q_index = 8

    if (CONDUCTION==0 and VISCOSITY==1):
        dp_index = 8

    if (CONDUCTION==1 and VISCOSITY==1):
        q_index  = 8
        dp_index = 9

for i in xrange(len(eigenvecs)):
    print("--------------------------")
    print("Eigenvalue   = ", eigenvecs[i][0])
    print(delta_rho,  " = ", eigenvecs[i][1][0][0])
    print(delta_u,    " = ", eigenvecs[i][1][0][1])
    print(delta_u1,   " = ", eigenvecs[i][1][0][2])
    print(delta_u2,   " = ", eigenvecs[i][1][0][3])
    print(delta_u3,   " = ", eigenvecs[i][1][0][4])

    if (EVOLVE_B_FIELDS):
        print(delta_B1, " = ", eigenvecs[i][1][0][5])
        print(delta_B2, " = ", eigenvecs[i][1][0][6])
        print(delta_B3, " = ", eigenvecs[i][1][0][7])

    if (CONDUCTION):
        print(delta_q, " = ", eigenvecs[i][1][0][q_index])

    if (VISCOSITY):
        print(delta_dp," = ", eigenvecs[i][1][0][dp_index])
︡6cfda601-82dc-4778-8bfd-b556a8bc6e2e︡{"html":"<div align='center'>Linearized system : </div>"}︡{"stdout":"\n\n"}︡{"html":"<div align='center'>$\\displaystyle \\left(\\begin{array}{r}\n\\delta_{\\rho_{\\mathit{dt}}} \\\\\n\\delta_{u_{\\mathit{dt}}} \\\\\n\\delta_{\\mathit{u1}_{\\mathit{dt}}} \\\\\n\\delta_{\\mathit{u2}_{\\mathit{dt}}} \\\\\n\\delta_{\\mathit{u3}_{\\mathit{dt}}} \\\\\n\\delta_{\\mathit{B1}_{\\mathit{dt}}} \\\\\n\\delta_{\\mathit{B2}_{\\mathit{dt}}} \\\\\n\\delta_{\\mathit{B3}_{\\mathit{dt}}} \\\\\n\\delta_{q_{\\mathit{dt}}} \\\\\n\\delta_{\\mathit{dp}_{\\mathit{dt}}}\n\\end{array}\\right)$  =  $\\displaystyle \\left(\\begin{array}{rrrrrrrrrr}\n0 &amp; 0 &amp; -i \\, k_{1} \\rho_{0} &amp; -i \\, k_{2} \\rho_{0} &amp; 0 &amp; 0 &amp; 0 &amp; 0 &amp; 0 &amp; 0 \\\\\n0 &amp; 0 &amp; -i \\, \\Gamma k_{1} u_{0} &amp; -i \\, \\Gamma k_{2} u_{0} &amp; 0 &amp; 0 &amp; 0 &amp; 0 &amp; \\frac{-i \\, B_{10} k_{1} - i \\, B_{20} k_{2}}{\\sqrt{B_{10}^{2} + B_{20}^{2}}} &amp; 0 \\\\\n-\\frac{{\\left(i \\, \\Gamma^{3} - 2 i \\, \\Gamma^{2} + i \\, \\Gamma\\right)} B_{10}^{2} k_{1} \\phi u_{0}^{2} + {\\left(i \\, \\Gamma^{3} - 2 i \\, \\Gamma^{2} + i \\, \\Gamma\\right)} B_{10} B_{20} k_{2} \\phi u_{0}^{2}}{{\\left(2 \\, \\Gamma \\rho_{0}^{2} u_{0} + {\\left(\\Gamma^{2} - {\\left(\\Gamma^{3} - 2 \\, \\Gamma^{2} + \\Gamma\\right)} \\phi\\right)} \\rho_{0} u_{0}^{2} + \\rho_{0}^{3}\\right)} B_{10}^{2} + {\\left(2 \\, \\Gamma \\rho_{0}^{2} u_{0} + {\\left(\\Gamma^{2} - {\\left(\\Gamma^{3} - 2 \\, \\Gamma^{2} + \\Gamma\\right)} \\phi\\right)} \\rho_{0} u_{0}^{2} + \\rho_{0}^{3}\\right)} B_{20}^{2}} &amp; \\frac{{\\left({\\left(-i \\, \\Gamma + i\\right)} k_{1} \\rho_{0} + {\\left({\\left(i \\, \\Gamma^{3} - 2 i \\, \\Gamma^{2} + i \\, \\Gamma\\right)} k_{1} \\phi + {\\left(-i \\, \\Gamma^{2} + i \\, \\Gamma\\right)} k_{1}\\right)} u_{0}\\right)} B_{10}^{4} + {\\left({\\left(-i \\, \\Gamma + i\\right)} k_{2} \\rho_{0} + {\\left({\\left(i \\, \\Gamma^{3} - 2 i \\, \\Gamma^{2} + i \\, \\Gamma\\right)} k_{2} \\phi + {\\left(-i \\, \\Gamma^{2} + i \\, \\Gamma\\right)} k_{2}\\right)} u_{0}\\right)} B_{10} B_{20}^{3} + {\\left({\\left(-i \\, \\Gamma + i\\right)} k_{1} \\rho_{0}^{2} + {\\left({\\left(i \\, \\Gamma^{3} - 2 i \\, \\Gamma^{2} + i \\, \\Gamma\\right)} k_{1} \\phi + {\\left(-2 i \\, \\Gamma^{2} + 2 i \\, \\Gamma\\right)} k_{1}\\right)} \\rho_{0} u_{0} + {\\left({\\left(i \\, \\Gamma^{4} - 2 i \\, \\Gamma^{3} + i \\, \\Gamma^{2}\\right)} k_{1} \\phi + {\\left(-i \\, \\Gamma^{3} + i \\, \\Gamma^{2}\\right)} k_{1}\\right)} u_{0}^{2}\\right)} B_{10}^{2} + {\\left({\\left(-i \\, \\Gamma + i\\right)} k_{1} \\rho_{0}^{2} + {\\left(-2 i \\, \\Gamma^{2} + 2 i \\, \\Gamma\\right)} k_{1} \\rho_{0} u_{0} + {\\left({\\left(-i \\, \\Gamma + i\\right)} k_{1} \\rho_{0} + {\\left({\\left(i \\, \\Gamma^{3} - 2 i \\, \\Gamma^{2} + i \\, \\Gamma\\right)} k_{1} \\phi + {\\left(-i \\, \\Gamma^{2} + i \\, \\Gamma\\right)} k_{1}\\right)} u_{0}\\right)} B_{10}^{2} + {\\left({\\left(i \\, \\Gamma^{4} - 3 i \\, \\Gamma^{3} + 3 i \\, \\Gamma^{2} - i \\, \\Gamma\\right)} k_{1} \\phi + {\\left(-i \\, \\Gamma^{3} + i \\, \\Gamma^{2}\\right)} k_{1}\\right)} u_{0}^{2}\\right)} B_{20}^{2} + {\\left({\\left({\\left(-i \\, \\Gamma + i\\right)} k_{2} \\rho_{0} + {\\left({\\left(i \\, \\Gamma^{3} - 2 i \\, \\Gamma^{2} + i \\, \\Gamma\\right)} k_{2} \\phi + {\\left(-i \\, \\Gamma^{2} + i \\, \\Gamma\\right)} k_{2}\\right)} u_{0}\\right)} B_{10}^{3} + {\\left({\\left(i \\, \\Gamma^{3} - 2 i \\, \\Gamma^{2} + i \\, \\Gamma\\right)} k_{2} \\phi \\rho_{0} u_{0} + {\\left(i \\, \\Gamma^{3} - 2 i \\, \\Gamma^{2} + i \\, \\Gamma\\right)} k_{2} \\phi u_{0}^{2}\\right)} B_{10}\\right)} B_{20}}{{\\left(2 \\, \\Gamma \\rho_{0} u_{0} + {\\left(\\Gamma^{2} - {\\left(\\Gamma^{3} - 2 \\, \\Gamma^{2} + \\Gamma\\right)} \\phi\\right)} u_{0}^{2} + \\rho_{0}^{2}\\right)} B_{10}^{4} + {\\left(2 \\, \\Gamma \\rho_{0} u_{0} + {\\left(\\Gamma^{2} - {\\left(\\Gamma^{3} - 2 \\, \\Gamma^{2} + \\Gamma\\right)} \\phi\\right)} u_{0}^{2} + \\rho_{0}^{2}\\right)} B_{20}^{4} + {\\left(3 \\, \\Gamma \\rho_{0}^{2} u_{0} + {\\left(3 \\, \\Gamma^{2} - {\\left(\\Gamma^{3} - 2 \\, \\Gamma^{2} + \\Gamma\\right)} \\phi\\right)} \\rho_{0} u_{0}^{2} + {\\left(\\Gamma^{3} - {\\left(\\Gamma^{4} - 2 \\, \\Gamma^{3} + \\Gamma^{2}\\right)} \\phi\\right)} u_{0}^{3} + \\rho_{0}^{3}\\right)} B_{10}^{2} + {\\left(3 \\, \\Gamma \\rho_{0}^{2} u_{0} + {\\left(3 \\, \\Gamma^{2} - {\\left(\\Gamma^{3} - 2 \\, \\Gamma^{2} + \\Gamma\\right)} \\phi\\right)} \\rho_{0} u_{0}^{2} + {\\left(\\Gamma^{3} - {\\left(\\Gamma^{4} - 2 \\, \\Gamma^{3} + \\Gamma^{2}\\right)} \\phi\\right)} u_{0}^{3} + 2 \\, {\\left(2 \\, \\Gamma \\rho_{0} u_{0} + {\\left(\\Gamma^{2} - {\\left(\\Gamma^{3} - 2 \\, \\Gamma^{2} + \\Gamma\\right)} \\phi\\right)} u_{0}^{2} + \\rho_{0}^{2}\\right)} B_{10}^{2} + \\rho_{0}^{3}\\right)} B_{20}^{2}} &amp; 0 &amp; 0 &amp; 0 &amp; \\frac{{\\left(i \\, \\Gamma k_{1} u_{0} + i \\, k_{1} \\rho_{0}\\right)} B_{10}^{3} + {\\left(i \\, \\Gamma k_{1} u_{0} + i \\, k_{1} \\rho_{0}\\right)} B_{10} B_{20}^{2} + {\\left(i \\, \\Gamma^{2} k_{1} u_{0}^{2} + 2 i \\, \\Gamma k_{1} \\rho_{0} u_{0} + i \\, k_{1} \\rho_{0}^{2}\\right)} B_{10} + {\\left(2 i \\, \\Gamma k_{2} \\rho_{0} u_{0} + i \\, k_{2} \\rho_{0}^{2} + {\\left(i \\, \\Gamma^{2} k_{2} + {\\left(-i \\, \\Gamma^{3} + 2 i \\, \\Gamma^{2} - i \\, \\Gamma\\right)} k_{2} \\phi\\right)} u_{0}^{2}\\right)} B_{20}}{3 \\, \\Gamma \\rho_{0}^{2} u_{0} + {\\left(3 \\, \\Gamma^{2} - {\\left(\\Gamma^{3} - 2 \\, \\Gamma^{2} + \\Gamma\\right)} \\phi\\right)} \\rho_{0} u_{0}^{2} + {\\left(\\Gamma^{3} - {\\left(\\Gamma^{4} - 2 \\, \\Gamma^{3} + \\Gamma^{2}\\right)} \\phi\\right)} u_{0}^{3} + {\\left(2 \\, \\Gamma \\rho_{0} u_{0} + {\\left(\\Gamma^{2} - {\\left(\\Gamma^{3} - 2 \\, \\Gamma^{2} + \\Gamma\\right)} \\phi\\right)} u_{0}^{2} + \\rho_{0}^{2}\\right)} B_{10}^{2} + {\\left(2 \\, \\Gamma \\rho_{0} u_{0} + {\\left(\\Gamma^{2} - {\\left(\\Gamma^{3} - 2 \\, \\Gamma^{2} + \\Gamma\\right)} \\phi\\right)} u_{0}^{2} + \\rho_{0}^{2}\\right)} B_{20}^{2} + \\rho_{0}^{3}} &amp; \\frac{{\\left(i \\, \\Gamma k_{2} u_{0} + i \\, k_{2} \\rho_{0}\\right)} B_{10}^{3} + {\\left(i \\, \\Gamma k_{2} u_{0} + i \\, k_{2} \\rho_{0}\\right)} B_{10} B_{20}^{2} + {\\left(i \\, \\Gamma^{2} k_{2} u_{0}^{2} + 2 i \\, \\Gamma k_{2} \\rho_{0} u_{0} + i \\, k_{2} \\rho_{0}^{2}\\right)} B_{10} + {\\left(-2 i \\, \\Gamma k_{1} \\rho_{0} u_{0} - i \\, k_{1} \\rho_{0}^{2} + {\\left(-i \\, \\Gamma^{2} k_{1} + {\\left(i \\, \\Gamma^{3} - 2 i \\, \\Gamma^{2} + i \\, \\Gamma\\right)} k_{1} \\phi\\right)} u_{0}^{2}\\right)} B_{20}}{3 \\, \\Gamma \\rho_{0}^{2} u_{0} + {\\left(3 \\, \\Gamma^{2} - {\\left(\\Gamma^{3} - 2 \\, \\Gamma^{2} + \\Gamma\\right)} \\phi\\right)} \\rho_{0} u_{0}^{2} + {\\left(\\Gamma^{3} - {\\left(\\Gamma^{4} - 2 \\, \\Gamma^{3} + \\Gamma^{2}\\right)} \\phi\\right)} u_{0}^{3} + {\\left(2 \\, \\Gamma \\rho_{0} u_{0} + {\\left(\\Gamma^{2} - {\\left(\\Gamma^{3} - 2 \\, \\Gamma^{2} + \\Gamma\\right)} \\phi\\right)} u_{0}^{2} + \\rho_{0}^{2}\\right)} B_{10}^{2} + {\\left(2 \\, \\Gamma \\rho_{0} u_{0} + {\\left(\\Gamma^{2} - {\\left(\\Gamma^{3} - 2 \\, \\Gamma^{2} + \\Gamma\\right)} \\phi\\right)} u_{0}^{2} + \\rho_{0}^{2}\\right)} B_{20}^{2} + \\rho_{0}^{3}} &amp; 0 &amp; \\frac{{\\left(\\Gamma u_{0} + \\rho_{0}\\right)} B_{10}}{{\\left(2 \\, \\Gamma \\rho_{0} \\tau u_{0} + {\\left(\\Gamma^{2} - {\\left(\\Gamma^{3} - 2 \\, \\Gamma^{2} + \\Gamma\\right)} \\phi\\right)} \\tau u_{0}^{2} + \\rho_{0}^{2} \\tau\\right)} \\sqrt{B_{10}^{2} + B_{20}^{2}}} &amp; -\\frac{{\\left(2 i \\, \\Gamma k_{1} u_{0} + 2 i \\, k_{1} \\rho_{0}\\right)} B_{10}^{4} + {\\left(2 i \\, \\Gamma k_{2} u_{0} + 2 i \\, k_{2} \\rho_{0}\\right)} B_{10} B_{20}^{3} + {\\left(2 i \\, \\Gamma^{2} k_{1} u_{0}^{2} + 4 i \\, \\Gamma k_{1} \\rho_{0} u_{0} + 2 i \\, k_{1} \\rho_{0}^{2}\\right)} B_{10}^{2} + {\\left(-2 i \\, \\Gamma k_{1} \\rho_{0} u_{0} + {\\left(2 i \\, \\Gamma k_{1} u_{0} + 2 i \\, k_{1} \\rho_{0}\\right)} B_{10}^{2} - i \\, k_{1} \\rho_{0}^{2} + {\\left(-i \\, \\Gamma^{2} k_{1} + {\\left(i \\, \\Gamma^{3} - 2 i \\, \\Gamma^{2} + i \\, \\Gamma\\right)} k_{1} \\phi\\right)} u_{0}^{2}\\right)} B_{20}^{2} + {\\left({\\left(2 i \\, \\Gamma k_{2} u_{0} + 2 i \\, k_{2} \\rho_{0}\\right)} B_{10}^{3} + {\\left(6 i \\, \\Gamma k_{2} \\rho_{0} u_{0} + 3 i \\, k_{2} \\rho_{0}^{2} + {\\left(3 i \\, \\Gamma^{2} k_{2} + {\\left(-i \\, \\Gamma^{3} + 2 i \\, \\Gamma^{2} - i \\, \\Gamma\\right)} k_{2} \\phi\\right)} u_{0}^{2}\\right)} B_{10}\\right)} B_{20}}{3 \\, {\\left({\\left(2 \\, \\Gamma \\rho_{0} u_{0} + {\\left(\\Gamma^{2} - {\\left(\\Gamma^{3} - 2 \\, \\Gamma^{2} + \\Gamma\\right)} \\phi\\right)} u_{0}^{2} + \\rho_{0}^{2}\\right)} B_{10}^{4} + {\\left(2 \\, \\Gamma \\rho_{0} u_{0} + {\\left(\\Gamma^{2} - {\\left(\\Gamma^{3} - 2 \\, \\Gamma^{2} + \\Gamma\\right)} \\phi\\right)} u_{0}^{2} + \\rho_{0}^{2}\\right)} B_{20}^{4} + {\\left(3 \\, \\Gamma \\rho_{0}^{2} u_{0} + {\\left(3 \\, \\Gamma^{2} - {\\left(\\Gamma^{3} - 2 \\, \\Gamma^{2} + \\Gamma\\right)} \\phi\\right)} \\rho_{0} u_{0}^{2} + {\\left(\\Gamma^{3} - {\\left(\\Gamma^{4} - 2 \\, \\Gamma^{3} + \\Gamma^{2}\\right)} \\phi\\right)} u_{0}^{3} + \\rho_{0}^{3}\\right)} B_{10}^{2} + {\\left(3 \\, \\Gamma \\rho_{0}^{2} u_{0} + {\\left(3 \\, \\Gamma^{2} - {\\left(\\Gamma^{3} - 2 \\, \\Gamma^{2} + \\Gamma\\right)} \\phi\\right)} \\rho_{0} u_{0}^{2} + {\\left(\\Gamma^{3} - {\\left(\\Gamma^{4} - 2 \\, \\Gamma^{3} + \\Gamma^{2}\\right)} \\phi\\right)} u_{0}^{3} + 2 \\, {\\left(2 \\, \\Gamma \\rho_{0} u_{0} + {\\left(\\Gamma^{2} - {\\left(\\Gamma^{3} - 2 \\, \\Gamma^{2} + \\Gamma\\right)} \\phi\\right)} u_{0}^{2} + \\rho_{0}^{2}\\right)} B_{10}^{2} + \\rho_{0}^{3}\\right)} B_{20}^{2}\\right)}} \\\\\n-\\frac{{\\left(i \\, \\Gamma^{3} - 2 i \\, \\Gamma^{2} + i \\, \\Gamma\\right)} B_{10} B_{20} k_{1} \\phi u_{0}^{2} + {\\left(i \\, \\Gamma^{3} - 2 i \\, \\Gamma^{2} + i \\, \\Gamma\\right)} B_{20}^{2} k_{2} \\phi u_{0}^{2}}{{\\left(2 \\, \\Gamma \\rho_{0}^{2} u_{0} + {\\left(\\Gamma^{2} - {\\left(\\Gamma^{3} - 2 \\, \\Gamma^{2} + \\Gamma\\right)} \\phi\\right)} \\rho_{0} u_{0}^{2} + \\rho_{0}^{3}\\right)} B_{10}^{2} + {\\left(2 \\, \\Gamma \\rho_{0}^{2} u_{0} + {\\left(\\Gamma^{2} - {\\left(\\Gamma^{3} - 2 \\, \\Gamma^{2} + \\Gamma\\right)} \\phi\\right)} \\rho_{0} u_{0}^{2} + \\rho_{0}^{3}\\right)} B_{20}^{2}} &amp; \\frac{{\\left({\\left(-i \\, \\Gamma + i\\right)} k_{1} \\rho_{0} + {\\left({\\left(i \\, \\Gamma^{3} - 2 i \\, \\Gamma^{2} + i \\, \\Gamma\\right)} k_{1} \\phi + {\\left(-i \\, \\Gamma^{2} + i \\, \\Gamma\\right)} k_{1}\\right)} u_{0}\\right)} B_{10} B_{20}^{3} + {\\left({\\left(-i \\, \\Gamma + i\\right)} k_{2} \\rho_{0} + {\\left({\\left(i \\, \\Gamma^{3} - 2 i \\, \\Gamma^{2} + i \\, \\Gamma\\right)} k_{2} \\phi + {\\left(-i \\, \\Gamma^{2} + i \\, \\Gamma\\right)} k_{2}\\right)} u_{0}\\right)} B_{20}^{4} + {\\left({\\left(-i \\, \\Gamma + i\\right)} k_{2} \\rho_{0}^{2} + {\\left(-2 i \\, \\Gamma^{2} + 2 i \\, \\Gamma\\right)} k_{2} \\rho_{0} u_{0} + {\\left({\\left(i \\, \\Gamma^{4} - 3 i \\, \\Gamma^{3} + 3 i \\, \\Gamma^{2} - i \\, \\Gamma\\right)} k_{2} \\phi + {\\left(-i \\, \\Gamma^{3} + i \\, \\Gamma^{2}\\right)} k_{2}\\right)} u_{0}^{2}\\right)} B_{10}^{2} + {\\left({\\left(-i \\, \\Gamma + i\\right)} k_{2} \\rho_{0}^{2} + {\\left({\\left(-i \\, \\Gamma + i\\right)} k_{2} \\rho_{0} + {\\left({\\left(i \\, \\Gamma^{3} - 2 i \\, \\Gamma^{2} + i \\, \\Gamma\\right)} k_{2} \\phi + {\\left(-i \\, \\Gamma^{2} + i \\, \\Gamma\\right)} k_{2}\\right)} u_{0}\\right)} B_{10}^{2} + {\\left({\\left(i \\, \\Gamma^{3} - 2 i \\, \\Gamma^{2} + i \\, \\Gamma\\right)} k_{2} \\phi + {\\left(-2 i \\, \\Gamma^{2} + 2 i \\, \\Gamma\\right)} k_{2}\\right)} \\rho_{0} u_{0} + {\\left({\\left(i \\, \\Gamma^{4} - 2 i \\, \\Gamma^{3} + i \\, \\Gamma^{2}\\right)} k_{2} \\phi + {\\left(-i \\, \\Gamma^{3} + i \\, \\Gamma^{2}\\right)} k_{2}\\right)} u_{0}^{2}\\right)} B_{20}^{2} + {\\left({\\left({\\left(-i \\, \\Gamma + i\\right)} k_{1} \\rho_{0} + {\\left({\\left(i \\, \\Gamma^{3} - 2 i \\, \\Gamma^{2} + i \\, \\Gamma\\right)} k_{1} \\phi + {\\left(-i \\, \\Gamma^{2} + i \\, \\Gamma\\right)} k_{1}\\right)} u_{0}\\right)} B_{10}^{3} + {\\left({\\left(i \\, \\Gamma^{3} - 2 i \\, \\Gamma^{2} + i \\, \\Gamma\\right)} k_{1} \\phi \\rho_{0} u_{0} + {\\left(i \\, \\Gamma^{3} - 2 i \\, \\Gamma^{2} + i \\, \\Gamma\\right)} k_{1} \\phi u_{0}^{2}\\right)} B_{10}\\right)} B_{20}}{{\\left(2 \\, \\Gamma \\rho_{0} u_{0} + {\\left(\\Gamma^{2} - {\\left(\\Gamma^{3} - 2 \\, \\Gamma^{2} + \\Gamma\\right)} \\phi\\right)} u_{0}^{2} + \\rho_{0}^{2}\\right)} B_{10}^{4} + {\\left(2 \\, \\Gamma \\rho_{0} u_{0} + {\\left(\\Gamma^{2} - {\\left(\\Gamma^{3} - 2 \\, \\Gamma^{2} + \\Gamma\\right)} \\phi\\right)} u_{0}^{2} + \\rho_{0}^{2}\\right)} B_{20}^{4} + {\\left(3 \\, \\Gamma \\rho_{0}^{2} u_{0} + {\\left(3 \\, \\Gamma^{2} - {\\left(\\Gamma^{3} - 2 \\, \\Gamma^{2} + \\Gamma\\right)} \\phi\\right)} \\rho_{0} u_{0}^{2} + {\\left(\\Gamma^{3} - {\\left(\\Gamma^{4} - 2 \\, \\Gamma^{3} + \\Gamma^{2}\\right)} \\phi\\right)} u_{0}^{3} + \\rho_{0}^{3}\\right)} B_{10}^{2} + {\\left(3 \\, \\Gamma \\rho_{0}^{2} u_{0} + {\\left(3 \\, \\Gamma^{2} - {\\left(\\Gamma^{3} - 2 \\, \\Gamma^{2} + \\Gamma\\right)} \\phi\\right)} \\rho_{0} u_{0}^{2} + {\\left(\\Gamma^{3} - {\\left(\\Gamma^{4} - 2 \\, \\Gamma^{3} + \\Gamma^{2}\\right)} \\phi\\right)} u_{0}^{3} + 2 \\, {\\left(2 \\, \\Gamma \\rho_{0} u_{0} + {\\left(\\Gamma^{2} - {\\left(\\Gamma^{3} - 2 \\, \\Gamma^{2} + \\Gamma\\right)} \\phi\\right)} u_{0}^{2} + \\rho_{0}^{2}\\right)} B_{10}^{2} + \\rho_{0}^{3}\\right)} B_{20}^{2}} &amp; 0 &amp; 0 &amp; 0 &amp; \\frac{{\\left(i \\, \\Gamma k_{1} u_{0} + i \\, k_{1} \\rho_{0}\\right)} B_{20}^{3} + {\\left(-2 i \\, \\Gamma k_{2} \\rho_{0} u_{0} - i \\, k_{2} \\rho_{0}^{2} + {\\left(-i \\, \\Gamma^{2} k_{2} + {\\left(i \\, \\Gamma^{3} - 2 i \\, \\Gamma^{2} + i \\, \\Gamma\\right)} k_{2} \\phi\\right)} u_{0}^{2}\\right)} B_{10} + {\\left(i \\, \\Gamma^{2} k_{1} u_{0}^{2} + 2 i \\, \\Gamma k_{1} \\rho_{0} u_{0} + {\\left(i \\, \\Gamma k_{1} u_{0} + i \\, k_{1} \\rho_{0}\\right)} B_{10}^{2} + i \\, k_{1} \\rho_{0}^{2}\\right)} B_{20}}{3 \\, \\Gamma \\rho_{0}^{2} u_{0} + {\\left(3 \\, \\Gamma^{2} - {\\left(\\Gamma^{3} - 2 \\, \\Gamma^{2} + \\Gamma\\right)} \\phi\\right)} \\rho_{0} u_{0}^{2} + {\\left(\\Gamma^{3} - {\\left(\\Gamma^{4} - 2 \\, \\Gamma^{3} + \\Gamma^{2}\\right)} \\phi\\right)} u_{0}^{3} + {\\left(2 \\, \\Gamma \\rho_{0} u_{0} + {\\left(\\Gamma^{2} - {\\left(\\Gamma^{3} - 2 \\, \\Gamma^{2} + \\Gamma\\right)} \\phi\\right)} u_{0}^{2} + \\rho_{0}^{2}\\right)} B_{10}^{2} + {\\left(2 \\, \\Gamma \\rho_{0} u_{0} + {\\left(\\Gamma^{2} - {\\left(\\Gamma^{3} - 2 \\, \\Gamma^{2} + \\Gamma\\right)} \\phi\\right)} u_{0}^{2} + \\rho_{0}^{2}\\right)} B_{20}^{2} + \\rho_{0}^{3}} &amp; \\frac{{\\left(i \\, \\Gamma k_{2} u_{0} + i \\, k_{2} \\rho_{0}\\right)} B_{20}^{3} + {\\left(2 i \\, \\Gamma k_{1} \\rho_{0} u_{0} + i \\, k_{1} \\rho_{0}^{2} + {\\left(i \\, \\Gamma^{2} k_{1} + {\\left(-i \\, \\Gamma^{3} + 2 i \\, \\Gamma^{2} - i \\, \\Gamma\\right)} k_{1} \\phi\\right)} u_{0}^{2}\\right)} B_{10} + {\\left(i \\, \\Gamma^{2} k_{2} u_{0}^{2} + 2 i \\, \\Gamma k_{2} \\rho_{0} u_{0} + {\\left(i \\, \\Gamma k_{2} u_{0} + i \\, k_{2} \\rho_{0}\\right)} B_{10}^{2} + i \\, k_{2} \\rho_{0}^{2}\\right)} B_{20}}{3 \\, \\Gamma \\rho_{0}^{2} u_{0} + {\\left(3 \\, \\Gamma^{2} - {\\left(\\Gamma^{3} - 2 \\, \\Gamma^{2} + \\Gamma\\right)} \\phi\\right)} \\rho_{0} u_{0}^{2} + {\\left(\\Gamma^{3} - {\\left(\\Gamma^{4} - 2 \\, \\Gamma^{3} + \\Gamma^{2}\\right)} \\phi\\right)} u_{0}^{3} + {\\left(2 \\, \\Gamma \\rho_{0} u_{0} + {\\left(\\Gamma^{2} - {\\left(\\Gamma^{3} - 2 \\, \\Gamma^{2} + \\Gamma\\right)} \\phi\\right)} u_{0}^{2} + \\rho_{0}^{2}\\right)} B_{10}^{2} + {\\left(2 \\, \\Gamma \\rho_{0} u_{0} + {\\left(\\Gamma^{2} - {\\left(\\Gamma^{3} - 2 \\, \\Gamma^{2} + \\Gamma\\right)} \\phi\\right)} u_{0}^{2} + \\rho_{0}^{2}\\right)} B_{20}^{2} + \\rho_{0}^{3}} &amp; 0 &amp; \\frac{{\\left(\\Gamma u_{0} + \\rho_{0}\\right)} B_{20}}{{\\left(2 \\, \\Gamma \\rho_{0} \\tau u_{0} + {\\left(\\Gamma^{2} - {\\left(\\Gamma^{3} - 2 \\, \\Gamma^{2} + \\Gamma\\right)} \\phi\\right)} \\tau u_{0}^{2} + \\rho_{0}^{2} \\tau\\right)} \\sqrt{B_{10}^{2} + B_{20}^{2}}} &amp; -\\frac{{\\left(2 i \\, \\Gamma k_{1} u_{0} + 2 i \\, k_{1} \\rho_{0}\\right)} B_{10} B_{20}^{3} + {\\left(2 i \\, \\Gamma k_{2} u_{0} + 2 i \\, k_{2} \\rho_{0}\\right)} B_{20}^{4} + {\\left(-2 i \\, \\Gamma k_{2} \\rho_{0} u_{0} - i \\, k_{2} \\rho_{0}^{2} + {\\left(-i \\, \\Gamma^{2} k_{2} + {\\left(i \\, \\Gamma^{3} - 2 i \\, \\Gamma^{2} + i \\, \\Gamma\\right)} k_{2} \\phi\\right)} u_{0}^{2}\\right)} B_{10}^{2} + {\\left(2 i \\, \\Gamma^{2} k_{2} u_{0}^{2} + 4 i \\, \\Gamma k_{2} \\rho_{0} u_{0} + {\\left(2 i \\, \\Gamma k_{2} u_{0} + 2 i \\, k_{2} \\rho_{0}\\right)} B_{10}^{2} + 2 i \\, k_{2} \\rho_{0}^{2}\\right)} B_{20}^{2} + {\\left({\\left(2 i \\, \\Gamma k_{1} u_{0} + 2 i \\, k_{1} \\rho_{0}\\right)} B_{10}^{3} + {\\left(6 i \\, \\Gamma k_{1} \\rho_{0} u_{0} + 3 i \\, k_{1} \\rho_{0}^{2} + {\\left(3 i \\, \\Gamma^{2} k_{1} + {\\left(-i \\, \\Gamma^{3} + 2 i \\, \\Gamma^{2} - i \\, \\Gamma\\right)} k_{1} \\phi\\right)} u_{0}^{2}\\right)} B_{10}\\right)} B_{20}}{3 \\, {\\left({\\left(2 \\, \\Gamma \\rho_{0} u_{0} + {\\left(\\Gamma^{2} - {\\left(\\Gamma^{3} - 2 \\, \\Gamma^{2} + \\Gamma\\right)} \\phi\\right)} u_{0}^{2} + \\rho_{0}^{2}\\right)} B_{10}^{4} + {\\left(2 \\, \\Gamma \\rho_{0} u_{0} + {\\left(\\Gamma^{2} - {\\left(\\Gamma^{3} - 2 \\, \\Gamma^{2} + \\Gamma\\right)} \\phi\\right)} u_{0}^{2} + \\rho_{0}^{2}\\right)} B_{20}^{4} + {\\left(3 \\, \\Gamma \\rho_{0}^{2} u_{0} + {\\left(3 \\, \\Gamma^{2} - {\\left(\\Gamma^{3} - 2 \\, \\Gamma^{2} + \\Gamma\\right)} \\phi\\right)} \\rho_{0} u_{0}^{2} + {\\left(\\Gamma^{3} - {\\left(\\Gamma^{4} - 2 \\, \\Gamma^{3} + \\Gamma^{2}\\right)} \\phi\\right)} u_{0}^{3} + \\rho_{0}^{3}\\right)} B_{10}^{2} + {\\left(3 \\, \\Gamma \\rho_{0}^{2} u_{0} + {\\left(3 \\, \\Gamma^{2} - {\\left(\\Gamma^{3} - 2 \\, \\Gamma^{2} + \\Gamma\\right)} \\phi\\right)} \\rho_{0} u_{0}^{2} + {\\left(\\Gamma^{3} - {\\left(\\Gamma^{4} - 2 \\, \\Gamma^{3} + \\Gamma^{2}\\right)} \\phi\\right)} u_{0}^{3} + 2 \\, {\\left(2 \\, \\Gamma \\rho_{0} u_{0} + {\\left(\\Gamma^{2} - {\\left(\\Gamma^{3} - 2 \\, \\Gamma^{2} + \\Gamma\\right)} \\phi\\right)} u_{0}^{2} + \\rho_{0}^{2}\\right)} B_{10}^{2} + \\rho_{0}^{3}\\right)} B_{20}^{2}\\right)}} \\\\\n0 &amp; 0 &amp; 0 &amp; 0 &amp; 0 &amp; 0 &amp; 0 &amp; \\frac{i \\, B_{10} k_{1} + i \\, B_{20} k_{2}}{B_{10}^{2} + B_{20}^{2} + \\Gamma u_{0} + \\rho_{0}} &amp; 0 &amp; 0 \\\\\n0 &amp; 0 &amp; i \\, B_{20} k_{2} &amp; -i \\, B_{10} k_{2} &amp; 0 &amp; 0 &amp; 0 &amp; 0 &amp; 0 &amp; 0 \\\\\n0 &amp; 0 &amp; -i \\, B_{20} k_{1} &amp; i \\, B_{10} k_{1} &amp; 0 &amp; 0 &amp; 0 &amp; 0 &amp; 0 &amp; 0 \\\\\n0 &amp; 0 &amp; 0 &amp; 0 &amp; i \\, B_{10} k_{1} + i \\, B_{20} k_{2} &amp; 0 &amp; 0 &amp; 0 &amp; 0 &amp; 0 \\\\\n\\frac{{\\left({\\left(i \\, \\Gamma^{3} - 2 i \\, \\Gamma^{2} + i \\, \\Gamma\\right)} k_{1} \\phi \\rho_{0} u_{0}^{2} + {\\left(i \\, \\Gamma^{4} - 2 i \\, \\Gamma^{3} + i \\, \\Gamma^{2}\\right)} k_{1} \\phi u_{0}^{3}\\right)} B_{10} + {\\left({\\left(i \\, \\Gamma^{3} - 2 i \\, \\Gamma^{2} + i \\, \\Gamma\\right)} k_{2} \\phi \\rho_{0} u_{0}^{2} + {\\left(i \\, \\Gamma^{4} - 2 i \\, \\Gamma^{3} + i \\, \\Gamma^{2}\\right)} k_{2} \\phi u_{0}^{3}\\right)} B_{20}}{{\\left(2 \\, \\Gamma \\rho_{0}^{2} u_{0} + {\\left(\\Gamma^{2} - {\\left(\\Gamma^{3} - 2 \\, \\Gamma^{2} + \\Gamma\\right)} \\phi\\right)} \\rho_{0} u_{0}^{2} + \\rho_{0}^{3}\\right)} \\sqrt{B_{10}^{2} + B_{20}^{2}}} &amp; -\\frac{{\\left({\\left(i \\, \\Gamma^{3} - 2 i \\, \\Gamma^{2} + i \\, \\Gamma\\right)} k_{1} \\phi \\rho_{0} u_{0} + {\\left(i \\, \\Gamma^{3} - 2 i \\, \\Gamma^{2} + i \\, \\Gamma\\right)} k_{1} \\phi u_{0}^{2}\\right)} B_{10} + {\\left({\\left(i \\, \\Gamma^{3} - 2 i \\, \\Gamma^{2} + i \\, \\Gamma\\right)} k_{2} \\phi \\rho_{0} u_{0} + {\\left(i \\, \\Gamma^{3} - 2 i \\, \\Gamma^{2} + i \\, \\Gamma\\right)} k_{2} \\phi u_{0}^{2}\\right)} B_{20}}{{\\left(2 \\, \\Gamma \\rho_{0} u_{0} + {\\left(\\Gamma^{2} - {\\left(\\Gamma^{3} - 2 \\, \\Gamma^{2} + \\Gamma\\right)} \\phi\\right)} u_{0}^{2} + \\rho_{0}^{2}\\right)} \\sqrt{B_{10}^{2} + B_{20}^{2}}} &amp; 0 &amp; 0 &amp; 0 &amp; -\\frac{{\\left(i \\, \\Gamma^{3} - 2 i \\, \\Gamma^{2} + i \\, \\Gamma\\right)} B_{10}^{2} k_{1} \\phi u_{0}^{2} + {\\left(i \\, \\Gamma^{3} - 2 i \\, \\Gamma^{2} + i \\, \\Gamma\\right)} B_{20}^{2} k_{1} \\phi u_{0}^{2}}{{\\left(2 \\, \\Gamma \\rho_{0} u_{0} + {\\left(\\Gamma^{2} - {\\left(\\Gamma^{3} - 2 \\, \\Gamma^{2} + \\Gamma\\right)} \\phi\\right)} u_{0}^{2} + \\rho_{0}^{2}\\right)} \\sqrt{B_{10}^{2} + B_{20}^{2}}} &amp; -\\frac{{\\left(i \\, \\Gamma^{3} - 2 i \\, \\Gamma^{2} + i \\, \\Gamma\\right)} B_{10}^{2} k_{2} \\phi u_{0}^{2} + {\\left(i \\, \\Gamma^{3} - 2 i \\, \\Gamma^{2} + i \\, \\Gamma\\right)} B_{20}^{2} k_{2} \\phi u_{0}^{2}}{{\\left(2 \\, \\Gamma \\rho_{0} u_{0} + {\\left(\\Gamma^{2} - {\\left(\\Gamma^{3} - 2 \\, \\Gamma^{2} + \\Gamma\\right)} \\phi\\right)} u_{0}^{2} + \\rho_{0}^{2}\\right)} \\sqrt{B_{10}^{2} + B_{20}^{2}}} &amp; 0 &amp; -\\frac{\\Gamma^{2} u_{0}^{2} + 2 \\, \\Gamma \\rho_{0} u_{0} + \\rho_{0}^{2}}{2 \\, \\Gamma \\rho_{0} \\tau u_{0} + {\\left(\\Gamma^{2} - {\\left(\\Gamma^{3} - 2 \\, \\Gamma^{2} + \\Gamma\\right)} \\phi\\right)} \\tau u_{0}^{2} + \\rho_{0}^{2} \\tau} &amp; \\frac{{\\left(2 i \\, \\Gamma^{3} - 4 i \\, \\Gamma^{2} + 2 i \\, \\Gamma\\right)} B_{10} k_{1} \\phi u_{0}^{2} + {\\left(2 i \\, \\Gamma^{3} - 4 i \\, \\Gamma^{2} + 2 i \\, \\Gamma\\right)} B_{20} k_{2} \\phi u_{0}^{2}}{3 \\, {\\left(2 \\, \\Gamma \\rho_{0} u_{0} + {\\left(\\Gamma^{2} - {\\left(\\Gamma^{3} - 2 \\, \\Gamma^{2} + \\Gamma\\right)} \\phi\\right)} u_{0}^{2} + \\rho_{0}^{2}\\right)} \\sqrt{B_{10}^{2} + B_{20}^{2}}} \\\\\n0 &amp; 0 &amp; \\frac{{\\left(-2 i \\, \\Gamma^{2} + 2 i \\, \\Gamma\\right)} B_{10}^{2} k_{1} \\psi \\rho_{0} u_{0} + {\\left(i \\, \\Gamma^{2} - i \\, \\Gamma\\right)} B_{20}^{2} k_{1} \\psi \\rho_{0} u_{0} + {\\left(-3 i \\, \\Gamma^{2} + 3 i \\, \\Gamma\\right)} B_{10} B_{20} k_{2} \\psi \\rho_{0} u_{0}}{{\\left(\\Gamma u_{0} + \\rho_{0}\\right)} B_{10}^{2} + {\\left(\\Gamma u_{0} + \\rho_{0}\\right)} B_{20}^{2}} &amp; -\\frac{{\\left(3 i \\, \\Gamma^{2} - 3 i \\, \\Gamma\\right)} B_{10} B_{20} k_{1} \\psi \\rho_{0} u_{0} + {\\left(-i \\, \\Gamma^{2} + i \\, \\Gamma\\right)} B_{10}^{2} k_{2} \\psi \\rho_{0} u_{0} + {\\left(2 i \\, \\Gamma^{2} - 2 i \\, \\Gamma\\right)} B_{20}^{2} k_{2} \\psi \\rho_{0} u_{0}}{{\\left(\\Gamma u_{0} + \\rho_{0}\\right)} B_{10}^{2} + {\\left(\\Gamma u_{0} + \\rho_{0}\\right)} B_{20}^{2}} &amp; 0 &amp; 0 &amp; 0 &amp; 0 &amp; 0 &amp; -\\frac{1}{\\tau}\n\\end{array}\\right)$ $\\displaystyle \\left(\\begin{array}{r}\n\\delta_{\\rho} \\\\\n\\delta_{u} \\\\\n\\delta_{u_{1}} \\\\\n\\delta_{u_{2}} \\\\\n\\delta_{u_{3}} \\\\\n\\delta_{B_{1}} \\\\\n\\delta_{B_{2}} \\\\\n\\delta_{B_{3}} \\\\\n\\delta_{q} \\\\\n\\delta_{\\mathit{dp}}\n\\end{array}\\right)$</div>"}︡{"stdout":"\n\n\n"}︡{"html":"<div align='center'>Analytic eigenvalues and eigenvectors in the $k_1, k_2 \\rightarrow 0$ limit : </div>"}︡{"tex":{"tex":"\\left[\\left(-\\frac{\\Gamma^{2} u_{0}^{2} + 2 \\, \\Gamma \\rho_{0} u_{0} + \\rho_{0}^{2}}{2 \\, \\Gamma \\rho_{0} \\tau u_{0} + {\\left(\\Gamma^{2} - {\\left(\\Gamma^{3} - 2 \\, \\Gamma^{2} + \\Gamma\\right)} \\phi\\right)} \\tau u_{0}^{2} + \\rho_{0}^{2} \\tau}, \\left[\\left(0,\\,0,\\,1,\\,\\frac{B_{20}}{B_{10}},\\,0,\\,0,\\,0,\\,0,\\,-\\frac{\\sqrt{B_{10}^{2} + B_{20}^{2}} {\\left(\\Gamma u_{0} + \\rho_{0}\\right)}}{B_{10}},\\,0\\right)\\right], 1\\right), \\left(-\\frac{1}{\\tau}, \\left[\\left(0,\\,0,\\,0,\\,0,\\,0,\\,0,\\,0,\\,0,\\,0,\\,1\\right)\\right], 1\\right), \\left(0, \\left[\\left(1,\\,0,\\,0,\\,0,\\,0,\\,0,\\,0,\\,0,\\,0,\\,0\\right), \\left(0,\\,1,\\,0,\\,0,\\,0,\\,0,\\,0,\\,0,\\,0,\\,0\\right), \\left(0,\\,0,\\,1,\\,0,\\,0,\\,0,\\,0,\\,0,\\,0,\\,0\\right), \\left(0,\\,0,\\,0,\\,1,\\,0,\\,0,\\,0,\\,0,\\,0,\\,0\\right), \\left(0,\\,0,\\,0,\\,0,\\,1,\\,0,\\,0,\\,0,\\,0,\\,0\\right), \\left(0,\\,0,\\,0,\\,0,\\,0,\\,1,\\,0,\\,0,\\,0,\\,0\\right), \\left(0,\\,0,\\,0,\\,0,\\,0,\\,0,\\,1,\\,0,\\,0,\\,0\\right), \\left(0,\\,0,\\,0,\\,0,\\,0,\\,0,\\,0,\\,1,\\,0,\\,0\\right)\\right], 8\\right)\\right]","display":true}}︡{"html":"<div align='center'>Numerical eigenvalues and eigenvectors for $k_1 = $ $\\displaystyle 2 \\, \\pi$  , $k_2$ =  $\\displaystyle 4 \\, \\pi$ :\n</div>"}︡{"stdout":"--------------------------\n('Eigenvalue   = ', -0.12805439673099883 + 8.014946591406833*I)\n(delta_rho, ' = ', 0.3393169459380907 + 0.0008079151965263154*I)\n(delta_u, ' = ', 0.9063481051620624)\n(delta_u1, ' = ', -0.07794513888657582 - 0.00039012826680087653*I)\n(delta_u2, ' = ', -0.1774386622175456 - 0.0037779544489260283*I)\n(delta_u3, ' = ', -5.143561037828635e-19 - 2.7330171274215138e-17*I)\n(delta_B1, ' = ', -0.008833474513140352 + 0.0005499642274790809*I)\n(delta_B2, ' = ', 0.004416737256570164 - 0.00027498211373948496*I)\n(delta_B3, ' = ', -3.0941861711314e-18 - 1.6892063388076385e-17*I)\n(delta_q, ' = ', -0.0008859147076993578 + 0.001227692746437654*I)\n(delta_dp, ' = ', 0.15906832627246215 + 0.02047390198972448*I)\n--------------------------\n('Eigenvalue   = ', -0.12805439673099872 - 8.014946591406837*I)\n(delta_rho, ' = ', 0.33931694593809103 - 0.0008079151965255937*I)\n(delta_u, ' = ', 0.9063481051620623)\n(delta_u1, ' = ', 0.07794513888657595 - 0.0003901282668009043*I)\n(delta_u2, ' = ', 0.17743866221754587 - 0.0037779544489259936*I)\n(delta_u3, ' = ', -7.921299858524113e-18 + 2.842784923556555e-17*I)\n(delta_B1, ' = ', -0.008833474513140399 - 0.0005499642274790781*I)\n(delta_B2, ' = ', 0.004416737256570194 + 0.00027498211373948615*I)\n(delta_B3, ' = ', 6.1428095550824016e-18 - 1.9624287390581625e-17*I)\n(delta_q, ' = ', 0.0008859147076991034 + 0.001227692746437668*I)\n(delta_dp, ' = ', 0.15906832627246198 - 0.020473901989724153*I)\n--------------------------\n('Eigenvalue   = ', -0.5533585207638139 + 3.626257128688842*I)\n(delta_rho, ' = ', -0.5185225240822474 + 0.1792647678001852*I)\n(delta_u, ' = ', 0.551617073639381)\n(delta_u1, ' = ', -0.008463122479547778 + 0.011862022608466305*I)\n(delta_u2, ' = ', 0.1617546637187075 - 0.0348280808236026*I)\n(delta_u3, ' = ', 9.128719756426342e-18 - 1.590628273090525e-17*I)\n(delta_B1, ' = ', -0.05973794979640739 + 0.033517075061509215*I)\n(delta_B2, ' = ', 0.029868974898203685 - 0.016758537530754604*I)\n(delta_B3, ' = ', 1.0183059831665314e-17 - 2.401795841471482e-17*I)\n(delta_q, ' = ', -0.5233486841539436 + 0.04767672501939452*I)\n(delta_dp, ' = ', -0.2909106062057661 + 0.021594520553364402*I)\n--------------------------\n('Eigenvalue   = ', -0.5533585207638108 - 3.626257128688849*I)\n(delta_rho, ' = ', -0.5185225240822464 - 0.1792647678001874*I)\n(delta_u, ' = ', 0.551617073639382)\n(delta_u1, ' = ', 0.008463122479547853 + 0.011862022608466373*I)\n(delta_u2, ' = ', -0.16175466371870748 - 0.034828080823603495*I)\n(delta_u3, ' = ', 1.1150888301102058e-17 + 2.3900040672376563e-17*I)\n(delta_B1, ' = ', -0.059737949796407556 - 0.0335170750615094*I)\n(delta_B2, ' = ', 0.029868974898203757 + 0.01675853753075467*I)\n(delta_B3, ' = ', -1.317599308756542e-17 - 3.254129323467699e-17*I)\n(delta_q, ' = ', 0.5233486841539429 + 0.04767672501939605*I)\n(delta_dp, ' = ', -0.2909106062057659 - 0.021594520553365606*I)\n--------------------------\n('Eigenvalue   = ', -0.032582265861863496 + 2.374578675580773*I)\n(delta_rho, ' = ', -0.1655542903438733 - 0.04080602067344037*I)\n(delta_u, ' = ', -0.2504999182918082 - 0.11191260373025272*I)\n(delta_u1, ' = ', 0.40459909465918575 + 0.007755126791439568*I)\n(delta_u2, ' = ', -0.171121719863687 + 0.004262514878919793*I)\n(delta_u3, ' = ', -9.645045787414414e-16 + 1.4828604899990338e-15*I)\n(delta_B1, ' = ', 0.7329049844550491)\n(delta_B2, ' = ', -0.3664524922275243 + 4.163336342344337e-17*I)\n(delta_B3, ' = ', -1.804783461299747e-15 + 2.7917305065258617e-15*I)\n(delta_q, ' = ', -0.0326128713301966 + 8.128139992545115e-05*I)\n(delta_dp, ' = ', 0.16598310644864028 + 0.05039326687801304*I)\n--------------------------\n('Eigenvalue   = ', -0.03258226586186431 - 2.3745786755807763*I)\n(delta_rho, ' = ', -0.16555429034387398 + 0.04080602067343827*I)\n(delta_u, ' = ', -0.2504999182918073 + 0.11191260373025616*I)\n(delta_u1, ' = ', -0.40459909465918553 + 0.007755126791439585*I)\n(delta_u2, ' = ', 0.1711217198636864 + 0.004262514878919727*I)\n(delta_u3, ' = ', -1.5054087107539244e-15 - 1.837173724820065e-15*I)\n(delta_B1, ' = ', 0.7329049844550489)\n(delta_B2, ' = ', -0.3664524922275243 - 2.949029909160572e-16*I)\n"}︡{"stdout":"(delta_B3, ' = ', 2.8068396606476423e-15 + 3.4311812625229985e-15*I)\n(delta_q, ' = ', 0.03261287133019803 + 8.128139992654576e-05*I)\n(delta_dp, ' = ', 0.16598310644863995 - 0.050393266878013786*I)\n--------------------------\n('Eigenvalue   = ', -0.6181191433730994 + 1.0433970309024512e-16*I)\n(delta_rho, ' = ', 0.4040113802128661 + 8.847089727481716e-16*I)\n(delta_u, ' = ', 0.8019719572160686)\n(delta_u1, ' = ', -7.714098457234364e-17 + 0.015408497913300819*I)\n(delta_u2, ' = ', 9.085614205428527e-17 - 0.027576905582327595*I)\n(delta_u3, ' = ', 1.084506875678059e-17 - 3.0538732062589785e-18*I)\n(delta_B1, ' = ', 0.15004037846992804 + 2.2985086056692694e-16*I)\n(delta_B2, ' = ', -0.07502018923496388 - 1.5851035761738075e-16*I)\n(delta_B3, ' = ', -9.97834766818435e-18 - 3.26741761286195e-17*I)\n(delta_q, ' = ', -3.217912047936977e-16 + 0.012238978432566371*I)\n(delta_dp, ' = ', -0.4053738279795086 - 3.5908775952719907e-16*I)\n--------------------------\n('Eigenvalue   = ', -4.440892098448468e-16 + 2.266205629205168*I)\n(delta_rho, ' = ', 6.246213773297642e-16 - 5.297411822987512e-16*I)\n(delta_u, ' = ', 1.0681167110669678e-15 - 6.277665219964732e-16*I)\n(delta_u1, ' = ', -7.939580692067271e-16 + 1.3821827675399122e-15*I)\n(delta_u2, ' = ', 2.503468520328191e-16 - 5.625194608100289e-16*I)\n(delta_u3, ' = ', 0.4580286124143448 - 2.498001805406602e-16*I)\n(delta_B1, ' = ', -1.0786127273034217e-15 + 2.2956684523441033e-15*I)\n(delta_B2, ' = ', 7.63222135775721e-16 - 1.3344795443535536e-15*I)\n(delta_B3, ' = ', 0.888937450110968)\n(delta_q, ' = ', -2.330240761856706e-17 - 3.810298757137183e-17*I)\n(delta_dp, ' = ', -5.038916700638819e-16 + 5.457522786131321e-16*I)\n--------------------------\n('Eigenvalue   = ', 9.27665875131909e-16 - 2.2662056292051673*I)\n(delta_rho, ' = ', -4.713117917886574e-16 - 4.183843249618291e-16*I)\n(delta_u, ' = ', -1.0831458726035688e-15 - 8.901442254634802e-16*I)\n(delta_u1, ' = ', -4.729873686298035e-16 - 1.3819627791887942e-15*I)\n(delta_u2, ' = ', 2.032172998273417e-16 + 5.797717646254965e-16*I)\n(delta_u3, ' = ', -0.45802861241434456 + 1.4224732503009818e-16*I)\n(delta_B1, ' = ', 6.070822727255268e-16 + 2.3389917136370904e-15*I)\n(delta_B2, ' = ', -5.134648960301807e-16 - 1.3437802188853436e-15*I)\n(delta_B3, ' = ', 0.8889374501109683)\n(delta_q, ' = ', -3.544040430521065e-17 - 4.557627989462774e-17*I)\n(delta_dp, ' = ', 4.1278407879801587e-16 + 6.359566503676294e-16*I)\n--------------------------\n('Eigenvalue   = ', 3.312381435919343e-16 + 1.77799477463703e-16*I)\n(delta_rho, ' = ', 0.3253956867279845 + 8.881784197001252e-16*I)\n(delta_u, ' = ', 0.6507913734559689)\n(delta_u1, ' = ', 1.3532792389002784e-16 + 2.0119301408311264e-16*I)\n(delta_u2, ' = ', 1.8551698744061074e-16 + 9.444958137599478e-17*I)\n(delta_u3, ' = ', -6.112978052977572e-17 - 8.235079128524141e-17*I)\n(delta_B1, ' = ', 0.39047482407358103 + 2.7755575615628914e-16*I)\n(delta_B2, ' = ', 0.5640191903285056 + 1.6653345369377348e-16*I)\n(delta_B3, ' = ', -1.4622897949435524e-16 - 9.081961660846472e-17*I)\n(delta_q, ' = ', -4.0351114185470864e-16 - 1.414148178047857e-16*I)\n(delta_dp, ' = ', -6.596007078225568e-17 - 1.5658553488433385e-16*I)\n"}︡
︠d39e6e67-8c79-4c5a-8606-03aa2a4281dc︠
eigenvecs_right = M.subs(k1=0, k2=0).eigenvectors_right()
︡0a9fda00-ebd4-4321-8077-7f04bcd70229︡
︠08755446-1f48-4d5e-8f64-a920ea7f0f73︠
solve(1/eigenvecs_right[0][0] == 0, phi)
︡99da33df-3840-46ac-808a-320c60e97e3a︡{"tex":{"tex":"\\left[\\phi = \\frac{\\Gamma^{2} u_{0}^{2} + 2 \\, \\Gamma \\rho_{0} u_{0} + \\rho_{0}^{2}}{{\\left(\\Gamma^{3} - 2 \\, \\Gamma^{2} + \\Gamma\\right)} u_{0}^{2}}\\right]","display":true}}︡
︠e4c88ae8-dd4b-4935-8fd5-8037e5091dac︠









