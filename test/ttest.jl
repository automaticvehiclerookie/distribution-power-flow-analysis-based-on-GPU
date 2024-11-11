using ModelingToolkit, NonlinearSolve
const baseMVA = 100.0
@variables θn_ac2,Vn_ac2,Vn_dc1,Vn_dc2
@parameters G11,B11,G12,B12,G21,B21,G22,B22
@parameters Pac2,Qac2
@parameters Gdc11,Gdc12,Gdc21,Gdc22,a,b,c,Pdc1,Pdc2
# Ps1=-60.0;Ps2=0.0;
# Qs1=-40.0;Qs2=0.0;

# Gt1=0.1193;Bt1=-8.919;
# Gt2=0.1193;Bt2=-8.919;
# bf1=0.0887;bf2=0.0887;
# Gc1=0.0037;Bc1=-6.0872;
# Gc2=0.0037;Bc2=-6.0872;
# Vn_ac1=1.06;θn_ac1=0.0;
# Vn_ac3=1.0;

# Gnm = [
#         6.25 -5 -1.25 0.0 0.0 ;
#         -5 10.8333 -1.66667 -1.66667 -2.5 ;
#         -1.25 -1.66667 12.9167 -10 0.0 ;
#         0.0 -1.66667 -10 12.9167 -1.25 ;
#         0.0 -2.5 0.0 -1.25 3.75 
#     ]
# Bnm = [
#         -18.695 15.0 3.75 0.0 0.0;
#         15.0 -32.415 5 5 7.5;
#         3.75 5 -38.695 30.0 0.0 ;
#         0.0 5 30.0 -38.695 3.75 ;
#         0.0 7.5 0.0 3.75 -11.21 
#     ]
# Pac2=0.2;Pac3=-0.45;Pac4=-0.4;Pac5=-0.6;
# Qac2=-0.1;Qac3=-0.15;Qac4=-0.05;Qac5=-0.1;
# Pc=Uc1^2 *Gc1-Uc1*Uf1*(Gc1*cos(θc1-θf1)+Bc1*sin(θc1-θf1))
# Qc=-Uc1^2 *Bc1-Uc1*Uf1*(Gc1*sin(θc1-θf1)-Bc1*cos(θc1-θf1))
# Ic=sqrt((Uc1^2 *Gc-Uc1*Uf1*(Gc*cos(θc1-θf1)+Bc*sin(θc1-θf1)))^2+(-Uc1^2 *Bc-Uc1*Uf1*(Gc*sin(θc1-θf1)-Bc*cos(θc1-θf1)))^2)/Vn_dc1
eqs4=[
    Vn_ac2*Vn_ac1*(G21*cos(θn_ac2-θn_ac1)+B21*sin(θn_ac2-θn_ac1))+Vn_ac2*Vn_ac2*(G22*cos(θn_ac2-θn_ac2)+B22*sin(θn_ac2-θn_ac2))+Pac2 ~ 0,
    Vn_ac2*Vn_ac1*(G21*sin(θn_ac2-θn_ac1)-B21*cos(θn_ac2-θn_ac1))+Vn_ac2*Vn_ac2*(G22*sin(θn_ac2-θn_ac2)+B22*cos(θn_ac2-θn_ac2))+Qac2 ~ 0,
]
# Gdc = [
#         19.2308 -19.2308;
#          -19.2308 19.2308
#     ]


eqs=vcat(eqs1,eqs2,eqs3,eqs4,eqs5)
parameters=[Gdc11,Gdc12,Gdc21,Gdc22,a,b,c,Pdc1,Pdc2,Gt1,Bt1,bf1,Gc1,Bc1,Vn_ac1,θn_ac1,G11,B11,G12,B12,G21,B21,G22,B22,Pac2,Qac2,Ps2,Qs2,Pdc1,Pdc2]
@mtkbuild ns = NonlinearSystem(eqs, [θn_ac2,Vn_ac2,Vn_dc1,Vn_dc2,θf1,Uf1,θc1,Uc1,M1], parameters)
u0=[
    θn_ac2 => 0.0,
    Vn_ac2 => 1.0,
    Vn_dc1 => 1.0,
    Vn_dc2 => 0.9,
    θf1 => 0.0,
    Uf1 => 1.0,
    θc1 => 0.0,
    Uc1 => 0.9,
    M1 => 1.0,
]

ps=[
    Ps2 => -40/baseMVA,
    Qs2 => -60/baseMVA,
    Gt1 => 0.1193,
    Bt1 => -8.919,
    bf1 => 0.0887,
    Gc1 => 0.0037,
    Bc1 => -6.0872,
    Vn_ac1 => 1.06,
    θn_ac1 => 0.0,
    Pac2 => 0.2,
    Qac2 => 0.1,
    Gdc11 => 19.2308,
    Gdc12 => -19.2308,
    Gdc21 => -19.2308,
    Gdc22 => 19.2308,
    G11 => 5.0,
    B11 => -14.97,
    G12 => -5.0,
    B12 => 15.0,
    G21 => -5.0,
    B21 => 15.0,
    G22 => 5.0,
    B22 => -14.97,
    a => 1.103,
    b => 0.887,
    c => 2.885,
    Pdc1 => 10/baseMVA,
    Pdc2 => 10/baseMVA,
]
prob = NonlinearProblem(ns, u0, ps,jac = true)
sol = solve(prob, NewtonRaphson())
