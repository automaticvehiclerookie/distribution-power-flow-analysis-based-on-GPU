using ModelingToolkit, NonlinearSolve
const baseMVA = 100.0
@variables θn_ac2,Vn_ac2,θn_ac4,Vn_ac4,θn_ac5,Vn_ac5,θn_ac3,Vn_dc1,Vn_dc2,θn_f1,Vn_f1,M1,Vn_c1,θn_c1,θn_f2,M2,Vn_c2,θn_c2
@parameters Vn_ac1,θn_ac1,Vn_ac3,Vn_f2
@parameters G11,B11,G12,B12,G13,B13,G14,B14,G15,B15
@parameters G21,B21,G22,B22,G23,B23,G24,B24,G25,B25
@parameters G31,B31,G32,B32,G33,B33,G34,B34,G35,B35
@parameters G41,B41,G42,B42,G43,B43,G44,B44,G45,B45
@parameters G51,B51,G52,B52,G53,B53,G54,B54,G55,B55
@parameters Pac2,Qac2,Pac3,Qac3,Pac4,Qac4,Pac5,Qac5,Pdc1,Pdc2
@parameters Gdc11,Gdc12,Gdc21,Gdc22,eta,tanphi,bf,Btf,Gtf,Ps1,Qs1,Gc,Bc,Ps2

eqs3=[
    -Vn_ac2^2*Gtf+Vn_ac2*Vn_f1*(Gtf*cos(θn_ac2-θn_f1)+Btf*sin(θn_ac2-θn_f1))-Ps1 ~ 0,
    Vn_ac2^2*Btf+Vn_ac2*Vn_f1*(Gtf*sin(θn_ac2-θn_f1)-Btf*cos(θn_ac2-θn_f1))-Qs1 ~ 0,
    -Vn_ac3^2*Gtf+Vn_ac3*Vn_f2*(Gtf*cos(θn_ac3-θn_f2)+Btf*sin(θn_ac3-θn_f2))-Ps2 ~ 0,
]
eqs4=[
    Vn_ac2*Vn_ac1*(G21*cos(θn_ac2-θn_ac1)+B21*sin(θn_ac2-θn_ac1))+Vn_ac2*Vn_ac2*(G22*cos(θn_ac2-θn_ac2)+B22*sin(θn_ac2-θn_ac2))+Vn_ac2*Vn_ac3*(G23*cos(θn_ac2-θn_ac3)+B23*sin(θn_ac2-θn_ac3))+Vn_ac2*Vn_ac4*(G24*cos(θn_ac2-θn_ac4)+B24*sin(θn_ac2-θn_ac4))+Vn_ac2*Vn_ac5*(G25*cos(θn_ac2-θn_ac5)+B25*sin(θn_ac2-θn_ac5))+Pac2-Ps1 ~ 0,
    Vn_ac2*Vn_ac1*(G21*sin(θn_ac2-θn_ac1)-B21*cos(θn_ac2-θn_ac1))+Vn_ac2*Vn_ac2*(G22*sin(θn_ac2-θn_ac2)-B22*cos(θn_ac2-θn_ac2))+Vn_ac2*Vn_ac3*(G23*sin(θn_ac2-θn_ac3)-B23*cos(θn_ac2-θn_ac3))+Vn_ac2*Vn_ac4*(G24*sin(θn_ac2-θn_ac4)-B24*cos(θn_ac2-θn_ac4))+Vn_ac2*Vn_ac5*(G25*sin(θn_ac2-θn_ac5)-B25*cos(θn_ac2-θn_ac5))+Qac2-Qs1 ~ 0,
    Vn_ac3*Vn_ac1*(G31*cos(θn_ac3-θn_ac1)+B31*sin(θn_ac3-θn_ac1))+Vn_ac3*Vn_ac2*(G32*cos(θn_ac3-θn_ac2)+B32*sin(θn_ac3-θn_ac2))+Vn_ac3*Vn_ac3*(G33*cos(θn_ac3-θn_ac3)+B33*sin(θn_ac3-θn_ac3))+Vn_ac3*Vn_ac4*(G34*cos(θn_ac3-θn_ac4)+B34*sin(θn_ac3-θn_ac4))+Vn_ac3*Vn_ac5*(G35*cos(θn_ac3-θn_ac5)+B35*sin(θn_ac3-θn_ac5))+Pac3-Ps2 ~ 0,
    Vn_ac4*Vn_ac1*(G41*cos(θn_ac4-θn_ac1)+B41*sin(θn_ac4-θn_ac1))+Vn_ac4*Vn_ac2*(G42*cos(θn_ac4-θn_ac2)+B42*sin(θn_ac4-θn_ac2))+Vn_ac4*Vn_ac3*(G43*cos(θn_ac4-θn_ac3)+B43*sin(θn_ac4-θn_ac3))+Vn_ac4*Vn_ac4*(G44*cos(θn_ac4-θn_ac4)+B44*sin(θn_ac4-θn_ac4))+Vn_ac4*Vn_ac5*(G45*cos(θn_ac4-θn_ac5)+B45*sin(θn_ac4-θn_ac5))+Pac4 ~ 0,
    Vn_ac4*Vn_ac1*(G41*sin(θn_ac4-θn_ac1)-B41*cos(θn_ac4-θn_ac1))+Vn_ac4*Vn_ac2*(G42*sin(θn_ac4-θn_ac2)-B42*cos(θn_ac4-θn_ac2))+Vn_ac4*Vn_ac3*(G43*sin(θn_ac4-θn_ac3)-B43*cos(θn_ac4-θn_ac3))+Vn_ac4*Vn_ac4*(G44*sin(θn_ac4-θn_ac4)-B44*cos(θn_ac4-θn_ac4))+Vn_ac4*Vn_ac5*(G45*sin(θn_ac4-θn_ac5)-B45*cos(θn_ac4-θn_ac5))+Qac4 ~ 0,
    Vn_ac5*Vn_ac1*(G51*cos(θn_ac5-θn_ac1)+B51*sin(θn_ac5-θn_ac1))+Vn_ac5*Vn_ac2*(G52*cos(θn_ac5-θn_ac2)+B52*sin(θn_ac5-θn_ac2))+Vn_ac5*Vn_ac3*(G53*cos(θn_ac5-θn_ac3)+B53*sin(θn_ac5-θn_ac3))+Vn_ac5*Vn_ac4*(G54*cos(θn_ac5-θn_ac4)+B54*sin(θn_ac5-θn_ac4))+Vn_ac5*Vn_ac5*(G55*cos(θn_ac5-θn_ac5)+B55*sin(θn_ac5-θn_ac5))+Pac5 ~ 0,
    Vn_ac5*Vn_ac1*(G51*sin(θn_ac5-θn_ac1)-B51*cos(θn_ac5-θn_ac1))+Vn_ac5*Vn_ac2*(G52*sin(θn_ac5-θn_ac2)-B52*cos(θn_ac5-θn_ac2))+Vn_ac5*Vn_ac3*(G53*sin(θn_ac5-θn_ac3)-B53*cos(θn_ac5-θn_ac3))+Vn_ac5*Vn_ac4*(G54*sin(θn_ac5-θn_ac4)-B54*cos(θn_ac5-θn_ac4))+Vn_ac5*Vn_ac5*(G55*sin(θn_ac5-θn_ac5)-B55*cos(θn_ac5-θn_ac5))+Qac5 ~ 0,
]

eqs5=[
    Vn_dc1*Vn_dc1*Gdc11+Vn_dc1*Vn_dc2*Gdc12+Pdc1+eta*(Vn_c1^2*Gc-Vn_c1*Vn_f1*(Gc*cos(θn_c1-θn_f1)+Bc*sin(θn_c1-θn_f1))) ~ 0,
    Vn_dc2*Vn_dc1*Gdc21+Vn_dc2*Vn_dc2*Gdc22+Pdc2+eta*(Vn_c2^2*Gc-Vn_c2*Vn_f2*(Gc*cos(θn_c2-θn_f2)+Bc*sin(θn_c2-θn_f2))) ~ 0,
    -Vn_f1^2*Btf-Vn_f1*Vn_ac2*(Gtf*sin(θn_f1-θn_ac2)-Btf*cos(θn_f1-θn_ac2))-Vn_f1^2*bf-Vn_f1^2*Bc-Vn_f1*Vn_c1*(Gc*sin(θn_f1-θn_c1)-Bc*cos(θn_f1-θn_c1)) ~ 0,
    Vn_f1^2*Gc-Vn_f1*Vn_c1*(Gc*cos(θn_f1-θn_c1)+Bc*sin(θn_f1-θn_c1))+Vn_f1^2*Gtf-Vn_f1*Vn_ac2*(Gtf*cos(θn_f1-θn_ac2)+Btf*sin(θn_f1-θn_ac2)) ~ 0,
    -Vn_f2^2*Btf-Vn_f2*Vn_ac3*(Gtf*sin(θn_f2-θn_ac3)-Btf*cos(θn_f2-θn_ac3))-Vn_f2^2*bf-Vn_f2^2*Bc-Vn_f2*Vn_c2*(Gc*sin(θn_f2-θn_c2)-Bc*cos(θn_f2-θn_c2)) ~ 0,
    Vn_f2^2*Gc-Vn_f2*Vn_c2*(Gc*cos(θn_f2-θn_c2)+Bc*sin(θn_f2-θn_c2))+Vn_f2^2*Gtf-Vn_f2*Vn_ac3*(Gtf*cos(θn_f2-θn_ac3)+Btf*sin(θn_f2-θn_ac3)) ~ 0,
]
eqs6=[
    M1-(sqrt(2)*Vn_c1)/(sqrt(3)*Vn_dc1) ~ 0,
    M2-(sqrt(2)*Vn_c2)/(sqrt(3)*Vn_dc2) ~ 0,
    
]
eqs=vcat(eqs3,eqs4,eqs5,eqs6)
parameters=[Vn_ac1,θn_ac1,Vn_ac3,G11,B11,G12,B12,G13,B13,G14,B14,G15,B15,G21,B21,G22,B22,G23,B23,G24,B24,G25,B25,G31,B31,G32,B32,G33,B33,G34,B34,G35,B35,G41,B41,G42,B42,G43,B43,G44,B44,G45,B45,G51,B51,G52,B52,G53,B53,G54,B54,G55,B55,Pac2,Qac2,Pac3,Qac3,Pac4,Qac4,Pac5,Qac5,Pdc1,Pdc2,Gdc11,Gdc12,Gdc21,Gdc22,eta,bf,Btf,Gtf,Ps1,Qs1,Gc,Bc,Vn_f2,Ps2]
@mtkbuild ns = NonlinearSystem(eqs, [θn_ac2,Vn_ac2,θn_ac4,Vn_ac4,θn_ac5,Vn_ac5,θn_ac3,Vn_dc1,Vn_dc2,Vn_f1,θn_f1,M1,Vn_c1,θn_c1,θn_f2,M2,Vn_c2,θn_c2], parameters)
u0=[
    θn_ac2 => 0.0,
    Vn_ac2 => 1.0,
    Vn_dc1 => 1.0,
    Vn_dc2 => 0.9,
    θn_ac4 => 0.0,
    Vn_ac4 => 1.0,
    θn_ac5 => 0.0,
    Vn_ac5 => 1.0,
    θn_ac3 => 0.0,
    Ps2 => 0.0,
    Vn_f1 => 1.0,
    θn_f1 => 0.0,
    M1 => 0.9,
    M2 => 0.9,
    Vn_c1 => 0.9,
    θn_c1 => 0.0,
    θn_f2 => 0.0,
    Vn_c2 => 0.9,
    θn_c2 => 0.0
]
ps=[
    Vn_ac1 => 1.06,
    θn_ac1 => 0.0,
    Vn_ac3 => 1.0,
    Pac2 => -0.2,
    Qac2 => 0.1,
    Pac3 => 0.45,
    Qac3 => 0.15,
    Pac4 => 0.4,
    Qac4 => 0.05,
    Pac5 => 0.6,
    Qac5 => 0.1,
    Gdc11 => 19.2308,
    Gdc12 => -19.2308,
    Gdc21 => -19.2308,
    Gdc22 => 19.2308,
    eta => 1,
    Pdc1 => 10.0/baseMVA,
    Pdc2 => 10.0/baseMVA,
    G11 =>6.25,
    B11 => -18.698,
    G12 => -5.0,
    B12 => 15.0,
    G13 => -1.25,
    B13 => 3.75,
    G14 => 0.0,
    B14 => 0.0,
    G15 => 0.0,
    B15 => 0.0,
    G21 => -5.0,
    B21 => 15.0,
    G22 => 10.8333,
    B22 => -32.415,
    G23 => -1.6667,
    B23 => 5.0,
    G24 => -1.6667,
    B24 => 5.0,
    G25 => -2.5,
    B25 => 7.5,
    G31 => -1.25,
    B31 => 3.75,
    G32 => -1.6667,
    B32 => 5.0,
    G33 => 12.9167,
    B33 => -38.695,
    G34 => -10.0,
    B34 => 30.0,
    G35 => -0.0,
    B35 => 0.0,
    G41 => 0.0,
    B41 => 0.0,
    G42 => -1.6667,
    B42 => 5.0,
    G43 => -10.0,
    B43 => 30.0,
    G44 => 12.9167,
    B44 => -38.695,
    G45 => -1.25,
    B45 => 3.75,
    G51 => 0.0,
    B51 => 0.0,
    G52 => -2.5,
    B52 => 7.5,
    G53 => -0.0,
    B53 => 0.0,
    G54 => -1.25,
    B54 => 3.75,
    G55 => 3.75,
    B55 => -11.21,
    Gtf => 0.1193,
    Btf => -8.919,
    bf => 0.0887,
    Ps1 => -0.6,
    Qs1 => -0.4,
    Bc => -6.08716599381522,
    Gc => 0.0037053603565955803,
    Vn_f2 => 0.9,
    Ps2 => 0.0
]
prob = NonlinearProblem(ns, u0, ps,jac = true)
sol = solve(prob, NewtonRaphson(),show_trace = Val(true))