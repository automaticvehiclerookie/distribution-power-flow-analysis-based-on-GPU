using ModelingToolkit, NonlinearSolve
const baseMVA = 100.0
@variables θn_ac2,Vn_dc1,Vn_dc2,θn_f,M,θn_c,Vn_c
@parameters Vn_ac1,θn_ac1,Vn_ac2,Vn_f
@parameters G11,B11,G12,B12,G21,B21,G22,B22
@parameters Pac2,Qac2,Pdc1,Pdc2
@parameters Gdc11,Gdc12,Gdc21,Gdc22,eta,Gtf,Btf,bf,Ps,Bc,Gc,Qs,a,b,c

eqs3=[
    -Vn_ac2^2*Gtf+Vn_ac2*Vn_f*(Gtf*cos(θn_ac2-θn_f)+Btf*sin(θn_ac2-θn_f))-Ps ~ 0,
    # Vn_ac2^2*Btf+Vn_ac2*Vn_f*(Gtf*sin(θn_ac2-θn_f)-Btf*cos(θn_ac2-θn_f))-Qs ~ 0,
]

eqs4=[
    Vn_ac2*Vn_ac1*(G21*cos(θn_ac2-θn_ac1)+B21*sin(θn_ac2-θn_ac1))+Vn_ac2*Vn_ac2*(G22*cos(θn_ac2-θn_ac2)+B22*sin(θn_ac2-θn_ac2))+Pac2+Vn_ac2^2*Gtf-Vn_ac2*Vn_f*(Gtf*cos(θn_ac2-θn_f)+Btf*sin(θn_ac2-θn_f))~ 0,
    # Vn_ac2*Vn_ac1*(G21*sin(θn_ac2-θn_ac1)-B21*cos(θn_ac2-θn_ac1))+Vn_ac2*Vn_ac2*(G22*sin(θn_ac2-θn_ac2)-B22*cos(θn_ac2-θn_ac2))+Qac2-Vn_ac2^2*Btf-Vn_ac2*Vn_f*(Gtf*sin(θn_ac2-θn_f)-Btf*cos(θn_ac2-θn_f))~ 0,
]

eqs5=[
    Vn_dc1*Vn_dc1*Gdc11+Vn_dc1*Vn_dc2*Gdc12+Pdc1+eta*(Vn_c^2*Gc-Vn_c*Vn_f*(Gc*cos(θn_c-θn_f)+Bc*sin(θn_c-θn_f)))+a*((Vn_c^2*Gc-Vn_c*Vn_f*(Gc*cos(θn_c-θn_f)+Bc*sin(θn_c-θn_f)))^2 +(-Vn_c^2*Bc-Vn_c*Vn_f*(Gc*sin(θn_c-θn_f)-Bc*cos(θn_c-θn_f)))^2)*inv(Vn_c)^2+b*(sqrt.(((Vn_c^2*Gc-Vn_c*Vn_f*(Gc*cos(θn_c-θn_f)+Bc*sin(θn_c-θn_f)))^2 +(-Vn_c^2*Bc-Vn_c*Vn_f*(Gc*sin(θn_c-θn_f)-Bc*cos(θn_c-θn_f)))^2)*inv(Vn_c)^2))+c ~ 0,
    Vn_dc2*Vn_dc1*Gdc21+Vn_dc2*Vn_dc2*Gdc22+Pdc2 ~ 0,
]
eqs6=[
    M-(sqrt(2)*Vn_c)/(sqrt(3)*Vn_dc1) ~ 0,
    -Vn_f^2-Vn_f*Vn_ac2*(Gtf*sin(θn_f-θn_ac2)-Btf*cos(θn_f-θn_ac2))-Vn_f^2 *bf-Vn_f^2*Bc-Vn_f*Vn_c*(Gc*sin(θn_f-θn_c)-Bc*cos(θn_f-θn_c)) ~ 0,
    Vn_f^2*Gc-Vn_f*Vn_c*(Gc*cos(θn_f-θn_c)+Bc*sin(θn_f-θn_c))+Vn_f^2*Gtf-Vn_f*Vn_ac2*(Gtf*cos(θn_f-θn_ac2)+Btf*sin(θn_f-θn_ac2)) ~ 0,
]
eqs=vcat(eqs3,eqs4,eqs5,eqs6)
parameters=[Vn_ac1,θn_ac1,G11,B11,G12,B12,G21,B21,G22,B22,Pac2,Qac2,Pdc1,Pdc2,Gdc11,Gdc12,Gdc21,Gdc22,eta,Gtf,Btf,bf,Ps,Bc,Gc,Qs,a,b,c,Vn_ac2,Vn_f]
@mtkbuild ns = NonlinearSystem(eqs, [θn_ac2,Vn_dc1,Vn_dc2,Vn_f,M,θn_c,Vn_c], parameters)
u0=[
    θn_ac2 => 0.0,
    Vn_dc1 => 1.0,
    Vn_dc2 => 0.9,
    θn_f => 0.0,
    M =>0.9,
    θn_c => 0.0,
    Vn_c => 0.9
]
ps=[
    Vn_ac1 => 1.06,
    θn_ac1 => 0.0,
    Pac2 => -0.4,
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
    eta => 1,
    Pdc1 => 10.0/baseMVA,
    Pdc2 => 10.0/baseMVA,
    Gtf => 0.1193,
    Btf => -8.919,
    bf => 0.0887,
    Ps => -0.6,
    Qs => -0.4,
    Bc => -6.08716599381522,
    Gc => 0.0037053603565955803,
    a => 0.0001,
    b => 0.0,
    c => 0.0,
    Vn_f => 0.973,
    Vn_ac2 => 1.016
]
prob = NonlinearProblem(ns, u0, ps,jac = true)
sol = solve(prob, NewtonRaphson())