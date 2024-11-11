using ModelingToolkit, NonlinearSolve

N = 7
@parameters Gnm[1:N, 1:N] Bnm[1:N, 1:N] Gdc[1:N, 1:N] U[1:N, 1:N] W[1:N] ηi[6] ηi[7] ηr[1] ηr[2] phi[1] phi[2]
@parameters Vn[1] θnm[1] Pn[2] Vn[2] Pn[3] Qn[3] Pn[4] Qn[4] Pn[5] Qn[5] Pn[6] Vn[7]
@variables θnm_un[2] Qn_un[2] Vn_un[3] θnm_un[3] Vn_un[4] θnm_un[4] Vn_un[5] θnm_un[5] Vn_un[6] Pn_un[7] P_ac[1] Q_ac[1] P_dc[6] P_ac[2] Q_ac[2] P_dc[7]
@variables s_ac1 s_dc1 s_ac2 s_dc2  # 引入松弛变量
equations = []
function apply_constraints()
    equations=[]
    push!(equations, 0 ~ P_ac[2]*P_dc[7])
    #miss
    push!(equations, 0 ~ -Q_ac[2] + P_ac[2] * tan(phi[2]))
    push!(equations, 0 ~ P_ac[2] - s_ac1^2)
    push!(equations, 0 ~ P_dc[7] - s_dc1^2)
    push!(equations, 0 ~ P_ac[1]*P_dc[6])
    #miss
    push!(equations, 0 ~ -Q_ac[1] + P_ac[1] * tan(phi[1]))
    push!(equations, 0 ~ P_ac[1] - s_ac2^2)
    push!(equations, 0 ~ P_dc[6] - s_dc2^2)
    #miss
    push!(equations, 0 ~ Vn[1] - 0.99*Vn_un[6])
    return equations
end
function main_equations()
    equations=[]
    #For the AC bus 
    a=0.0
    b=0.0
    a= -Pn[3] + U[3,1] * ((1-W[3])*(1-W[1]) * (Vn_un[3]^2 * Gnm[3,1] - Vn_un[3] * Vn[1] * (Gnm[3,1] * cos(θnm_un[3]-θnm[1]) + Bnm[3,1] * sin(θnm_un[3]-θnm[1]))))+ U[3,2] * ((1-W[3])*(1-W[2]) * (Vn_un[3]^2 * Gnm[3,2] - Vn_un[3] * Vn[2] * (Gnm[3,2] * cos(θnm_un[3]-θnm_un[2]) + Bnm[3,2] * sin(θnm_un[3]-θnm_un[2]))))+ U[3,4] * ((1-W[3])*(1-W[4]) * (Vn_un[3]^2 * Gnm[3,4] - Vn_un[3] * Vn_un[4] * (Gnm[3,4] * cos(θnm_un[3]-θnm_un[4]) + Bnm[3,4] * sin(θnm_un[3]-θnm_un[4]))))+ U[3,5] * ((1-W[3])*(1-W[5]) * (Vn_un[3]^2 * Gnm[3,5] - Vn_un[3] * Vn_un[5] * (Gnm[3,5] * cos(θnm_un[3]-θnm_un[5]) + Bnm[3,5] * sin(θnm_un[3]-θnm_un[5]))))
    push!(equations, 0 ~ a)
    b= -Qn[3] + U[3,1] * ((1-W[3])*(1-W[1]) * (-Vn_un[3]^2 * Bnm[3,1] - Vn_un[3] * Vn[1] * (Gnm[3,1] * sin(θnm_un[3]-θnm[1]) - Bnm[3,1] * cos(θnm_un[3]-θnm[1]))))+ U[3,2] * ((1-W[3])*(1-W[2]) * (-Vn_un[3]^2 * Bnm[3,2] - Vn_un[3] * Vn[2] * (Gnm[3,2] * sin(θnm_un[3]-θnm_un[2]) - Bnm[3,2] * cos(θnm_un[3]-θnm_un[2]))))+ U[3,4] * ((1-W[3])*(1-W[4]) * (-Vn_un[3]^2 * Bnm[3,4] - Vn_un[3] * Vn_un[4] * (Gnm[3,4] * sin(θnm_un[3]-θnm_un[4]) - Bnm[3,4] * cos(θnm_un[3]-θnm_un[4]))))+ U[3,5] * ((1-W[3])*(1-W[5]) * (-Vn_un[3]^2 * Bnm[3,5] - Vn_un[3] * Vn_un[5] * (Gnm[3,5] * sin(θnm_un[3]-θnm_un[5]) - Bnm[3,5] * cos(θnm_un[3]-θnm_un[5]))))
    push!(equations, 0 ~ b)

    a=0.0
    b=0.0
    a= -Pn[4] + U[4,1] * ((1-W[4])*(1-W[1]) * (Vn_un[4]^2 * Gnm[4,1] - Vn_un[4] * Vn[1] * (Gnm[4,1] * cos(θnm_un[4]-θnm[1]) + Bnm[4,1] * sin(θnm_un[4]-θnm[1]))))+ U[4,2] * ((1-W[4])*(1-W[2]) * (Vn_un[4]^2 * Gnm[4,2] - Vn_un[4] * Vn[2] * (Gnm[4,2] * cos(θnm_un[4]-θnm_un[2]) + Bnm[4,2] * sin(θnm_un[4]-θnm_un[2]))))+ U[4,3] * ((1-W[4])*(1-W[3]) * (Vn_un[4]^2 * Gnm[4,3] - Vn_un[4] * Vn_un[3] * (Gnm[4,3] * cos(θnm_un[4]-θnm_un[3]) + Bnm[4,3] * sin(θnm_un[4]-θnm_un[3]))))+ U[4,5] * ((1-W[4])*(1-W[5]) * (Vn_un[4]^2 * Gnm[4,5] - Vn_un[4] * Vn_un[5] * (Gnm[4,5] * cos(θnm_un[4]-θnm_un[5]) + Bnm[4,5] * sin(θnm_un[4]-θnm_un[5]))))
    push!(equations, 0 ~ a)
    b= -Qn[4] + U[4,1] * ((1-W[4])*(1-W[1]) * (-Vn_un[4]^2 * Bnm[4,1] - Vn_un[4] * Vn[1] * (Gnm[4,1] * sin(θnm_un[4]-θnm[1]) - Bnm[4,1] * cos(θnm_un[4]-θnm[1]))))+ U[4,2] * ((1-W[4])*(1-W[2]) * (-Vn_un[4]^2 * Bnm[4,2] - Vn_un[4] * Vn[2] * (Gnm[4,2] * sin(θnm_un[4]-θnm_un[2]) - Bnm[4,2] * cos(θnm_un[4]-θnm_un[2]))))+ U[4,3] * ((1-W[4])*(1-W[3]) * (-Vn_un[4]^2 * Bnm[4,3] - Vn_un[4] * Vn_un[3] * (Gnm[4,3] * sin(θnm_un[4]-θnm_un[3]) - Bnm[4,3] * cos(θnm_un[4]-θnm_un[3]))))+ U[4,5] * ((1-W[4])*(1-W[5]) * (-Vn_un[4]^2 * Bnm[4,5] - Vn_un[4] * Vn_un[5] * (Gnm[4,5] * sin(θnm_un[4]-θnm_un[5]) - Bnm[4,5] * cos(θnm_un[4]-θnm_un[5]))))
    push!(equations, 0 ~ b)

    a=0.0
    b=0.0
    a= -Pn[5] + U[5,1] * ((1-W[5])*(1-W[1]) * (Vn_un[5]^2 * Gnm[5,1] - Vn_un[5] * Vn[1] * (Gnm[5,1] * cos(θnm_un[5]-θnm[1]) + Bnm[5,1] * sin(θnm_un[5]-θnm[1]))))+ U[5,2] * ((1-W[5])*(1-W[2]) * (Vn_un[5]^2 * Gnm[5,2] - Vn_un[5] * Vn[2] * (Gnm[5,2] * cos(θnm_un[5]-θnm_un[2]) + Bnm[5,2] * sin(θnm_un[5]-θnm_un[2]))))+ U[5,3] * ((1-W[5])*(1-W[3]) * (Vn_un[5]^2 * Gnm[5,3] - Vn_un[5] * Vn_un[3] * (Gnm[5,3] * cos(θnm_un[5]-θnm_un[3]) + Bnm[5,3] * sin(θnm_un[5]-θnm_un[3]))))+ U[5,4] * ((1-W[5])*(1-W[4]) * (Vn_un[5]^2 * Gnm[5,4] - Vn_un[5] * Vn_un[4] * (Gnm[5,4] * cos(θnm_un[5]-θnm_un[4]) + Bnm[5,4] * sin(θnm_un[5]-θnm_un[4]))))
    push!(equations, 0 ~ a)
    b= -Qn[5] + U[5,1] * ((1-W[5])*(1-W[1]) * (-Vn_un[5]^2 * Bnm[5,1] - Vn_un[5] * Vn[1] * (Gnm[5,1] * sin(θnm_un[5]-θnm[1]) - Bnm[5,1] * cos(θnm_un[5]-θnm[1]))))+ U[5,2] * ((1-W[5])*(1-W[2]) * (-Vn_un[5]^2 * Bnm[5,2] - Vn_un[5] * Vn[2] * (Gnm[5,2] * sin(θnm_un[5]-θnm_un[2]) - Bnm[5,2] * cos(θnm_un[5]-θnm_un[2]))))+ U[5,3] * ((1-W[5])*(1-W[3]) * (-Vn_un[5]^2 * Bnm[5,3] - Vn_un[5] * Vn_un[3] * (Gnm[5,3] * sin(θnm_un[5]-θnm_un[3]) - Bnm[5,3] * cos(θnm_un[5]-θnm_un[3]))))+ U[5,4] * ((1-W[5])*(1-W[4]) * (-Vn_un[5]^2 * Bnm[5,4] - Vn_un[5] * Vn_un[4] * (Gnm[5,4] * sin(θnm_un[5]-θnm_un[4]) - Bnm[5,4] * cos(θnm_un[5]-θnm_un[4]))))
    push!(equations, 0 ~ b)

    #VSC AC bus
    a = 0.0
    a= P_ac[1] - ηi[6]*P_dc[6]+U[1,2] * ((1-W[1])*(1-W[2]) * (Vn[1]^2 * Gnm[1,2] - Vn[1] * Vn[2] * (Gnm[1,2] * cos(θnm[1]-θnm_un[2]) + Bnm[1,2] * sin(θnm[1]-θnm_un[2]))))+ U[1,3] * ((1-W[1])*(1-W[3]) * (Vn[1]^2 * Gnm[1,3] - Vn[1] * Vn_un[3] * (Gnm[1,3] * cos(θnm[1]-θnm_un[3]) + Bnm[1,3] * sin(θnm[1]-θnm_un[3]))))+ U[1,4] * ((1-W[1])*(1-W[4]) * (Vn[1]^2 * Gnm[1,4] - Vn[1] * Vn_un[4] * (Gnm[1,4] * cos(θnm[1]-θnm_un[4]) + Bnm[1,4] * sin(θnm[1]-θnm_un[4]))))+ U[1,5] * ((1-W[1])*(1-W[5]) * (Vn[1]^2 * Gnm[1,5] - Vn[1] * Vn_un[5] * (Gnm[1,5] * cos(θnm[1]-θnm_un[5]) + Bnm[1,5] * sin(θnm[1]-θnm_un[5]))))
    push!(equations, 0 ~ a)

    a = 0.0
    b = 0.0
    #miss
    a = -Pn[2] + P_ac[2]-ηi[7]*P_dc[7] +U[2,1] * ((1-W[2])*(1-W[1]) * (Vn[2]^2 * Gnm[2,1] - Vn[2] * Vn[1] * (Gnm[2,1] * cos(θnm_un[2]-θnm[1]) + Bnm[2,1] * sin(θnm_un[2]-θnm[1]))))+ U[2,3] * ((1-W[2])*(1-W[3]) * (Vn[2]^2 * Gnm[2,3] - Vn[2] * Vn_un[3] * (Gnm[2,3] * cos(θnm_un[2]-θnm_un[3]) + Bnm[2,3] * sin(θnm_un[2]-θnm_un[3]))))+ U[2,4] * ((1-W[2])*(1-W[4]) * (Vn[2]^2 * Gnm[2,4] - Vn[2] * Vn_un[4] * (Gnm[2,4] * cos(θnm_un[2]-θnm_un[4]) + Bnm[2,4] * sin(θnm_un[2]-θnm_un[4]))))+ U[2,5] * ((1-W[2])*(1-W[5]) * (Vn[2]^2 * Gnm[2,5] - Vn[2] * Vn_un[5] * (Gnm[2,5] * cos(θnm_un[2]-θnm_un[5]) + Bnm[2,5] * sin(θnm_un[2]-θnm_un[5]))))
    b=-Qn_un[2]+ Q_ac[2]+U[2,1] * ((1-W[2])*(1-W[1]) * (-Vn[2]^2 * Bnm[2,1] - Vn[2] * Vn[1] * (Gnm[2,1] * sin(θnm_un[2]-θnm[1]) - Bnm[2,1] * cos(θnm_un[2]-θnm[1]))))+ U[2,3] * ((1-W[2])*(1-W[3]) * (-Vn[2]^2 * Bnm[2,3] - Vn[2] * Vn_un[3] * (Gnm[2,3] * sin(θnm_un[2]-θnm_un[3]) - Bnm[2,3] * cos(θnm_un[2]-θnm_un[3]))))+ U[2,4] * ((1-W[2])*(1-W[4]) * (-Vn[2]^2 * Bnm[2,4] - Vn[2] * Vn_un[4] * (Gnm[2,4] * sin(θnm_un[2]-θnm_un[4]) - Bnm[2,4] * cos(θnm_un[2]-θnm_un[4]))))+ U[2,5] * ((1-W[2])*(1-W[5]) * (-Vn[2]^2 * Bnm[2,5] - Vn[2] * Vn_un[5] * (Gnm[2,5] * sin(θnm_un[2]-θnm_un[5]) - Bnm[2,5] * cos(θnm_un[2]-θnm_un[5]))))
    push!(equations, 0 ~ a)
    push!(equations, 0 ~ b)

    #VSC DC bus
    a = 0.0
    a=-Pn[6] + P_dc[6] - ηr[1]*P_ac[1]+U[6,7] * (W[6]*W[7] * (Gdc[6,7] * (Vn_un[6]^2 - Vn_un[6] * Vn[7])))
    #miss
    push!(equations, 0 ~ a)

    a=0.0
    a=-Pn_un[7]+ P_dc[7]- ηr[2]*P_ac[2]+U[7,6] * (W[7]*W[6] * (Gdc[7,6] * (Vn[7]^2 - Vn[7] * Vn_un[6])))
    #miss
    push!(equations, 0 ~ a)  
    return equations
end
variables=[]
variables=vcat(θnm_un[2], Qn_un[2], Vn_un[3], θnm_un[3], Vn_un[4], θnm_un[4], Vn_un[5], θnm_un[5], Vn_un[6], Pn_un[7], P_ac[1], Q_ac[1], P_dc[6], P_ac[2], Q_ac[2], P_dc[7], s_ac1, s_dc1, s_ac2, s_dc2)
parameters = []
for i in 1:N
    for j in 1:N
        parameters=vcat(parameters, Gnm[i,j], Bnm[i,j], Gdc[i,j], U[i,j], W[i])
    end
end
parameters=vcat(parameters, ηi[6], ηi[7], ηr[1], ηr[2], phi[1], phi[2], Vn[1], θnm[1], Pn[2], Vn[2], Pn[3], Qn[3], Pn[4], Qn[4], Pn[5], Qn[5], Pn[6], Vn[7])
eq1=main_equations()
eq2=apply_constraints()
equations=vcat(equations, eq1, eq2)
sys = NonlinearSystem(equations, variables, parameters, name=:ACDCPowerFlow)
# Define initial conditions and parameters for the modeltest.jl system
u0_modeltest = Dict()
ps_modeltest = Dict()

#initialize the variables
u0_modeltest[θnm_un[2]] = 0.0
u0_modeltest[Qn_un[2]] = 10.0
u0_modeltest[Vn_un[3]] = 1.0
u0_modeltest[θnm_un[3]] = 0.0
u0_modeltest[Vn_un[4]] = 1.0
u0_modeltest[θnm_un[4]] = 0.0
u0_modeltest[Vn_un[5]] = 1.0
u0_modeltest[θnm_un[5]] = 0.0
u0_modeltest[Vn_un[6]] = 1.0707070707070707
u0_modeltest[Pn_un[7]] = 0.0
u0_modeltest[P_ac[1]] = 60.0
u0_modeltest[Q_ac[1]] = 103.92
u0_modeltest[P_dc[6]] = 0.0
u0_modeltest[P_ac[2]] = 0.0
u0_modeltest[Q_ac[2]] = 0.0
u0_modeltest[P_dc[7]] = 0.0
u0_modeltest[s_ac1] = 7.746
u0_modeltest[s_dc1] = 0.0
u0_modeltest[s_ac2] = 0.0
u0_modeltest[s_dc2] =0.0

#initialize the parameters
G=[
    6.25 -5 -1.25 0.0 0.0 0.0 0.0;
    -5 10.8333 -1.6667 -1.6667 -2.5 0.0 0.0;
    -1.25 -1.6667 12.9167 -10 0.0 0.0 0.0;
    0.0 -1.6667 -10 12.6667 -1.25 0.0 0.0;
    0.0 -2.5 0.0 -1.25 3.75 0.0 0.0;
    0.0 0.0 0.0 0.0 0.0 0.0 0.0;
    0.0 0.0 0.0 0.0 0.0 0.0 0.0
]
for i in 1:7
    for j in 1:7
        ps_modeltest[Gnm[i,j]] = G[i,j]
    end
end
B=[
    -18.695 15.0 -1.25 0.0 0.0 0.0 0.0;
    15.0 -32.415 5 5 7.5 0.0 0.0;
    3.75 5 -38.695 30.0 0.0 0.0 0.0;
    0.0 5 30.0 -38.695 3.75 0.0 0.0;
    0.0 7.5 0.0 3.75 -11.21 0.0 0.0;
    0.0 0.0 0.0 0.0 0.0 0.0 0.0;
    0.0 0.0 0.0 0.0 0.0 0.0 0.0 
]
for i in 1:7
    for j in 1:7
        ps_modeltest[Bnm[i,j]] = B[i,j]
    end 
end
Wa=[0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0]
for i in 1:7
    ps_modeltest[W[i]] = Wa[i]
end
Ua=[
    1.0 1.0 1.0 0.0 0.0 1.0 0.0;
    1.0 1.0 1.0 1.0 1.0 0.0 1.0;
    1.0 1.0 1.0 1.0 0.0 0.0 0.0;
    0.0 1.0 1.0 1.0 1.0 0.0 0.0;
    0.0 1.0 0.0 1.0 1.0 0.0 0.0;
    1.0 0.0 0.0 0.0 0.0 1.0 1.0;
    0.0 1.0 0.0 0.0 0.0 1.0 1.0
]
for i in 1:7
    for j in 1:7
        ps_modeltest[U[i,j]] = Ua[i,j]
    end
end
Gdca=[
    0.0 0.0 0.0 0.0 0.0 0.0 0.0;
    0.0 0.0 0.0 0.0 0.0 0.0 0.0;
    0.0 0.0 0.0 0.0 0.0 0.0 0.0;
    0.0 0.0 0.0 0.0 0.0 0.0 0.0;
    0.0 0.0 0.0 0.0 0.0 0.0 0.0;
    0.0 0.0 0.0 0.0 0.0 0.0 19.231;
    0.0 0.0 0.0 0.0 0.0 19.231 0.0
]
for i in 1:7
    for j in 1:7
        ps_modeltest[Gdc[i,j]] = Gdca[i,j]
    end
end
    ps_modeltest[ηi[6]] = 0.887
    ps_modeltest[ηi[7]] = 0.887


    ps_modeltest[ηr[1]] = 0.887
    ps_modeltest[ηr[2]] = 0.887
    ps_modeltest[phi[1]] = 1.047
    ps_modeltest[phi[2]] = 1.047
   
    ps_modeltest[θnm[1]] = 0.0
    ps_modeltest[Pn[2]] = 20.0
    ps_modeltest[Pn[3]] = 45.0
    ps_modeltest[Pn[4]] = 40.0
    ps_modeltest[Pn[5]] = 60.0
    ps_modeltest[Pn[6]] = 0.0
    
    ps_modeltest[Qn[3]] = 15.0
    ps_modeltest[Qn[4]] = 5.0
    ps_modeltest[Qn[5]] = 10.0
    
    ps_modeltest[Vn[1]] = 1.06
    ps_modeltest[Vn[2]] = 1.0
    ps_modeltest[Vn[7]] = 1.0

    simplified_sys = structural_simplify(sys)
    prob_modeltest = NonlinearProblem(simplified_sys, u0_modeltest, ps_modeltest)
    sol_modeltest = solve(prob_modeltest, NewtonRaphson())
