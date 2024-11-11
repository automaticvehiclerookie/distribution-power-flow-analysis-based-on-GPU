#Input : baseMVA, bus, branch, pdc, vdc
#Output: Ybus:contains the DC buses and AC buses,In the power flow analysis,
#we need to delete the DC admittance in Ybus diagonal elements
#Bbus_dc: the DC bus admittance matrix without the AC buses
function makeYacdcbus(baseMVA, bus, branch,ref,pv,pq,pdc,vdc)
    ## define named indices into bus, branch matrices
    (PQ, PV, REF,Pdc,Vdc,NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA,
    VM,VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN) = PowerFlow.idx_bus();
    (F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C,TAP, SHIFT, BR_STATUS, 
 PF, QF, PT, QT, MU_SF, MU_ST, BRANCHMODE, PC, QC, ETACR,ETACI, PHI, M) = PowerFlow.idx_brch();

    # constants
    nb = size(bus, 1)          # number of buses
    nl = size(branch, 1)       # number of lines
    # check that bus numbers are equal to indices to bus (one set of bus numbers)

    if any(bus[:, BUS_I] .!= (1:nb))
        error("makeYbus: buses must be numbered consecutively in bus matrix; use ext2int() to convert to internal ordering")
    end
    # for each branch, compute the elements of the branch admittance matrix where
    stat = branch[:, BR_STATUS]                    # ones at in-service branches
    Ys = stat ./ (branch[:, BR_R] .+ 1im * branch[:, BR_X])  # series admittance
    Bc = stat .* branch[:, BR_B]                           # line charging susceptance
    tap = ones(nl,1)                              # default tap ratio = 1
    i = findall(branch[:, TAP] .!= 0)                       # indices of non-zero tap ratios
    tap[i] = branch[i, TAP]                        # assign non-zero tap ratios
    tap = tap .* exp.(1im*pi/180 * branch[:, SHIFT]) # add phase shifters
    Ytt = Ys + 1im*Bc/2
    Yff = Ytt ./ (tap .* conj.(tap))
    Yft = - Ys ./ conj.(tap)
    Ytf = - Ys ./ tap
     # compute shunt admittance
     Ysh = (bus[:, GS] .+ 1im * bus[:, BS]) / baseMVA  # vector of shunt admittances
     # bus indices
    f = branch[:, F_BUS]                           # list of "from" buses
    t = branch[:, T_BUS]                           # list of "to" buses
    ## build connection matrices
    Cf = sparse(1:nl, f, vec(ones(nl, 1)), nl, nb);      ## connection matrix for line & from buses
    Ct = sparse(1:nl, t, vec(ones(nl, 1)), nl, nb);      ## connection matrix for line & to buses
    ## build Yf and Yt such that Yf * V is the vector of complex branch currents injected
    ## at each branch's "from" bus, and Yt is the same for the "to" bus end
    Yf = sparse(1:nl, 1:nl, vec(Yff), nl, nl) * Cf + sparse(1:nl, 1:nl, vec(Yft), nl, nl) * Ct
    Yt = sparse(1:nl, 1:nl, vec(Ytf), nl, nl) * Cf + sparse(1:nl, 1:nl, vec(Ytt), nl, nl) * Ct
    ## build Ybus
    Ybus = Cf' * Yf + Ct' * Yt +             ## branch admittances
        sparse(1:nb, 1:nb, Ysh, nb, nb);    ## shunt admittance
    Cy=collect(1:size(Ybus,1))
    Ybus=Ybus+Diagonal(vec(sum(Ybus[:,vcat(pdc,vdc)],dims=2)))
    Ybus=Ybus[vcat(ref,pv,pq),vcat(ref,pv,pq)]
    ## build Bbus for DC network
    bus_dc=bus[vcat(pdc,vdc),:]
    Bbus_dc=zeros(size(bus_dc,1),size(bus_dc,1))
    branch_dc=branch[findall(branch[:, BRANCHMODE] .== 2),:]
    Cdc_row=findall(branch_dc[:,F_BUS] .== bus_dc[:,BUS_I])
    Cdc_col=findall(branch_dc[:,T_BUS] .== bus_dc[:,BUS_I])
    Bbus_dc[Cdc_row,Cdc_col] .= -1.0./branch_dc[:,BR_R]
    Bbus_dc[Cdc_col,Cdc_row] .= -1.0./branch_dc[:,BR_R]
    Cdc=collect(1:size(bus_dc,1))
    #TODO:计算自导纳
    Bbus_dc+=Diagonal(vec(-sum(Bbus_dc,dims=2)))
    # end

    return Ybus, Bbus_dc,Yf, Yt
end