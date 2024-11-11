function makeacSbus(baseMVA, bus, gen, Vm, Sg=nothing, nargout=1)
    # Define named indices into bus, gen matrices
    (F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, TAP, SHIFT, BR_STATUS,
  PF, QF, PT, QT, MU_SF, MU_ST, MU_ANGMIN, MU_ANGMAX) = PowerFlow.idx_ACbrch();
  (GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, 
  MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, 
  QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF) = PowerFlow.idx_gen();
    nb = size(bus, 1)
    Sd = makeSdaczip(baseMVA, bus)
    
    if nargout == 2
        Sbus = []
        if isempty(Vm)
            dSbus_dVm = spzeros(nb, nb)
        else
            dSbus_dVm = -(spdiagm(0 => Sd.i + 2 .* Vm .* Sd.z))
        end
        return dSbus_dVm
    else
    # Compute per-bus generation in p.u.
        on = findall(gen[:, GEN_STATUS] .> 0)  # which generators are on?
        gbus = gen[on, GEN_BUS]  # what buses are they at?
        ngon = length(on)
        Cg = sparse(gbus, 1:ngon, 1, nb, ngon)  # connection matrix
        if Sg != nothing
            Sbusg = Cg * Sg[on]
        else
            Sbusg = Cg * (gen[on, PG] .+ 1im * gen[on, QG]) / baseMVA
        end
        # element i, j is 1 if gen on(j) at bus i is ON
        Sbusd = Sd.p .+ Sd.i .* Vm .+ Sd.z .* Vm.^2
        Sbus = Sbusg - Sbusd
        return Sbus
    end
end