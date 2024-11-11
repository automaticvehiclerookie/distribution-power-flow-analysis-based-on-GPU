mutable struct Sd
    z::Vector{ComplexF64}
    i::Vector{ComplexF64}
    p::Vector{ComplexF64}
end

function makeSdaczip(baseMVA, bus)
    (PQ, PV, REF,NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA,
     VM,VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN) =PowerFlow.idx_ACbus();
     (Pdc, Vdc, isolated, BUSDC_I, BUSAC_I, GRID, PDC, VDC, BASEKVDC, VDCMAX, VDCMIN, CDC) =PowerFlow.idx_DCbus();
        pw = [1 0 0]
        qw = pw
    z = (bus[:, PD] * pw[3]  + 1im * bus[:, QD] * qw[3]) / baseMVA
    i = (bus[:, PD] * pw[2]  + 1im * bus[:, QD] * qw[2]) / baseMVA
    p = (bus[:, PD] * pw[1]  + 1im * bus[:, QD] * qw[1]) / baseMVA
    sd = Sd(z, i, p)
    return sd
end