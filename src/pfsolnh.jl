function pfsolnh(baseMVA, busAC,busDC, gen, branchAC, branchDC, Ybus, Gbus, Yf, Yt, V_ac,Vm_dc, ref, pv, pq, pdc, vdc)
    (PQ, PV, REF,NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA,
     VM,VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN) =PowerFlow.idx_ACbus();
    (Pdc, Vdc, isolated, BUSDC_I, TYPE_DC, GRID, PDC, VDC, BASEKVDC, VDCMAX, VDCMIN, CDC) =PowerFlow.idx_DCbus();
 (F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, TAP, SHIFT, BR_STATUS,
  PF, QF, PT, QT, MU_SF, MU_ST, MU_ANGMIN, MU_ANGMAX) = PowerFlow.idx_ACbrch();
    (F_BUSDC, T_BUSDC, BR_RDC, LDC, CDC, RATE_A, RATE_B, RATE_C, STATUS) = PowerFlow.idx_DCbrch();
    (BUSDC_I, BUSAC_I, P_G, Q_G, VTAR, RTF, XTF, BF, RC, XC, BASEKVAC, VMMAX, VMMIN,
     IMAX, STATUS, LOSSA, LOSSB, LOSSCREC, LOSSCINV)=PowerFlow.idx_conv()
(GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, 
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, 
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF) = PowerFlow.idx_gen();

    ##----- update bus voltages -----
    busAC[:, VM] = abs.(V_ac)
    busAC[:, VA] = angle.(V_ac) * 180 / π
    busDC[:, VDC] = abs.(Vm_dc)
    #precision control
    busAC = round.(busAC, digits=6)
    busDC = round.(busDC, digits=6)
    # ##----- update Qg for gens at PV/slack buses and Pg for slack bus(es) -----
    # ## generator info
    # on=findall((gen[:, GEN_STATUS].>0).*(busAC[Int.(gen[:, GEN_BUS]), BUS_TYPE] .!= PQ))   # Which generators are on and not at PQ buses?
    # off = findall(x -> x <= 0, gen[:, GEN_STATUS])  # Which generators are off?
    # Sbus = V_ac[gbus] .* conj.(Ybus[gbus, :] * V_ac)
    return busAC,busDC, gen, branchAC, branchDC
end