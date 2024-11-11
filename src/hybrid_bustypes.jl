function hybrid_bustypes(busAC::Matrix{Float64}, busDC::Matrix{Float64},converter::Matrix{Float64})
    # constants
    (PQ, PV, REF,NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA,
     VM,VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN) =PowerFlow.idx_ACbus();
     (Pdc, Vdc, isolated, BUSDC_I, TYPE_DC, GRID, PDC, VDC, BASEKVDC, VDCMAX, VDCMIN, CDC) =PowerFlow.idx_DCbus();
     (BUSDC_I, BUSAC_I, P_G, Q_G, VTAR, RTF, XTF, BF, RC, XC, BASEKVAC, VMMAX, VMMIN,
     IMAX, STATUS, LOSSA, LOSSB, LOSSCREC, LOSSCINV)=PowerFlow.idx_conv()
    # get the bustype from the bus matrix

    ref = findall(busAC[:, BUS_TYPE] .== REF)
    pv  = findall(busAC[:, BUS_TYPE] .== PV)
    pq  = findall(busAC[:, BUS_TYPE] .== PQ)
    pdc = findall(busDC[:, TYPE_DC] .== Pdc)
    vdc = findall(busDC[:, TYPE_DC] .== Vdc)
    return ref, pv, pq, pdc, vdc
end