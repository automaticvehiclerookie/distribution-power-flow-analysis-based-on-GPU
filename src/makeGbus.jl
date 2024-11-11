function makeGbus(baseMVA, busDC, branchDC)
    (F_BUSDC, T_BUSDC, BR_RDC, LDC, CDC, RATE_A, RATE_B, RATE_C, STATUS) = PowerFlow.idx_DCbrch()
    (Pdc, Vdc, isolated, BUSDC_I, BUSAC_I, GRID, PDC, VDC, BASEKVDC, VDCMAX, VDCMIN, CDC) =PowerFlow.idx_DCbus();
    Gbus_dc=zeros(size(busDC,1),size(busDC,1))
    Cdc_row=findall(branchDC[:,F_BUSDC] .== busDC[:,BUSDC_I])
    Cdc_col=findall(branchDC[:,T_BUSDC] .== busDC[:,BUSDC_I])
    Gbus_dc[Cdc_row,Cdc_col] .= -1.0./branchDC[:,BR_RDC]
    Gbus_dc[Cdc_col,Cdc_row] .= -1.0./branchDC[:,BR_RDC]
    Gbus_dc+=Diagonal(vec(-sum(Gbus_dc,dims=2)))
end