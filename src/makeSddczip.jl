function makeSddczip(baseMVA, bus)
     (Pdc, Vdc, isolated, BUSDC_I, BUSAC_I, GRID, PDC, VDC, BASEKVDC, VDCMAX, VDCMIN, CDC) =PowerFlow.idx_DCbus();
        pw = [1 0 0]
        qw = pw
    z = (bus[:, PDC] * pw[3]  ) / baseMVA
    i = (bus[:, PDC] * pw[2]  ) / baseMVA
    p = (bus[:, PDC] * pw[1] ) / baseMVA
    sd = Sd(z, i, p)
    return sd
end