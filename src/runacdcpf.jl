"""
   Main function to call the AC/DC hybrid power flow function
    Input: case file
    Output: results of the power flow as a dictionary
    Example:
    bus, gen, branch = run_hybridpf(casefile)
"""
function runacdcpf(mpc::Dict)
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

baseMVA = mpc["baseMVA"];
busAC = mpc["busAC"];
busDC = mpc["busDC"];
gen = mpc["genAC"];
branchAC = mpc["branchAC"];
branchDC = mpc["branchDC"];
converter = mpc["converter"];
alg="NR";

(busAC, genAC, branchAC) = PowerFlow.ext2int(busAC, gen, branchAC);
#get the bus type from bus matrix
(ref, pv, pq, pdc, vdc) =PowerFlow.hybrid_bustypes(busAC, busDC,converter);

## generator info
on = findall(genAC[:, GEN_STATUS] .> 0)  # which generators are on?
gbus = genAC[on, GEN_BUS]  # what buses are they at?
##-----  run the power flow  ----- 
its = 0; 
V0  = busAC[:, VM] .* exp.(1im * pi/180 * busAC[:, VA])
vcb = ones(size(V0));           ## create mask of voltage-controlled buses
vcb[pq] .= 0;                    ## exclude PQ buses
gbus=Int.(gbus);
k = findall(Bool.(vcb[gbus]));            ## in-service gens at v-c buses
V0[gbus[k]] = gen[on[k], VG] ./ abs.(V0[gbus[k]]).* V0[gbus[k]];
#make Ybus for AC/DC hybrid system
#make Gbus for DC network
#TODO: 交流网连通，平衡节点位于交流网络中
Ybus,Yf,Yt = PowerFlow.makeYbus(baseMVA, busAC, branchAC);
Gbus=PowerFlow.makeGbus(baseMVA, busDC, branchDC);
V0dc = busDC[:, VDC] ;
#newton power flow
repeat=1;
while (repeat>0)
    Sbusac = Vm -> PowerFlow.makeacSbus(baseMVA, busAC, gen, Vm);
    Sbusdc = Vmdc -> PowerFlow.makedcSbus(baseMVA, busDC, gen, V0dc);
    if alg == "NR"
        V_ac,Vm_dc,converged,iterations= PowerFlow.acdcnewtonpf(Ybus,Gbus,Sbusac,Sbusdc,pv,pq,pdc,vdc,busAC,busDC,converter,baseMVA,V0dc,V0)
    end
    its += iterations;
    #TODO:修改pfsoln函数，使其能够处理直流网络和转换器
    busAC,busDC, gen, branchAC,branchDC = pfsolnh(baseMVA, busAC,busDC, gen, branchAC, branchDC, Ybus, Gbus, Yf, Yt, V_ac,Vm_dc, ref, pv, pq, pdc, vdc);
    repeat=0;
end
#return the results
mpc["busAC"] = busAC;
mpc["busDC"] = busDC;
mpc["branchAC"] = branchAC;
mpc["branchDC"] = branchDC;
return mpc
end