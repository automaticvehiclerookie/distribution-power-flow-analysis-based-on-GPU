function case7h_new()
    mpc = Dict{String, Any}();
    mpc["version"] = "2";
    mpc["baseMVA"] = 100;
    #TODO:modify the data to matpower format
    ## bus data
    # bus_i type Pd Qd Gs Bs area Vm Va baseKV zone Vmax Vmin
    #type: 1 PQ bus
    #      2 PV bus
    #      3 ref bus
    #      4 P bus
    #      5 Q bus
    #      6 isolated bus
    mpc["bus"] = [
        1   3   0   0   0   0   1   1.06    0   345 1   1.1 0.9;
        2   1   20  10  0   0   1   1   0   345 1   1.1 0.9;
        3   2   45  15  0   0   1   1   0   345 1   1.1 0.9;
        4   1   40  5   0   0   1   1   0   345 1   1.1 0.9;
        5   1   60  10  0   0   1   1   0   345 1   1.1 0.9;
        6   4   0   0   0   0   1   1   0   345 1   1.1 0.9;
        7   4   47   0   0   0   1   1   0   345 1   1.1 0.9
        ];

    ## generator data
    # bus Pg Qg Qmax Qmin Vg mBase status Pmax Pmin Pc1 Pc2 Qc1min Qc1max Qc2min Qc2max ramp_agc ramp_10 ramp_30 ramp_q_apf apf
    mpc["genAC"] = [
        1   0   0   500 -500    1.06    100 1   250 10;
  	    2   40  0   300 -300    1   100 1   300 10
    ];

    ## branch data
    #Branchtype: 1 for AC, 2 for DC,3 for converter
    # fbus tbus r x b rateA rateB rateC ratio angle status Branchtype P_conv Q_conv etaCrec etaCinv tanphi M
    mpc["branch"] = [
        1   2   0.02    0.06    0.06    100 100 100 0   0   1   1   0   0;
        1   3   0.08    0.24    0.05    100 100 100 0   0   1   1   0   0;
        2   3   0.06    0.18    0.04    100 100 100 0   0   1   1   0   0;
        2   4   0.06    0.18    0.04    100 100 100 0   0   1   1   0   0;
        2   5   0.04    0.12    0.03    100 100 100 0   0   1   1   0   0;
        3   4   0.01    0.03    0.02    100 100 100 0   0   1   1   0   0;
        4   5   0.08    0.24    0.05    100 100 100 0   0   1   1   0   0;
        6   7   0.052   0.0     0.0     100 100 100 0   0   1   2   0   0;
        6   2   0.00158  0.275  0.0     100 100 100 0   0   1   3   -60  -40;
        7   3   0.00158  0.275  0.0     100 100 100 0   0   1   3   0   0;
    ];
    return mpc
end