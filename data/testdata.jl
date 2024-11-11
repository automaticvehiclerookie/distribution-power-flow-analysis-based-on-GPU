function testdata()
    mpc = Dict{String, Any}();
    mpc["version"] = "2";
    mpc["baseMVA"] = 100;
    #TODO:
    ## bus data
    # bus_i type Pd Qd Gs Bs area Vm Va baseKV zone Vmax Vmin
    mpc["busAC"] = [
        1   3   0   0   0   0   1   1.06    0   345 1   1.1 0.9;
        2   1   20  10  0   0   1   1   0   345 1   1.1 0.9;
        ];
    # busdc_i  typedc grid    Pdc     Vdc     basekVdc    Vdcmax  Vdcmin  Cdc 
    mpc["busDC"] = [
        1   1   1   10   1   345 1.1 0.9 0;
        2   1   1   10   0.9   345 1.1 0.9 0;
        ]

    ## generator data
    # bus Pg Qg Qmax Qmin Vg mBase status Pmax Pmin Pc1 Pc2 Qc1min Qc1max Qc2min Qc2max ramp_agc ramp_10 ramp_30 ramp_q_apf apf
    mpc["genAC"] = [
        1   0   0   500 -500    1.06    100 1   250 10;
  	    2   0  0   300 -300    1   100 1   300 10
    ];

    ## branch data
    # fbus tbus r x b rateA rateB rateC ratio angle status angmin angmax 
    mpc["branchAC"] = [
        1   2   0.02    0.06    0.06    100 100 100 0   0   1;
    ];
    # fbusdc  tbusdc  r      l    c   rateA   rateB   rateC   status
    mpc["branchDC"] = [
        1   2   0.052   0   0   100 100 100 1
    ];
    # busdc_i busac_i P_g   Q_g   Vtar    rtf     xtf     bf     rc     xc     basekVac    Vmmax   Vmmin   Imax    status   LossA LossB  LossCrec LossCinv 
    mpc["converter"]=[
        1   2   -60 -40 1   0.0015  0.1121  0.0887  0.0001  0.16428 345 1.1 0.9 1.2 1   0.001   0.0   0.0    4.371;
    ]
    return mpc
end
