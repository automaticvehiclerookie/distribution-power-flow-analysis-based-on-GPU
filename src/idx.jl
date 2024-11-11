function idx_ACbus()
    # define AC bus types
    PQ      = 1;
    PV      = 2;
    REF     = 3;
    NONE    = 4;

   
    # define the indices
    BUS_I       = 1;    ## bus number (1 to 29997)
    BUS_TYPE    = 2;    ## bus type (1 - PQ bus, 2 - PV bus, 3 - reference bus, 4 - isolated bus)
    PD          = 3;    ## Pd, real power demand (MW)
    QD          = 4;    ## Qd, reactive power demand (MVAr)
    GS          = 5;    ## Gs, shunt conductance (MW at V = 1.0 p.u.)
    BS          = 6;    ## Bs, shunt susceptance (MVAr at V = 1.0 p.u.)
    BUS_AREA    = 7;    ## area number, 1-100
    VM          = 8;    ## Vm, voltage magnitude (p.u.)
    VA          = 9;    ## Va, voltage angle (degrees)
    BASE_KV     = 10;   ## baseKV, base voltage (kV)
    ZONE        = 11;   ## zone, loss zone (1-999)
    VMAX        = 12;   ## maxVm, maximum voltage magnitude (p.u.)      (not in PTI format)
    VMIN        = 13;   ## minVm, minimum voltage magnitude (p.u.)      (not in PTI format)

    ## included in opf solution, not necessarily in input
    ## assume objective function has units, u
    LAM_P       = 14;   ## Lagrange multiplier on real power mismatch (u/MW)
    LAM_Q       = 15;   ## Lagrange multiplier on reactive power mismatch (u/MVAr)
    MU_VMAX     = 16;   ## Kuhn-Tucker multiplier on upper voltage limit (u/p.u.)
    MU_VMIN     = 17;   ## Kuhn-Tucker multiplier on lower voltage limit (u/p.u.)

return PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN
end
#-------------------------------------------------
function idx_DCbus()
    # define DC bus types
    Pdc = 1;
    Vdc = 2;
    isolated=3;
    # define the indices
    BUSDC_I=1; # bus number (1 to 29997)
    TYPE_DC=2; # type_dc, converter type
    GRID=3; # grid number
    PDC=4; # Pdc, DC power injection (MW)
    VDC=5; # Vdc, DC voltage (p.u.)
    BASEKVDC=6; # basekVdc, base DC voltage (kV)
    VDCMAX=7; # Vdcmax, maximum DC voltage (p.u.)
    VDCMIN=8; # Vdcmin, minimum DC voltage (p.u.)
    CDC=9; # cdc, voltage source converter loss coefficient (p.u. MW^2/MVA^2)
    return Pdc, Vdc, isolated, BUSDC_I, TYPE_DC, GRID, PDC, VDC, BASEKVDC, VDCMAX, VDCMIN, CDC
end
#-------------------------------------------------
## define the indices
function idx_ACbrch()
    F_BUS       = 1;    ## f, from bus number
    T_BUS       = 2;    ## t, to bus number
    BR_R        = 3;    ## r, resistance (p.u.)
    BR_X        = 4;    ## x, reactance (p.u.)
    BR_B        = 5;    ## b, total line charging susceptance (p.u.)
    RATE_A      = 6;    ## rateA, MVA rating A (long term rating)
    RATE_B      = 7;    ## rateB, MVA rating B (short term rating)
    RATE_C      = 8;    ## rateC, MVA rating C (emergency rating)
    TAP         = 9;    ## ratio, transformer off nominal turns ratio
    SHIFT       = 10;   ## angle, transformer phase shift angle (degrees)
    BR_STATUS   = 11;   ## initial branch status, 1 - in service, 0 - out of service

    ## included in power flow solution, not necessarily in input
    PF          = 12;   ## real power injected at "from" bus end (MW)       (not in PTI format)
    QF          = 13;   ## reactive power injected at "from" bus end (MVAr) (not in PTI format)
    PT          = 14;   ## real power injected at "to" bus end (MW)         (not in PTI format)
    QT          = 15;   ## reactive power injected at "to" bus end (MVAr)   (not in PTI format)

    ## included in opf solution, not necessarily in input
    ## assume objective function has units, u
    MU_SF       = 16;   ## Kuhn-Tucker multiplier on MVA limit at "from" bus (u/MVA)
    MU_ST       = 17;   ## Kuhn-Tucker multiplier on MVA limit at "to" bus (u/MVA)
    MU_ANGMIN   = 18;   ## Kuhn-Tucker multiplier lower angle difference limit (u/degree)
    MU_ANGMAX   = 19;   ## Kuhn-Tucker multiplier upper angle difference limit (u/degree)
    return F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, MU_ANGMIN, MU_ANGMAX
end
#-------------------------------------------------
function idx_DCbrch()
    F_BUSDC=1; # fbusdc, from bus number (1 to 29997)
    T_BUSDC=2; # tbusdc, to bus number (1 to 29997)
    R=3; # r, resistance (p.u.)
    L=4; # l, inductance (p.u.)
    C=5; # c, capacitance (p.u.)
    RATEA=6; # rateA, MVA rating A (long term rating)
    RATEB=7; # rateB, MVA rating B (short term rating)
    RATEC=8; # rateC, MVA rating C (emergency rating)
    STATUS=9; # initial branch status, 1 - in service, 0 - out of service
    return F_BUSDC, T_BUSDC, R, L, C, RATEA, RATEB, RATEC, STATUS
end
#-------------------------------------------------
function idx_conv()
    BUSDC_I=1; # bus number (1 to 29997)
    BUSAC_I=2; # bus number (1 to 29997)
    P_G=3; # P_g, DC power injection (MW)
    Q_G=4; # Q_g, DC power injection (MVAr)
    VTAR=5; # Vtar, target voltage (p.u.)
    RTF=6; # rtf, transformer off nominal turns ratio
    XTF=7; # xtf, transformer reactance (p.u.)
    BF=8; # bf, transformer susceptance (p.u.)
    RC=9; # rc, converter ac side resistance (p.u.)
    XC=10; # xc, converter ac side reactance (p.u.)
    BASEKVAC=11; # basekVac, base AC voltage (kV)
    VMMAX=12; # Vmmax, maximum converter voltage magnitude (p.u.)
    VMMIN=13; # Vmmin, minimum converter voltage magnitude (p.u.)
    IMAX=14; # Imax, maximum converter current (p.u.)
    STATUS=15; # status, converter status
    LOSSA=16; # LossA, constant loss coefficient
    LOSSB=17; # LossB, linear loss coefficient
    LOSSCREC=18; # LossCrec, recovery loss coefficient
    LOSSCINV=19; # LossCinv, inverter loss coefficient
    return BUSDC_I, BUSAC_I, P_G, Q_G, VTAR, RTF, XTF, BF, RC, XC, BASEKVAC, VMMAX, VMMIN, IMAX, STATUS, LOSSA, LOSSB, LOSSCREC, LOSSCINV
end
#-------------------------------------------------
function idx_gen()
    ## define the indices
    GEN_BUS     = 1;    ## bus number
    PG          = 2;    ## Pg, real power output (MW)
    QG          = 3;    ## Qg, reactive power output (MVAr)
    QMAX        = 4;    ## Qmax, maximum reactive power output at Pmin (MVAr)
    QMIN        = 5;    ## Qmin, minimum reactive power output at Pmin (MVAr)
    VG          = 6;    ## Vg, voltage magnitude setpoint (p.u.)
    MBASE       = 7;    ## mBase, total MVA base of this machine, defaults to baseMVA
    GEN_STATUS  = 8;    ## status, 1 - machine in service, 0 - machine out of service
    PMAX        = 9;    ## Pmax, maximum real power output (MW)
    PMIN        = 10;   ## Pmin, minimum real power output (MW)
    PC1         = 11;   ## Pc1, lower real power output of PQ capability curve (MW)
    PC2         = 12;   ## Pc2, upper real power output of PQ capability curve (MW)
    QC1MIN      = 13;   ## Qc1min, minimum reactive power output at Pc1 (MVAr)
    QC1MAX      = 14;   ## Qc1max, maximum reactive power output at Pc1 (MVAr)
    QC2MIN      = 15;   ## Qc2min, minimum reactive power output at Pc2 (MVAr)
    QC2MAX      = 16;   ## Qc2max, maximum reactive power output at Pc2 (MVAr)
    RAMP_AGC    = 17;   ## ramp rate for load following/AGC (MW/min)
    RAMP_10     = 18;   ## ramp rate for 10 minute reserves (MW)
    RAMP_30     = 19;   ## ramp rate for 30 minute reserves (MW)
    RAMP_Q      = 20;   ## ramp rate for reactive power (2 sec timescale) (MVAr/min)
    APF         = 21;   ## area participation factor

    ## included in opf solution, not necessarily in input
    ## assume objective function has units, u
    MU_PMAX     = 22;   ## Kuhn-Tucker multiplier on upper Pg limit (u/MW)
    MU_PMIN     = 23;   ## Kuhn-Tucker multiplier on lower Pg limit (u/MW)
    MU_QMAX     = 24;   ## Kuhn-Tucker multiplier on upper Qg limit (u/MVAr)
    MU_QMIN     = 25;   ## Kuhn-Tucker multiplier on lower Qg limit (u/MVAr)

    ## Note: When a generator's PQ capability curve is not simply a box and the
    ## upper Qg limit is binding, the multiplier on this constraint is split into
    ## it's P and Q components and combined with the appropriate MU_Pxxx and
    ## MU_Qxxx values. Likewise for the lower Q limits.
    return GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, PC1, PC2, QC1MIN, QC1MAX, QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF, MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN
end