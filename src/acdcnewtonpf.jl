function acdcnewtonpf(Ybus,Gbus,Sbusac,Sbusdc,pv,pq,pdc,vdc,busAC,busDC,converter,baseMVA,V0dc,V0,alg="NR",tol0=1e-8,max_it0=100)
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

    tol = tol0
    max_it = max_it0
    lin_solver = Char[]
    # Initialize
    converged = false
    i = 0
    V_ac = V0
    #AC part initialize
    Va_ac = angle.(V_ac)
    Vm_ac = abs.(V_ac)
    Vm_dc=V0dc
    Vn_f=ones(size(converter,1))
    Va_f=zeros(size(converter,1))
    Vn_c=0.9.*ones(size(converter,1))
    Va_c=zeros(size(converter,1))
    M=0.9.*ones(size(converter,1))
     # Set up indexing for updating variables
    npv = length(pv)
    npq = length(pq)
    npdc = length(pdc)
    nvdc = length(vdc)
    a=collect(1:npq);
    dict_pq = Dict(zip(pq, a))
    a=collect(1:npv);
    dict_pv = Dict(zip(pv, a))
    a=collect(1:npdc);
    dict_pdc = Dict(zip(pdc, a))
    index_converterpq=Int.(converter[findall(busAC[Int.(converter[:,BUSAC_I]),BUS_TYPE].==PQ),BUSAC_I])
    index_converterpv=Int.(converter[findall(busAC[Int.(converter[:,BUSAC_I]),BUS_TYPE].==PV),BUSAC_I])
    index_converterpdc=Int.(converter[findall(busDC[Int.(converter[:,BUSDC_I]),TYPE_DC].==Pdc),BUSDC_I])
    # a=collect(1:size(findall(busAC[Int.(converter[:,BUSAC_I]),BUS_TYPE].==PQ),1));
    # dict_converterpq=Dict(zip(converter[findall(busAC[Int.(converter[:,BUSAC_I]),BUS_TYPE].==PQ),BUSAC_I], a))
    # a=collect(1:size(findall(busAC[Int.(converter[:,BUSAC_I]),BUS_TYPE].==PV),1));
    # dict_converterpv=Dict(zip(converter[findall(busAC[Int.(converter[:,BUSAC_I]),BUS_TYPE].==PV),BUSAC_I], a))
    j1 = 1; j2 = npq; # j1:j2 - V angle of pq buses
    j3 = j2 + 1; j4 = j2 + npq; # j3:j4 - V angle of pq buses
    j5 = j4 + 1; j6 = j4 + npv; # j5:j6 - V mag of pq buses
    j7 = j6 + 1; j8 = j6 + npdc; # j7:j8 - Pdc of dc buses
    j9 = j8 + 1; j10 = j8 + size(converter,1); # j9:j10 - Vaf of converter
    # t11 = j10 + 1; t12 = j10 + size(converter,1); # j11:j12 - Vnf of converter
    t11=j10.+findall(busAC[Int.(converter[:,BUSAC_I]),BUS_TYPE].==PQ)#Vnf of converter
    # t13 = t12 + 1; t14 = t12 + size(converter,1); # j13:j14 - Vac of converter
    # t15 = t14 + 1; t16 = t14 + size(converter,1); # j15:j16 - Vnc of converter
    # t17 = t16 + 1; t18 = t16 + size(converter,1); # j17:j18 - M of converter
    t12=j10+length(findall(busAC[Int.(converter[:,BUSAC_I]),BUS_TYPE].==PQ))+1;t13=j10+length(findall(busAC[Int.(converter[:,BUSAC_I]),BUS_TYPE].==PQ))+ size(converter,1)#Vac of converter
    t14=t13+1;t15=t13+size(converter,1)#Vnc of converter
    t16=t15+1;t17=t15+size(converter,1)#M of converter


    #extract the resistance and reactance of the converter
    rtf=converter[:,RTF]
    xtf=converter[:,XTF]
    bf=converter[:,BF]
    rc=converter[:,RC]
    xc=converter[:,XC]
    Gtf,Btf,Gc,Bc = PowerFlow.makereactorvsc(rtf,xtf,rc,xc)
    #preoperation 
    Cc_ac=zeros(size(busAC,1),size(converter,1))
    for i in 1:size(converter,1)
        Cc_ac[Int(converter[i,BUSAC_I]),Int(converter[i,BUSDC_I])]=1
    end
    Cc_dc=zeros(size(busDC,1),size(converter,1))
    for i in 1:size(converter,1)
        Cc_dc[Int(converter[i,BUSDC_I]),Int(converter[i,BUSDC_I])]=1
    end
    a=converter[:,LOSSA]
    b=converter[:,LOSSB]
    c=converter[:,LOSSCREC]
    #define count
    count=size(converter,1)
    index=zeros(count,5)
    stat=1;
    for i in 1:count
        index[i,1]=stat
        index[i,2]=stat+count
        index[i,3]=stat+2*count
        index[i,4]=stat+3*count
        index[i,5]=stat+4*count
        stat+=1
    end
    # Evaluate F(x0)
    #calculate VSC mismatch vector
    Vns=Vm_ac[Int.(converter[:,BUSAC_I])]
    Vas=Va_ac[Int.(converter[:,BUSAC_I])]
    Vndc=Vm_dc[Int.(converter[:,BUSDC_I])]
    Ps=-Vns.^2 .*Gtf .+Vns .*Vn_f.*(Gtf.*cos.(Vas.-Va_f).+Btf.*sin.(Vas.-Va_f))
    Qs=Vns.^2 .*Btf .+Vns .*Vn_f.*(Gtf.*sin.(Vas.-Va_f).-Btf.*cos.(Vas.-Va_f))
    misPs=Ps-converter[:,P_G]/baseMVA
    misQs=Qs-converter[:,Q_G]/baseMVA
    # Pf=Vn_f.*Vns.*(-Gtf.*cos.(Va_f.-Vas).-Btf.*sin.(Va_f.-Vas)) .+Vn_f.*Vn_c.*(-Gc.*cos.(Va_f.-Va_c).-Bc.*sin.(Va_f.-Va_c))
    Pf=Vn_f.^2 .*Gtf .-Vn_f.*Vns.*(Gtf.*cos.(Va_f.-Vas).+Btf.*sin.(Va_f.-Vas)) +Vn_f.^2 .*Gc.+Vn_f.*Vn_c.*(Gc.*cos.(Va_f.-Va_c)+Bc.*sin.(Va_f.-Va_c))
    # Qf=Vn_f.*Vns.*(-Gtf.*sin.(Va_f.-Vas).+Btf.*cos.(Va_f.-Vas)) .+Vn_f.*Vn_c.*(-Gc.*sin.(Va_f.-Va_c)+Bc.*cos.(Va_f.-Va_c)).-Vn_f.^2 .*bf
    Qf=-Vn_f.^2 .*Btf .-Vn_f.*Vns.*(Gtf.*sin.(Va_f.-Vas).-Btf.*cos.(Va_f.-Vas))-Vn_f.^2 .*Bc.-Vn_f.*Vn_c.*(Gc.*sin.(Va_f.-Va_c)-Bc.*cos.(Va_f.-Va_c))-Vn_f.^2 .*bf
    misPf=Pf
    misQf=Qf
    misM=sqrt(2).*Vn_c .* inv.(sqrt(3).*Vndc)-M
    Pc=Vn_c.^2 .*Gc .-Vn_c .*Vn_f.*(Gc.*cos.(Va_c.-Va_f).+Bc.*sin.(Va_c.-Va_f))
    Qc=-Vn_c.^2 .*Bc .-Vn_c .*Vn_f.*(Gc.*sin.(Va_c.-Va_f).-Bc.*cos.(Va_c.-Va_f))
    Ic=sqrt.(Pc.^2 .+Qc.^2).*inv.(Vn_c)
    # Ic=(Vn_c.^2 +Vn_f.^2) .*(Gc.^2 .+Bc.^2)-2 .*Vn_c .*Vn_f.*cos.(Va_c-Va_f).*(Gc.^2 .+Bc.^2)
    Pcdc=-Pc-(a.*Ic.^2 .+b.*Ic .+c)
    #calculate AC mismatch vector
    mis_ac = V_ac .* conj.(Ybus * V_ac) - Sbusac(Vm_ac) -Cc_ac*(Ps+1im*Qs)
    mis_dc = Vm_dc .* conj.(Gbus * Vm_dc) - Sbusdc(Vm_dc) -Cc_dc*Pcdc
    F=[real(mis_ac[pq]);imag(mis_ac[pq]);real(mis_ac[pv]);real(mis_dc[pdc]);misPs;misQs;misPf;misQf[findall(busAC[Int.(converter[:,BUSAC_I]),BUS_TYPE].==PQ)];misM]
    normF = PowerFlow.norm(F, Inf)
    if normF < tol
        converged = true
    end
    while (!converged && i < max_it)
        # Update iteration counter
        i += 1
        dSbus_dVa_ac, dSbus_dVm_ac = PowerFlow.dSbus_dV(Ybus, V_ac)
        j11 = real(dSbus_dVa_ac[pq, pq])
        j12 = real(dSbus_dVm_ac[pq, pq])
        j13 = real(dSbus_dVa_ac[pq, pv])
        j21 = imag(dSbus_dVa_ac[pq, pq])
        j22 = imag(dSbus_dVm_ac[pq, pq])
        j23 = imag(dSbus_dVa_ac[pq, pv])
        j31 = real(dSbus_dVa_ac[pv, pq])
        j32 = real(dSbus_dVm_ac[pv, pq])
        j33 = real(dSbus_dVa_ac[pv, pv])
        J=[j11 j12 j13;j21 j22 j23;j31 j32 j33]
        dSbus_dVm_dc = PowerFlow.dSbus_dVdc(Gbus, Vm_dc)
        P=dSbus_dVm_dc[pdc,pdc]
        dP_dVm_matrix,dP_dVvsc_ac_matrix,dPs_dVas,dPs_dVns,dQs_dVas,dQs_dVns,dPs_dVa_f,dPs_dVm_f,dQs_dVa_f,dQs_dVm_f=PowerFlow.makevsc_acjacob(Vns,Vas,Vn_f,Va_f,Gtf,Btf,busAC,converter,npq,npv,dict_pq,dict_pv,index,index_converterpq,index_converterpv)
        Jacobianac=[J+dP_dVm_matrix zeros(npq*2+npv,npdc) dP_dVvsc_ac_matrix]
        dP_dVvsc_dc_matrix=PowerFlow.makevsc_dcjacob(Vn_c,Va_c,Vn_f,Va_f,Gc,Bc,busAC,busDC,converter,Pc,a,b,Ic,pdc,dict_pdc,index,index_converterpdc)
        Jacobiandc=[zeros(npdc,npq*2+npv) P dP_dVvsc_dc_matrix]
        Jacobianvsc=PowerFlow.makevsc_vscjacobi(Vn_f,Va_f,Vns,Vas,Vn_c,Va_c,Gtf,Btf,Gc,Bc,bf,busAC,busDC,converter,Vndc,npq,npv,npdc,dPs_dVas,dPs_dVns,dQs_dVas,dQs_dVns,dPs_dVa_f,dPs_dVm_f,dQs_dVa_f,dQs_dVm_f,dict_pq,dict_pv,dict_pdc,index,index_converterpq,index_converterpv,index_converterpdc)
        Jacobian=sparse([Jacobianac;Jacobiandc;Jacobianvsc])
        @time dx=junlinsolve(Jacobian,-F)
        # @time dx=-Jacobian\F
        # Update voltage
        Va_ac[pq] += dx[j1:j2]
        Vm_ac[pq] += dx[j3:j4]
        if npv > 0
            Va_ac[pv] += dx[j5:j6]
        end
        Vm_dc[pdc] += dx[j7:j8]
        Va_f += dx[j9:j10]
        Vn_f[findall(busAC[Int.(converter[:,BUSAC_I]),BUS_TYPE].==PQ)] += dx[t11]
        Va_c += dx[t12:t13]
        Vn_c += dx[t14:t15]
        M += dx[t16:t17]

        V_ac = Vm_ac .* exp.(1im * Va_ac)
        Vns=Vm_ac[Int.(converter[:,BUSAC_I])]
        Vas=Va_ac[Int.(converter[:,BUSAC_I])]
        Vndc=Vm_dc[Int.(converter[:,BUSDC_I])]
        Ps=-Vns.^2 .*Gtf .+Vns .*Vn_f.*(Gtf.*cos.(Vas.-Va_f).+Btf.*sin.(Vas.-Va_f))
        Qs=Vns.^2 .*Btf .+Vns .*Vn_f.*(Gtf.*sin.(Vas.-Va_f).-Btf.*cos.(Vas.-Va_f))
        misPs=Ps-converter[:,P_G]/baseMVA
        misQs=Qs-converter[:,Q_G]/baseMVA
        # Pf=Vn_f.*Vns.*(-Gtf.*cos.(Va_f.-Vas).-Btf.*sin.(Va_f.-Vas)) .+Vn_f.*Vn_c.*(-Gc.*cos.(Va_f.-Va_c).-Bc.*sin.(Va_f.-Va_c))
        Pf=Vn_f.^2 .*Gtf .-Vn_f.*Vns.*(Gtf.*cos.(Va_f.-Vas).+Btf.*sin.(Va_f.-Vas)) +Vn_f.^2 .*Gc.+Vn_f.*Vn_c.*(Gc.*cos.(Va_f.-Va_c)+Bc.*sin.(Va_f.-Va_c))
        # Qf=Vn_f.*Vns.*(-Gtf.*sin.(Va_f.-Vas).+Btf.*cos.(Va_f.-Vas)) .+Vn_f.*Vn_c.*(-Gc.*sin.(Va_f.-Va_c)+Bc.*cos.(Va_f.-Va_c)).-Vn_f.^2 .*bf
        Qf=-Vn_f.^2 .*Btf .-Vn_f.*Vns.*(Gtf.*sin.(Va_f.-Vas).-Btf.*cos.(Va_f.-Vas))-Vn_f.^2 .*Bc.-Vn_f.*Vn_c.*(Gc.*sin.(Va_f.-Va_c)-Bc.*cos.(Va_f.-Va_c))-Vn_f.^2 .*bf
        misPf=Pf
        misQf=Qf
        misM=sqrt(2).*Vn_c .* inv.(sqrt(3).*Vndc)-M
        Pc=Vn_c.^2 .*Gc .-Vn_c .*Vn_f.*(Gc.*cos.(Va_c.-Va_f).+Bc.*sin.(Va_c.-Va_f))
        Qc=-Vn_c.^2 .*Bc .-Vn_c .*Vn_f.*(Gc.*sin.(Va_c.-Va_f).-Bc.*cos.(Va_c.-Va_f))
        Ic=sqrt.(Pc.^2 .+Qc.^2).*inv.(Vn_c)
        # Ic=(Vn_c.^2 +Vn_f.^2) .*(Gc.^2 .+Bc.^2)-2 .*Vn_c .*Vn_f.*cos.(Va_c-Va_f).*(Gc.^2 .+Bc.^2)
        Pcdc=-Pc-(a.*Ic.^2 .+b.*Ic .+c)
        #calculate AC mismatch vector
        mis_ac = V_ac .* conj.(Ybus * V_ac) - Sbusac(Vm_ac) -Cc_ac*(Ps+1im*Qs)
        mis_dc = Vm_dc .* conj.(Gbus * Vm_dc) - Sbusdc(Vm_dc) -Cc_dc*Pcdc
        F=[real(mis_ac[pq]);imag(mis_ac[pq]);real(mis_ac[pv]);real(mis_dc[pdc]);misPs;misQs;misPf;misQf[findall(busAC[Int.(converter[:,BUSAC_I]),BUS_TYPE].==PQ)];misM]
        normF = PowerFlow.norm(F, Inf)
        if normF < tol
            converged = true
        end
    end
    return V_ac,Vm_dc,converged,i
end