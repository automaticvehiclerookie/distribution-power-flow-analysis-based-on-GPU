function makevsc_vscjacobi(Vn_f,Va_f,Vns,Vas,Vn_c,Va_c,Gtf,Btf,Gc,Bc,bf,busAC,busDC,converter,Vndc,npq,npv,npdc,dPs_dVas,dPs_dVns,dQs_dVas,dQs_dVns,dPs_dVa_f,dPs_dVm_f,dQs_dVa_f,dQs_dVm_f,dict_pq,dict_pv,dict_pdc,index,index_converterpq,index_converterpv,index_converterpdc)
    (BUSDC_I, BUSAC_I, P_G, Q_G, VTAR, RTF, XTF, BF, RC, XC, BASEKVAC, VMMAX, VMMIN,
     IMAX, STATUS, LOSSA, LOSSB, LOSSCREC, LOSSCINV)=PowerFlow.idx_conv()
     (PQ, PV, REF,NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA,
     VM,VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN) =PowerFlow.idx_ACbus();
     (Pdc, Vdc, isolated, BUSDC_I, TYPE_DC, GRID, PDC, VDC, BASEKVDC, VDCMAX, VDCMIN, CDC) =PowerFlow.idx_DCbus();
     nc=size(converter,1)
    dPf_dVa_f=-Vn_f.*Vns.*(-Gtf.*sin.(Va_f.-Vas).+Btf.*cos.(Va_f.-Vas)) .-Vn_f.*Vn_c.*(-Gc.*sin.(Va_f.-Va_c).+Bc.*cos.(Va_f.-Va_c))
    dPf_dVm_f=-Vns.*(Gtf.*cos.(Va_f.-Vas).+Btf.*sin.(Va_f.-Vas)) .-Vn_c.*(Gc.*cos.(Va_f.-Va_c).+Bc.*sin.(Va_f.-Va_c))+2*Vn_f.*(Gtf.+Gc)
    dPf_dVa_c=-Vn_f.*Vn_c.*(Gc.*sin.(Va_f.-Va_c).-Bc.*cos.(Va_f.-Va_f))
    dPf_dVm_c=Vn_f.*(-Gc.*cos.(Va_f.-Va_c).-Bc.*sin.(Va_f.-Va_c))
    dPf_dVms=-Vn_f.*(Gtf.*cos.(Va_f.-Vas).+Btf.*sin.(Va_f.-Vas))
    dPf_dVas=-Vn_f.*Vns.*(Gtf.*sin.(Va_f.-Vas).-Btf.*cos.(Va_f.-Vas))
    dQf_dVa_f=-Vn_f.*Vns.*(Gtf.*cos.(Va_f.-Vas).+Btf.*sin.(Va_f.-Vas)) .-Vn_f.*Vn_c.*(Gc.*cos.(Va_f.-Va_c).+Bc.*sin.(Va_f.-Va_c))
    dQf_dVm_f=-Vns.*(Gtf.*sin.(Va_f.-Vas).-Btf.*cos.(Va_f.-Vas)) .-Vn_c.*(Gc.*sin.(Va_f.-Va_c).-Bc.*cos.(Va_f.-Va_c)).-2*Vn_f.*(bf.+Btf.+Bc)
    dQf_dVa_c=-Vn_f.*Vn_c.*(-Gc.*cos.(Va_f.-Va_c).-Bc.*sin.(Va_f.-Va_c))
    dQf_dVm_c=-Vn_f.*(Gc.*sin.(Va_f.-Va_c).-Bc.*cos.(Va_f.-Va_c))
    dQf_dVms=-Vn_f.*(Gtf.*sin.(Va_f.-Vas).-Btf.*cos.(Va_f.-Vas))
    dQf_dVas=-Vn_f.*Vns.*(-Gtf.*cos.(Va_f.-Vas).-Btf.*sin.(Va_f.-Vas))

    dM_dVc=sqrt(2) .* inv.(sqrt(3).*Vndc)
    dM_dVdc=-sqrt(2).*Vn_c .* inv.(sqrt(3).*Vndc.^2)
    #TODO:将其转化为矩阵形式
    index_PQ = map(k -> dict_pq[k],index_converterpq )
    index_PV = map(k -> dict_pv[k],index_converterpv )
    index_PDC = map(k -> dict_pdc[k],index_converterpdc )

    dP_dVvsc_vsc_matrix=zeros(5*size(converter,1),2*npq+npv+npdc+5*size(converter,1))
    #dPs/dVas and dPs/dVns in pq buses
    dP_dVvsc_vsc_matrix[CartesianIndex.(Int.(index[findall(busAC[Int.(converter[:,BUSAC_I]),BUS_TYPE].==PQ),1]),index_PQ)]=dPs_dVas[findall(busAC[Int.(converter[:,BUSAC_I]),BUS_TYPE].==PQ)]
    dP_dVvsc_vsc_matrix[CartesianIndex.(Int.(index[findall(busAC[Int.(converter[:,BUSAC_I]),BUS_TYPE].==PQ),1]),index_PQ.+npq)]=dPs_dVns[findall(busAC[Int.(converter[:,BUSAC_I]),BUS_TYPE].==PQ)]
    dP_dVvsc_vsc_matrix[CartesianIndex.(Int.(index[findall(busAC[Int.(converter[:,BUSAC_I]),BUS_TYPE].==PQ),2]),index_PQ)]=dQs_dVas[findall(busAC[Int.(converter[:,BUSAC_I]),BUS_TYPE].==PQ)]
    dP_dVvsc_vsc_matrix[CartesianIndex.(Int.(index[findall(busAC[Int.(converter[:,BUSAC_I]),BUS_TYPE].==PQ),2]),index_PQ.+npq)]=dQs_dVns[findall(busAC[Int.(converter[:,BUSAC_I]),BUS_TYPE].==PQ)]

    #dPs/dVas and dPs/dVas in pv buses
    dP_dVvsc_vsc_matrix[CartesianIndex.(Int.(index[findall(busAC[Int.(converter[:,BUSAC_I]),BUS_TYPE].==PV),1]),index_PV.+2*npq)]=dPs_dVas[findall(busAC[Int.(converter[:,BUSAC_I]),BUS_TYPE].==PV)]

    #vsc part 
    dP_dVvsc_vsc_matrix[CartesianIndex.(Int.(index[findall(busDC[Int.(converter[:,BUSDC_I]),BUS_TYPE].==Pdc),1]),npq*2+npv+npdc.+Int.(index[findall(busDC[Int.(converter[:,BUSDC_I]),BUS_TYPE].==Pdc),1]))]=dPs_dVa_f[findall(busDC[Int.(converter[:,BUSDC_I]),BUS_TYPE].==Pdc)]
    dP_dVvsc_vsc_matrix[CartesianIndex.(Int.(index[findall(busDC[Int.(converter[:,BUSDC_I]),BUS_TYPE].==Pdc),1]),npq*2+npv+npdc.+Int.(index[findall(busDC[Int.(converter[:,BUSDC_I]),BUS_TYPE].==Pdc),2]))]=dPs_dVm_f[findall(busDC[Int.(converter[:,BUSDC_I]),BUS_TYPE].==Pdc)]
    dP_dVvsc_vsc_matrix[CartesianIndex.(Int.(index[findall(busDC[Int.(converter[:,BUSDC_I]),BUS_TYPE].==Pdc),2]),npq*2+npv+npdc.+Int.(index[findall(busDC[Int.(converter[:,BUSDC_I]),BUS_TYPE].==Pdc),1]))]=dQs_dVa_f[findall(busDC[Int.(converter[:,BUSDC_I]),BUS_TYPE].==Pdc)]
    dP_dVvsc_vsc_matrix[CartesianIndex.(Int.(index[findall(busDC[Int.(converter[:,BUSDC_I]),BUS_TYPE].==Pdc),2]),npq*2+npv+npdc.+Int.(index[findall(busDC[Int.(converter[:,BUSDC_I]),BUS_TYPE].==Pdc),2]))]=dQs_dVm_f[findall(busDC[Int.(converter[:,BUSDC_I]),BUS_TYPE].==Pdc)]

    dP_dVvsc_vsc_matrix[CartesianIndex.(Int.(index[findall(busAC[Int.(converter[:,BUSAC_I]),BUS_TYPE].==PQ),3]),index_PQ)]=dPf_dVas[findall(busAC[Int.(converter[:,BUSAC_I]),BUS_TYPE].==PQ)]
    dP_dVvsc_vsc_matrix[CartesianIndex.(Int.(index[findall(busAC[Int.(converter[:,BUSAC_I]),BUS_TYPE].==PQ),3]),index_PQ.+npq)]=dPf_dVms[findall(busAC[Int.(converter[:,BUSAC_I]),BUS_TYPE].==PQ)]
    dP_dVvsc_vsc_matrix[CartesianIndex.(Int.(index[findall(busAC[Int.(converter[:,BUSAC_I]),BUS_TYPE].==PV),3]),index_PV.+2*npq)]=dPf_dVas[findall(busAC[Int.(converter[:,BUSAC_I]),BUS_TYPE].==PV)]
    dP_dVvsc_vsc_matrix[CartesianIndex.(Int.(index[findall(busDC[Int.(converter[:,BUSDC_I]),BUS_TYPE].==Pdc),3]),npq*2+npv+npdc.+Int.(index[findall(busDC[Int.(converter[:,BUSDC_I]),BUS_TYPE].==Pdc),1]))]=dPf_dVa_f[findall(busDC[Int.(converter[:,BUSDC_I]),BUS_TYPE].==Pdc)]
    dP_dVvsc_vsc_matrix[CartesianIndex.(Int.(index[findall(busDC[Int.(converter[:,BUSDC_I]),BUS_TYPE].==Pdc),3]),npq*2+npv+npdc.+Int.(index[findall(busDC[Int.(converter[:,BUSDC_I]),BUS_TYPE].==Pdc),2]))]=dPf_dVm_f[findall(busDC[Int.(converter[:,BUSDC_I]),BUS_TYPE].==Pdc)]
    dP_dVvsc_vsc_matrix[CartesianIndex.(Int.(index[findall(busDC[Int.(converter[:,BUSDC_I]),BUS_TYPE].==Pdc),3]),npq*2+npv+npdc.+Int.(index[findall(busDC[Int.(converter[:,BUSDC_I]),BUS_TYPE].==Pdc),3]))]=dPf_dVa_c[findall(busDC[Int.(converter[:,BUSDC_I]),BUS_TYPE].==Pdc)]
    dP_dVvsc_vsc_matrix[CartesianIndex.(Int.(index[findall(busDC[Int.(converter[:,BUSDC_I]),BUS_TYPE].==Pdc),3]),npq*2+npv+npdc.+Int.(index[findall(busDC[Int.(converter[:,BUSDC_I]),BUS_TYPE].==Pdc),4]))]=dPf_dVm_c[findall(busDC[Int.(converter[:,BUSDC_I]),BUS_TYPE].==Pdc)]
    dP_dVvsc_vsc_matrix[CartesianIndex.(Int.(index[findall(busAC[Int.(converter[:,BUSAC_I]),BUS_TYPE].==PQ),4]),index_PQ)]=dQf_dVas[findall(busAC[Int.(converter[:,BUSAC_I]),BUS_TYPE].==PQ)]
    dP_dVvsc_vsc_matrix[CartesianIndex.(Int.(index[findall(busAC[Int.(converter[:,BUSAC_I]),BUS_TYPE].==PQ),4]),index_PQ.+npq)]=dQf_dVms[findall(busAC[Int.(converter[:,BUSAC_I]),BUS_TYPE].==PQ)]
    dP_dVvsc_vsc_matrix[CartesianIndex.(Int.(index[findall(busAC[Int.(converter[:,BUSAC_I]),BUS_TYPE].==PV),4]),index_PV.+2*npq)]=dQf_dVas[findall(busAC[Int.(converter[:,BUSAC_I]),BUS_TYPE].==PV)]
    dP_dVvsc_vsc_matrix[CartesianIndex.(Int.(index[findall(busDC[Int.(converter[:,BUSDC_I]),BUS_TYPE].==Pdc),4]),npq*2+npv+npdc.+Int.(index[findall(busDC[Int.(converter[:,BUSDC_I]),BUS_TYPE].==Pdc),1]))]=dQf_dVa_f[findall(busDC[Int.(converter[:,BUSDC_I]),BUS_TYPE].==Pdc)]
    dP_dVvsc_vsc_matrix[CartesianIndex.(Int.(index[findall(busDC[Int.(converter[:,BUSDC_I]),BUS_TYPE].==Pdc),4]),npq*2+npv+npdc.+Int.(index[findall(busDC[Int.(converter[:,BUSDC_I]),BUS_TYPE].==Pdc),2]))]=dQf_dVm_f[findall(busDC[Int.(converter[:,BUSDC_I]),BUS_TYPE].==Pdc)]
    dP_dVvsc_vsc_matrix[CartesianIndex.(Int.(index[findall(busDC[Int.(converter[:,BUSDC_I]),BUS_TYPE].==Pdc),4]),npq*2+npv+npdc.+Int.(index[findall(busDC[Int.(converter[:,BUSDC_I]),BUS_TYPE].==Pdc),3]))]=dQf_dVa_c[findall(busDC[Int.(converter[:,BUSDC_I]),BUS_TYPE].==Pdc)]
    dP_dVvsc_vsc_matrix[CartesianIndex.(Int.(index[findall(busDC[Int.(converter[:,BUSDC_I]),BUS_TYPE].==Pdc),4]),npq*2+npv+npdc.+Int.(index[findall(busDC[Int.(converter[:,BUSDC_I]),BUS_TYPE].==Pdc),4]))]=dQf_dVm_c[findall(busDC[Int.(converter[:,BUSDC_I]),BUS_TYPE].==Pdc)]
    dP_dVvsc_vsc_matrix[CartesianIndex.(Int.(index[findall(busDC[Int.(converter[:,BUSDC_I]),BUS_TYPE].==Pdc),5]),npq*2+npv.+index_PDC[Int.(converter[:,BUSDC_I])])]=dM_dVdc[findall(busDC[Int.(converter[:,BUSDC_I]),BUS_TYPE].==Pdc)]
    dP_dVvsc_vsc_matrix[CartesianIndex.(Int.(index[findall(busDC[Int.(converter[:,BUSDC_I]),BUS_TYPE].==Pdc),5]),npq*2+npv+npdc.+Int.(index[findall(busDC[Int.(converter[:,BUSDC_I]),BUS_TYPE].==Pdc),4]))]=dM_dVc[findall(busDC[Int.(converter[:,BUSDC_I]),BUS_TYPE].==Pdc)]
    dP_dVvsc_vsc_matrix[CartesianIndex.(Int.(index[findall(busDC[Int.(converter[:,BUSDC_I]),BUS_TYPE].==Pdc),5]),npq*2+npv+npdc.+Int.(index[findall(busDC[Int.(converter[:,BUSDC_I]),BUS_TYPE].==Pdc),5]))].=-1
    # dP_dVvsc_vsc_matrix[Int.(index[:,1]),index_PV.+2*npq]=copy(dPs_dVas)
    if isempty(findall(busAC[Int.(converter[:,BUSAC_I]),BUS_TYPE].==PV))==false
        delete_cols=npq*2+npv+npdc+nc .+findall(busAC[Int.(converter[:,BUSAC_I]),BUS_TYPE].==PV)
        dP_dVvsc_vsc_matrix = dP_dVvsc_vsc_matrix[:,Not( delete_cols)]
        delete_rows=nc .+findall(busAC[Int.(converter[:,BUSAC_I]),BUS_TYPE].==PV)
        dP_dVvsc_vsc_matrix = dP_dVvsc_vsc_matrix[Not( delete_rows),:]
    end
    return dP_dVvsc_vsc_matrix
end