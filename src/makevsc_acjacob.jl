function makevsc_acjacob(Vns,Vas,Vn_f,Va_f,Gtf,Btf,busAC,converter,npq,npv,dict_pq,dict_pv,index,index_converterpq,index_converterpv)
    (BUSDC_I, BUSAC_I, P_G, Q_G, VTAR, RTF, XTF, BF, RC, XC, BASEKVAC, VMMAX, VMMIN,
     IMAX, STATUS, LOSSA, LOSSB, LOSSCREC, LOSSCINV)=PowerFlow.idx_conv()
     (PQ, PV, REF,NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA,
     VM,VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN) =PowerFlow.idx_ACbus()
    nc=size(converter,1)
    dPs_dVa_f=Vns.*Vn_f.*(Gtf.*sin.(Vas.-Va_f).-Btf.*cos.(Vas.-Va_f))
    dPs_dVm_f=Vns.*(Gtf.*cos.(Vas.-Va_f).+Btf.*sin.(Vas.-Va_f))
    dPs_dVns=-2*Vns.*Gtf .+Vn_f.*(Gtf.*cos.(Vas.-Va_f).+Btf.*sin.(Vas.-Va_f))
    dPs_dVas=Vns.*Vn_f.*(-Gtf.*sin.(Vas.-Va_f).+Btf.*cos.(Vas.-Va_f))
    dQs_dVa_f=Vns.*Vn_f.*(-Gtf.*cos.(Vas.-Va_f).-Btf.*sin.(Vas.-Va_f))
    dQs_dVm_f=Vns.*(Gtf.*sin.(Vas.-Va_f).-Btf.*cos.(Vas.-Va_f))
    dQs_dVns=2*Vns.*Btf .+Vn_f.*(Gtf.*sin.(Vas.-Va_f).-Btf.*cos.(Vas.-Va_f))
    dQs_dVas=Vns.*Vn_f.*(Gtf.*cos.(Vas.-Va_f).+Btf.*sin.(Vas.-Va_f))
    #TODO:将其转化为矩阵形式
    dP_dVm_matrix=zeros(npq*2+npv,npq*2+npv)
    index_Ps_Vas = map(k -> dict_pq[k],index_converterpq)
    index_Ps_Vas_pv = map(k -> dict_pv[k],index_converterpv)
    dP_dVm_matrix[CartesianIndex.(index_Ps_Vas,index_Ps_Vas)]=-dPs_dVas[findall(busAC[Int.(converter[:,BUSAC_I]),BUS_TYPE].==PQ)]
    dP_dVm_matrix[CartesianIndex.(index_Ps_Vas,index_Ps_Vas.+npq)]=-dPs_dVns[findall(busAC[Int.(converter[:,BUSAC_I]),BUS_TYPE].==PQ)]
    dP_dVm_matrix[CartesianIndex.(index_Ps_Vas.+npq,index_Ps_Vas)]=-dQs_dVas[findall(busAC[Int.(converter[:,BUSAC_I]),BUS_TYPE].==PQ)]
    dP_dVm_matrix[CartesianIndex.(index_Ps_Vas.+npq,index_Ps_Vas.+npq)]=-dQs_dVns[findall(busAC[Int.(converter[:,BUSAC_I]),BUS_TYPE].==PQ)]
    dP_dVm_matrix[CartesianIndex.(index_Ps_Vas_pv.+npq*2,index_Ps_Vas_pv.+npq*2)]=-dPs_dVas[findall(busAC[Int.(converter[:,BUSAC_I]),BUS_TYPE].==PV)]
    
    dP_dVvsc_ac_matrix=zeros(npq*2+npv,5*size(converter,1))
    dP_dVvsc_ac_matrix[CartesianIndex.(index_Ps_Vas,Int.(index[findall(busAC[Int.(converter[:,BUSAC_I]),BUS_TYPE].==PQ),1]))]=copy(-dPs_dVa_f[findall(busAC[Int.(converter[:,BUSAC_I]),BUS_TYPE].==PQ)])
    dP_dVvsc_ac_matrix[CartesianIndex.(index_Ps_Vas,Int.(index[findall(busAC[Int.(converter[:,BUSAC_I]),BUS_TYPE].==PQ),2]))]=copy(-dPs_dVm_f[findall(busAC[Int.(converter[:,BUSAC_I]),BUS_TYPE].==PQ)])
    dP_dVvsc_ac_matrix[CartesianIndex.(index_Ps_Vas.+npq,Int.(index[findall(busAC[Int.(converter[:,BUSAC_I]),BUS_TYPE].==PQ),1]))]=copy(-dQs_dVa_f[findall(busAC[Int.(converter[:,BUSAC_I]),BUS_TYPE].==PQ)])
    dP_dVvsc_ac_matrix[CartesianIndex.(index_Ps_Vas.+npq,Int.(index[findall(busAC[Int.(converter[:,BUSAC_I]),BUS_TYPE].==PQ),2]))]=copy(-dQs_dVm_f[findall(busAC[Int.(converter[:,BUSAC_I]),BUS_TYPE].==PQ)])
    dP_dVvsc_ac_matrix[CartesianIndex.(index_Ps_Vas_pv.+npq*2,Int.(index[findall(busAC[Int.(converter[:,BUSAC_I]),BUS_TYPE].==PV),1]))]=copy(-dPs_dVa_f[findall(busAC[Int.(converter[:,BUSAC_I]),BUS_TYPE].==PV)])
    dP_dVvsc_ac_matrix[CartesianIndex.(index_Ps_Vas_pv.+npq*2,Int.(index[findall(busAC[Int.(converter[:,BUSAC_I]),BUS_TYPE].==PV),2]))]=copy(-dPs_dVm_f[findall(busAC[Int.(converter[:,BUSAC_I]),BUS_TYPE].==PV)])
    if isempty(findall(busAC[Int.(converter[:,BUSAC_I]),BUS_TYPE].==PV))==false
        delete_cols=nc .+findall(busAC[Int.(converter[:,BUSAC_I]),BUS_TYPE].==PV)
        dP_dVvsc_ac_matrix = dP_dVvsc_ac_matrix[:,Not( delete_cols)]
    end
    return dP_dVm_matrix,dP_dVvsc_ac_matrix,dPs_dVas,dPs_dVns,dQs_dVas,dQs_dVns,dPs_dVa_f,dPs_dVm_f,dQs_dVa_f,dQs_dVm_f
end