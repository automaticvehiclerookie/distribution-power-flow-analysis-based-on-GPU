function makevsc_dcjacob(Vn_c,Va_c,Vn_f,Va_f,Gc,Bc,busAC,busDC,converter,Pc,a,b,Ic,pdc,dict_pdc,index,index_converterpdc)
    (BUSDC_I, BUSAC_I, P_G, Q_G, VTAR, RTF, XTF, BF, RC, XC, BASEKVAC, VMMAX, VMMIN,
     IMAX, STATUS, LOSSA, LOSSB, LOSSCREC, LOSSCINV)=PowerFlow.idx_conv()
     (PQ, PV, REF,NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA,
     VM,VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN) =PowerFlow.idx_ACbus()
    nc=size(converter,1)
    dPc_dVa_f=-Vn_c.*Vn_f.*(-Gc.*sin.(Va_c.-Va_f).-Bc.*cos.(Va_c.-Va_f))
    dPc_dVm_f=-Vn_c.*(Gc.*cos.(Va_c.-Va_f).+Bc.*sin.(Va_c.-Va_f))
    dPc_dVm_c=2 .*Vn_c.*Gc .-Vn_f.*(Gc.*cos.(Va_c.-Va_f).+Bc.*sin.(Va_c.-Va_f))
    dPc_dVa_c=-Vn_c.*Vn_f.*(-Gc.*sin.(Va_c.-Va_f).+Bc.*cos.(Va_c.-Va_f))
    dIc2_dVa_f=-2 .*Vn_c.*Vn_f.*sin.(Va_c.-Va_f).*(Gc.^2+Bc.^2)
    dIc2_dVm_f=2 .*Vn_f.*(Gc.^2+Bc.^2)-2 .*Vn_c.*cos.(Va_c.-Va_f).*(Gc.^2+Bc.^2)
    dIc2_dVm_c=2 .*Vn_c.*(Gc.^2+Bc.^2)-2 .*Vn_f.*cos.(Va_c.-Va_f).*(Gc.^2+Bc.^2)
    # dIc_dVa_c=-inv.(Vn_c.*(cos.(Va_c)).^2).*(dPc_dVa_c.*cos.(Va_c)+Pc.*sin.(Va_c))
    # dIc_dVm_c=inv.(Vn_c.^2 .*cos.(Va_c)).*(dPc_dVm_c.*Vn_c-Pc)
    
    
    # dPcdc_dVa_f=-dPc_dVa_f-2 .*a .*(Pc.*dPc_dVa_f).*inv.(Vn_c.^2 .*cos.(Va_c).^2)-b.*dPc_dVa_f.*inv.(Vn_c .*cos.(Va_c))
    # dPcdc_dVm_f=-dPc_dVm_f-2 .*a .*(Pc.*dPc_dVm_f).*inv.(Vn_c.^2 .*cos.(Va_c).^2)-b.*dPc_dVm_f.*inv.(Vn_c .*cos.(Va_c))

    # dPcdc_dVa_c=-dPcdc_dVa_f
    # dPcdc_dVm_c=-dPc_dVm_c-a.*(2 .*Pc.*dPc_dVm_c.*Vn_c.^2 -2 .*Vn_c.*Pc.^2).*inv.(cos.(Va_c).^2 .*Vn_c.^4)-b.*(dPc_dVm_c.*Vn_c-Pc).*inv.(Vn_c.^2 .*cos.(Va_c))
    dPcdc_dVa_f=-dPc_dVa_f-(a+b.*inv.(2 .*Ic)).*dIc2_dVa_f
    dPcdc_dVm_f=-dPc_dVm_f-(a+b.*inv.(2 .*Ic)).*dIc2_dVm_f#have been modified
    dPcdc_dVa_c=-dPcdc_dVa_f
    dPcdc_dVm_c=-dPc_dVm_c-(a+b.*inv.(2 .*Ic)).*dIc2_dVm_c#have been modified
    index_Pdc_vsc = map(k -> dict_pdc[k],index_converterpdc )
    dP_dVvsc_dc_matrix=zeros(size(busDC,1),5*size(converter,1))
    dP_dVvsc_dc_matrix[CartesianIndex.(index_Pdc_vsc[Int.(converter[:,BUSDC_I])], Int.(index[:,1]))].=dPcdc_dVa_f
    dP_dVvsc_dc_matrix[CartesianIndex.(index_Pdc_vsc[Int.(converter[:,BUSDC_I])],Int.(index[:,2]))].=dPcdc_dVm_f
    dP_dVvsc_dc_matrix[CartesianIndex.(index_Pdc_vsc[Int.(converter[:,BUSDC_I])],Int.(index[:,3]))].=dPcdc_dVa_c
    dP_dVvsc_dc_matrix[CartesianIndex.(index_Pdc_vsc[Int.(converter[:,BUSDC_I])],Int.(index[:,4]))].=dPcdc_dVm_c
    if isempty(findall(busAC[Int.(converter[:,BUSAC_I]),BUS_TYPE].==PV))==false
        delete_cols=nc .+findall(busAC[Int.(converter[:,BUSAC_I]),BUS_TYPE].==PV)
        dP_dVvsc_dc_matrix = dP_dVvsc_dc_matrix[:,Not( delete_cols)]
    end
    return dP_dVvsc_dc_matrix
end