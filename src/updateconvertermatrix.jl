function updateconvertermatrix(converter_matrixac,converter_matrixdc,converter,npv,npq,np,nv,pq,pv,pdc,vdc,d,c,slackac,slackdc,dict_dc,dict_pq,dict_pv)
    (PQ, PV, REF,Pdc,Vdc,NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA,
    VM,VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN) =PowerFlow.idx_bus();
    (F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C,TAP, SHIFT, BR_STATUS, 
    PF, QF, PT, QT, MU_SF, MU_ST, BRANCHMODE, PC, QC, ETACR,ETACI, PHI, M) = PowerFlow.idx_brch();
    static_matrixac=[0.0 0.0 0.0;1.0 0.0 0.0]
    static_matrixdc=[0.0 0.0 0.0]
    a=collect(1:size(vdc,1))
    dict_vdc = Dict(vdc[i] => a[i] for i in 1:length(a))
    for i in 1:size(converter,1)
        static_matrixac[1,2]=2*slackac[i]
        static_matrixac[1,3]=-2*slackdc[i]*converter[i,ETACI]
        static_matrixdc[1,2]=2*slackdc[i]
        static_matrixdc[1,3]=-2*slackac[i]*converter[i,ETACR]
        converter_matrixac[i]=copy(static_matrixac)
        converter_matrixdc[i]=copy(static_matrixdc)
    end
    jacobian_converter_matrix=sparse(zeros(2*npv+2*npq+np+3*size(converter,1),2*npv+2*npq+np+3*size(converter,1)))
    for i in 1:length(c)
        if(c[i] in pv)
            jacobian_converter_matrix[2*npq+getindex(dict_pv,c[i]),2*npv+2*npq+np+3*(i-1)+1:2*npv+2*npq+np+3*i]=converter_matrixac[i][1,:]
        else
            jacobian_converter_matrix[getindex(dict_pq,c[i]),2*npv+2*npq+np+3*(i-1)+1:2*npv+2*npq+np+3*i]=converter_matrixac[i][1,:]
            jacobian_converter_matrix[getindex(dict_pq,c[i])+npq,2*npv+2*npq+np+3*(i-1)+1:2*npv+2*npq+np+3*i]=converter_matrixac[i][2,:]
        end
    end
    # extracted_rows = jacobian_converter_matrix[pv, :]
    # remaining_rows = jacobian_converter_matrix[setdiff(1:size(jacobian_converter_matrix, 1), pv), :]
    # jacobian_converter_matrix= vcat(extracted_rows, remaining_rows)
    for i in 1:length(d)
        if(d[i] in pdc)
            jacobian_converter_matrix[2*npv+2*npq+getindex.(Ref(dict_dc), d[i]),2*npv+2*npq+np+3*(i-1)+1:2*npv+2*npq+np+3*i]=converter_matrixdc[i][1,:]
        else
            jacobian_converter_matrix[2*npv+2*npq+np+getindex.(Ref(dict_vdc), d[i]),2*npv+2*npq+np+3*(i-1)+1:2*npv+2*npq+np+3*i]=converter_matrixdc[i][1,:]
        end
    end
    #construct the dynamic part of Jacobian matrix

    dynamic_matrix=[1.0 0.0 0.0;0.0 0.0 0.0]
    dynamic_matrixpv=[0.0 0.0 0.0]
    for i in 1:size(converter,1)
        dynamic_matrix[1,2]=-2*slackac[i]*converter[i,PHI]
        dynamic_matrix[1,3]=2*slackdc[i]*converter[i,PHI]*converter[i,ETACI]
        dynamic_matrix[2,2]=slackdc[i]
        dynamic_matrix[2,3]=slackac[i]
        jacobian_converter_matrix[2*npv+npq*2+np+3*(i-1)+1:2*npv+npq*2+np+3*i-1,2*npv+2*npq+np+3*(i-1)+1:2*npv+2*npq+np+3*i]=dynamic_matrix
    end
    Vmatrix_converter=sparse(zeros(2*npv+2*npq+np+3*size(converter,1),2*npv+2*npq+np+3*size(converter,1)))
    for i in 1:size(converter,1)
        if(c[i] in pq)
            Vmatrix_converter[2*npv+2*npq+np+3*i,npq+getindex(dict_pq,c[i])]=-converter[i,M]
        else
            Vmatrix_converter[2*npv+2*npq+np+3*i,2*npq+npv+getindex(dict_pv,c[i])]=-converter[i,M]
        end
        if(d[i] in pdc)
            Vmatrix_converter[2*npv+2*npq+np+3*i,2*npv+npq*2+getindex.(Ref(dict_dc), d[i])]=1
        else
            Vmatrix_converter[2*npv+2*npq+np+3*i,2*npv+npq*2+np+getindex.(Ref(dict_vdc), d[i])]=1
        end
    end
    return jacobian_converter_matrix,Vmatrix_converter
end