function bustypes(bus::Matrix{Float64}, gen::Matrix{Float64})
    # constants
    (PQ, PV, REF,Pdc,Vdc,NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA,
    VM,VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN) = idx_bus();
    (GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, 
        MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, 
        QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF) = idx_gen();

    # get generator status
    nb = size(bus, 1)
    ng = size(gen, 1)
    Cg = sparse(gen[:, GEN_BUS], 1:ng, gen[:, GEN_STATUS] .> 0, nb, ng)  # gen connection matrix
    bus_gen_status = Cg * ones(ng, 1)  # number of generators at each bus that are ON
    
    # form index lists for slack, PV, and PQ buses
    bus_gen_status = vec(bus_gen_status)
    map!(x -> x != 0 ? true : x, bus_gen_status, bus_gen_status)
    ref = findall(bus[:, BUS_TYPE] .== REF .* bus_gen_status )  # reference bus index
    #pv  = findall(bus[:, BUS_TYPE] .== PV  .* bus_gen_status )  # PV bus indices
    pv  = findall(bus[:, BUS_TYPE] .== PV)  # PV bus indices
    m = bus[:, BUS_TYPE] .== PQ;
    n = bus_gen_status;
    x = ones(size(n));
    y = x-n;
    s = xor(m,y)
    s = Bool.(s)
    s = vec(s)
    #if use ac powerfolw,change this
    #pq  = findall(s)  # PQ bus indices
    pq= findall(bus[:, BUS_TYPE] .== PQ)
    # pick a new reference bus if for some reason there is none (may have been shut down)
    if isempty(ref)
        ref = pv[1]  # use the first PV bus
        pv = pv[2:end]  # delete it from PV list
    end

    return ref, pv, pq
end

function xor(m,n)
    l = length(m);
    s = ones(1,l);
    for i = 1:l
        if(m[i]==0&&n[i]==0)
            s[i]=0;
        else
            s[i]=1;
        end
    end
    return s
end