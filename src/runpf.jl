"""
   Main function to call the AC power flow function
    Input: case file
    Output: results of the power flow as a dictionary
    Example:
    bus, gen, branch = runpf(casefile)
"""

# Step 1: Define the runpf function
function runpf(mpc, opt::Dict{String})
    # Step 2.1: Define the data structures
    # Define named indices into bus, gen, branch matrices
    (PQ, PV, REF,Pdc,Vdc,NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA,
    VM,VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN) = PowerFlow.idx_bus();
    (F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C,TAP, SHIFT,
       BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST,ANGMIN, ANGMAX, MVSC1,
        MVSC2, BRANCHMODE, ETACR,ETACI, PHI, MU_ANGMIN, MU_ANGMAX) = PowerFlow.idx_brch();
    (GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, 
        MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, 
        QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF) = PowerFlow.idx_gen();
    
    # options
    qlim = opt["PF"]["ENFORCE_Q_LIMS"];         ## enforce Q limits on gens?
    dc = opt["PF"]["DC"];  ## use DC formulation?  

    # Step 2.2: Convert the data into the required format
    baseMVA = mpc["baseMVA"];
    bus =  mpc["bus"];
    gen = mpc["gen"];
    branch = mpc["branch"];
    # convert the external data to internal data 
    (bus, gen, branch) = PowerFlow.ext2int(bus, gen, branch);
    ## get bus index lists of each type of bus
    (ref, pv, pq) = PowerFlow.bustypes(bus, gen);
    ## generator info
    on = findall(gen[:, GEN_STATUS] .> 0)  # which generators are on?
    gbus = gen[on, GEN_BUS]  # what buses are they at?
    # Step 2.3: Run the power flow
    ##-----  run the power flow  ----- 
    alg = opt["PF"]["PF_ALG"];
    its = 0;            ## total iterations
    if(dc==1)
        # @printf(" -- DC Power Flow (%s)\n", solver);
        #initial state
        Va0 = bus[:, VA] * (pi/180);
        #build B matrices and phase shift injections
        (B,Bf,Pbusinj,Pfinj) = makeBdc(baseMVA, bus, branch);
        # compute complex bus power injections (generation - load)
        # adjusted for phase shifters and real shunts
        Pbus = real(makeSbus(baseMVA, bus, gen, zeros(size(bus, 1)))) - Pbusinj - bus[:, GS] / baseMVA
        #"run" the power flow
        Va = dcpf(B, Pbus, Va0, ref,pv,pq);
        its = 1;
        #update data matrices with solution
        (rows, cols) = size(branch)

        # Determine the number of columns to add
        cols_to_add = 17 - cols
        # If cols_to_add is greater than 0, add more columns
        if cols_to_add > 0
            # Add columns filled with zeros
            branch = [branch zeros(rows, cols_to_add)]
        end
            branch[:, [QF, QT]] = zeros(size(branch, 1), 2);
            Va = Va[1]  # Extract the vector from the tuple
            branch[:, PF] = (Bf * Va + Pfinj) * baseMVA
            branch[:, PT] = -branch[:, PF];
            bus[:, VM] = ones(size(bus, 1), 1);
            bus[:, VA] = Va * (180/pi);
            #update Pg for slack generator (1st gen at ref bus)
            # (note: other gens at ref bus are accounted for in Pbus)
            #      Pg = Pinj + Pload + Gs
            #      newPg = oldPg + newPinj - oldPinj
            refgen = zeros(size(ref));
            for k = 1:length(ref)
                temp = findall(x -> x == ref[k], gbus)
                refgen[k] = on[temp[1]]
            end
            gen[Int.(refgen), PG] = gen[Int.(refgen), PG] + (B[Int.(ref), :] * Va - Pbus[Int.(ref)]) * baseMVA
    else
        # @printf(" -- AC Power Flow (%s)\n", solver);
        ## initialize
        V0  = bus[:, VM] .* exp.(1im * pi/180 * bus[:, VA])
        vcb = ones(size(V0));           ## create mask of voltage-controlled buses
        vcb[pq] .= 0;                    ## exclude PQ buses
        gbus=Int.(gbus);
        k = findall(Bool.(vcb[gbus]));            ## in-service gens at v-c buses
        V0[gbus[k]] = gen[on[k], VG] ./ abs.(V0[gbus[k]]).* V0[gbus[k]];
        if(qlim>0)
            ref0 = ref;                         ## save index and angle of
            Varef0 = bus[ref0, VA];             ##   original reference bus(es)
            limited=[];                         ## list of indices of gens @ Q lims
            fixedQg = zeros(size(gen, 1), 1);    ## Qg of gens at Q limits
        end
        ## build admittance matrices
        (Ybus, Yf, Yt) = PowerFlow.makeYbus(baseMVA, bus, branch);
        repeat=1;
        while (repeat>0)
            ## function for computing V dependent complex bus power injections
                ## (generation - load)
            Sbus = Vm -> PowerFlow.makeSbus(baseMVA, bus, gen, Vm);
            ## run the power flow is NR is selected
            if alg == "NR"
                V, success, iterations = newtonpf(Ybus, Sbus, V0, ref, pv, pq, opt["PF"]["PF_TOL"], opt["PF"]["PF_MAX_IT"], opt["PF"]["NR_ALG"]);
                its += iterations;
            end
            ## update data matrices with solution
            bus, gen, branch = pfsoln(baseMVA, bus, gen, branch, Ybus, Yf, Yt, V, ref, pv, pq);
            if(qlim>0&&success)
                #find gens with violated Q constraints
                mx = findall((gen[:, GEN_STATUS] .> 0) .& (gen[:, QG] .> gen[:, QMAX] .+ opt["OPF"]["OPF_VIOLATION"]))
                mn = findall((gen[:, GEN_STATUS] .> 0) .& (gen[:, QG] .< gen[:, QMIN] .- opt["OPF"]["OPF_VIOLATION"]))
                if ~isempty(mx) || ~isempty(mn) #we have some Q limit violations
                    #first check for INFEASIBILITY
                    infeas = union(mx, mn)
                    #union of scalars is a row vector
                    remaining = findall((gen[:, GEN_STATUS] .> 0) .& ((bus[Int.(gen[:, GEN_BUS]), BUS_TYPE] .== PV) .| (bus[Int.(gen[:, GEN_BUS]), BUS_TYPE] .== REF)))
                    if length(infeas) == length(remaining) && all(infeas == remaining) && (isempty(mx) || isempty(mn))
                        #all remaining PV/REF gens are violating AND all are
                        #violating same limit (all violating Qmin or all Qmax)
                        success = 0
                        break;
                    end   
                    #one at a time?
                    if qlim == 2    # fix largest violation, ignore the rest
                        violations = vcat(gen[mx, QG] .- gen[mx, QMAX], gen[mn, QMIN] .- gen[mn, QG])
                        k = argmax(violations)[1]
                        if k > length(mx)
                            mn = mn[k-length(mx)]
                            mx = Int[]
                        else
                            mx = mx[k]
                            mn = Int[]
                        end
                    end
                    #save corresponding limit values
                    fixedQg[mx,1] = gen[mx, QMAX];
                    fixedQg[mn,1] = gen[mn, QMIN]; 
                    mx = [mx;mn];

                    #convert to PQ bus
                    gen[mx, QG] = fixedQg[mx]; #set Qg to binding limit
                    bus[Int.(gen[mx, GEN_BUS]), BUS_TYPE] .= PQ
                    #update bus index lists of each type of bus
                    ref_temp = ref;
                    (ref, pv, pq) = bustypes(bus, gen);
                    # previous line can modify lists to select new REF bus
                    # if there was none, so we should update bus with these
                    # just to keep them consistent
                    if ref != ref_temp
                        if isa(bus[ref, BUS_TYPE],Vector)
                            bus[ref, BUS_TYPE][1] = REF;
                        else
                            bus[ref, BUS_TYPE] = REF;
                        end
                        bus[pv, BUS_TYPE] .= PV
                    end
                    limited = [limited; mx];
                else
                    repeat = 0;
                end
            else
                repeat = 0;
            end
        end
        if qlim>0 && !isempty(limited)
            if ref != ref0
            ## adjust voltage angles to make original ref bus correct
                bus[:, VA] .= bus[:, VA] .- bus[ref0, VA] .+ Varef0
            end
        end
    end
        ## Step 2.4: -----  output results  -----
        ## convert back to original bus numbering & print results
        mpc["bus"] = bus
        mpc["gen"] = gen
        mpc["branch"] = branch
        if(dc==0)
            mpc["success"] = success
        end
        mpc["iterations"] = its
        #comment: print the final voltage magnitudes
        # print(mpc["bus"],its);
    
    return mpc
end


mutable struct MPC
    version::String
    baseMVA::Int
    bus::Array{Float64,2}
    gen::Array{Float64,2}
    branch::Array{Float64,2}
    gencost::Array{Float64,2}
end

function runpf(mpc::MPC)
    # Step 2.1: Define the data structures
    # Define named indices into bus, gen, branch matrices
    (PQ, PV, REF, IAC,IDC,Pdc,Vdc,NONE, BUS_I, BUS_TYPE, PD,
     QD, GS, BS, BUS_AREA, VM,VA, BASE_KV, ZONE, VMAX, VMIN,
      LAM_P, LAM_Q, MU_VMAX, MU_VMIN) = idx_bus();
    (F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, 
        TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, 
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX) = idx_brch();
    (GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, 
        MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, 
        QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF) = idx_gen();
    # Step 2.2: Convert the data into the required format
    baseMVA = mpc["baseMVA"];
    bus =  mpc["bus"];
    gen = mpc["gen"];
    branch = mpc["branch"];
    #add zero columns to branch for flows if needed
    if size(branch, 2) < QT
        branch = [branch zeros(size(branch, 1), QT-size(branch, 2))]
    end
    # convert the external data to internal data 
    (bus, gen, branch) = ext2int(bus, gen, branch);
    ## get bus index lists of each type of bus
    (ref, pv, pq) = bustypes(bus, gen);
    ## generator info
    on = findall(gen[:, GEN_STATUS] .> 0)  # which generators are on?
    gbus = gen[on, GEN_BUS]  # what buses are they at?
    # Step 2.3: Run the power flow
    ##-----  run the power flow  ----- 
    alg="NR";
    solver = "Newton";
    its = 0;            ## total iterations
    # @printf(" -- AC Power Flow (%s)\n", solver);
    ## initialize
    V0  = bus[:, VM] .* exp.(1im * pi/180 * bus[:, VA])
    vcb = ones(size(V0));           ## create mask of voltage-controlled buses
    vcb[pq] .= 0;                    ## exclude PQ buses
    gbus=Int.(gbus);
    k = findall(Bool.(vcb[gbus]));            ## in-service gens at v-c buses
    V0[gbus[k]] = gen[on[k], VG] ./ abs.(V0[gbus[k]]).* V0[gbus[k]];
    ref0 = ref;                         ## save index and angle of
    Varef0 = bus[ref0, VA];             ##   original reference bus(es)
    ## build admittance matrices
    (Ybus, Yf, Yt) = makeYbus(baseMVA, bus, branch);
    ## function for computing V dependent complex bus power injections
        ## (generation - load)
    Sbus = Vm -> makeSbus(baseMVA, bus, gen, Vm);
    ## run the power flow
    V, success, iterations = newtonpf(Ybus, Sbus, V0, ref, pv, pq, 1.0000e-08, 10);
    its += iterations;
    ## update data matrices with solution
    bus, gen, branch = pfsoln(baseMVA, bus, gen, branch, Ybus, Yf, Yt, V, ref, pv, pq);
    if ref != ref0
    ## adjust voltage angles to make original ref bus correct
        bus[:, VA] = bus[:, VA] - bus[ref0, VA] + Varef0;
    end
    (bus, gen, branch) = ext2int(bus, gen, branch);
    ## Step 2.4: -----  output results  -----
    ## convert back to original bus numbering & print results
    mpc["bus"] = bus
    mpc["gen"] = gen
    mpc["branch"] = branch
    mpc["success"] = success
    mpc["iterations"] = its
    #comment: print the final voltage magnitudes
    # print(mpc["bus"],its);
    return mpc
end