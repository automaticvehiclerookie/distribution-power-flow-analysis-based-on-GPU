"""
    Define the power flow module with different functions
"""

module PowerFlow
    using Printf
    using SparseArrays
    using LinearAlgebra
    using PrettyTables
    using AMD
    using SuiteSparse
    using IterativeSolvers
    using IncompleteLU
    using KrylovKit
    using Krylov
    using LinearOperators
    using CUDA.CUSPARSE
    using CUDA
    using AlgebraicMultigrid
    using JuMP
    using Ipopt
    using NonlinearSolve
    using InvertedIndices
    # using KrylovPreconditioners
    # using different packages based on the operating system
    if Sys.iswindows()
        using CUDA
    # else
    #     using Metal
    end
    include("idx.jl")
    include("bustypes.jl")
    include("ext2int.jl")
    include("makeYbus.jl")
    include("newtonpf.jl")
    include("makeSbus.jl")
    include("makeSdaczip.jl")
    include("makeSddczip.jl")
    include("mplinsolve.jl")
    include("total_load.jl")
    include("pfsoln.jl")
    include("dSbus_dV.jl")
    include("runpf.jl")
    include("settings.jl")
    include("rundcpf.jl")
    include("makeBdc.jl")
    include("dcpf.jl")
    include("makereactorvsc.jl")
    include("hybrid_bustypes.jl")
    include("makehYbus.jl")
    include("runmodel.jl")
    include("hybridnewtonpf.jl")
    include("makeGbus.jl")
    include("makeacSbus.jl")
    include("makedcSbus.jl")
    include("dSbus_dVac.jl")
    include("dSbus_dVdc.jl")
    include("acdcnewtonpf.jl")
    include("junlinsolve.jl")
    include("updateconvertermatrix.jl")
    include("makevsc_acjacob.jl")
    include("makevsc_dcjacob.jl")
    include("makevsc_vscjacobi.jl")
    include("pfsolnh.jl")
    include("runacdcpf.jl")
    # include("data_structure.jl")
    # export idx_bus, idx_brch, idx_gen, bustypes, makeYbus, newtonpf, makeSbus, makeSdzip, mplinsolve, total_load, pfsoln, dSbus_dV, MPC
end

export PowerFlow