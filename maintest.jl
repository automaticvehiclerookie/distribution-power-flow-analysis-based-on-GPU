"""
    Main function for the AC power flow
    testfile:case7h.jl and testdata.jl
"""

# Detect the current working operating system
if Sys.iswindows()
    # Add the path to the data folder
    push!(LOAD_PATH, pwd()*"\\src\\")
    include(pwd()*"\\data\\testdata.jl")
else
    # Add the path to the data folder
    push!(LOAD_PATH, pwd()*"/src/")
    include(pwd()*"/data/case118.jl")
    using AppleAccelerate
end
# push!(LOAD_PATH, pwd()*"\\data\\");
using PowerFlow
# using MATLAB
# using Plots
# Start a MATLAB engine
# mat"addpath('C:\\Codes\\matpower')"
# mpc = case118();
#opt = PowerFlow.options() # The initial settings 
# opt["PF"]["NR_ALG"] = "bicgstab";
# opt["PF"]["ENFORCE_Q_LIMS"]=0
mpc=testdata()
#run power flow
 @time mpc = PowerFlow.runacdcpf(mpc)
 println(mpc["busAC"])
 println(mpc["busDC"])