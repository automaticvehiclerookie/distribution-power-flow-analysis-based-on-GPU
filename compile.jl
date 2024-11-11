using PackageCompiler
using Pkg

# 激活项目环境并安装依赖项
Pkg.activate("C:/Users/13733/Desktop/DistributionPowerFlow")
Pkg.instantiate()

# 编译为可执行文件
create_app("C:/Users/13733/Desktop/DistributionPowerFlow", 
           "C:/Users/13733/Desktop/output"; 
           precompile_execution_file="C:/Users/13733/Desktop/DistributionPowerFlow/maintest.jl")