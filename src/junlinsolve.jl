function junlinsolve(A, b)
    n = size(A,1)
        F = ilu(A, τ = 0.05)
        opM = LinearOperator(Float64, n, n, false, false, (y, v) -> forward_substitution!(y, F, v))
        opN = LinearOperator(Float64, n, n, false, false, (y, v) -> backward_substitution!(y, F, v))
        x , stats = bicgstab(A, b, history=false, M=opM, N=opN)
        return x
    
end