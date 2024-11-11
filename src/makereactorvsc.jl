function makereactorvsc(rtf,xtf,rc,xc)
    ztf=rtf.+im*xtf
    ytf = inv.(ztf)
    Gtf=real(ytf)
    Btf=imag(ytf)
    zc=rc+im*xc
    yc=inv.(zc)
    Gc=real(yc)
    Bc=imag(yc)
    return Gtf,Btf,Gc,Bc
end