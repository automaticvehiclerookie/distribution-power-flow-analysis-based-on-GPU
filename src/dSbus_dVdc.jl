function dSbus_dVdc(Gbus_dc, V_dc, vcart=0)
    n = length(V_dc)
    Ibus = Gbus_dc * V_dc+Diagonal(Gbus_dc)*V_dc
    diagIbus = Diagonal(Ibus)
    dSbus_dV_dc=Gbus_dc.*V_dc
    dSbus_dV_dc = dSbus_dV_dc - Diagonal(diag(dSbus_dV_dc))
    dSbus_dV_dc = dSbus_dV_dc+diagIbus
    return dSbus_dV_dc
end