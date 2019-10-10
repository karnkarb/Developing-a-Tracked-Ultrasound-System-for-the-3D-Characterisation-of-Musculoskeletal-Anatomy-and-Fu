function [FRE_mean,FRE_max,FRE_SD] = FRE_Feedback(Z_cross_wire_plot, XH_ac, tr_reg_probe_col, tr_reg_phantom_col)

for p = 3:3:size(Z_cross_wire_plot,1)
    
    q = p/3;
    Reg_XT(p-2,:) = tr_reg_probe_col{1, q}.b * Z_cross_wire_plot(p-2,:) * tr_reg_probe_col{1, q}.T + tr_reg_probe_col{1, q}.c(1,:);
    Reg_XT(p-1,:) = tr_reg_probe_col{1, q}.b * Z_cross_wire_plot(p-1,:) * tr_reg_probe_col{1, q}.T + tr_reg_probe_col{1, q}.c(1,:);
    Reg_XT(p,:) = tr_reg_probe_col{1, q}.b * Z_cross_wire_plot(p,:) * tr_reg_probe_col{1, q}.T + tr_reg_probe_col{1, q}.c(1,:);
    Reg_XH(p-2,:) = tr_reg_phantom_col{1, q}.b * Reg_XT(p-2,:) * tr_reg_phantom_col{1, q}.T + tr_reg_phantom_col{1, q}.c(1,:);
    Reg_XH(p-1,:) = tr_reg_phantom_col{1, q}.b * Reg_XT(p-1,:) * tr_reg_phantom_col{1, q}.T + tr_reg_phantom_col{1, q}.c(1,:);
    Reg_XH(p,:) = tr_reg_phantom_col{1, q}.b * Reg_XT(p,:) * tr_reg_phantom_col{1, q}.T + tr_reg_phantom_col{1, q}.c(1,:);

end

FRE_diff = Reg_XH - XH_ac;
FRE_norm = vecnorm(FRE_diff,2,2);
FRE_mean = mean(FRE_norm);
FRE_max = max(FRE_norm);
FRE_SD = std(FRE_norm);

end