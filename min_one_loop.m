function [New_Shape_polynum,New_Boundary_Conditions,New_Total_energy] = min_one_loop(Shape_polynum,Boundary_Conditions,Physical_Conditions,minimaztion_res,Coeff_n,Configuration)

New_Shape_polynum=Shape_polynum;
New_Boundary_Conditions=Boundary_Conditions;

[New_Boundary_Conditions,~,~] = minimize_z(New_Shape_polynum,New_Boundary_Conditions,Physical_Conditions,minimaztion_res.length_res,Configuration);

if strcmp(Configuration.Topology,'Pore')
    [New_Boundary_Conditions,~,~] = minimize_pore_radius(New_Shape_polynum,New_Boundary_Conditions,Physical_Conditions,minimaztion_res);
    [New_Shape_polynum,New_Total_energy]=Minimize_MP_polynom_PORE(New_Shape_polynum,New_Boundary_Conditions,Physical_Conditions,minimaztion_res.MP,minimaztion_res.MP_steps,Coeff_n);
end


if strcmp(Configuration.Topology,'Stalk')
    [New_Shape_polynum,~] = Minimize_N_order_one_iter("MP",New_Shape_polynum,New_Boundary_Conditions,Physical_Conditions,minimaztion_res.MP,minimaztion_res.MP_steps,Coeff_n);
    [New_Shape_polynum,~] = Minimize_N_order_one_iter("proxy",New_Shape_polynum,New_Boundary_Conditions,Physical_Conditions,minimaztion_res.proxy,minimaztion_res.Proxy_steps,Coeff_n);
    [New_Shape_polynum,New_Total_energy] = Minimize_N_order_one_iter("distal",New_Shape_polynum,New_Boundary_Conditions,Physical_Conditions,minimaztion_res.distal,minimaztion_res.Distal_steps,Coeff_n);

end



end

