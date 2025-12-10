function [New_Shape_polynum,New_Boundary_Conditions,Energy_end] =...NEW_minimaztion_res
    min_Multi_Coeff_loop(Shape_polynum,Boundary_Conditions,Physical_Conditions,minimaztion_res,Configuration)

minimaztion_res5=minimaztion_res;
minimaztion_res6=minimaztion_res;
minimaztion_res7=minimaztion_res;
minimaztion_res8=minimaztion_res;
minimaztion_res9=minimaztion_res;
minimaztion_res10=minimaztion_res;

minimaztion_res5.MP=minimaztion_res.MP5;
minimaztion_res6.MP=minimaztion_res.MP6;
minimaztion_res7.MP=minimaztion_res.MP7;
minimaztion_res8.MP=minimaztion_res.MP8;
minimaztion_res9.MP=minimaztion_res.MP9;
minimaztion_res10.MP=minimaztion_res.MP10;

minimaztion_res5.proxy=minimaztion_res.proxy5;
minimaztion_res6.proxy=minimaztion_res.proxy6;
minimaztion_res7.proxy=minimaztion_res.proxy7;
minimaztion_res8.proxy=minimaztion_res.proxy8;
minimaztion_res9.proxy=minimaztion_res.proxy9;
minimaztion_res10.proxy=minimaztion_res.proxy10;

minimaztion_res5.distal=minimaztion_res.distal5;
minimaztion_res6.distal=minimaztion_res.distal6;
minimaztion_res7.distal=minimaztion_res.distal7;
minimaztion_res8.distal=minimaztion_res.distal8;
minimaztion_res9.distal=minimaztion_res.distal9;
minimaztion_res10.distal=minimaztion_res.distal10;



if strcmp(Configuration.Topology,'Stalk') 
    energy_strat = Stalk_energy(Shape_polynum,Boundary_Conditions,Physical_Conditions);
elseif strcmp(Configuration.Topology,'Pore') 
    energy_strat= Pore_energy(Shape_polynum,Boundary_Conditions,Physical_Conditions);
end
fprintf('Start Energy: %.2f\n',energy_strat);

Energy_end=0;
New_Shape_polynum=Shape_polynum;
New_Boundary_Conditions=Boundary_Conditions;
int=1;
while (energy_strat-Energy_end)>minimaztion_res.Enrgy_MIN_RES && int<=minimaztion_res.MAX_LOOPs
    if strcmp(Configuration.Topology,'Stalk') 
        energy_strat = Stalk_energy(New_Shape_polynum,New_Boundary_Conditions,Physical_Conditions);
    elseif strcmp(Configuration.Topology,'Pore') 
        energy_strat= Pore_energy(New_Shape_polynum,New_Boundary_Conditions,Physical_Conditions);
    end
    [New_Shape_polynum,New_Boundary_Conditions,New_Total_energy]...
        = min_Multi_loop(New_Shape_polynum,New_Boundary_Conditions,Physical_Conditions,minimaztion_res5,5,Configuration);

    [New_Shape_polynum,New_Boundary_Conditions,New_Total_energy]...
        = min_Multi_loop(New_Shape_polynum,New_Boundary_Conditions,Physical_Conditions,minimaztion_res6,6,Configuration);

    [New_Shape_polynum,New_Boundary_Conditions,New_Total_energy]...
        = min_Multi_loop(New_Shape_polynum,New_Boundary_Conditions,Physical_Conditions,minimaztion_res7,7,Configuration);

    [New_Shape_polynum,New_Boundary_Conditions,New_Total_energy]...
        = min_Multi_loop(New_Shape_polynum,New_Boundary_Conditions,Physical_Conditions,minimaztion_res8,8,Configuration);

    [New_Shape_polynum,New_Boundary_Conditions,New_Total_energy]...
        = min_Multi_loop(New_Shape_polynum,New_Boundary_Conditions,Physical_Conditions,minimaztion_res9,9,Configuration);

    [New_Shape_polynum,New_Boundary_Conditions,Energy_end]...
        = min_Multi_loop(New_Shape_polynum,New_Boundary_Conditions,Physical_Conditions,minimaztion_res10,10,Configuration);


    fprintf('Energy:%f , Delta energy=%f\n',Energy_end,energy_strat-Energy_end);
    

    int=int+1;
end



end

