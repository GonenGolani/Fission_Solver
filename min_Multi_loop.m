function [New_Shape_polynum,New_Boundary_Conditions,New_Total_energy] = min_Multi_loop(Shape_polynum,Boundary_Conditions,Physical_Conditions,minimaztion_res,min_coeff,Configuration)
    int=1;
    count=0;
    if strcmp(Configuration.Topology,'Stalk') 
        Old_energy = Stalk_energy(Shape_polynum,Boundary_Conditions,Physical_Conditions);
    elseif strcmp(Configuration.Topology,'Pore') 
        Old_energy= Pore_energy(Shape_polynum,Boundary_Conditions,Physical_Conditions);
    end
    start_energy=Old_energy;

    Old_Shape_polynum=Shape_polynum;
    Old_Boundary_Conditions=Boundary_Conditions;
    start_Shape_polynum=Shape_polynum;
    start_Boundary_Conditions=Boundary_Conditions;
    energy_change=100;
    while int<=minimaztion_res.MAX_LOOPs && abs(energy_change)>minimaztion_res.Enrgy_MIN_RES
       
        [New_Shape_polynum,New_Boundary_Conditions,New_Total_energy] = min_one_loop(Old_Shape_polynum,Old_Boundary_Conditions,Physical_Conditions,minimaztion_res,min_coeff,Configuration);
        energy_change=New_Total_energy-Old_energy;

        if New_Total_energy<Old_energy
            count=0;
            Old_Shape_polynum=New_Shape_polynum;
            Old_energy=New_Total_energy;
            Old_Boundary_Conditions=New_Boundary_Conditions;
        end
        if abs(energy_change)<minimaztion_res.Enrgy_MIN_RES || New_Total_energy>Old_energy
            if minimaztion_res.MP>minimaztion_res.min_res_limit
                minimaztion_res.MP=minimaztion_res.MP/2;
            end
            if minimaztion_res.proxy>minimaztion_res.min_res_limit
                minimaztion_res.proxy=minimaztion_res.proxy/2;
            end
            if minimaztion_res.distal>minimaztion_res.min_res_limit
                minimaztion_res.distal=minimaztion_res.distal/2;
            end
            if minimaztion_res.length_res>minimaztion_res.min_res_limit
                minimaztion_res.length_res=minimaztion_res.length_res/2;
            end
            count=count+1;
        end
       

        if count==100
            break;
        end
        int=int+1;
    
    end

    if New_Total_energy>start_energy
        New_Shape_polynum=start_Shape_polynum;
        New_Boundary_Conditions=start_Boundary_Conditions;
        New_Total_energy=start_energy;
    end

end

   
