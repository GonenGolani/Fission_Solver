function [out_shape_polynum,out_Boundary_Conditions,Energy_out] = Scan_parameters(Bud_Radius,theta_struct,edge_angle_stuct,Shape_polynum,Physical_Conditions,Boundary_Conditions,minimaztion_res,exclude_points,location)
    
    % Theta is the angel between the z axis and vector (between (0,0) to (Z,Rad))
    % Energy_out = out vector (theta, edge angle, energy)
    % exclude_points - matrix keeping data points that are not calcualted
    
    strat_minimaztion_res=minimaztion_res;
    Edge_angle=edge_angle_stuct.max;
    int_tot=1;
    while Edge_angle>edge_angle_stuct.min
        theta=theta_struct.min;
        minimaztion_res=strat_minimaztion_res;
        last_shape_polynum=Shape_polynum;
        while theta<theta_struct.max
            
            %do not calc at excluded points
            if(or(min(abs(exclude_points(:,:)-Edge_angle))>0.001,min(abs(exclude_points(:,:)-theta))>0.001))
                % reconfigure boundary conditions
                Boundary_Conditions.Rad=sin(Edge_angle*pi/180)*Bud_Radius;
                Boundary_Conditions.dydx1=tan(Edge_angle*pi/180);
                Boundary_Conditions.Theta_p1=-pi/2+Edge_angle*pi/180;
                Boundary_Conditions.Theta_d1=pi/2+Edge_angle*pi/180;
                Boundary_Conditions.Theta_d0=90*pi/180;
                Boundary_Conditions.dydx0=tan(pi/4);

                %scan angle up
                [New_Shape_polynum,new_Total_energy,minimaztion_res_opt] = min_Multi_Coeff_loop(last_shape_polynum,Boundary_Conditions,Physical_Conditions,minimaztion_res);
                minimaztion_res=minimaztion_res_opt;
                if and(new_Total_energy<100000, new_Total_energy>-100)
                    last_shape_polynum=New_Shape_polynum;
                    Energy_out(int_tot,1)=theta;
                    Energy_out(int_tot,2)=Edge_angle;
                    Energy_out(int_tot,3)=new_Total_energy;
                    out_shape_polynum(int_tot)=New_Shape_polynum;
                    out_Boundary_Conditions(int_tot)=Boundary_Conditions;
                    fprintf('Theta:%.1f , Edge Angle=%f , Energy %f\n',theta,Edge_angle,new_Total_energy);
                    save_data(Energy_out,out_shape_polynum,out_Boundary_Conditions,location);

                    int_tot=int_tot+1;
                end
            end
            theta=theta+theta_struct.res;


        end
        Edge_angle=Edge_angle-edge_angle_stuct.res;
    end
    

end

