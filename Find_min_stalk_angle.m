function [min_Shape_polynum,Min_energy,min_angle] = Find_min_stalk_angle(Shape_polynum,Boundary_Conditions,Physical_Conditions,minimaztion_res,Stalk_angle_struct)




number_of_loops=round((Stalk_angle_struct.Max-Stalk_angle_struct.start_angle)/Stalk_angle_struct.res);
Total_energy=zeros(1,number_of_loops);
%Shape_polynum_vec=zeros(1,number_of_loops);
Angle_vec=zeros(1,number_of_loops);
First_up_int=round((Stalk_angle_struct.Max-Stalk_angle_struct.start_angle)/Stalk_angle_struct.res);
First_down_int=First_up_int-1;


int=First_up_int;

% scan up
Stalk_angle=Stalk_angle_struct.start_angle;
while Stalk_angle<=Stalk_angle_struct.Max
    Boundary_Conditions.dydx0=tan(Stalk_angle);
    [Shape_polynum,Total_energy(int)] = min_Multi_Coeff_loop(Shape_polynum,Boundary_Conditions,Physical_Conditions,minimaztion_res);
    Shape_polynum_vec(int)=Shape_polynum;
    Angle_vec(int)=Stalk_angle;
    int=int+1;
    Stalk_angle=Stalk_angle+Stalk_angle_struct.res;
end
    

% scan down
Stalk_angle=Stalk_angle_struct.start_angle-Stalk_angle_struct.res;
int=First_down_int;
Shape_polynum=Shape_polynum_vec(First_up_int);
while Stalk_angle>=Stalk_angle_struct.Min-0.00001
    Boundary_Conditions.dydx0=tan(Stalk_angle);
    [Shape_polynum,Total_energy(int)] = min_Multi_Coeff_loop(Shape_polynum,Boundary_Conditions,Physical_Conditions,minimaztion_res);
    Shape_polynum_vec(int)=Shape_polynum;
    Angle_vec(int)=Stalk_angle;
    int=int-1;
    Stalk_angle=Stalk_angle-Stalk_angle_struct.res;
end

    
    [Min_energy,Index]=min(Total_energy);
    min_Shape_polynum=Shape_polynum_vec(Index);
    min_angle=Angle_vec(Index);
    Min_energy;
    plot(Angle_vec*180/pi,Total_energy)
    
end

