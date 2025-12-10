function [New_Boundary_Conditions,min_energy,min_z1] = minimize_z(Shape_polynum,Old_Boundary_Conditions,Physical_Conditions,minimaztion_res,Configuration)

z_vector=linspace(Old_Boundary_Conditions.z1-minimaztion_res*100,Old_Boundary_Conditions.z1+minimaztion_res*100,100);
if Old_Boundary_Conditions.z1-minimaztion_res*100<1
    z_vector=linspace(1,Old_Boundary_Conditions.z1+minimaztion_res*100,100);
end

int=1;
vector_length=length(z_vector);
New_Boundary_Conditions=Old_Boundary_Conditions;
while int<=vector_length
    Old_Boundary_Conditions.z1=z_vector(int);
    if strcmp(Configuration.Topology,'Stalk') 
        Energy_vector(int) = Stalk_energy(Shape_polynum,Old_Boundary_Conditions,Physical_Conditions);
    elseif strcmp(Configuration.Topology,'Pore') 
        Energy_vector(int)= Pore_energy(Shape_polynum,Old_Boundary_Conditions,Physical_Conditions);
    end
    int=int+1;
end
    Energy_vector(Energy_vector < -100) = NaN;
    [min_energy,min_index]=min(Energy_vector);
    min_z1=z_vector(min_index);
    New_Boundary_Conditions.z1=min_z1;

end