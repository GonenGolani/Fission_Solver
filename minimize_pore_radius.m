function [New_Boundary_Conditions,min_energy,min_pore_radius] = minimize_pore_radius(Shape_polynum,Old_Boundary_Conditions,Physical_Conditions,minimaztion_res)

pore_radius_vector=linspace(Old_Boundary_Conditions.Rad0-minimaztion_res.length_res*10,Old_Boundary_Conditions.Rad0+minimaztion_res.length_res*10,20);

if Old_Boundary_Conditions.Rad0-minimaztion_res.length_res*10<Physical_Conditions.lipid_length0
    pore_radius_vector=linspace(Physical_Conditions.lipid_length0+0.01,Old_Boundary_Conditions.Rad0+minimaztion_res.length_res*10,20);
end


int=1;
vector_length=length(pore_radius_vector);
New_Boundary_Conditions=Old_Boundary_Conditions;

while int<=vector_length
    New_Boundary_Conditions.Rad0=pore_radius_vector(int);
    Energy_vector(int)= Pore_energy(Shape_polynum,New_Boundary_Conditions,Physical_Conditions);
    int=int+1;    
end
    Energy_vector(Energy_vector < -100) = NaN;
    [min_energy,min_index]=min(Energy_vector);
    min_pore_radius=pore_radius_vector(min_index);
    New_Boundary_Conditions.Rad0=min_pore_radius;

end