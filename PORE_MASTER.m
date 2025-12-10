%master script - stalk energy minimaztion
clear all;
clc
close all;

location=pwd;
[Physical_Conditions,Boundary_Conditions,minimaztion_res,Stalk_nucleus_Rad,theta_struct,edge_angle_stuct,Rad0_stuct,exclude_points,Shape_polynum] = read_files(location);

%load boundary conditions to vectors:

Boundary_Conditions.Theta_d0=180*pi/180;
Boundary_Conditions.dydx0=0;
[out_shape_polynum,out_Boundary_Conditions,Energy_out,out_Energy_Rad0] = Scan_parameters(Stalk_nucleus_Rad,theta_struct,edge_angle_stuct,Rad0_stuct,Shape_polynum,Physical_Conditions,Boundary_Conditions,minimaztion_res,exclude_points,location);
%save_data(Energy_out,out_shape_polynum,out_Boundary_Conditions,location)
%min energy configuration
[min_energy,index]=min(Energy_out(:,3));


plot_Pore(out_shape_polynum(index),out_Boundary_Conditions(index),Physical_Conditions.lipid_length0);

h =  findobj('type','figure'); NUMBER_OF_FIGURES = length(h);

figure(NUMBER_OF_FIGURES+1)
scatter3(Energy_out(:,1),Energy_out(:,2),Energy_out(:,3));
xlabel('Theta');
ylabel('Phi');
zlabel('Energy');

plot_Energy_Rad0(out_Energy_Rad0)



prompt = 'End Of file';
input(prompt)
exit;
