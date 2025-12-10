%master script - stalk energy minimaztion
clear all;
clc
close all;
do_plots=0;
location=pwd;
[Physical_Conditions,Boundary_Conditions,minimaztion_res,Shape_polynum,Configuration] = read_files(location);

Boundary_Conditions.Rad=sin(Configuration.Edge_angle*pi/180)*Configuration.Bud_radius;

%load boundary conditions to vectors:
if strcmp(Configuration.Topology,'Stalk') 
    Boundary_Conditions.Theta_d0=90*pi/180;
    Boundary_Conditions.dydx0=tan(pi/4);
    Boundary_Conditions.Rad0=0.01;
    Boundary_Conditions.z1=Boundary_Conditions.Rad/2;
elseif strcmp(Configuration.Topology,'Pore') 
    Configuration.Edge_angle=90-Configuration.Edge_angle;
    Boundary_Conditions.Theta_d0=180*pi/180;
    Boundary_Conditions.dydx0=0;
    Boundary_Conditions.Rad0=2;
    Boundary_Conditions.z1=Boundary_Conditions.Rad/2;
else
    fprintf('Error: no Topology loaded!\n')
end

Boundary_Conditions.dydx1=tan(Configuration.Edge_angle*pi/180);
Boundary_Conditions.Theta_p1=-pi/2+Configuration.Edge_angle*pi/180;
Boundary_Conditions.Theta_d1=pi/2+Configuration.Edge_angle*pi/180;

[Min_Shape_polynum,Min_Boundary_Conditions,new_Total_energy] = ...
    min_Multi_Coeff_loop(Shape_polynum,Boundary_Conditions,Physical_Conditions,minimaztion_res,Configuration);

cd('Simulation Output')

if strcmp(Configuration.Topology,'Stalk') 
    new_Total_energy = Stalk_energy(Min_Shape_polynum,Min_Boundary_Conditions,Physical_Conditions);
    plot_stalk(Min_Shape_polynum,Min_Boundary_Conditions,Physical_Conditions.lipid_length0,Configuration,do_plots);

elseif strcmp(Configuration.Topology,'Pore') 
    new_Total_energy= Pore_energy(Min_Shape_polynum,Min_Boundary_Conditions,Physical_Conditions);
    plot_Pore(Min_Shape_polynum,Min_Boundary_Conditions,Physical_Conditions.lipid_length0,Configuration,do_plots);
    pore_water_gap=2*(Min_Boundary_Conditions.Rad0-Physical_Conditions.lipid_length0);
    fid=fopen('pore_water_gap_diameter.txt','wt'); fprintf (fid,'%f',pore_water_gap); fclose (fid);

end
    %plot min configuartion

fid=fopen('Final_energy.txt','wt'); fprintf (fid,'%f',new_Total_energy); fclose (fid);


save('Min_Shape_polynum.mat','Min_Shape_polynum');            
save('Min_Boundary_Conditions.mat','Min_Boundary_Conditions');
save('Minimized_configuration')

cd(location);

exit;
