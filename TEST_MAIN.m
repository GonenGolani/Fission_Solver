%master script - stalk energy minimaztion
clear all;
clc
close all;
do_plots=0;
location=pwd;
[Physical_Conditions,~,minimaztion_res,~,Configuration] = read_files(location);

load('Min_Shape_polynum.mat')
load('Min_Boundary_Conditions.mat');

Min_Boundary_Conditions.Rad=sin(Configuration.Edge_angle*pi/180)*Configuration.Bud_radius;



if strcmp(Configuration.Topology,'Stalk') 
    new_Total_energy = Stalk_energy(Min_Shape_polynum,Min_Boundary_Conditions,Physical_Conditions);
    plot_stalk(Min_Shape_polynum,Min_Boundary_Conditions,Physical_Conditions.lipid_length0,Configuration,do_plots);

elseif strcmp(Configuration.Topology,'Pore') 
    new_Total_energy= Pore_energy(Min_Shape_polynum,Min_Boundary_Conditions,Physical_Conditions);
    plot_Pore(Min_Shape_polynum,Min_Boundary_Conditions,Physical_Conditions.lipid_length0,Configuration,do_plots);
end
    %plot min configuartion

disp(new_Total_energy); 
