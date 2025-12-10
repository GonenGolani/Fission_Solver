function [Physical_Conditions,Boundary_Conditions,minimaztion_res,Shape_polynum,Configuration] = read_files(location)
% out =1 all OK
% out=0 not good
cd(location); 
cd('physical conditions');
Physical_Conditions.lipid_length0=importdata('lipid_length.txt');
Physical_Conditions.J0_proxy=importdata('J0_proxy.txt');
Physical_Conditions.J0_dist=importdata('J0_dist.txt');
Physical_Conditions.gamma=importdata('gamma.txt');
Physical_Conditions.kappa_m=importdata('kappa_m.txt');
Physical_Conditions.kappa_bar_m=importdata('kappa_bar_m.txt');
Physical_Conditions.kappa_t=importdata('kappa_t.txt');


cd(location); cd('Boundary conditions');
Boundary_Conditions.MAX_INT=importdata('MAX_INT.txt');

cd(location); cd('Boundary conditions'); cd('MP') 
Boundary_Conditions.z0=importdata('z0.txt');
Boundary_Conditions.z1=importdata('z1.txt');
Boundary_Conditions.dydx0=importdata('dydx0.txt');
Boundary_Conditions.dydx1=importdata('dydx1.txt');
Boundary_Conditions.Rad=importdata('Rad.txt');
Boundary_Conditions.Rad0=importdata('Rad0.txt');


cd(location); cd('Boundary conditions'); cd('Proximal');
Boundary_Conditions.Theta_p0=importdata('Theta_p0.txt');
Boundary_Conditions.Theta_p1=importdata('Theta_p1.txt');
Boundary_Conditions.dTheta_p0ds=importdata('dTheta_p0ds.txt');
Boundary_Conditions.dTheta_p1ds=importdata('dTheta_p1ds.txt');

cd(location); cd('Boundary conditions'); cd('Distal');
Boundary_Conditions.Theta_d0=importdata('Theta_d0.txt');
Boundary_Conditions.Theta_d1=importdata('Theta_d1.txt');
Boundary_Conditions.dTheta_d0ds=importdata('dTheta_d0ds.txt');
Boundary_Conditions.dTheta_d1ds=importdata('dTheta_d1ds.txt');

cd(location); cd('Minimaztion'); cd('General');
minimaztion_res.MAX_LOOPs=importdata('MAX_LOOPs.txt');
minimaztion_res.Enrgy_MIN_RES=importdata('Enrgy_MIN_RES.txt');
minimaztion_res.MP_steps=importdata('MP_steps.txt');
minimaztion_res.Proxy_steps=importdata('Proxy_steps.txt');
minimaztion_res.Distal_steps=importdata('Distal_steps.txt');
minimaztion_res.min_res_limit=importdata('min_res_limit.txt');
minimaztion_res.Energy_diverge_up=importdata('Energy_diverge_up.txt');
minimaztion_res.Energy_diverge_down=importdata('Energy_diverge_down.txt');
minimaztion_res.length_res=importdata('length_res.txt');


cd(location); cd('Minimaztion'); cd('Coeff 5');
minimaztion_res.distal5=importdata('distal5.txt');
minimaztion_res.MP5=importdata('MP5.txt');
minimaztion_res.proxy5=importdata('proxy5.txt');

cd(location); cd('Minimaztion'); cd('Coeff 6');
minimaztion_res.distal6=importdata('distal6.txt');
minimaztion_res.MP6=importdata('MP6.txt');
minimaztion_res.proxy6=importdata('proxy6.txt');

cd(location); cd('Minimaztion'); cd('Coeff 7');
minimaztion_res.distal7=importdata('distal7.txt');
minimaztion_res.MP7=importdata('MP7.txt');
minimaztion_res.proxy7=importdata('proxy7.txt');

cd(location); cd('Minimaztion'); cd('Coeff 8');
minimaztion_res.distal8=importdata('distal8.txt');
minimaztion_res.MP8=importdata('MP8.txt');
minimaztion_res.proxy8=importdata('proxy8.txt');

cd(location); cd('Minimaztion'); cd('Coeff 9');
minimaztion_res.distal9=importdata('distal9.txt');
minimaztion_res.MP9=importdata('MP9.txt');
minimaztion_res.proxy9=importdata('proxy9.txt');

cd(location); cd('Minimaztion'); cd('Coeff 10');
minimaztion_res.distal10=importdata('distal10.txt');
minimaztion_res.MP10=importdata('MP10.txt');
minimaztion_res.proxy10=importdata('proxy10.txt');

cd(location); cd('computed Parameters'); 
Configuration.Bud_radius=importdata('Bud_radius.txt');
Configuration.Edge_angle=importdata('Edge_angle.txt');
Configuration.Topology=importdata('Topology.txt');

if exist('initial_ploy.mat')==2
    Shape_polynum=importdata('initial_ploy.mat');
end
if exist('initial_ploy.mat')==0
    initial_poly=[0,0,0,0,0,0,0,0,0,0,0];
    Shape_polynum.MP=initial_poly;
    Shape_polynum.proxy_angle=initial_poly;
    Shape_polynum.distal_angle=initial_poly;
end



cd(location);
end

