%array size
clear all;
clc
close all;
N=10;
MAX_INT=1000;
lipid_length0=1.5;
J0_proxy=-0.1;
J0_dist=-0.1;
gamma=0.25;
kappa_m=10;
kappa_bar_m=-5;
kappa_t=0.25;
Rad=10;


%midplane boundary conditions
z0=0;
z1=4;
dydx0=tan(pi/4);
dydx1=0;
MP_an=[0,0,0,0,0,0,0,0,0,0];


%midplane first gauss shape
mid_plane_first_poly=poly_force_boundary(z0,z1,dydx0,dydx1,Rad,MP_an);



%proximal lipid director
Theta_p0=0;
Theta_p1=-pi/2;
dTheta_p0ds=0;
dTheta_p1ds=0;
Proxy_angle_an=[0,0,0,0,0,0,0,0,0,0];

proximal_angle_poly=poly_force_boundary(Theta_p0,Theta_p1,dTheta_p0ds,dTheta_p1ds,Rad,Proxy_angle_an);

%distal lipid director
Theta_d0=pi/2;
Theta_d1=pi/2;
dTheta_d0ds=0;
dTheta_d1ds=0;
distal_angle_an=[0,0,0,0,0,0,0,0,0,0];

Distal_angle_poly=poly_force_boundary(Theta_d0,Theta_d1,dTheta_d0ds,dTheta_d1ds,Rad,distal_angle_an);



Rad_vec=linspace(0,Rad,MAX_INT);
dr=Rad/MAX_INT;
% mid plane position and tangent angle
mid_plane_z=mypolyval(mid_plane_first_poly,Rad_vec);
tangent_angle=atan(mypolyder(mid_plane_first_poly,Rad_vec));


% proximal monolayer tilt and position
proxy_lipid_director_angle=mypolyval(proximal_angle_poly,Rad_vec);
proxy_tilt_angle=(pi/2+proxy_lipid_director_angle-tangent_angle);
proxy_tilt=tan(pi/2+proxy_lipid_director_angle-tangent_angle);
proxy_lipid_length=lipid_length0*(1+proxy_tilt.^2).^0.5;

proxy_plane_r=Rad_vec+cos(proxy_lipid_director_angle).*proxy_lipid_length;
proxy_plane_z=mid_plane_z+sin(proxy_lipid_director_angle).*proxy_lipid_length;

%proxy splays
proxy_para_splay=para_splay(proxy_plane_r,tangent_angle,proxy_tilt_angle,dr);
proxy_meridian_splay=Meridian_splay(dr,tangent_angle,proxy_tilt_angle,MAX_INT);

% distal monolayer tilt and position
dist_lipid_director_angle=mypolyval(Distal_angle_poly,Rad_vec);
dist_tilt_angle=pi/2-dist_lipid_director_angle+tangent_angle;
dist_tilt=tan(pi/2-dist_lipid_director_angle+tangent_angle);
dist_lipid_length=lipid_length0*(1+dist_tilt.^2).^0.5;

dist_plane_r=Rad_vec+cos(dist_lipid_director_angle).*dist_lipid_length;
dist_plane_z=mid_plane_z+sin(dist_lipid_director_angle).*dist_lipid_length;

%distal splays
dist_para_splay=para_splay(dist_plane_r,tangent_angle,-dist_tilt_angle,dr);
dist_meridian_splay=Meridian_splay(dr,tangent_angle,-dist_tilt_angle,MAX_INT);


%total and Gaussian and energy splay

proxy_total_splay=proxy_para_splay+proxy_meridian_splay;
proxy_Gauss_splay=proxy_para_splay.*proxy_meridian_splay;
proxy_area_vec=Area_of_revolution_vector(proxy_plane_r,proxy_plane_z,MAX_INT);
proxy_energy_density=(0.5*kappa_m*(proxy_total_splay-J0_proxy).^2+kappa_bar_m*proxy_Gauss_splay+0.5*kappa_t*proxy_tilt.^2+gamma);
proxy_energy=sum(proxy_energy_density.*proxy_area_vec);


dist_total_splay=dist_para_splay+dist_meridian_splay;
dist_Gauss_splay=dist_para_splay.*dist_meridian_splay;
dist_area_vec=Area_of_revolution_vector(dist_plane_r,dist_plane_z,MAX_INT);
dist_energy_density=(0.5*kappa_m*(dist_total_splay-J0_dist).^2+kappa_bar_m*dist_Gauss_splay+0.5*kappa_t*dist_tilt.^2+gamma);
dist_energy=sum(dist_energy_density.*dist_area_vec);

Energy_total=(proxy_energy+dist_energy-pi*Rad^2*(0.5*kappa_m*(J0_dist^2+J0_proxy^2)+gamma))*2

%plotting
%creat lipid director lines
X_lipids_distal(:,1)=Rad_vec';
X_lipids_distal(:,2)=dist_plane_r';
Y_lipids_distal(:,1)=mid_plane_z';
Y_lipids_distal(:,2)=dist_plane_z';

X_lipids_proxy(:,1)=Rad_vec';
X_lipids_proxy(:,2)=proxy_plane_r';
Y_lipids_proxy(:,1)=mid_plane_z';
Y_lipids_proxy(:,2)=proxy_plane_z';

figure(1)
hold on;
plot(Rad_vec,mid_plane_z,'b');
plot(dist_plane_r,dist_plane_z);
plot(proxy_plane_r,proxy_plane_z);
line(X_lipids_distal',Y_lipids_distal');
line(X_lipids_proxy',Y_lipids_proxy');
xlabel('radial length from stalk');
ylabel('hight');

figure(2)
hold on;
plot(Rad_vec,tangent_angle*180/pi);
plot(Rad_vec,proxy_lipid_director_angle*180/pi);
plot(Rad_vec,dist_lipid_director_angle*180/pi);

xlabel('radial length from stalk');
ylabel('Angle');
legend('tangent angle','proximal lipid angle','distal lipid angle');

figure(3)
hold on
plot(Rad_vec,proxy_tilt);
plot(Rad_vec,dist_tilt);
xlabel('radial length from stalk');
ylabel('tilt');
legend('proxy','distal');

figure(4)
hold on
plot(Rad_vec,(tangent_angle+proxy_tilt_angle)*180/pi);
plot(Rad_vec,(tangent_angle-dist_tilt_angle)*180/pi);
xlabel('radial length from stalk');
ylabel('tangent+tilt angle');
legend('proxy','distal');


figure(5)
hold on
plot(Rad_vec,proxy_meridian_splay);
plot(Rad_vec,proxy_para_splay);
plot(Rad_vec,proxy_total_splay);
xlabel('radial length from stalk');
ylabel('splay Proxy');
legend('meridian','parallel','total');

figure(6)
hold on
plot(Rad_vec,dist_meridian_splay);
plot(Rad_vec,dist_para_splay);
plot(Rad_vec,dist_total_splay);
xlabel('radial length from stalk');
ylabel('splay Distal');
legend('meridian','parallel','total');

figure(7)
hold on
plot(Rad_vec,0.5*kappa_m*(proxy_total_splay-J0_proxy).^2);
plot(Rad_vec,kappa_bar_m*proxy_Gauss_splay);
plot(Rad_vec,0.5*kappa_t*proxy_tilt.^2);
plot(Rad_vec,gamma);
xlabel('radial length from stalk');
ylabel('Energy density proxy');
legend('(delta Splay)^2','saddle splay','tilt','tension');

figure(8)
hold on
plot(Rad_vec,0.5*kappa_m*(dist_total_splay-J0_proxy).^2);
plot(Rad_vec,kappa_bar_m*dist_Gauss_splay);
plot(Rad_vec,0.5*kappa_t*dist_tilt.^2);
plot(Rad_vec,gamma);
xlabel('radial length from stalk');
ylabel('Energy density distal');
legend('(delta Splay)^2','saddle splay','tilt','tension');
