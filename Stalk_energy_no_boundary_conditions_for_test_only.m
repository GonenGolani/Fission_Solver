function [Total_energy] = Stalk_energy_no_boundary_conditions_for_test_only(Shape_polynum,Boundary_Conditions,Physical_Conditions,MAX_INT)
    % finds energy of stalk given initla lipids diredctor and mid-plane polynom
    % first 4 cells in polynom are zero, but they do not matter
    % Shape_polynum:
    %   MP - mid-plane polynom
    %   proxy_angle - proximal monolayer lipid director polynom
    %   distal_angle - distal monolayer lipid director polynom
    %
    % Boundary_Conditions: structure keeping all 12 boundary conditions

    

    %physical conditions:
    kappa_m=Physical_Conditions.kappa_m;
    kappa_bar_m=Physical_Conditions.kappa_bar_m;
    kappa_t=Physical_Conditions.kappa_t;
    J0_proxy=Physical_Conditions.J0_proxy;
    J0_dist=Physical_Conditions.J0_dist;
    gamma=Physical_Conditions.gamma;
    lipid_length0=Physical_Conditions.lipid_length0;
    

    %midplane boundary conditions
    z0=Boundary_Conditions.z0;
    z1=Boundary_Conditions.z1;
    dydx0=Boundary_Conditions.dydx0;
    dydx1=Boundary_Conditions.dydx1;
    Rad=Boundary_Conditions.Rad;
    
    %midplane first gauss shape
    mid_plane_first_poly=Shape_polynum.MP;


    %proximal lipid director
    Theta_p0=Boundary_Conditions.Theta_p0;
    Theta_p1=Boundary_Conditions.Theta_p1;
    dTheta_p0ds=Boundary_Conditions.dTheta_p0ds;
    dTheta_p1ds=Boundary_Conditions.dTheta_p1ds;

    proximal_angle_poly=Shape_polynum.proxy_angle;

    %distal lipid director
    Theta_d0=Boundary_Conditions.Theta_d0;
    Theta_d1=Boundary_Conditions.Theta_d1;
    dTheta_d0ds=Boundary_Conditions.dTheta_d0ds;
    dTheta_d1ds=Boundary_Conditions.dTheta_d1ds;

    distal_angle_poly=Shape_polynum.distal_angle;
    
    
    Rad_vec=linspace(0,Rad,MAX_INT);
    dr=Rad/MAX_INT;
    % mid plane position and tangent angle
    mid_plane_z=mypolyval(mid_plane_first_poly,Rad_vec);
    tangent_angle=atan(mypolyder(mid_plane_first_poly,Rad_vec));


    % proximal monolayer tilt and position
    % from here each degreen of freedom is in diffrent dimension
    % 1 - Radial direction of all terms
    % 2 - MD polynum coefficients
    % 3 - proxy polynum coefficients
    % 4 - distal polynum coefficients
    
    proxy_lipid_director_angle=mypolyval(proximal_angle_poly,Rad_vec);
    proxy_tilt_angle=(pi/2+proxy_lipid_director_angle-tangent_angle);
    proxy_tilt=tan(pi/2+proxy_lipid_director_angle-tangent_angle);
    proxy_lipid_length=lipid_length0*(1+proxy_tilt.^2).^0.5;

    proxy_plane_r=Rad_vec+cos(proxy_lipid_director_angle).*proxy_lipid_length;
    proxy_plane_z=mid_plane_z+sin(proxy_lipid_director_angle).*proxy_lipid_length;

    %proxy splays
    proxy_para_splay=para_splay(proxy_plane_r,tangent_angle,proxy_tilt_angle*0,dr);
    proxy_meridian_splay=Meridian_splay(dr,tangent_angle,proxy_tilt_angle*0,MAX_INT);

    % distal monolayer tilt and position
    dist_lipid_director_angle=mypolyval(distal_angle_poly,Rad_vec);
    dist_tilt_angle=pi/2-dist_lipid_director_angle+tangent_angle;
    dist_tilt=tan(pi/2-dist_lipid_director_angle+tangent_angle);
    dist_lipid_length=lipid_length0*(1+dist_tilt.^2).^0.5;

    dist_plane_r=Rad_vec+cos(dist_lipid_director_angle).*dist_lipid_length;
    dist_plane_z=mid_plane_z+sin(dist_lipid_director_angle).*dist_lipid_length;

    %distal splays
    dist_para_splay=para_splay(dist_plane_r,-tangent_angle,dist_tilt_angle,dr);
    dist_meridian_splay=Meridian_splay(dr,-tangent_angle,dist_tilt_angle,MAX_INT);


    %total and Gaussian and energy splay

    proxy_total_splay=proxy_para_splay+proxy_meridian_splay;
    proxy_Gauss_splay=proxy_para_splay.*proxy_meridian_splay;
    proxy_area_vec=Area_of_revolution_vector(proxy_plane_r,proxy_plane_z,MAX_INT);
    proxy_energy_density=(0.5*kappa_m*(proxy_total_splay-J0_proxy).^2+kappa_bar_m*proxy_Gauss_splay+0.5*kappa_t*proxy_tilt.^2+gamma);
    proxy_energy=sum(proxy_energy_density.*proxy_area_vec,2);


    dist_total_splay=dist_para_splay+dist_meridian_splay;
    dist_Gauss_splay=dist_para_splay.*dist_meridian_splay;
    dist_area_vec=Area_of_revolution_vector(dist_plane_r,dist_plane_z,MAX_INT);
    dist_energy_density=(0.5*kappa_m*(dist_total_splay-J0_dist).^2+kappa_bar_m*dist_Gauss_splay+0.5*kappa_t*dist_tilt.^2+gamma);
    dist_energy=sum(dist_energy_density.*dist_area_vec,2);

    Total_energy=(proxy_energy+dist_energy-pi*Rad^2*(0.5*kappa_m*(J0_dist^2+J0_proxy^2)+2*gamma));

end

