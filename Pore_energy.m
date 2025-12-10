function [Total_energy] = Pore_energy(Shape_polynum,Boundary_Conditions,Physical_Conditions)
    
    % finds energy of Pore given initla lipids diredctor and mid-plane polynom
    % the ploynum is inverse than stalk x(y)=sum(an*y*(n-1))
    % first 4 cells in polynom are zero, but they do not matter
    % Shape_polynum:
    %   MP - mid-plane polynom
    %   proxy_angle - proximal monolayer lipid director polynom
    %   distal_angle - distal monolayer lipid director polynom
    %
    % Boundary_Conditions: structure keeping all 12 boundary conditions
    MAX_INT=Boundary_Conditions.MAX_INT;

    %physical conditions:
    kappa_m=Physical_Conditions.kappa_m;
    kappa_bar_m=Physical_Conditions.kappa_bar_m;
    kappa_t=Physical_Conditions.kappa_t;
    J0_proxy=Physical_Conditions.J0_proxy;
    J0_dist=Physical_Conditions.J0_dist;
    gamma=Physical_Conditions.gamma;
    lipid_length0=Physical_Conditions.lipid_length0;
    

    z0=Boundary_Conditions.z0;
    z1=Boundary_Conditions.z1;
    dydx0=Boundary_Conditions.dydx0;
    dydx1=Boundary_Conditions.dydx1;
    Rad1=Boundary_Conditions.Rad;
    Rad0=Boundary_Conditions.Rad0;
    %midplane first gauss shape
    mid_plane_first_poly=poly_force_boundary(Rad0,Rad1,dydx0,dydx1,z1,z0,Shape_polynum.MP);



    mid_plane_z=linspace(z0,z1,MAX_INT);


    % mid plane position and tangent angle
    Rad_vec=mypolyval(mid_plane_first_poly,mid_plane_z);
    dr=diff(Rad_vec,1,2); dr(:,MAX_INT)=dr(:,MAX_INT-1);
    tangent_angle=pi/2-atan2(mypolyder(mid_plane_first_poly,mid_plane_z),1);


    
    % proximal monolayer tilt and position
    proxy_lipid_director_angle=-pi/2+tangent_angle;
    proxy_tilt_angle=(pi/2+proxy_lipid_director_angle-tangent_angle);
    proxy_tilt=tan(pi/2+proxy_lipid_director_angle-tangent_angle);
    proxy_lipid_length=lipid_length0*(1+proxy_tilt.^2).^0.5;

    proxy_plane_r=Rad_vec+cos(proxy_lipid_director_angle).*proxy_lipid_length;
    proxy_plane_z=mid_plane_z+sin(proxy_lipid_director_angle).*proxy_lipid_length;

    %proxy splays
    proxy_para_splay=para_splay(abs(proxy_plane_r),tangent_angle,proxy_tilt_angle,dr);
    proxy_meridian_splay=Meridian_splay(dr,tangent_angle,proxy_tilt_angle,MAX_INT);

    % distal monolayer tilt and position
    dist_lipid_director_angle=pi/2+tangent_angle;
    dist_tilt_angle=pi/2-dist_lipid_director_angle+tangent_angle;
    dist_tilt=tan(pi/2-dist_lipid_director_angle+tangent_angle);
    dist_lipid_length=lipid_length0*(1+dist_tilt.^2).^0.5;

    dist_plane_r=Rad_vec+cos(dist_lipid_director_angle).*dist_lipid_length;
    dist_plane_z=mid_plane_z+sin(dist_lipid_director_angle).*dist_lipid_length;

    %distal splays
    dist_para_splay=para_splay(abs(dist_plane_r),-tangent_angle,dist_tilt_angle,dr);
    dist_meridian_splay=Meridian_splay(dr,-tangent_angle,dist_tilt_angle,MAX_INT);


    %total and Gaussian and energy splay

    proxy_total_splay=proxy_para_splay+proxy_meridian_splay;
    proxy_Gauss_splay=proxy_para_splay.*proxy_meridian_splay;
    proxy_area_vec=Area_of_revolution_vector(abs(proxy_plane_r),proxy_plane_z,MAX_INT);


    dist_total_splay=dist_para_splay+dist_meridian_splay;
    dist_Gauss_splay=dist_para_splay.*dist_meridian_splay;
    dist_area_vec=Area_of_revolution_vector(abs(dist_plane_r),dist_plane_z,MAX_INT);
    

    
    proxy_energy_density=0.5*kappa_m*(proxy_total_splay-J0_proxy).^2+kappa_bar_m*proxy_Gauss_splay+0.5*kappa_t*proxy_tilt.^2+gamma-0.5*kappa_m*J0_proxy^2;
    proxy_energy=sum(proxy_energy_density.*proxy_area_vec,2);
    
    dist_energy_density=0.5*kappa_m*(dist_total_splay-J0_dist).^2+kappa_bar_m*dist_Gauss_splay+0.5*kappa_t*dist_tilt.^2+gamma-0.5*kappa_m*J0_dist^2;
    dist_energy=sum(dist_energy_density.*dist_area_vec,2);
    

    %add panaelty for any negative value
    negative_energy_penalty = (any(dist_plane_r < 0,2)+ any(dist_plane_z < 0,2)+any(proxy_plane_r < 0,2)+ any(proxy_plane_z < 0,2))*10^6;
 
    Total_energy=proxy_energy+dist_energy+negative_energy_penalty;

end

