function plot_Pore(Shape_polynum,Boundary_Conditions,lipid_length0,Configuration,do_plots)
    % finds energy of stalk given initla lipids diredctor and mid-plane polynom
    % first 4 cells in polynom are zero, but they do not matter
    % Shape_polynum:
    %   MP - mid-plane polynom
    %   proxy_angle - proximal monolayer lipid director polynom
    %   distal_angle - distal monolayer lipid director polynom
    %
    % Boundary_Conditions: structure keeping all 12 boundary conditions
    MAX_INT=Boundary_Conditions.MAX_INT;
    
    close all;

    %midplane boundary conditions
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
    dr=diff(Rad_vec); dr(MAX_INT)=dr(MAX_INT-1);
    tangent_angle=pi/2-atan2(mypolyder(mid_plane_first_poly,mid_plane_z),1);


    
    % proximal monolayer tilt and position
    proxy_lipid_director_angle=-pi/2+tangent_angle;
    proxy_tilt_angle=(pi/2+proxy_lipid_director_angle-tangent_angle);
    proxy_tilt=tan(pi/2+proxy_lipid_director_angle-tangent_angle);
    proxy_lipid_length=lipid_length0*(1+proxy_tilt.^2).^0.5;

    proxy_plane_r=Rad_vec+cos(proxy_lipid_director_angle).*proxy_lipid_length;
    proxy_plane_z=mid_plane_z+sin(proxy_lipid_director_angle).*proxy_lipid_length;

    %proxy splays
    proxy_para_splay=para_splay(proxy_plane_r,tangent_angle,proxy_tilt_angle,dr);
    proxy_meridian_splay=Meridian_splay(dr,tangent_angle,proxy_tilt_angle,MAX_INT);

    % distal monolayer tilt and position
    dist_lipid_director_angle=pi/2+tangent_angle;
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


    dist_total_splay=dist_para_splay+dist_meridian_splay;
    dist_Gauss_splay=dist_para_splay.*dist_meridian_splay;
    dist_area_vec=Area_of_revolution_vector(dist_plane_r,dist_plane_z,MAX_INT); 



    %creat lipid director lines
    n=10;
  
    X_lipids_distal(:,1)=Rad_vec([1:n:MAX_INT MAX_INT])';
    X_lipids_distal(:,2)=dist_plane_r([1:n:MAX_INT MAX_INT])';
    Y_lipids_distal(:,1)=mid_plane_z([1:n:MAX_INT MAX_INT])';
    Y_lipids_distal(:,2)=dist_plane_z([1:n:MAX_INT MAX_INT])';

    X_lipids_proxy(:,1)=Rad_vec([1:n:MAX_INT MAX_INT])';
    X_lipids_proxy(:,2)=proxy_plane_r([1:n:MAX_INT MAX_INT])';
    Y_lipids_proxy(:,1)=mid_plane_z([1:n:MAX_INT MAX_INT])';
    Y_lipids_proxy(:,2)=proxy_plane_z([1:n:MAX_INT MAX_INT])';

   %phi vector
   edge_angle=atan(dydx1^-1);
   phi_vector(1,:)=linspace(Rad1,Rad1+lipid_length0*cos(edge_angle));
   phi_vector(2,:)= z1+(phi_vector(1,:)-Rad1)*dydx1^-1;
   
   % bud vector plot 
   angle_theta=linspace(-90+edge_angle/pi*180,90,200);
   arc_theta_vector(1,:)=Configuration.Bud_radius*cos(angle_theta*pi/180);
   arc_theta_vector(2,:)=Configuration.Bud_radius*sin(angle_theta*pi/180)+z1+Configuration.Bud_radius*cos(edge_angle);
   

     
    figure(1)
    hold on;
    plot(Rad_vec,mid_plane_z,'b');
    plot(dist_plane_r,dist_plane_z);
    plot(proxy_plane_r,proxy_plane_z);
    
    %phi and theta angle plot
    p_phi1=plot(phi_vector(1,:),phi_vector(2,:));
    p_phi2=plot(linspace(Rad1,Rad1+lipid_length0),z1*ones(1,100));
    p_theta2=plot(arc_theta_vector(1,:),arc_theta_vector(2,:));
    p_phi1.LineStyle='--';   p_phi1.Color='k';     p_phi2.LineStyle='--';   p_phi2.Color='k'; p_theta2.LineStyle='--'; p_theta2.Color='k'; 
    
    
    %phi_txt = ['\phi=' num2str(edge_angle*180/pi) '^o'];
    %t1=text(Rad1+lipid_length0/4,z1+lipid_length0/4,phi_txt);
    %t2=text(Rad1/2,(Rad1^2+z1^2-Rad1/2.^2).^0.5-lipid_length0/4,theta_txt); 
    
    line(X_lipids_distal',Y_lipids_distal','Color','magenta');
    line(X_lipids_proxy',Y_lipids_proxy','Color','red');
    xL=xlabel('radial length from stalk [nm]');
    yL=ylabel('z [nm]');
    xlim([0 (Rad1^2+z1^2)^0.5+1.2])
    ylim([0 (Rad1^2+z1^2)^0.5+1.2])
    set(gca,'FontSize',18);
    set(gca,'FontWeight','bold');
    xL.FontSize=18;
    xL.FontWeight='bold';
    yL.FontSize=18;
    yL.FontWeight='bold';
    set(gcf,'color','white')
    saveas(gcf,'Pore shape','fig');
    %saveas(gcf,'Pore shape','emf');
    %plotting
    
 if do_plots   
    figure(2)
    hold on;
    plot(Rad_vec,tangent_angle*180/pi);
    plot(Rad_vec,proxy_lipid_director_angle*180/pi);
    plot(Rad_vec,dist_lipid_director_angle*180/pi);

    xL=xlabel('radial length from stalk [nm]');
    yL=ylabel('Angle');
    legend('tangent angle','proximal lipid angle','distal lipid angle');
    set(gca,'FontSize',18);
    set(gca,'FontWeight','bold');
    xL.FontSize=18;
    xL.FontWeight='bold';
    yL.FontSize=18;
    yL.FontWeight='bold';


    figure(4)
    hold on
    plot(Rad_vec,(tangent_angle+proxy_tilt_angle)*180/pi);
    plot(Rad_vec,(tangent_angle-dist_tilt_angle)*180/pi);
    plot(Rad_vec,(tangent_angle)*180/pi);
    xL=xlabel('radial length from stalk');
    yL=ylabel('tangent+tilt angle');
    legend('proxy','distal','tangent');
    set(gca,'FontSize',18);
    set(gca,'FontWeight','bold');
    xL.FontSize=18;
    xL.FontWeight='bold';
    yL.FontSize=18;
    yL.FontWeight='bold';

    


    figure(5)
    hold on
    plot(Rad_vec,proxy_meridian_splay);
    plot(Rad_vec,proxy_para_splay);
    plot(Rad_vec,proxy_total_splay);
    xL=xlabel('radial length from stalk');
    yL=ylabel('splay Proxy');
    legend('meridian','parallel','total');
    set(gca,'FontSize',18);
    set(gca,'FontWeight','bold');
    xL.FontSize=18;
    xL.FontWeight='bold';
    yL.FontSize=18;
    yL.FontWeight='bold';

    

    figure(6)
    hold on
    plot(Rad_vec,dist_meridian_splay);
    plot(Rad_vec,dist_para_splay);
    plot(Rad_vec,dist_total_splay);
    xL=xlabel('radial length from stalk');
    yL=ylabel('splay Distal');
    legend('meridian','parallel','total');
    set(gca,'FontSize',18);
    set(gca,'FontWeight','bold');
    xL.FontSize=18;
    xL.FontWeight='bold';
    yL.FontSize=18;
    yL.FontWeight='bold';

    

    figure(7)
    hold on
    plot(Rad_vec,0.5*(proxy_total_splay).^2);
    plot(Rad_vec,proxy_Gauss_splay);
    plot(Rad_vec,0.5*proxy_tilt.^2);
    xL=xlabel('radial length from stalk');
    yL=ylabel('Energy density proxy');
    legend('(Splay)^2/kappa_bar','saddle splay/kappa_bar_m','tilt^2/kappa_t');
    set(gca,'FontSize',18);
    set(gca,'FontWeight','bold');
    xL.FontSize=18;
    xL.FontWeight='bold';
    yL.FontSize=18;
    yL.FontWeight='bold';

    

    figure(8)
    hold on
    plot(Rad_vec,0.5*(dist_total_splay).^2);
    plot(Rad_vec,dist_Gauss_splay);
    plot(Rad_vec,0.5*dist_tilt.^2);
    xL=xlabel('radial length from stalk');
    yL=ylabel('Energy density distal');
    legend('(Splay)^2/kappa_bar','saddle splay/kappa_bar_m','tilt^2/kappa_t');
    set(gca,'FontSize',18);
    set(gca,'FontWeight','bold');
    xL.FontSize=18;
    xL.FontWeight='bold';
    yL.FontSize=18;
    yL.FontWeight='bold';

    
   
    figure(9)
    hold on
    plot(Rad_vec(1:MAX_INT-1),diff(tangent_angle)*lipid_length0); %-1 becuse of throwing out
    xL=xlabel('radial length from stalk');
    yL=ylabel('Delta tangent');
    set(gca,'FontSize',18);
    set(gca,'FontWeight','bold');
    xL.FontSize=18;
    xL.FontWeight='bold';
    yL.FontSize=18;
    yL.FontWeight='bold';
 end
    
    
    

   
    
end

