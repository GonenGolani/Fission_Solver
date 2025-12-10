function [New_Shape_polynum,New__Total_energy] = Minimize_N_order_one_iter(min_dim,Shape_polynum,Boundary_Conditions,Physical_Conditions,minimaztion_res,search_range,Coeff_n)
    % Coeff_n holds the coeficieant order that is minimized
    % Shape_polynum:
    %   MP - mid-plane polynom
    %   proxy_angle - proximal monolayer lipid director polynom
    %   distal_angle - distal monolayer lipid director polynom
    %   min_dim - minize the coeff of given dimension:"MP", "proxy" or "distal" 
    
    %find last iteration coeffcient value
    MP_last_an_value=Shape_polynum.MP(Coeff_n);
    proxy_last_an_value=Shape_polynum.proxy_angle(Coeff_n);
    distal_last_an_value=Shape_polynum.distal_angle(Coeff_n);
    New_Shape_polynum=Shape_polynum;
    temp_shape_polynum=Shape_polynum;
    %serach range
    coeff_length=minimaztion_res*search_range;
    
    if min_dim=="MP"
        MP_coeff_vector=linspace(MP_last_an_value-coeff_length/2,MP_last_an_value+coeff_length/2,search_range); 
        MP_coeff_vector(search_range+1)=MP_last_an_value;
        temp_shape_polynum.MP=Shape_polynum.MP.*ones(1,search_range+1)';
        temp_shape_polynum.MP(:,Coeff_n)=MP_coeff_vector;
        

    end

    
    if min_dim=="proxy"
        proxy_coeff_vector=linspace(proxy_last_an_value-coeff_length/2,proxy_last_an_value+coeff_length/2,search_range);   
        proxy_coeff_vector(search_range+1)=proxy_last_an_value;
        temp_shape_polynum.proxy_angle=Shape_polynum.proxy_angle.*ones(1,search_range+1)';
        temp_shape_polynum.proxy_angle(:,Coeff_n)=proxy_coeff_vector;
    end
    
    
    if min_dim=="distal"
        distal_coeff_vector=linspace(distal_last_an_value-coeff_length/2,distal_last_an_value+coeff_length/2,search_range); 
        distal_coeff_vector(search_range+1)=distal_last_an_value;
        temp_shape_polynum.distal_angle=Shape_polynum.distal_angle.*ones(1,search_range+1)';
        temp_shape_polynum.distal_angle(:,Coeff_n)=distal_coeff_vector;
    end
    
    
    Energy_matrix=Stalk_energy(temp_shape_polynum,Boundary_Conditions,Physical_Conditions);
    
    Energy_matrix(Energy_matrix < -100) = NaN;


    [New__Total_energy,Index]=min(Energy_matrix);
    
    
    if min_dim=="MP"
        New_Shape_polynum.MP(Coeff_n)=MP_coeff_vector(Index);        
    end
    if min_dim=="proxy"
        New_Shape_polynum.proxy_angle(Coeff_n)=proxy_coeff_vector(Index);        
    end
    if min_dim=="distal"
        New_Shape_polynum.distal_angle(Coeff_n)=distal_coeff_vector(Index);        
    end
    
end

