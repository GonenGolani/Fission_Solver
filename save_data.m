function [nothing] = save_data(new_energy,new_shape_polynum,new_Boundary_Conditions,location)
    cd(location);

        save('Theta_Phi_Energy.txt','new_energy','-ascii');
    
        save('shape_polynum.mat','new_shape_polynum');            

        save('Boundary_Conditions.mat','new_Boundary_Conditions');            
    
    nothing=0;
end

