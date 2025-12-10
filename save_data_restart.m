function [nothing] = save_data_restart(new_energy,new_shape_polynum,new_Boundary_Conditions,location)
    cd(location);
    if exist('out_energy.txt')==2
        old_energy=importdata('out_energy.txt');
        saved_data=[old_energy;new_energy];
        save('out_energy.txt','saved_data','-ascii');       
    end
    
    if exist('shape_polynum.mat')==2
        old_shape_polynum=importdata('shape_polynum.mat');
        saved_shape_polynum=[old_shape_polynum;new_shape_polynum'];  
        save('shape_polynum.mat','saved_shape_polynum'); 
    end
    
    if exist('Boundary_Conditions.mat')==2
        old_Boundary_Conditions=importdata('Boundary_Conditions.mat');
        saved_Boundary_Conditions=[old_Boundary_Conditions;new_Boundary_Conditions'];  
        save('Boundary_Conditions.mat','saved_Boundary_Conditions'); 
    end
    
    if exist('out_energy.txt')==0
        save('out_energy.txt','new_energy','-ascii');
    end
    
    if exist('shape_polynum.mat')==0
        save('shape_polynum.mat','new_shape_polynum');            
    end

    if exist('Boundary_Conditions.mat')==0
        save('Boundary_Conditions.mat','new_Boundary_Conditions');            
    end
    
    nothing=0;
end

