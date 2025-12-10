function out = Area_of_revolution_vector(r_vector,z_vector,MAX_INT)

    dr_vec=diff(r_vector,1,2);
    dz_vec=diff(z_vector,1,2);
    ds_vec=(dr_vec.^2+dz_vec.^2).^0.5;
    ds_vec(:,MAX_INT)=ds_vec(:,MAX_INT-1);
    out=(2*pi*ds_vec.*r_vector);

end

