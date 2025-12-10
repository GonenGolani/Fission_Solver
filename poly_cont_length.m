function out = Area_of_revolution(r_vector,z_vector,MAX_INT)

    dr_vec=diff(r_vector);
    dz_vec=diff(z_vector);
    ds_vec=(dr_vec.^2+dz_vec.^2).^0.5;
    ds_vec(MAX_INT)=ds_vec(MAX_INT-1);
    out=sum(2*pi*ds_vec.*r_vector);

end

