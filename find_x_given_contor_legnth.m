function out = find_x_given_contor_legnth(angle_vec,S_max,MAX_INT)

    ds=S_max/MAX_INT;
    dx = ds*cos(angle_vec);
    out=sum(ds);

end

