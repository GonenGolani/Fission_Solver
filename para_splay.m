function out = para_splay(radial_vector,Tangent_Angle,tilt_angle,dr)
    out=sin(Tangent_Angle+tilt_angle)./(radial_vector+dr);
end
