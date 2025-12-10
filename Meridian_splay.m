function out = Meridian_splay(dr,Tangent_Angle,tilt_angle,MAX_INT)
     delta_angle=diff(tilt_angle+Tangent_Angle,1,2);

     delta_angle(:,MAX_INT)=0;

     out=cos(Tangent_Angle).*cos(tilt_angle).*delta_angle./dr;
end

