function [a_n_out] = Midplane_poly_assign(Ang_0,Ang_R,z0,z1,Rad,a_n_in)
% Input= vector a_n_in with ploynum coeffients with length N.
% first 4 does not matter becuse found from boundary conditions
% boundary condtions:
% outer edge radius: Rad
% angle at rho=0: Ang_0 at 0 in deg
% angle at=rho=R: Ang_R at 0 in deg
% hight at rho=0: z0 
% hight at rho=Rad: z1 
% output= assign to 4 first places values using. give vector a_n_out 

an_length=length(a_n_in);

a_n_out=a_n_in;

a_n_out(1)=z0;
a_n_out(2)=tan(Ang_0);

% solve equation for n=3
a3_vector=(4*ones(1,an_length)-linspace(1,an_length,an_length)).*a_n_out;
a_n_out(3)=(3*z1-Rad*tan(Ang_R)-mypolyval_ex_n(a3_vector,Rad,3))/Rad^2;

% solve equation for n=4
a4_vector=(3*ones(1,an_length)-linspace(1,an_length,an_length)).*a_n_out;
a_n_out(4)=-(2*z1-Rad*tan(Ang_R)-mypolyval_ex_n(a4_vector,Rad,4))/Rad^3;



end

