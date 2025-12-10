function [a_n_out] = poly_force_boundary(y0,y1,dydx0,dydx1,Rad1,Rad0,a_n_in)

% Input= vector a_n_in with ploynum coeffients with length N.
% first 4 does not matter becuse found from boundary conditions
% boundary condtions:
% outer edge radius: Rad
% angle at rho=Rad0: Ang_0 at 0 in deg
% angle at=rho=Rad: Ang_R at 0 in deg
% hight at rho=Rad0: z0 
% hight at rho=Rad: z1 
% output= assign to 4 first places values using. give vector a_n_out 
% solved by inversing the matrix



A_mat=[1 Rad0 Rad0^2 Rad0^3 ; 1 Rad1 Rad1^2 Rad1^3 ; 0 1 2*Rad0 3*Rad0^2 ; 0 1 2*Rad1 3*Rad1^2];
A_inv=inv(A_mat);


arr_size=size(a_n_in);
number_of_lines=arr_size(1);
an_length=arr_size(2);

thorw_ot_first_4_vector=ones(1,an_length);
thorw_ot_first_4_vector(1:4)=0;
% RHS side solution vector
RHS_vector(:,1)=y0-mypolyval(a_n_in.*thorw_ot_first_4_vector,Rad0);
RHS_vector(:,2)=y1-mypolyval(a_n_in.*thorw_ot_first_4_vector,Rad1);
RHS_vector(:,3)=dydx0-mypolyder(a_n_in.*thorw_ot_first_4_vector,Rad0);
RHS_vector(:,4)=dydx1-mypolyder(a_n_in.*thorw_ot_first_4_vector,Rad1);

%out an vector 1-4
a_1_to_4=A_inv*RHS_vector';

a_n_out=a_n_in;
a_n_out(:,1:4)=a_1_to_4';


end

