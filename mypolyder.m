function out = mypolyder(an,x)
% in: p(x)= a1+a2*x+a3*x^2+...+aN*x^(N-1)
%out: dpdx = 0 + a2 +2*a3*x+...+aN*(N-1)*x^(N-2)
    arr_size=size(an);
    number_of_lines=arr_size(1);
    number_of_N=arr_size(2);
    x_vector(1,1,:)=x;
    depth_vec(1,1,:)=ones(1,length(x));
    
    one_mat=ones(number_of_N);
    one_mat_bigger=ones(number_of_N+1);
    diag_min_1_matrix=diag(one_mat(1,:),-1);
    x_matrix=x_vector.*(tril(one_mat_bigger)-diag_min_1_matrix-diag(one_mat_bigger(1,:)));
    x_matrix=x_matrix+diag_min_1_matrix+tril(one_mat_bigger)';
       
    x_matrix(number_of_N+1,:,:)=[];
    x_matrix(:,number_of_N+1,:)=[];
    x_matrix(1,1,:)=0;
    x_pow_n_vector=permute(prod(x_matrix,2),[1 3 2]);
    n_matrix=diag(linspace(1,number_of_N,number_of_N)-1);
    out=(an*n_matrix)*x_pow_n_vector;
    
end




