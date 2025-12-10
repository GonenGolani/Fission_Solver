function out = mypolyval(an,x)
% p(x)= a1+a2*x+a3*x^2+...+aN*x^N
    arr_size=size(an);
    number_of_lines=arr_size(1);
    number_of_N=arr_size(2);
    x_vector(1,1,:)=x;
    depth_vec(1,1,:)=ones(1,length(x));
    
    one_mat=ones(number_of_N);
    x_matrix=x_vector.*(tril(one_mat)-diag(one_mat(1,:)))+tril(one_mat)';
    x_matrix(1,1,:)=1;
    x_pow_n_vector=permute(prod(x_matrix,2),[1 3 2]);
    out=(an)*x_pow_n_vector;
end

