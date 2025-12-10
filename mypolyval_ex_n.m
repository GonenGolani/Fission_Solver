function out = mypolyval_ex_n(an,x,n)
% p(x)= a1+a2*x+a3*x^2+...+aN*x^
% take out ac n=n ->0

    arr_size=size(an);
    number_of_lines=arr_size(1);
    number_of_N=arr_size(2);
    
    one_mat=ones(number_of_N);
    x_matrix=x*(tril(one_mat)-diag(one_mat(1,:)))+tril(one_mat)';
    x_matrix(1,1)=1;
    x_pow_n_vector=prod(x_matrix');

    an(1:arr_size(1),n)=0;
    out=an*x_pow_n_vector';
end
