%IT method
%input: a matrix A;
%output: a matrix L and S s.t. A=L+S;
function [L,S]=IT(A,iter_times)
    [m,n]=size(A);
    lambda=10/sqrt(m*n);
    tau=0.9*norm(A,2);
    sigema=0.5;
    Y=zeros(m,n);
    for i=1:iter_times
        L=svto(Y,tau);
        S=sto(Y,lambda*tau);
        Y=Y+sigema*(A-L-S);
    end
end