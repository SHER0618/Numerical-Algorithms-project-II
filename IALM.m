%IALM method
%input:a matrix A
%output:a matrix L and S s.t. A=L+S;
function [L,S]=IALM(A,iter_times)
    [m,n]=size(A);
    lambda=1/sqrt(m*n);
    mu=10*lambda;
    S=zeros(m,n);
    Y=zeros(m,n);
    for i=1:iter_times
        L=svto(A-S+Y/mu,1/mu);
        S=sto(A-L+Y/mu,lambda/mu);
        Y=Y+mu*(A-L-S);
    end
end
