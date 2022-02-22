%APG method
%input:a matrix A;
%output: a matrix L and S s.t. A=L+S;
function [L,S]=APG(A,iter_times)
    [m,n]=size(A);
    lambda=10/sqrt(m*n);
    Lf=2;
    ita=0.9;
    miu=ita*norm(A,2);
    miu_min=1e-9*miu;
    t_last=1;
    t_new=1;
    L_last=zeros(m,n);
    S_last=zeros(m,n);
    L_new=A;
    S_new=S_last;
    for i=1:iter_times
        YL=L_new+(t_last-1)*(L_new-L_last)/t_new;
        YS=S_new+(t_last-1)*(S_new-S_last)/t_new;
        t_last=t_new;
        L_last=L_new;
        S_last=S_new;
        L_new=svto(YL+(A-YL-YS)/Lf,miu/2);
        S_new=sto(YS+(A-YL-YS)/Lf,lambda*miu/2);
        t_new=(1+sqrt(1+4*t_last^2))/2;
        miu=max(ita*miu,miu_min);
    end
    L=L_new;
    S=S_new;
end