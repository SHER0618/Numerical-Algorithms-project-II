%Singular value threshold operator
%input: a matrix Q and lambda
%output: a matrix B
function B=svto(Q,lambda)
[S,V,D]=svd(Q);
B=S*sto(V,lambda)*D';
end