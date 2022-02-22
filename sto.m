%soft threshold operator
%input: a matrix Q and lambda
%output: a matrix B
function B=sto(Q,lambda)
B=sign(Q).*max(abs(Q)-lambda,0);
end