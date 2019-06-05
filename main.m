% A = [ 1.000 0.000 0.000; 5.000 6.000 7.000; 0.000 1.000 0.000];
A = spconvert(load('../dataset/orsirr_2.mtx'));

[dim, ~] = size(A);
M = eye(dim);
t=5;
maxiter=50;
use_par = true;
% use_par = false;
err_thresh = 0.3;

debug = false;


[Mfinal] = spai(A, M, t, use_par, err_thresh, maxiter, debug);

b = A*ones(dim, 1);

[x_star, flag, relres, iter] = bicgstab(A*Mfinal, b, 1e-8);

Mfinal * x_star
flag
relres
iter

