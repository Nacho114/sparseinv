A = [ 1.000 0.000 0.000; 5.000 6.000 7.000; 0.000 1.000 0.000];
% A = spconvert(load('../dataset/sherman1.mtx'));

% A = 3*eye(3)
% A(1, 2) = 2

[dim, ~] = size(A);

% intialize M
M = eye(dim);

% param
num_workers = 0;
max_iter = 1;
eps = 1e-8;
debug = false;
outer_max_iter = 4;
lfil = -1;
use_par = false;

% set up parallel
if use_par
    num_workers = 2;
    parpool(num_workers);
else
    num_workers = 0;
end

Mfinal = M;
M_ = M;
A_ = A;
for s = 1:outer_max_iter
    A_ = A_*Mfinal;
    M_ = M_ * Mfinal;
    [Mfinal] = mr(A_, eye(dim), num_workers, max_iter, eps, lfil, debug);
    
    display(norm(A*M_ - eye(dim), 'fro'))
end

% x_true = ones(dim, 1);
% b = A*x_true;
% [x,flag,relres,iter] = gmres(A, b);

poolobj = gcp('nocreate');
delete(poolobj);