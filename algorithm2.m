A = [ 1.000 0.000 0.000; 5.000 6.000 7.000; 0.000 1.000 0.000];
A = spconvert(load('../dataset/sherman1.mtx'));

% A = 3*eye(3)
% A(1, 2) = 2

[dim, ~] = size(A);

% intialize M
M = eye(dim);

% param
num_workers = 0;
max_iter = 224;
eps = 1e-8;
debug = false;
outer_max_iter = 1;
lfil = -1;
use_par = true;

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
    [Mfinal] = minres(A_, Mfinal, num_workers, max_iter, eps, lfil, debug);
    
    display(norm(A*Mfinal - eye(dim), 'fro'))
end


poolobj = gcp('nocreate');
delete(poolobj);