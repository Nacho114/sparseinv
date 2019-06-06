A = [ 1.000 0.000 0.000; 5.000 6.000 7.000; 0.000 1.000 0.000];
% A = spconvert(load('../dataset/sherman1.mtx'));

% A = 3*eye(3)
% A(1, 2) = 2

[dim, ~] = size(A);


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

% intialize M_final
Mfinal = eye(dim);
A_ = A;

sprintf('before preconditioning, error is %.5f\n', s, norm(A*Mfinal - eye(dim), 'fro'))

for s = 1:outer_max_iter
    % get a preconditioner for current A_ from MR
    [M] = mr(A_, eye(dim), num_workers, max_iter, eps, lfil, debug);
    
    % update A and Mfinal for restarts
    A_ = A_*M;
    Mfinal = Mfinal * M;
    
    sprintf('after outer iteration %d, error is %.5f\n', s, norm(A*Mfinal - eye(dim), 'fro'))
end

sprintf('in the end, after preconditioning, error is %.5f\n', s, norm(A*Mfinal - eye(dim), 'fro'))

% x_true = ones(dim, 1);
% b = A*x_true;
% [x,flag,relres,iter] = gmres(A, b);

poolobj = gcp('nocreate');
delete(poolobj);