% A = [ 1.000 0.000 0.000; 5.000 6.000 7.000; 0.000 1.000 0.000];
A = spconvert(load('../dataset/sherman1.mtx'));

% A = 3*eye(3)
% A(1, 2) = 2

[dim, ~] = size(A);


% param
inner_iter_mr = 1;
normalize=true;
% below this eps, we don't run an inner iteration 
% prevents nan
eps = 1e-8;
debug = false;
outer_max_iter = 1;
lfil = 10;
use_par = true;
alpha=1;

% set up parallel
if use_par
    num_workers = 2;
    parpool(num_workers);
else
    num_workers = 0;
end

if normalize
    display(cond(A))
    col_norms = sqrt(diag(A'*A))' + 1e-8;
    A = A ./ col_norms;
    A_tmp = A;
    display(cond(A))
end


A_ = A;
% intialize M_final
% if normalize
%     Mfinal = A_';
% else
Mfinal = eye(dim);
% end

if normalize
    Minit = alpha * A_'; 
else
    Minit = eye(dim);    
end
    
sprintf('before preconditioning, error is %.5f', norm(A*Mfinal - eye(dim), 'fro'))

for s = 1:outer_max_iter
    % get a preconditioner for current A_ from MR
    
    
    [M] = mr(A_, Minit, num_workers, inner_iter_mr, eps, lfil, debug);
    
    % update A and Mfinal for restarts
%     A_ = A_*M;
%     Mfinal = Mfinal * M;
    Minit = M;
    Mfinal  = M;
    sprintf('after outer iteration %d, error is %.5f', s, norm(A*Mfinal - eye(dim), 'fro'))
end

sprintf('in the end, after preconditioning, error is %.5f', norm(A*Mfinal - eye(dim), 'fro'))

x_true = ones(dim, 1);
b = A*x_true;
% run gmres(20)
gmres_max_iter = 500;
gmres_tol = 1e-5;
gmres_restarts = 20;
[x_star, flag, relres, iter] = gmres(A*Mfinal, b, gmres_restarts, gmres_tol, gmres_max_iter);
x = Mfinal * x_star;

display(iter)
display(relres)

poolobj = gcp('nocreate');
delete(poolobj);