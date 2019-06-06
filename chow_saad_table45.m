
A = spconvert(load('../dataset/sherman1.mtx'));

[dim, ~] = size(A);


max_iter = 1;
eps = 1e-8;
debug = false;
outer_max_iter_list = [1, 2, 3, 4, 5];
lfil = 10;
use_par = true;

% set up parallel
if use_par
    num_workers = 2;
    parpool(num_workers);
else
    num_workers = 0;
end

% info
relres_list = [];
iter_list = [];


for max_outer_iter = outer_max_iter_list
    M = eye(dim);
    Mfinal = M;
    M_ = M;
    A_ = A;

    for s = 1:max_outer_iter
        A_ = A_*Mfinal;
        M_ = M_ * Mfinal;
        [Mfinal] = mr(A_, E, num_workers, max_iter, eps, lfil, debug);

        display(norm(A*M_ - eye(dim), 'fro'))
    end
    
    x_true = ones(dim, 1);
    b = A*x_true;
    [x,flag,relres,iter] = gmres(A, b);
    relres_list = [relres_list relres];
    iter_list = [iter_list iter];

end



poolobj = gcp('nocreate');
delete(poolobj);