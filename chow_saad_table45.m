% Code to generate table 4.5 of Chow Saad

% seed
% rng(1)

path = "../dataset/";
datasets = ["sherman1" "sherman3" "saylr3"];
ext = ".mtx";

% general param
eps = 1e-8; % to prevent nan, stabalize normalisation
outer_max_iters = [1 2 3 4 5];
lfil = 10;
alpha = 1;

inner_iter_mr = 1;
normalize=true;

use_par = true;
debug = false;

% set up parallel
if use_par
    num_workers = 2;
    parpool(num_workers);
else
    num_workers = 0;
end

all_iter = zeros(3, 5);
all_relres = zeros(3, 5);
idx = 1;
for dataset = datasets
    
    A = spconvert(load(path+dataset+ext));
    [dim, ~] = size(A);

    if normalize
        display(cond(A))
        col_norms = sqrt(diag(A'*A))' + 1e-8;
        A = A ./ col_norms;
        A_tmp = A;
        display(cond(A))
    end

    for outer_max_iter = outer_max_iters


        A_ = A;
        Minit = alpha * A_'; 

        for s = 1:outer_max_iter
            % get a preconditioner for current A_ from MR

            [M] = mr(A_, Minit, num_workers, inner_iter_mr, eps, lfil, debug);

            Minit = M;
            Mfinal = M;
            
        end

        x_true = ones(dim, 1);
        b = A*x_true;
        % run gmres(20)
        gmres_max_iter = 500;
        gmres_tol = 1e-5;
        gmres_restarts = 20;
        [x_star, flag, relres, iter] = gmres(A*Mfinal, b, gmres_restarts, gmres_tol, gmres_max_iter);
        x = Mfinal * x_star;
        
        all_iter(idx, outer_max_iter) = iter(1) * iter(2);
        all_relres(idx, outer_max_iter) = relres;
        
    end
    idx = idx + 1;
end

fileID = fopen('chow_saad_table45_result.txt', 'w');

fprintf(fileID, 'iter\n');

for idx = 1:3
    curr_iter = all_iter(idx, :);
    fprintf(fileID, '%s ', datasets(idx));
    if idx == 3
        fprintf(fileID, '  ');
    end
    fprintf(fileID, '%d ', curr_iter');
    fprintf(fileID, '\n');
end

fprintf(fileID, 'relres\n');

for idx = 1:3
    curr_iter = all_relres(idx, :);
    fprintf(fileID, '%s ', datasets(idx));
    if idx == 3
        fprintf(fileID, '  ');
    end
    fprintf(fileID, '%d ', curr_iter');
    fprintf(fileID, '\n');
end

poolobj = gcp('nocreate');
delete(poolobj);