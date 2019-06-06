A = spconvert(load('../dataset/memplus.mtx'));

[dim, ~] = size(A);

err_thresh_list = [0.5 0.3];
[~, nb_err_thresh] = size(err_thresh_list);

t=10;
maxiter=5;
use_par = true;
% use_par = false;
debug = false;

% set up parallel
if use_par
    num_workers = 2;
    parpool(num_workers);
else
    num_workers = 0;
end

is_correct = [];
bic_relres = [];
bic_conv_iter = [];

for i = 1:nb_err_thresh

    M = eye(dim);
    err_thresh = err_thresh_list(i);
    [Mfinal] = spai(A, M, t, num_workers, err_thresh, maxiter, debug);

    b = A*ones(dim, 1);

    bic_thresh = 1e-6;
    bic_default_iter = 500;
    [x_star, flag, curr_relres, curr_iter] = bicgstab(A*Mfinal, b, bic_thresh, bic_default_iter);
    x = Mfinal * x_star;
    
    
    curr_is_corr = all(abs(x - ones(dim, 1)) > 1e-6);
    is_correct = [is_correct curr_is_corr];
    bic_relres = [bic_relres curr_relres];
    bic_conv_iter = [bic_conv_iter curr_iter];
    

end

poolobj = gcp('nocreate');
delete(poolobj);