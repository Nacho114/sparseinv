A = spconvert(load('../dataset/memplus.mtx'));

[dim, ~] = size(A);

err_thresh_list = [0.5 0.3];
bic_thresh_list = [1e-8 1e-12 1e-16];

[~, nb_err_thresh] = size(err_thresh_list);
[~, nb_bic_thresh] = size(bic_thresh_list);

t=10;
maxiter=5;
use_par = true;

debug = false;

% set up parallel
if use_par
    num_workers = 2;
    parpool(num_workers);
else
    num_workers = 0;
end

b = A*ones(dim, 1);

fileID = fopen('../results/results_benzi_tuma_table4.txt', 'w');
fprintf(fileID, 'Results at error_threshold (eps) for various BICGSTAB tolerance are as follows\n\n\n');

for i = 1:nb_err_thresh
    
    M = eye(dim);
    err_thresh = err_thresh_list(i);
    fprintf(fileID, 'At error threshold: %.1f\n', err_thresh);
    [Mfinal] = spai(A, M, t, num_workers, err_thresh, maxiter, debug);
    
    is_correct = [];
    bic_relres = [];
    bic_conv_iter = [];

    for h = 1:nb_bic_thresh
        bic_thresh = bic_thresh_list(h);
        bic_default_iter = 500;
        [x_star, flag, curr_relres, curr_iter] = bicgstab(A*Mfinal, b, bic_thresh, bic_default_iter);
        x = Mfinal * x_star;


        curr_is_corr = all(abs(x - ones(dim, 1)) > 1e-6);
        is_correct = [is_correct curr_is_corr];
        bic_relres = [bic_relres curr_relres];
        bic_conv_iter = [bic_conv_iter curr_iter];
    
    end
    
    fprintf(fileID, '%16s %16s %16s\r\n\n','tolerance ','iterations', 'relres');
    fprintf(fileID, '%16.e %16.1f %16.e\r\n', [bic_thresh_list', bic_conv_iter', bic_relres']');
    fprintf(fileID', '\n\n');
end

fclose(fileID);
poolobj = gcp('nocreate');
delete(poolobj);