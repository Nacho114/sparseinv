A = spconvert(load('../dataset/orsirr_2.mtx'));

[dim, ~] = size(A);

err_thresh_list = [0.6, 0.5, 0.4, 0.3, 0.2];
[~, nb_err_thresh] = size(err_thresh_list);

t=5;
maxiter=50;
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

criteria_values = zeros(5, 5);
criteria = ["eps", "F", "L2", "L1", "cond2", "nz(M)/nz(a)"];

Id = eye(dim);
for i = 1:nb_err_thresh

    M = Id;
    err_thresh = err_thresh_list(i);
    [Mfinal] = spai(A, M, t, num_workers, err_thresh, maxiter, debug);
    
    % Frobenius norm 
    criteria_values(1, i) = norm(A*Mfinal - Id, 'fro');
    
    % 2 norm 
    criteria_values(2, i) = norm(A*Mfinal - Id, 2);
    
    % 1 norm 
    criteria_values(3, i) = norm(A*Mfinal - Id, 1);
    
    % cond2
    criteria_values(4, i) = cond(A*Mfinal);
    
    % nz(M)/nz(a)
    criteria_values(5, i) = nnz(Mfinal)/nnz(A);

end

fileID = fopen('../results/grotehuckle_table2.txt', 'w');

fprintf(fileID, '%6s %6s %6s %6s %6s %13s \n', criteria);
fprintf(fileID, '%6.2f %6.2f %6.2f  %6.2f  %6.2f  %6.2f\n', [err_thresh_list' criteria_values(1,:)' criteria_values(2, :)' criteria_values(3,:)'  criteria_values(4,:)' criteria_values(5,:)' ]');  


poolobj = gcp('nocreate');
delete(poolobj);