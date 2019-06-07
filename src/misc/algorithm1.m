% Parameters
maxiter=50;
t=5;
debug = false;
% matrix_type options include ex1, ex2
matrix_type = 'orsirr2' 
use_par = true;
% matrix_type = 'ex1'
% use_par = false;
err_thresh = 0.6;
withf=true;
rem_below_avg = true;
time_method = false;

% intialize A
if strcmp(matrix_type, 'ex1')
    A = [ 1.000 0.000 0.000; 5.000 6.000 7.000; 0.000 1.000 0.000];
elseif strcmp(matrix_type, 'ex2')
    A=[1 8 0 5 7; 8 10 1 2 0; 0 0 0 1 0; 1 2 0 0 3; 2 1 -1 0 0.2];
elseif strcmp(matrix_type, 'orsirr2')
    B = load('../../dataset/orsirr_2.mtx')
    A = spconvert(load('../../dataset/orsirr_2.mtx'));
end

display(size(A))
[dim, ~] = size(A);

% intialize M
M = eye(dim);
Id = eye(dim);

% set up parallel
if use_par
    num_workers = 2;
    parpool(num_workers)
else
    num_workers = 0;
end

tic

parfor (k = 1:dim, num_workers)
    sprintf('------------- Column %d -------------', k)
    
    J = [k];
    m_final = M(:, k);
    
    for x = 1:maxiter
        
        % compute error
        col_err= norm(A*m_final - Id(:, k));
        sprintf('Before iter %d: column error is %.5f \n', x, col_err)
        
        [m_hat, r] = iterSingleColumn(A, J, k, debug);
        
        if debug
            sprintf('debug M')
            display(m_hat)
        end
        
        m_final(J) = m_hat;
        
        % break if norm is below the 'eps' error threshold
        if norm(r) < err_thresh
            sprintf('exited at iteration %d for column %d, norm is %.5f ', x, k, norm(r))
            break
        end
        
        if time_method
            f = @() updateJ(A, r, debug, t);
            timeit(f)
        end
        
        %J_star should be a row vector as J is 
        J_star = updateJ(A, J, r, debug, t, withf, rem_below_avg)';
        
        if debug
            display(J)
        end
        
        J = union(J, J_star);
                
    end
    
    col_err= norm(A*m_final - Id(:, k));
    sprintf('Column %d : At the end of iterations, error is %.5f \n', x, col_err)
    
    M(:, k) = m_final;
    
    sprintf('-------------------------------')
    
end

total_error = norm(A*M - eye(dim), 'fro')
toc
poolobj = gcp('nocreate');
delete(poolobj);
