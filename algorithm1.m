%A = [ 1.000 0.000 0.000; 5.000 6.000 7.000; 0.000 1.000 0.000];

% A = [ 1.000 0.000 0.000; 5.000 9.000 1000.000; 1.000 2.000 0.000];
A = spconvert(load('../dataset/orsirr_2.mtx'))
% A=[1 8 0 5 7; 8 10 1 2 0; 0 0 0 1 0; 1 2 0 0 3; 2 1 -1 0 0.2];
display(size(A))
[dim, ~] = size(A);
M = eye(dim);
% iter = 5;
maxiter=5;
t=1;
debug = false;
parpool(2);
err_thresh = -1;
Id = eye(dim);
withf=true;

tic
parfor k = 1:dim
% for k = 1:dim
    sprintf('------------- Column %d -------------', k)
    % hardcoded part starts here
    % k = 2;
    % J = [2 3];
    J = [k];
    m_final = M(:, k);
    % M(3, 2) = 1;
    % hardcoded part ends here
    
    for x = 1:maxiter
        
        % compute error
        col_err= norm(A*m_final - Id(:, k));

        sprintf('Before iter %d: column error is %.5f \n', x, col_err)
        
        [m_hat, r] = iterSingleColumn(A, J, k, debug);
        
        if norm(r) < err_thresh
            sprintf('exited at iteration %d for column %d, norm is %.5f ', x, k, norm(r))
            break
        end
        
        if debug
            sprintf('debug M')
%             display(M)
            display(m_hat)
        end
        
        m_final(J) = m_hat;
        
        if debug
%             display(M)
        end
%         f = @() updateJ(A, r, debug, t);
%         timeit(f)
        %J_star should be a row vector as J is 
        J_star = updateJ(A, J, r, debug, t, withf)';
        
        if debug
%             display(M)
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
