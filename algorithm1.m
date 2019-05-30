%A = [ 1.000 0.000 0.000; 5.000 6.000 7.000; 0.000 1.000 0.000];

% A = [ 1.000 0.000 0.000; 5.000 9.000 1000.000; 1.000 2.000 0.000];
A = spconvert(load('../dataset/orsirr_2.mtx'))
display(size(A))
[dim, ~] = size(A);
M = eye(dim);
iter = 5;
t=1;
debug = false;
for k = 1:dim
    sprintf('------------- Column %d -------------', k)
    % hardcoded part starts here
    % k = 2;
    % J = [2 3];
    J = [k];

    % M(3, 2) = 1;
    % hardcoded part ends here
    for x = 1:iter
        
        % compute error
        rem = A*M - eye(dim);
        col_err = norm(rem(:,k), 'fro');
        total_err = norm(rem, 'fro');
        sprintf('Before iter %d: \n total error is %.5f, column error is %.5f', x, total_err, col_err)
        
        [m_hat, r] = iterSingleColumn(A, M, J, k, debug);
        
        if debug
            sprintf('debug M')
            display(M)
            display(m_hat)
        end
        
        M(J, k) = m_hat;
        
        if debug
            display(M)
        end
        
        %J_star should be a row vector as J is 
        J_star = updateJ(A, r, debug, t)';
        
        if debug
            display(M)
            display(J)
        end
        
        J = union(J, J_star)
        
    end
    
    sprintf('-------------------------------')
    
end
