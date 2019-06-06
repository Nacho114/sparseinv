function [Mfinal] = spai(A, M, t, num_workers, err_thresh, maxiter, debug)

% Parameters

withf=true;
rem_below_avg = true;

[dim, ~] = size(A);


parfor (k = 1:dim, num_workers)

    J = [k];
    m_final = M(:, k);
    
    for x = 1:maxiter
        
        [m_hat, r] = iterSingleColumn(A, J, k, debug);
        
        m_final(J) = m_hat;
        
        % break if norm is below the 'eps' error threshold
        if norm(r) < err_thresh
            break
        end
        
  
        %J_star should be a row vector as J is 
        J_star = updateJ(A, J, r, debug, t, withf, rem_below_avg)';

        J = union(J, J_star);
                
    end
    
    M(:, k) = m_final;
    
    
end

Mfinal = M;
