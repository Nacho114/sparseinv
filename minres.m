function [Mfinal] = minres(A, M, num_workers, maxiter, eps, lfil, debug)    

[dim, ~] = size(A);
E = eye(dim);
Mfinal = M;

    % for each column
    parfor (j = 1:dim, num_workers)
%     for j = 1:dim
        ej = E(:, j);
        mj = M(:, j);
        
        % for max inner iter ni
        for i = 1:maxiter
            % compute residual
            rj = ej - A*mj;
            
            if debug
                sprintf('residual column %d iter %d res %.5f \n', j, i, norm(rj))
            end
            
            if norm(rj) < eps
                break
            end
            
            
            % TODO, use SAXPY? When A is sparse, see chow saad
            

            Arj = A*rj;
            
            alphaj = dot(rj, Arj) / dot(Arj, Arj);
            % minimise residual norm
            mj = mj + alphaj*rj;
            
            
            % apply numerical dropping
            if lfil >= 1 

                % get score
                pj = -2 * mj .* (A'* rj) + mj.^2;
                % sort by score
                [~,I] = sort(pj, 'descend');
                % keep topk items, and set rest to zero
                mj_ = zeros(size(mj));
                topk = I(1:lfil);
                mj_(topk) = mj(topk);
                mj = mj_;
                
            end
           
        end
        
        Mfinal(:,j) = mj;

    end
    
end