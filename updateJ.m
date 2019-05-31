function [J_star] = updateJ(A, J, r, debug, t, withf, rem_below_avg)
    
%     display(r)
    norm_r_sq = norm(r)^2;
    eps=1e-10;
    inf = 1e16;
    
    [dimA, ~] = size(A);
    pj_sqrs = ones(dimA, 1)*inf;
    
    L = find(abs(r)>eps);
    if debug
        display('J_tilda is ')
        display(A(L, :))
        display(abs(A(L, :)))
    end
    
    if withf
        % for example, J_tilda = [2, 3, 4], J=[2,3]
        J_tilda = find(sum(abs(A(L, :))~=0, 1))';
        % remove J from J_tilda
        J_tilda = setdiff(J_tilda, J);
        
    else
        J_tilda = [1:dimA]';
    end
    
    [dimJT, ~] = size(J_tilda); 
    
    % compute decrement by taking more things    
    sum_pj_sq = 0
    for idx = 1:dimJT
        j = J_tilda(idx);
        pj_sqrs(j) = norm_r_sq - (r'*A(:, j))^2 / norm(A(:, j))^2;
        sum_pj_sq = sum_pj_sq + pj_sqrs(j);
    end
    
    if rem_below_avg
        avg_pj_sq = sum_pj_sq/dimJT;
        pj_sqrs(pj_sqrs>avg_pj_sq) = inf;
    end
    
    if t == 1
        [~, J_star] = min(pj_sqrs);
    else
        t = min(t, sum(pj_sqrs < inf));
        display(t)
        [xs, idx] = sort(pj_sqrs);
        J_star = idx(1:t);
        %[~, J_star] = mink(pj_sqrs, t);
    end
    