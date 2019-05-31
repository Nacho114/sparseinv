function [J_star] = updateJ(A, r, debug, t)
    
%     display(r)
    norm_r = norm(r)^2;
    eps=1e-10;
    
    L = find(abs(r)>eps);
    if debug
        display('J_tilda is ')
        display(A(L, :))
        display(abs(A(L, :)))
    end
    
%     J_tilda = find(sum(abs(A(L, :))~=0, 1))';
    
    [dim, ~] = size(A);
    J_tilda = [1:dim]';
    
    %TODO: remove J from J_tilda
    
    % compute decrement by taking more things
%     [dim, ~] = size(J_tilda);
    pj_sqrs = zeros(dim, 1);
    
    for idx = 1:dim
        j = J_tilda(idx);
        pj_sqrs(j) = norm_r - (r'*A(:, j))^2 / norm(A(:, j))^2;

    end
    
    if t == 1
        [~, J_star] = min(pj_sqrs);
    else
        [xs, idx] = sort(pj_sqrs);
        J_star = idx(1:t);
%         [~, J_star] = mink(pj_sqrs, t);
    end
    