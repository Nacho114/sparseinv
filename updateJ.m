function [J_star] = updateJ(A, r)
    
    
    norm_r = norm(r)^2;
    eps=1e-10
    display(eps)
    
    L = find(r>eps)
    J_tilda = find(sum(A(L, :)~=0, 1))'
    %TODO: remove J from J_tilda
    
    % compute decrement by taking more things
    [dim, ~] = size(J_tilda);
    pj_sqrs = zeros(dim, 1);
    
    for idx = 1:dim
        j = J_tilda(idx)
        pj_sqrs(j) = norm_r - (r'*A(:, j))^2 / norm(A(:, j))^2;

    end
    
    [~, J_star] = min(pj_sqrs)
    