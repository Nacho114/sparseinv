function [I] = getIfromJ(A, J, debug)
% Get the I index set from A and J

A_hat = A(:, J);
I = find(sum(A_hat~=0, 2));

end