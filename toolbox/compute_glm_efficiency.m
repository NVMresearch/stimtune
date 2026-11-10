function eff = compute_glm_efficiency(X,c)
% Calculates the efficiency of a contrast
% 
% X - design matrix
% c - contrast vector or a matrix with one contrast per column

% If a vector is provided then ensure it's a column rather than a row
if isvector(c)
    c = c(:);
end

% Pad contrast if needed
if size(c,1) < size(X,2)
    npad = size(X,2) - size(c,1);
    c = padarray(c,[npad 0],'post');
end

% Efficiency from Eq. 3 of Friston et al. (1999), NeuroImage
eff = 1 / trace(c' * spm_pinv(X' * X) * c);

end