%% Caution:: order matters
% compute the relative error given two vectors
function out = rele(z,x)
out = norm(z - x)/norm(x);
