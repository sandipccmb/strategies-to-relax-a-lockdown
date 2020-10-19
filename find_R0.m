function R0 = find_R0(M, S0, infectious_inds, s_E)

% M  is model matrix
% S0 is initial number susceptible
% infectious_inds is the set of indices for the setting to be assessed
% s_E is the set of indices of E

Mlin   = M.lin;
Mlam   = M.lam;

% Adjust for the number susceptible
Mlam = Mlam.*repmat(S0,1,size(Mlam,2));

% Construct F matrix
mat = zeros(size(Mlin));
mat(s_E,:) = Mlam;
F = mat(infectious_inds, infectious_inds);

% Construct V matrix
m = -(Mlin - diag(M.mortvec));
V = m(infectious_inds, infectious_inds);

% Get R0
R0 = max(eigs(F*inv(V)));