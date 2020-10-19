function out = goveqs_basis3(t, in, M, i, s, p, r, agg, sel)

invec = in(1:i.nstates);
out   = zeros(length(in),1);

% --- Conctruct the linear terms
out(1:i.nstates) = M.lin*invec;

% Get new infections
lam = M.lam*invec;
newinfs = lam.*invec(s.S);

out(s.S) = out(s.S) - newinfs(:);
out(s.E) = out(s.E) + newinfs(:);

% Also reinfections
newinfs2 = p.imm*lam.*invec(s.R);
out(s.R) = out(s.R) - newinfs2(:);
out(s.E) = out(s.E) + newinfs2(:);

% Implement deaths
morts = M.mortvec.*invec;
out(1:i.nstates) = out(1:i.nstates) - morts;

% Get the auxiliaries
out(i.aux.inc)   = agg.inc*(sel.inc.*M.lin)*invec;
out(i.aux.hosp)  = agg.hosp*(sel.hosp.*M.lin)*invec;
out(i.aux.mort)  = agg.mor*morts; 
out(i.aux.hosp2) = sum((sel.hosp2.*M.lin)*invec);