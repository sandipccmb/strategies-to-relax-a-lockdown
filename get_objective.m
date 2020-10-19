function [out, aux] = get_objective(x, r, p, i, s, deaths_sm)

R0             = x(1);
cbeta          = 1-x(2);
p_death_report = x(3);
t_init         = x(4);

% --- Set up the models ---------------------------------------------------

% Without lockdown
r0 = r; p0 = p;
r0.beta = r.beta*R0;
M0 = make_model3(p0, r0, i, s, gps, prm);

% With lockdown
r1 = r0; p1 = p;
r1.beta = r.beta*R0*cbeta;
M1 = make_model3(p1, r1, i, s, gps, prm);


% --- Run the model -------------------------------------------------------
init = zeros(1,i.nx); seed = 10;
tmp = prm.N'; tmp2 = tmp(:)';
init(s.S) = tmp2;
init(i.IS1.(gps.geo{1}).ad) = seed; init(i.S.(gps.geo{1}).ad) = init(i.S.(gps.geo{1}).ad) - seed;

% --- Pre-lockdown
geq0 = @(t,in) goveqs_basis3(t, in, M0, i, s, p0, r0, agg, sel);
[t0,soln0] = ode15s(geq0, [t_init t_lockdown], init, odeset('Nonnegative',1:i.nx,'Refine',32));

% --- Post-lockdown
geq1 = @(t,in) goveqs_basis3(t, in, M1, i, s, p1, r1, agg, sel);
[t1,soln1] = ode15s(geq1, [t_lockdown t_end], soln0(end,:), odeset('Nonnegative',1:i.nx,'Refine',32));

ta    = [t0; t1(2:end)];
solna = [soln0; soln1(2:end,:)];

t     = ceil(t_init):t_end;
soln  = interp1(ta, solna, t);

dsol  = diff(soln,[],1);
% Find daily cases
sim_cases = sum(dsol(:,i.aux.inc),2);
% Find daily deaths
sim_morts = sum(dsol(:,i.aux.mort),2);


% Pad the data, or solution
comp = length(sim_morts) - length(deaths_sm);
if comp > 0
    deaths_sm = [zeros(comp,1); deaths_sm];
elseif comp < 0
    deaths_sm(1:(-comp)) = [];
end

N = sim_morts; 
k = deaths_sm;

if ~isempty(find(N<k))
    out = -Inf;
else
    vec = k*log(p_death_report) + (N-k)*log(1-p_death_report) + gammaln(N+1) - gammaln(k+1) - gammaln(N-k+1);    
    out = sum(vec);
end

aux.sim_cases = sim_cases;
aux.sim_morts = sim_morts;


