load data\EX_succ_e.mat;
load data\candidates_iML1515-EX_succ_e.mat;
load data\reduced\reduced_iML1515-EX_succ_e.mat 

% Parameter settings 
M = 300;
F = 1;
maxknock = 3; % Maximum total operating costs
% maxreg_all = 1; % The total cost ceiling for upward and downward adjustments
r = 500;
min_growth_scale = 0.05; 

solWT=optimizeCbModel(model);
w_lb = model.lb;
% Set the lower growth limit for the mutant to a percentage of the wild-type’s optimal value
w_lb(find(model.c==1)) = (floor(abs(solWT.f*1e6))*1e-6);    
m_lb = model.lb;
m_lb(find(model.c==1)) = (floor(abs(solWT.f*1e6))*1e-6)*min_growth_scale;

nrxn = size(model.S,2);
nvar = nrxn*2;
nmetab = size(model.S,1);

model.g = zeros(nrxn,1);
target_str = 'EX_succ_e';  
target_idx = find(strcmp(reducedModel.rxns, target_str));
if(length(target_idx)~=1)
    return
end
model.g(target_idx) = 1;

e_plus = [];
e_minus = [];
s_plus = [];
s_minus = [];

for i = 1:length(model.rxns)
    if contains(model.rxns{i}, 'plus_dem') 
        e_plus(end + 1) = i;
    elseif contains(model.rxns{i}, 'minus_dem')
        e_minus(end + 1) = i;
    end
end

s_plus = e_plus + 1; 
s_minus = e_minus + 1; 
num_es = length(s_plus);

% -----------------------------Constructing the c-mapping matrix: start------------------------------------------
c_plus = sparse(nrxn, nrxn);
c_minus = sparse(nrxn, nrxn);
for i = 1:length(e_plus)
    idx = e_plus(i);
    c_plus(idx, idx) = 1;
end

for i = 1:length(e_minus)
    idx = e_minus(i);
    c_minus(idx, idx) = 1; 
end

% -----------------------------Constructing the k-mapping matrix: start------------------------------------------
k_plus = sparse(nrxn, nrxn);
k_minus = sparse(nrxn, nrxn); 
for i = 1:length(s_plus)
    idx = s_plus(i);
    k_plus(idx, idx) = 1;
end

for i = 1:length(s_minus)
    idx = s_minus(i);
    k_minus(idx, idx) = 1;
end

% -----------------------------build mapping matrix-----------------------------------------
% pair_x: each row stores the (forward, backward) indices of regulated reversible reactions
% enforce y_f - y_b = 0 for regulated reactions
pair_x = model.pair_x;
n_pair_x = size(pair_x,1);
new_pair_x = sparse(nrxn,nrxn);
for i = 1:n_pair_x
    new_pair_x(pair_x(i,1),pair_x(i,1)) = 1;
    new_pair_x(pair_x(i,1),pair_x(i,2)) = -1;
end
new_pair_x = [new_pair_x,sparse(nrxn,nrxn);sparse(nrxn,nrxn),new_pair_x];

% new_pair_s enforces that the total selection over f/b and their pseudo-reactions is ≤ 1
new_pair_s = sparse(n_pair_x,nrxn);
for i = 1:size(pair_x,1)
    f_rxn = model.rxns(pair_x(i, 1)); 
    b_rxn = model.rxns(pair_x(i, 2));
    s_f_plusId = find(strcmp(model.rxns, sprintf('%s_plus_slack', f_rxn{1})));
    s_f_minusId = find(strcmp(model.rxns, sprintf('%s_minus_slack', f_rxn{1})));
    s_b_plusId = find(strcmp(model.rxns, sprintf('%s_plus_slack', b_rxn{1})));
    s_b_minusId = find(strcmp(model.rxns, sprintf('%s_minus_slack', b_rxn{1})));
    new_pair_s(i,s_f_plusId) = 1;
    new_pair_s(i,s_f_minusId) = 1;
    new_pair_s(i,s_b_plusId) = 1;
    new_pair_s(i,s_b_minusId) = 1;
    new_pair_s(i,pair_x(i,1)) = 1; 
end

new_pair_s = [new_pair_s,sparse(n_pair_x,nrxn);sparse(n_pair_x,nrxn),new_pair_s];

% For each regulated reaction, the up-regulation, down-regulation, and knockout options are mutually exclusive
% This constraint applies to all regulated reactions
new_pair_reg = sparse(num_es,nrxn);
for i = 1:num_es
    new_pair_reg(i,s_plus(i)) = 1;
    new_pair_reg(i,s_minus(i)) = 1;
    new_pair_reg(i,model.reg_ind(i)) = 1;
end
new_pair_reg = [new_pair_reg,sparse(num_es,nrxn);sparse(num_es,nrxn),new_pair_reg];

c_p = [zeros(size(model.g,1),1);model.g];
% The wild-type model is extended to match the dimensionality of the mutant model
% This extension is omitted in the paper for clarity, as it is mathematically equivalent and does not affect the results
hat_ub = [model.ub;model.ub]; 
hat_lb = [w_lb;m_lb];

hat_k1 = [k_plus, sparse(nrxn,nrxn);sparse(nrxn,nrxn),k_plus];
hat_k2 = [k_minus, sparse(nrxn,nrxn);sparse(nrxn,nrxn),k_minus];
hat_S = [model.S, sparse(nmetab,nrxn); sparse(nmetab,nrxn), model.S];
k_ko = sparse(nrxn,nrxn);
koIdx = model.reg_ind; 
for i = 1 : length(koIdx)
    k_ko(koIdx(i),koIdx(i)) = 1;
end
hat_k3 = [k_ko,sparse(nrxn,nrxn);sparse(nrxn,nrxn),k_ko];
diag_ub = [diag(model.ub),sparse(nrxn,nrxn);sparse(nrxn,nrxn),diag(model.ub)];

R1 = sparse(nrxn,nrxn);
R2 = sparse(nrxn,nrxn);
R3 = sparse(nrxn,nrxn);

for i = 1:length(s_plus)
    R1(s_plus(i),s_plus(i)) = r;
    R2(s_minus(i),s_minus(i)) = r;
end

for i = 1:length(koIdx)
    R3(koIdx(i),koIdx(i)) = r;
end

R11 = [R1,sparse(nrxn,nrxn);sparse(nrxn,nrxn),R1];
R22 = [R2,sparse(nrxn,nrxn);sparse(nrxn,nrxn),R2];
R33 = [R3,sparse(nrxn,nrxn);sparse(nrxn,nrxn),R3];

knock = ones(1,nrxn);
b_indices = find(endsWith(model.rxns, '_b'));
if ~isempty(b_indices)
    b_begin = b_indices(1); 
    b_end = b_indices(end);  
else
    b_begin = []; 
    b_end = [];
end
knock(b_begin:b_end) = 0;

y_plus_ub = zeros(nrxn,1);
y_minus_ub = zeros(nrxn,1);
y_x_ub = zeros(nrxn,1);
y_plus_ub(s_plus) = 1;
y_minus_ub(s_minus) = 1;
y_x_ub(koIdx) = 1;

y_plus_ub = [y_plus_ub;y_plus_ub];
y_minus_ub = [y_minus_ub;y_minus_ub];
y_x_ub = [zeros(nrxn,1);y_x_ub];


timestamp = datestr(now, 'yyyymmdd_HHMMSS');
filename = sprintf('experiment_%s_%s.xlsx', target_str, timestamp);
while(maxknock<12)
        hat_c1 = [-F * c_plus,c_plus];  
        hat_c2 = [c_minus, -F * c_minus]; 
        % x = [ z epsilon xi phi alpha beta pi rho sigma y_plus y_minus y_x];
        c = [ 
              zeros(nvar,1);
              zeros(2*nmetab,1);
              zeros(nrxn,1);
              zeros(nrxn,1);
              -hat_ub;
              -hat_ub;
              -hat_ub;
              hat_lb;
              -hat_ub;
              zeros(nvar,1);
              zeros(nvar,1);
              zeros(nvar,1);
              ];
    
            a = [
                 hat_S sparse(2*nmetab,2*nmetab) sparse(2*nmetab,nrxn) sparse(2*nmetab,nrxn) sparse(2*nmetab,nvar) sparse(2*nmetab,nvar) sparse(2*nmetab,nvar)  sparse(2*nmetab,nvar) sparse(2*nmetab,nvar) sparse(2*nmetab,nvar) sparse(2*nmetab,nvar)  sparse(2*nmetab,nvar);
                 hat_c1 sparse(nrxn,2*nmetab) sparse(nrxn,nrxn) sparse(nrxn,nrxn) sparse(nrxn,nvar) sparse(nrxn,nvar) sparse(nrxn,nvar) sparse(nrxn,nvar) sparse(nrxn,nvar) sparse(nrxn,nvar) sparse(nrxn,nvar) sparse(nrxn,nvar); 
                 hat_c2 sparse(nrxn,2*nmetab) sparse(nrxn,nrxn) sparse(nrxn,nrxn) sparse(nrxn,nvar) sparse(nrxn,nvar) sparse(nrxn,nvar) sparse(nrxn,nvar) sparse(nrxn,nvar) sparse(nrxn,nvar) sparse(nrxn,nvar) sparse(nrxn,nvar);
     
                 hat_k1 sparse(nvar,2*nmetab) sparse(nvar,nrxn) sparse(nvar,nrxn) sparse(nvar,nvar) sparse(nvar,nvar) sparse(nvar,nvar) sparse(nvar,nvar) sparse(nvar,nvar) diag_ub sparse(nvar,nvar) sparse(nvar,nvar);
                 hat_k1 sparse(nvar,2*nmetab) sparse(nvar,nrxn) sparse(nvar,nrxn) sparse(nvar,nvar) sparse(nvar,nvar) sparse(nvar,nvar) sparse(nvar,nvar) sparse(nvar,nvar) sparse(nvar,nvar) sparse(nvar,nvar) sparse(nvar,nvar);
                 hat_k2 sparse(nvar,2*nmetab) sparse(nvar,nrxn) sparse(nvar,nrxn) sparse(nvar,nvar) sparse(nvar,nvar) sparse(nvar,nvar) sparse(nvar,nvar) sparse(nvar,nvar) sparse(nvar,nvar) diag_ub sparse(nvar,nvar);
                 hat_k2 sparse(nvar,2*nmetab) sparse(nvar,nrxn) sparse(nvar,nrxn) sparse(nvar,nvar) sparse(nvar,nvar) sparse(nvar,nvar) sparse(nvar,nvar) sparse(nvar,nvar) sparse(nvar,nvar) sparse(nvar,nvar) sparse(nvar,nvar);
     
                 hat_k3 sparse(nvar,2*nmetab) sparse(nvar,nrxn) sparse(nvar,nrxn) sparse(nvar,nvar) sparse(nvar,nvar) sparse(nvar,nvar) sparse(nvar,nvar) sparse(nvar,nvar) sparse(nvar,nvar) sparse(nvar,nvar) sparse(nvar,nvar);
                 hat_k3 sparse(nvar,2*nmetab) sparse(nvar,nrxn) sparse(nvar,nrxn) sparse(nvar,nvar) sparse(nvar,nvar) sparse(nvar,nvar) sparse(nvar,nvar) sparse(nvar,nvar) sparse(nvar,nvar) sparse(nvar,nvar) diag_ub;
                 speye(nvar) sparse(nvar,2*nmetab) sparse(nvar,nrxn) sparse(nvar,nrxn) sparse(nvar,nvar) sparse(nvar,nvar) sparse(nvar,nvar) sparse(nvar,nvar) sparse(nvar,nvar) sparse(nvar,nvar) sparse(nvar,nvar) sparse(nvar,nvar);
                 -speye(nvar) sparse(nvar,2*nmetab) sparse(nvar,nrxn) sparse(nvar,nrxn) sparse(nvar,nvar) sparse(nvar,nvar) sparse(nvar,nvar) sparse(nvar,nvar) sparse(nvar,nvar) sparse(nvar,nvar) sparse(nvar,nvar) sparse(nvar,nvar);   
     
                 sparse(nvar,nvar) hat_S' hat_c1' hat_c2' -hat_k1' -hat_k2' -hat_k3' speye(nvar) -speye(nvar) -R11 -R22 -R33;
                 % The upper and lower limits (excluding knockout) are equal in u and v
                 sparse(nrxn,nvar) sparse(nrxn,2*nmetab) sparse(nrxn,nrxn) sparse(nrxn,nrxn) sparse(nrxn,nvar) sparse(nrxn,nvar) sparse(nrxn,nvar) sparse(nrxn,nvar) sparse(nrxn,nvar) [speye(nrxn),-speye(nrxn)] sparse(nrxn,nvar) sparse(nrxn,nvar);
                 sparse(nrxn,nvar) sparse(nrxn,2*nmetab) sparse(nrxn,nrxn) sparse(nrxn,nrxn) sparse(nrxn,nvar) sparse(nrxn,nvar) sparse(nrxn,nvar) sparse(nrxn,nvar) sparse(nrxn,nvar) sparse(nrxn,nvar) [speye(nrxn),-speye(nrxn)] sparse(nrxn,nvar); 
                 
                 sparse(nvar,nvar) sparse(nvar,2*nmetab) sparse(nvar,nrxn) sparse(nvar,nrxn) sparse(nvar,nvar) sparse(nvar,nvar) sparse(nvar,nvar) sparse(nvar,nvar) sparse(nvar,nvar) sparse(nvar,nvar) sparse(nvar,nvar) new_pair_x;
     
                 sparse(n_pair_x*2,nvar) sparse(n_pair_x*2,2*nmetab) sparse(n_pair_x*2,nrxn) sparse(n_pair_x*2,nrxn) sparse(n_pair_x*2,nvar) sparse(n_pair_x*2,nvar) sparse(n_pair_x*2,nvar) sparse(n_pair_x*2,nvar) sparse(n_pair_x*2,nvar) new_pair_s new_pair_s new_pair_s;
                 sparse(num_es*2,nvar) sparse(num_es*2,2*nmetab) sparse(num_es*2,nrxn) sparse(num_es*2,nrxn) sparse(num_es*2,nvar) sparse(num_es*2,nvar) sparse(num_es*2,nvar) sparse(num_es*2,nvar) sparse(num_es*2,nvar) new_pair_reg new_pair_reg new_pair_reg;
                 sparse(1,nvar) sparse(1,2*nmetab) sparse(1,nrxn) sparse(1,nrxn) sparse(1,nvar) sparse(1,nvar) sparse(1,nvar) sparse(1,nvar) sparse(1,nvar) [zeros(1,nrxn),knock] [zeros(1,nrxn),knock] [zeros(1,nrxn),knock]; 
                 % sparse(1,nvar) sparse(1,2*nmetab) sparse(1,nrxn) sparse(1,nrxn) sparse(1,nvar) sparse(1,nvar) sparse(1,nvar) sparse(1,nvar) sparse(1,nvar) [zeros(1,nrxn),knock] [zeros(1,nrxn),knock] sparse(1,nvar); 

                ];
    
        b = [ 
              zeros(2*nmetab, 1);
              zeros(nrxn,1);
              zeros(nrxn,1);
     
              hat_ub;
              zeros(nvar,1);
              hat_ub;
              zeros(nvar,1);
     
              zeros(nvar,1);
              hat_ub;
              hat_lb;
              -hat_ub;
     
              c_p;
              zeros(nrxn,1);
              zeros(nrxn,1);
    
              zeros(nvar,1);
    
              ones(n_pair_x*2,1);
              ones(num_es*2,1);
              maxknock;
              % maxreg_all;
            ];
    
        ctype = char([ ...
                      'E' * ones(2*nmetab,1);
                      'G' * ones(nrxn,1);
                      'G' * ones(nrxn,1);
     
                      'L' * ones(nvar,1);
                      'G' * ones(nvar,1);
                      'L' * ones(nvar,1);
                      'G' * ones(nvar,1);
                                 
                      'G' * ones(nvar,1);
                      'L' * ones(nvar,1);
                      'G' * ones(nvar,1);
                      'G' * ones(nvar,1);
     
                      'L' * ones(nvar,1);
                      'E' * ones(nrxn,1);
                      'E' * ones(nrxn,1);
    
                      'E' * ones(nvar,1);
     
                      'L' * ones(n_pair_x*2,1);
                      'L' * ones(num_es*2,1);
                      'L' * ones(1,1);
                      % 'L' * ones(1,1);
            ]);
    
    
        lb = [hat_lb;
              -Inf*ones(2*nmetab,1);
              zeros(nrxn,1);
              zeros(nrxn,1);
              zeros(nvar,1);
              zeros(nvar,1);
              zeros(nvar,1);
              zeros(nvar,1);
              zeros(nvar,1);
              zeros(nvar,1);
              zeros(nvar,1);
              zeros(nvar,1);
              ];
        ub = [hat_ub;
              Inf*ones(2*nmetab,1);
              Inf*ones(nrxn,1);
              Inf*ones(nrxn,1);
              Inf*ones(nvar,1);
              Inf*ones(nvar,1);
              Inf*ones(nvar,1);
              Inf*ones(nvar,1);
              Inf*ones(nvar,1);
              y_plus_ub;
              y_minus_ub;
              y_x_ub;
              ];
    
        vartype = char([ ...
              'C'*ones(nvar,1);
              'C'*ones(2*nmetab,1);
              'C' * ones(nrxn,1);
              'C' * ones(nrxn,1);
              'C'*ones(nvar,1);
              'C'*ones(nvar,1);
              'C'*ones(nvar,1);
              'C'*ones(nvar,1);                        
              'C'*ones(nvar,1);
              'B'*ones(nvar,1);
              'B'*ones(nvar,1);
              'B'*ones(nvar,1);
            ]);
    
        bilevelMILPproblem = struct();
        bilevelMILPproblem.A = a;
        bilevelMILPproblem.b = b;
        bilevelMILPproblem.c = c;
        bilevelMILPproblem.csense = ctype;
        bilevelMILPproblem.lb = lb;
        bilevelMILPproblem.ub = ub;
        bilevelMILPproblem.vartype = vartype;
        bilevelMILPproblem.osense = -1;

        solverParams.MIPFocus = 2;
        solverParams.Cuts = 1;
        solverParams.Method = 2;
        solverParams.FeasibilityTol=1e-9; 
        solverParams.TimeLimit = 3600*3;
        solverParams.OutputFlag = 1;
        solverParams.Heuristics = 1;
        solverParams.PoolSolutions = 30;
        solverParams.PoolSearchMode = 1;

        % solverParams.LogFile = sprintf('gurobi_log_%s_%d%d_%s.txt',target_str,maxknock,maxreg_all,timestamp);
        optKnockSol = solveCobraMILP(bilevelMILPproblem, solverParams);
        
        % Determine whether there are multiple solutions
        if isfield(optKnockSol, 'pool') && ~isempty(optKnockSol.pool)
            num_solutions = numel(optKnockSol.pool);
        else
            num_solutions = 1;
        end
        
        result_rows = {}; 

        for sol_idx = 1:num_solutions
            if num_solutions == 1
                x_full = optKnockSol.full;
                fmaxt = optKnockSol.obj;
            else
                x_full = optKnockSol.pool(sol_idx).xn;
                fmaxt = optKnockSol.pool(sol_idx).objval;
            end
        
            % === Decode variable ===
            z = x_full(1:nvar);
        
            y_plus = x_full(7*nvar+2*nmetab+nrxn+1:8*nvar+2*nmetab);
            y_minus = x_full(8*nvar+2*nmetab+nrxn+1:9*nvar+2*nmetab);
            y_x = x_full(9*nvar+2*nmetab+nrxn+1:end);
        
            yt = y_plus + y_minus + y_x;
            indices = find(yt == 1);
            delete_rxns = model.rxns(indices);
        
            % === Print to the console (tabular format) ===
            fprintf('\n============= Solution #%d =============\n', sol_idx);
            indices = indices(:);
            delete_rxns = delete_rxns(:);
            T = table(indices, delete_rxns, 'VariableNames', {'Index', 'Reaction'});
            disp(T);
            fprintf('Objective function value: %.4f\n', fmaxt);
            fprintf('Target product rate: %.4f\n', z(target_idx + nrxn));
            fprintf('Solution time: %.2fs\n', optKnockSol.time);
            disp('-----------------------------------');
        
            headers = {'F', 'maxknock', 'M', 'r', 'maxreg_all','times', ...
                       'Objective function value', 'Target product rate', 'Delete reaction number', ...
                       'Delete reaction name', 'Minimum growth rate ratio', 'Decode'};

            del_ids_str = strjoin(string(indices), ', ');
            del_names_str = strjoin(delete_rxns, ', ');
        
            data_row = {F, maxknock, M, r, maxreg_all, optKnockSol.time, ...
                        fmaxt, z(target_idx+nrxn), del_ids_str, ...
                        del_names_str, min_growth_scale, sol_idx};

            result_rows = [result_rows; data_row];  
        end
        % ======= Write everything at once once the loop has finished =======
        if ~isfile(filename)
            writecell(headers, filename, 'Sheet', 1, 'Range', 'A1');
            next_row = 2;
        else
            raw = readcell(filename);
            
            next_row = size(raw,1) + 1;
        end

        writecell(result_rows, filename, 'Sheet', 1, 'Range', ['A', num2str(next_row)]);
        maxknock = maxknock+1;
end