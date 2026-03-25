function reg = add_pseudo_reactions(sys, sub_reg_ind)

reg = sys;

% Direction flags for regulatory extensions
add_plus  = true;   % Up-regulation (>=)
add_minus = true;   % Down-regulation (<=)

% Each direction introduces two pseudo-reactions (demand and slack)
num_pm = (add_plus + add_minus) * length(sub_reg_ind);
num_pr = (add_plus + add_minus) * length(sub_reg_ind) * 2;

% Extend the model structure
reg.S(end+num_pm, end+num_pr) = 0;
reg.rxns(end+num_pr) = {''};
reg.lb(end+num_pr) = 0;
reg.ub(end+num_pr) = 0;
reg.c = zeros(size(reg.S,2),1);
reg.c(find(sys.c==1))=1;
reg.present = ones(size(reg.S,2),1);
reg.b = zeros(size(reg.S,1),1)

% Record slack reactions
reg.slackRxnInd = [];

% Starting indices of newly added pseudo-metabolites and pseudo-reactions
pm_ind = size(sys.S, 1) + 1;  
pr_ind = size(sys.S, 2) + 1;  

for k = 1:length(sub_reg_ind)

    sub_ind = sub_reg_ind(k);

    %% ---------------------------------------------------------
    %  (1) Up-regulation branch
    %% ---------------------------------------------------------
    if add_plus

        % Add pseudo-metabolite for the plus branch
        reg.S(pm_ind, :) = 0;  
        reg.S(pm_ind, sub_ind) = 1;   

        % Demand reaction in the plus branch
        plus_dem = pr_ind;
        reg.S(pm_ind, plus_dem) = -1;  
        reg.lb(plus_dem) = reg.lb(sub_ind);
        reg.ub(plus_dem) = reg.ub(sub_ind);
        reg.rxns{plus_dem} = sprintf('%s_plus_dem', reg.rxns{sub_ind});

        pr_ind = pr_ind + 1;

        % Slack reaction in the plus branch
        plus_slack = pr_ind;
        reg.S(pm_ind, plus_slack) = +1; 
        reg.lb(plus_slack) = reg.lb(sub_ind);
        reg.ub(plus_slack) = reg.ub(sub_ind);
        reg.rxns{plus_slack} = sprintf('%s_plus_slack', reg.rxns{sub_ind});

        reg.slackRxnInd(end+1) = plus_slack; 

        pr_ind = pr_ind + 1;
        pm_ind = pm_ind + 1; 
    end



    %% ---------------------------------------------------------
    %  (2) % Down-regulation branch
    %% ---------------------------------------------------------
    if add_minus

        % Add pseudo-metabolite for the minus branch
        reg.S(pm_ind, :) = 0;
        reg.S(pm_ind, sub_ind) = 1;  % 原反应 → pm_minus 产生

        % Demand reaction in the minus branch
        minus_dem = pr_ind;
        reg.S(pm_ind, minus_dem) = -1;  % 消耗 pm
        reg.lb(minus_dem) = reg.lb(sub_ind);
        reg.ub(minus_dem) = reg.ub(sub_ind);
        reg.rxns{minus_dem} = sprintf('%s_minus_dem', reg.rxns{sub_ind});

        pr_ind = pr_ind + 1;

        % Slack reaction in the minus branch
        minus_slack = pr_ind;
        reg.S(pm_ind, minus_slack) = -1;  % ★ 原作者结构：下调松弛是 -1
        reg.lb(minus_slack) = reg.lb(sub_ind);
        reg.ub(minus_slack) = reg.ub(sub_ind);
        reg.rxns{minus_slack} = sprintf('%s_minus_slack', reg.rxns{sub_ind});

        reg.slackRxnInd(end+1) = minus_slack;

        pr_ind = pr_ind + 1;
        pm_ind = pm_ind + 1;
    end

end
