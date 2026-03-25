load data\candidates_iML1515-EX_succ_e.mat;  
load data\reduced\reduced_iML1515-EX_succ_e.mat;

regulation.reg_ind = find(ismember(reducedModel.rxns, candidates.regList));
regulation.reg_ind = sort(regulation.reg_ind); 

% Record whether each regulated reaction is reversible (1: reversible, 0: irreversible)
raw_direction_reg_ind = [];
for i =  1:length(regulation.reg_ind)
    if(reducedModel.lb(regulation.reg_ind(i))<0 && reducedModel.ub(regulation.reg_ind(i))>0)
        direction = 1;
        raw_direction_reg_ind = [raw_direction_reg_ind,direction];
    else
        direction = 0;
        raw_direction_reg_ind = [raw_direction_reg_ind,direction];
    end
end

% Store the original number of regulated reactions for later forward/backward pairing
raw_nreg_ind = length(regulation.reg_ind);

% Add reversible reaction names to specificReactions for subsequent irreversible conversion
specificReactions = {};
num = 0;
for i = 1:length(regulation.reg_ind)
    index = regulation.reg_ind(i);
    if reducedModel.lb(index) < 0
        reaction = reducedModel.rxns{index};
        specificReactions{end + 1} = reaction;
        if reducedModel.ub(index) > 0
            num = num+1;
        end
    end
end

% Convert reactions to irreversible form
newModel = convertToIrreversible(reducedModel, 'sRxns', specificReactions,'OrderReactions', false,'flipOrientation', true);
% Append the newly created backward reactions of reversible regulated reactions
for i = length(newModel.rxns) - num + 1:length(newModel.rxns)
    regulation.reg_ind(end + 1) = i;
end

regulation.reg_down_up=[true(1,numel(regulation.reg_ind));true(1,numel(regulation.reg_ind))]; % 确定是上调还是下调
regulation.reg_bounds=[]; 

% Convert the remaining reactions to irreversible form
% Final reaction order:
% original reactions + backward part of regulated reactions + backward part of residual reactions + pseudo-reactions
residualReactions = {};
pair_residual_f = [];
for i = 1:length(newModel.rxns)
    if newModel.lb(i) < 0
        reaction = newModel.rxns{i};
        residualReactions{end + 1} = reaction;
        if newModel.ub(i) > 0
            pair_residual_f = [pair_residual_f;i];
        end
    end
end
newModel = convertToIrreversible(newModel, 'sRxns', residualReactions,'OrderReactions', false,'flipOrientation', true);

residual_num = length(pair_residual_f); 
pair_residual_b = [];
for i = length(newModel.rxns) - residual_num + 1:length(newModel.rxns)
    pair_residual_b = [pair_residual_b;i];
end

pair_residual = [pair_residual_f, pair_residual_b];
regModel=add_pseudo_reactions(newModel,regulation.reg_ind);

% Build the forward/backward mapping for regulated reversible reactions
pair_x = [];
k = 1;
for i = 1:raw_nreg_ind
    if raw_direction_reg_ind(i) == 1
        newRow_x = [regulation.reg_ind(i),regulation.reg_ind(raw_nreg_ind+k);];
        pair_x = [pair_x; newRow_x];
        k = k+1;
    end
end

regModel.pair_x = pair_x;
regModel.pair_residual = pair_residual;
regModel.reg_ind = regulation.reg_ind;

model = struct();
model.S = regModel.S;
model.rxns = regModel.rxns;
model.lb = regModel.lb;
model.ub = regModel.ub;
model.c = regModel.c;
model.b = regModel.b;
model.present = regModel.present;
model.reg_ind = regModel.reg_ind;
model.pair_residual = regModel.pair_residual;
model.pair_x = regModel.pair_x;

save(fullfile('data','EX_succ_e.mat'), 'model');

return ;