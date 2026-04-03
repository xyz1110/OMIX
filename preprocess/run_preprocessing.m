warning('off','all')

% Preprossing
clear;
close all;

targetRxn='EX_succ_e'; 

loadedData  = load('..\data\iML1515.mat');
model = loadedData.model;

% get biomass reaction
if isfield(model,'csense')
    if ~iscolumn(model.csense)
        model.csense=columnVector(model.csense);
    end
else
    model.csense=char('E'*ones(length(model.b),1));
end
biomassRxn=model.rxns{model.c==1};

% set uptake rate of oxygen and glucose
oxygenRxn='EX_o2_e';
substrate='EX_glc__D_e';

% change carbon source
model=changeRxnBounds(model,substrate, -10, 'l'); 

% limit reaction rate in realistic range
model.lb(model.lb<-100)=-100;
model.ub(model.ub>100)=100;
orimodel=model;

% when the model size is big, it is better to compress the model so that
% the linear reactions can be reduced
if size(model.S,1)>100
    % compress the model and get compressed candidate reactions
    [model,candidates]=preprocessing(orimodel,substrate,oxygenRxn,biomassRxn,targetRxn);
    regList = candidates;
    koList = candidates;
else
    transRxns=findTransRxns(model);
    regList=setdiff(model.rxns(~contains(model.rxns,'EX_')), transRxns); % exclude transport reactions
    koList=regList;
end

% make sure target reaction is in the reaction list of compressed model
if ~strcmp(model.rxns, targetRxn)
    newTargetRxn=findLumpRxns(model, targetRxn);
else
    newTargetRxn=targetRxn;
end

% remove unreasonable reactions from candidate knockout set
regList=setdiff(regList, {'ATPM', biomassRxn, newTargetRxn});
koList=setdiff(koList, {'ATPM', biomassRxn, newTargetRxn});

candidates = struct();
candidates.regList = regList;
candidates.koList = koList;

disp(['Preprocessing step completed, regulation candidates: ', num2str(length(regList)), ...
    '; knockout candidates:', num2str(length(koList))]);

save(fullfile('..','data',['candidates_' model.description '.mat']), 'candidates');
