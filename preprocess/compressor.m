function [modelReduced,biomassRxn,targetRxn, options]= compressor(model,options)

global CBT_MILP_SOLVER
solver = CBT_MILP_SOLVER;
bfRelation=true;

tic

%% Check if the model was already reduced or not

if isfield(model, 'reductionStatus')
    modelReduced = model;
    warning('The model was already reduced')
    return
end

%% Prepare variables

if ~exist('options','var') || isempty(options)
    options.verbFlag = false;
    options.loadFVAFlux = false;
else
    if ~isfield(options,'verbFlag'), options.verbFlag = false; end
    if ~isfield(options,'loadFVAFlux'), options.loadFVAFlux = true; end
    if ~isfield(options,'minBioPecentage'), options.minBioPecentage = 5; end
end

[~,n] = size(model.S);
biomassRxnID = findRxnIDs(model,options.biomassRxn);
if biomassRxnID == 0, error([num2str(options.biomassRxn) ' not found.']);end
targetRxnID = findRxnIDs(model,options.targetRxn);
if targetRxnID == 0, error([num2str(options.targetRxn) ' not found.']);end
oxygenRxnID = findRxnIDs(model,options.oxygenRxn);
if oxygenRxnID == 0, error([num2str(options.oxygenRxn) ' not found.']);end

reductionStatus = 1;


%% Calculate theoretical maximum production rate of target metabolite
modelMaxTarget = changeObjective(model,options.targetRxn,1);

solMaxTarget = optimizeCbModel(modelMaxTarget,'max');
solMaxTarget.x = round(solMaxTarget.x*10^6)/10^6;
theorMaxProdRate = solMaxTarget.x(targetRxnID);
if theorMaxProdRate == 0
    warning('Target metabolite cannot be produced in the original model.')
end

%% Flux Variability analysis based on biomass

if options.verbFlag == true
    disp('Model reduction start.')
end

fvaFile = fullfile('..', 'data', 'fluxFVA');
if ~exist(fvaFile, 'dir')
    mkdir(fvaFile);
end

fvaFile = fullfile(fvaFile, ['fluxFVA_' model.description '.mat']);

if options.loadFVAFlux == true && exist(fvaFile,'file')
    load(fvaFile);
else
    if options.verbFlag == true
        disp('Flux variability anaysis')
    end
    
    [minflux,maxflux]= fluxVariability(model,options.minBioPecentage); % minimal 5% of the wild type.
    fluxFVA = [minflux,maxflux];
    save(fvaFile,'fluxFVA')
end
WTsol=optimizeCbModel(model);
model.mins=fluxFVA(:,1); model.maxs=fluxFVA(:,2);
model.lb(model.mins<0 & model.mins>1e-3)=-1e-1;
model.ub(model.maxs>0 & model.maxs<1e-3)=1e-1;
WTfvasol=optimizeCbModel(model);

if (abs(WTsol.f-WTfvasol.f)/WTsol.f>5e-2)
    error('the model cannot be compressed because theoretical maximum growth is not guaranteed. Please increase flux bounds')
end

%% Remove reactions which cannot carry any fluxes in the given condition

if options.verbFlag == true
    disp('Remove blocked reactions from model')
end
removeRxnIDs = sum(abs([fluxFVA(:,1),fluxFVA(:,2)]),2)<=0;
modelRem = removeRxns(model,model.rxns(removeRxnIDs));

nonRemovedGene = sum(modelRem.rxnGeneMat,1)>0;
modelRem.rxnGeneMat = modelRem.rxnGeneMat(:,nonRemovedGene);
modelRem.genes = modelRem.genes(nonRemovedGene);

if ~isfield(modelRem, 'rules') || any(cellfun(@isempty,modelRem.rules))
    modelRem = generateRules(modelRem);
    modelRem.rules=columnVector(modelRem.rules);
end
modelRem.rules=updateGeneRule(modelRem.rules, length(modelRem.genes)); % writen by myself.
if length(modelRem.csense)>length(modelRem.b), modelRem.csense(length(modelRem.b)+1:end)=[];end

clear nonRemovedGene

if options.verbFlag == true
    disp([num2str(length(modelRem.rxns)) ' reactions' ...
        ', ' num2str(length(modelRem.mets)) ' metabolites, ' num2str(length(modelRem.genes)) ' genes'])
end


%% Combine adjacent reactions without branching

if options.verbFlag == true
    disp('Combine adjacent reactions without branching')
end

% Identify exchange reactions and its flux bounds
[isExRxns, ~] = findExcRxns(modelRem);

% Combine linear reactions (except reactions with metabolites in biomassRxn)
combineRest = 1;
rxnSets = modelRem.rxns;
preservedRxnIDs = findRxnIDs(modelRem,{options.targetRxn,options.biomassRxn,options.substrateRxn});
if ismember(0, preservedRxnIDs)
    disp('Warning: some preserved reactions do not exist in the reduced model.')
end
preservedRxnIDs=setdiff(preservedRxnIDs,0); % remove unfound rxns.
preservedIDs = find(all(modelRem.S(:,preservedRxnIDs)==0,2));
while combineRest == 1
    combineRest = 0;
    
    for tempMetID = preservedIDs'
        metAssocRxnIDs = find(modelRem.S(tempMetID,:));
        
        if size(metAssocRxnIDs ,2) == 2
            combineRest = 1;
            if isExRxns(metAssocRxnIDs(2)) >= 1
                metAssocRxnIDs = metAssocRxnIDs([2,1]);
            end
            
            %Convert stoichiometric matrix
            coeffRatio = modelRem.S(tempMetID,metAssocRxnIDs(1))/modelRem.S(tempMetID,metAssocRxnIDs(2));
            modelRem.S(:,metAssocRxnIDs(2)) = modelRem.S(:,metAssocRxnIDs(2)) * (-coeffRatio);
            modelRem.S(:,metAssocRxnIDs(1)) = ...
                modelRem.S(:,metAssocRxnIDs(1)) + modelRem.S(:,metAssocRxnIDs(2));
            modelRem.S(:,metAssocRxnIDs(2)) = 0 ;
            
            %Convert objective function
            modelRem.c(metAssocRxnIDs(2)) = modelRem.c(metAssocRxnIDs(2)) * (-coeffRatio);
            modelRem.c(metAssocRxnIDs(1)) = ...
                modelRem.c(metAssocRxnIDs(1)) + modelRem.c(metAssocRxnIDs(2));
            modelRem.c(metAssocRxnIDs(2)) = 0 ;
            
            % Update lower and upper bounds
            lb=modelRem.lb(metAssocRxnIDs);
            ub=modelRem.ub(metAssocRxnIDs);
            mins=modelRem.mins(metAssocRxnIDs);
            maxs=modelRem.maxs(metAssocRxnIDs);
            
            lb(2)=modelRem.lb(metAssocRxnIDs(2))/(-coeffRatio);
            ub(2)=modelRem.ub(metAssocRxnIDs(2))/(-coeffRatio);
            mins(2)=modelRem.mins(metAssocRxnIDs(2))/(-coeffRatio);
            maxs(2)=modelRem.maxs(metAssocRxnIDs(2))/(-coeffRatio);
            if coeffRatio > 0 % one reaction produces and the other consumes the same metabolite in opposite direction.
                [lb(2),ub(2)]=deal(ub(2),lb(2));
                [mins(2),maxs(2)]=deal(maxs(2),mins(2));
            end
            modelRem.lb(metAssocRxnIDs(1))=max(lb);
            modelRem.ub(metAssocRxnIDs(1))=min(ub);
            
            modelRem.mins(metAssocRxnIDs(1))=max(mins);
            modelRem.maxs(metAssocRxnIDs(1))=min(maxs);
            
            %testSol = optimizeCbModel(modelRem,'max')
            
            %Convert reaction and metabolite names
            modelRem.rxns(metAssocRxnIDs(1)) = ...
                {[modelRem.rxns{metAssocRxnIDs(1)} , '/' , modelRem.rxns{metAssocRxnIDs(2)}]};
            if isfield(modelRem,'rxnNames')
                modelRem.rxnNames(metAssocRxnIDs(1)) = ...
                    {[modelRem.rxnNames{metAssocRxnIDs(1)} , '/' , modelRem.rxnNames{metAssocRxnIDs(2)}]};
                modelRem.rxnNames(metAssocRxnIDs(2)) = {'none'};
            end
            rxnSets(metAssocRxnIDs(1)) = ...
                {[rxnSets{metAssocRxnIDs(1)} , '/' , rxnSets{metAssocRxnIDs(2)}]};
            rxnSets(metAssocRxnIDs(2)) = {'none'};
            modelRem.mets(tempMetID,1) = {'none'};
            
            %Convert genes
            modelRem.rxnGeneMat(metAssocRxnIDs(1),:) = ...
                modelRem.rxnGeneMat(metAssocRxnIDs(1),:) + modelRem.rxnGeneMat(metAssocRxnIDs(2),:);
            modelRem.rxnGeneMat(metAssocRxnIDs(2),:) = 0;
            if ~isempty(modelRem.grRules{metAssocRxnIDs(1)})
                modelRem.grRules(metAssocRxnIDs(1),1) = {modelRem.grRules{metAssocRxnIDs(1)}};
                if ~isempty(modelRem.grRules{metAssocRxnIDs(2)})
                    modelRem.grRules(metAssocRxnIDs(1),1)={[modelRem.grRules{metAssocRxnIDs(1)},' and ', ...
                        modelRem.grRules{metAssocRxnIDs(2)}]};
                end
            else
                modelRem.grRules(metAssocRxnIDs(1),1) = {modelRem.grRules{metAssocRxnIDs(2)}};
            end
            modelRem.grRules(metAssocRxnIDs(2),1) = {''};
            
            % Adjust gene rules
            if ~isempty(modelRem.rules{metAssocRxnIDs(1)})
                modelRem.rules(metAssocRxnIDs(1),1) = {modelRem.rules{metAssocRxnIDs(1)}};
                if ~isempty(modelRem.rules{metAssocRxnIDs(2)})
                    modelRem.rules(metAssocRxnIDs(1),1)={[modelRem.rules{metAssocRxnIDs(1)},' & ', ...
                        modelRem.rules{metAssocRxnIDs(2)}]};
                end
            else
                modelRem.rules(metAssocRxnIDs(1),1) = {modelRem.rules{metAssocRxnIDs(2)}};
            end
            modelRem.rules(metAssocRxnIDs(2),1)={''};
            
            %Change Subsystems
            if isempty(modelRem.subSystems{metAssocRxnIDs(1)})
                modelRem.subSystems{metAssocRxnIDs(1)}=modelRem.subSystems{metAssocRxnIDs(2)};
            elseif ~isempty(modelRem.subSystems{metAssocRxnIDs(2)})
                modelRem.subSystems{metAssocRxnIDs(1)}=[modelRem.subSystems{metAssocRxnIDs(1)}, '/' ,modelRem.subSystems{metAssocRxnIDs(2)}];
            end
                
            %modelRem.subSystems{metAssocRxnIDs(1)}=[modelRem.subSystems{metAssocRxnIDs(1)}, '/' ,modelRem.subSystems{metAssocRxnIDs(2)}];
            %subSel2=ismember(modelRem.subSystems{metAssocRxnIDs(2)},options.selSubSystems);
            %if (~subSel2)
             %   modelRem.subSystems{metAssocRxnIDs(1)}=modelRem.subSystems{metAssocRxnIDs(2)};
            %end
            
        end
        clear metAssocRxnsNum coefRatio
    end
end
clear combine tempBiomassRxnNum nonBiomassMets coeffRatio

combineRxnIDs = find((sum(modelRem.S.^2) == 0));
modelMove = removeRxns(modelRem,modelRem.rxns(combineRxnIDs));

modelMove.csense(size(modelMove.S,1)+1:end)=[];

modelMove.rxnGeneMat(modelMove.rxnGeneMat>=1) = 1;
nonCombinedGene = find(sum(modelMove.rxnGeneMat,1)>0);
modelMove.rxnGeneMat = modelMove.rxnGeneMat(:,nonCombinedGene);
modelMove.genes = modelMove.genes(nonCombinedGene);
clear nonCombinedGene

if options.verbFlag == true
    disp([num2str(length(modelMove.rxns)) ' reactions' ...
        ', ' num2str(length(modelMove.mets)) ' metabolites, ' num2str(length(modelMove.genes)) ' genes'])
end

% testSol = optimizeCbModel(modelCom,'max')
biomassRxnID = findRxnIDs(modelMove,options.biomassRxn);
targetRxnID = findRxnIDs(modelMove,options.targetRxn);
% oxygenRxnID = findRxnIDs(modelMove,oxygenRxn);

%% Identify reaction correlations of models between before and after reduction
if (bfRelation)
    rxnAssociations = cell(length(modelMove.rxns),100);
    rxnAssocMat = zeros(length(modelMove.rxns),length(model.rxns));
    combineNum = zeros(length(modelMove.rxns),1);
    for i = 1 : length(modelMove.rxns)
        slashLocations = find(modelMove.rxns{i,1} == '/');
        combineNum(i) = length(slashLocations)+1;
        switch length(slashLocations)
            case 0
                rxnAssociations(i,1)=modelMove.rxns(i,1);
                rxnAssocMat(i,strcmp(modelMove.rxns(i),model.rxns))= 1;
            case 1
                rxnAssociations(i,1) = {modelMove.rxns{i}(1:slashLocations(1)-1)};
                rxnAssocMat(i,strcmp(modelMove.rxns{i}(1:slashLocations(1)-1),model.rxns))= 1;
                rxnAssociations(i,2) = {modelMove.rxns{i}(slashLocations(1)+1:end)};
                rxnAssocMat(i,strcmp(modelMove.rxns{i}(slashLocations(1)+1:end),model.rxns))= 1;
            otherwise
                rxnAssociations(i,1) = {modelMove.rxns{i}(1:slashLocations(1)-1)};
                rxnAssocMat(i,strcmp(modelMove.rxns{i}(1:slashLocations(1)-1),model.rxns))= 1;
                for j = 2 : length(slashLocations)
                    rxnAssociations(i,j) = {modelMove.rxns{i}(slashLocations(j-1)+1:slashLocations(j)-1)};
                    rxnAssocMat(i,strcmp(modelMove.rxns{i}(slashLocations(j-1)+1:slashLocations(j)-1),model.rxns))= 1;
                end
                rxnAssociations(i,j+1) = {modelMove.rxns{i}(slashLocations(j)+1:end)};
                rxnAssocMat(i,strcmp(modelMove.rxns{i}(slashLocations(j)+1:end),model.rxns))= 1;
        end
    end
    rxnAssociations = rxnAssociations(:,1:max(combineNum));
    rxnAssocMat = sparse(logical(rxnAssocMat));
    clear slashLocations
    
    
    %% Identify correlations between gene sets and reactions in reduced model
    %     [geneSets,IA,~]=unique(modelMove.grRules);
    %     geneSetsRules=modelMove.rules(IA);
    
    [unqGeneSetsMat] = unique(modelMove.rxnGeneMat,'rows');
    geneSets = cell(size(unqGeneSetsMat,1),1);
    for i = 1 : size(unqGeneSetsMat,1)
        geneNum = find(unqGeneSetsMat(i,:));
        switch length(geneNum)
            case 0
                geneSets{i} = 'None or Unknown';
            case 1
                geneSets{i} = modelMove.genes{geneNum};
            case 2
                geneSets{i} = [modelMove.genes{geneNum(1)} '_' modelMove.genes{geneNum(2)}];
            otherwise
                tempGeneSets = modelMove.genes{geneNum(1)};
                for j = 1 : length(geneNum)-1
                    tempGeneSets = [tempGeneSets '_' modelMove.genes{geneNum(j+1)}];
                end
                geneSets{i} = tempGeneSets;
        end
    end
    
    [~, rxnGeneSets] = ismember(modelMove.rxnGeneMat,unqGeneSetsMat,'rows');
    geneSetRxnMat = zeros(length(geneSets),length(modelMove.rxns));
    for i = 1 : length(modelMove.rxns)
        geneSetRxnMat(rxnGeneSets(i),i) = 1;
    end
    geneSetRxnMat = sparse(geneSetRxnMat);
    clear tempGeneSets
    
    
    %% Identify reactions associated with gene sets
    
    geneSetAssocRxns = cell(length(geneSets),1);
    for i = 1 : length(geneSets)
        associatedRxns = modelMove.rxns(logical(geneSetRxnMat(i,:)));
        tempAssocRxnsName = associatedRxns{1};
        if length(associatedRxns) >= 2
            for j = 2 : length(associatedRxns)
                tempAssocRxnsName(end+1:end+length(associatedRxns{j})+1) =[',' associatedRxns{j}];
            end
        end
        geneSetAssocRxns(i) = {tempAssocRxnsName};
    end
end

if (bfRelation)
    tempFluxSol1 = cell(1,length(geneSets));
    tempFluxSol2 = cell(1,length(geneSets));
    if isempty(gcp('nocreate'))
        for i = 1 : length(geneSets)
            modelDel = modelMove;
            modelDel.lb(logical(geneSetRxnMat(i,:))) = 0;
            modelDel.ub(logical(geneSetRxnMat(i,:))) = 0;
            delSol = optimizeCbModel(modelDel,'max');
            tempFluxSol1{i} = delSol.x;
            modelDel = changeObjective(modelDel,modelDel.rxns(1:2),[0,1]);
            delSol = optimizeCbModel(modelDel,'max');
            tempFluxSol2{i} = delSol.x;
        end
    else
        parfor i = 1 : length(geneSets)
            changeCobraSolver(solver, 'LP', 0);
            modelDel = modelMove;
            modelDel.lb(logical(geneSetRxnMat(i,:))) = 0;
            modelDel.ub(logical(geneSetRxnMat(i,:))) = 0;
            delSol = optimizeCbModel(modelDel,'max');
            tempFluxSol1{i} = delSol.x;
            modelDel = changeObjective(modelDel,modelDel.rxns(1:2),[0,1]);
            delSol = optimizeCbModel(modelDel,'max');
            tempFluxSol2{i} = delSol.x;
        end
    end
    emptyIDs = find(cellfun('isempty',tempFluxSol1));
    for i = 1: length(emptyIDs)
        tempFluxSol1{emptyIDs(i)} = zeros(length(modelMove.rxns),1);
    end
    emptyIDs = find(cellfun('isempty',tempFluxSol2));
    for i = 1: length(emptyIDs)
        tempFluxSol2{emptyIDs(i)} = zeros(length(modelMove.rxns),1);
    end
    SKOfluxes1 = cat(2,tempFluxSol1{:});
    SKOfluxes1 = round(SKOfluxes1*10^6)/10^6;
    SKOfluxes2 = cat(2,tempFluxSol2{:});
    SKOfluxes2 = round(SKOfluxes2*10^6)/10^6;
    essentialGeneSets = (SKOfluxes1(1,:) <= modelMove.lb(biomassRxnID) | SKOfluxes2(2,:) == 0)';
    clear tempFluxSolution delRxnSets
end

%% Arrange results of model reduction

modelReduced = modelMove;
modelReduced.oriRxns = model.rxns;
modelReduced.reductionStatus = reductionStatus;
if (bfRelation)
    modelReduced.rxnAssociations = rxnAssociations;
    modelReduced.rxnAssocMat = rxnAssocMat;
    modelReduced.unqGeneSetsMat = unqGeneSetsMat;
    modelReduced.geneSets = geneSets;
    modelReduced.geneSetRxnMat = geneSetRxnMat;
    modelReduced.geneSetAssocRxns = geneSetAssocRxns;
    modelReduced.essentialGeneSets = essentialGeneSets;
    % modelReduced.rxnFormulas = rxnFormulas;
    essentialRxnSets=geneSetAssocRxns(essentialGeneSets);
else
    essentialRxnSets={};
end

modelReduced.rules=updateGeneRule(modelReduced.rules, length(modelReduced.genes));
%% find essential genes & reactions
grRatio = singleGeneDeletion(modelReduced);
cmpEssentialGenesId=(grRatio<options.minBioPecentage/100); % less than 5%
modelReduced.essentialGenes= modelReduced.genes(cmpEssentialGenesId);

biomassRxn = modelReduced.rxns{biomassRxnID};
targetRxn = modelReduced.rxns{targetRxnID};
% oxygenRxn = modelReduced.rxns{oxygenRxnID};

if options.verbFlag == true
    if modelReduced.reductionStatus == 1
        disp('Model reduction succeeds.')
    end
    toc
end
end

function rulesOut=updateGeneRule(rules,nGenes)
rulesOut=rules;
nRules=length(rules);
for ri=1:nRules
    % find all gene indices
    infeasibleGenes=regexp(rulesOut{ri},'\d*','Match');
    for pj=infeasibleGenes
        if (str2num(pj{1})>nGenes) % assign true to the rule of genes that are not included in the model
            rulesOut{ri}=strrep(rulesOut{ri},['x(',pj{1},')'], 'true');
        end
    end
    if ~isempty(rulesOut{ri})
        rulesOut{ri}=['(',rulesOut{ri},')']; % nest rules
    end
end
end
