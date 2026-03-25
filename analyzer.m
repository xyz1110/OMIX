
load data\candidates_iML1515-EX_succ_e.mat;
load data\reduced\reduced_iML1515-EX_succ_e.mat

model = reducedModel;
irreSign=1; % the reactions in the solution are from an irreversible model
global scale;
scale = 0.05;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if irreSign
    %Flip all pure backward reactions
    backReacs = model.lb < 0 & model.ub <= 0;

    model.S(:,backReacs) = -model.S(:,backReacs);
    templbs = -model.ub(backReacs);
    model.ub(backReacs) = -model.lb(backReacs);
    model.lb(backReacs) = templbs;
    model.c(backReacs) = - model.c(backReacs); %Also flip the objective coefficient, as otherwise the target changes.

    % Convert to irreversible rxns
    [imodel,matchRev,~,irrev2rev] = convertToIrreversible(model,'OrderReactions',false);

    model=imodel;
end

%%==========================================
solution={
    {'EX_ac_e/ACtex(3)'}
    {'SUCDi(3)'}
    {'POR5_f(3)'}
    {'POR5_b(3)'}
    {'PSERT/PSP_L(2)'}
    {'PDH(1)'}
};


targetRxn='EX_succ_e';

F = [];
for i = 1 : length(solution)
    F = [F,0]; 
end

%%=========================================


LP=testFV(model, solution, targetRxn,  F);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function LP=testFV(model, solution, targetRxn,  F)
biomassRxn=model.rxns{model.c==1};

solution0={
    };
F0 = [];
LP0 = buildLPRegFold(model, solution0, targetRxn, F0);
[~, ~, h(1)]=productionEnvelopeForRegulation(model, LP0, 'k', targetRxn, biomassRxn);
hold on;

LP=buildLPRegFold(model, solution, targetRxn, F);

if LP.stat
    disp('Plot production envelopes with F:');
    disp(F);
    [~, ~, h(end+1)]=productionEnvelopeForRegulation(model,LP, 'b',targetRxn, biomassRxn, 'fold'); % able to plot PE of design strategies containing up/down-regulations by F fold
end

end



function LP=buildLPRegFold(model, regList, targetRxn, F)

assert(length(regList)==length(F))

% regList should follow the format {'FUM(1)','ATPS4rpp(2)','DXPS(3)'}
% 1: up-regulation, 2: down-regulation, 3: knockout

global scale;

% wild typle flux bounds
[w_lb, w_ub] =deal(model.lb, model.ub);
solWT=optimizeCbModel(model);
w_lb(model.c==1)=floor(solWT.f*1e6)*1e-6; 
% Set the mutant minimum growth rate as a fraction of the wild-type optimum
m_lb = model.lb;
m_lb(model.c==1)=solWT.f*scale;

% mutant flux bounds
[m_lb, m_ub]=deal(m_lb, model.ub);

[~, nRxns]=size(model.S);

upDown=zeros(nRxns,1);
ko_vector=zeros(1, nRxns);

rxnIDs=[];

Fi=zeros(length(regList),1);
for i=1:length(regList)
    strs=split(regList{i}, '(');
    [rxn, action]=deal(strs{1},str2double(strs{2}(1)));
    rxnID=findRxnIDs(model, rxn);
    rxnIDs(end+1)=rxnID;
    if action == 1
        upDown(rxnID)=1; % F>
        Fi(i)=1;
    elseif action == 2
        upDown(rxnID)=-1; 
        Fi(i)=-1; % F>0.5
    elseif action==3
       ko_vector(rxnID)=1; % F <tol
    end
end

[~, idx]=sort(rxnIDs,'ascend');
F=columnVector(F(idx));
Fi=Fi(idx);

m_ub(ko_vector==1)=F(Fi==0); % use F value to limit knockout flux

% set objective
c = zeros(size(model.c));
obj=c;
c(findRxnIDs(model, targetRxn))=1;
obj=[obj;c];

% create LP constraints: Ax=b, Cx<=d
A=[blkdiag(model.S, model.S); % su=0, sv=0
    %sparse(sum(ko_vector),nRxns), selMatrix(ko_vector) %v=0 for knockouts
    ]; 

b=zeros(size(A,1),1); %b=0

C=[];
if any(Fi==1)
    C=[C;(1+F(Fi==1)).*selMatrix(upDown==1),-selMatrix(upDown==1)]; % up: (1+F)*u-v<=0
end

if any(Fi==-1)
    C=[C; -selMatrix(upDown==-1), (1+F(Fi==-1)).*selMatrix(upDown==-1)]; % down: -u+(1+F)*v<=0
end
d=zeros(size(C,1),1);

csense=char('E'*ones(size(A,1),1)); 
dsense=char('L'*ones(size(d)));

lb=[w_lb; m_lb]; 
ub=[w_ub; m_ub];

[LP.c, LP.S, LP.dxdt, LP.C, LP.d, LP.csense, LP.dsense, LP.lb, LP.ub] = deal(obj, A, b, C, d, csense, dsense, lb, ub);

[solution] = optimizeCbModel(LP, 'min');
for i = 1:size(model.S,2)
    fprintf('%d\t \t%8.3f\t \t%8.3f\n', i, solution.x(i),solution.x(i+size(model.S,2)));  % 打印格式：编号、反应名称、速率（带3位小数）
end
if solution.stat==1
    disp('Model is feasible')
    LP.stat=1;
else
    disp('Model is infeasible')
    LP.stat=0;
end
end