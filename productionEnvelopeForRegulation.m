function [biomassValues, targetValues, lineHandle] = productionEnvelopeForRegulation(model, LP, lineColor, targetRxn, biomassRxn, regType, nPts)

if (nargin < 3)
    lineColor = 'k';
end
if (nargin < 4)
    % Target flux
    targetRxn = '{ 1 EX_etoh(e)(24) }';
end
if (nargin < 5)
    % Biomass flux
    biomassRxn = '{ 1 Biomass_Ecoli_core_w_GAM(13) }';
end
if (nargin < 6)
    regType = 'fold'; % fold change of flux 
end
if (nargin < 7)
    nPts = 20;
end


% vector
vecBiomass=ismember(model.rxns,biomassRxn);
vecTarget=ismember(model.rxns,targetRxn);

% Run FBA to get upper bound for biomass
if strcmp(regType,'fold')
    LP.c=[zeros(length(vecBiomass),1); vecBiomass];
else
    LP.c=[vecBiomass;vecBiomass];
end
solMax = optimizeCbModel(LP,'max');
solMin = optimizeCbModel(LP,'min');

% Create biomass range vector
biomassValues = linspace(solMin.f,solMax.f,nPts);

% Max/min for target production
LP.C(end+1,:)=-LP.c';
LP.csense(end+1)=char('L');
LP.d(end+1,:)=0;

if strcmp(regType,'fold')
    LP.c=[zeros(length(vecTarget),1); vecTarget];
else
    LP.c=[vecTarget;vecTarget];
end
for i = 1:length(biomassValues)
    LP.d(end,:)=-floor(biomassValues(i)*1e6)/1e6; % change growth value
    sol = optimizeCbModel(LP,'max');
    if (sol.stat > 0)
        targetUpperBound(i) = sol.f;
    else
        targetUpperBound(i) = NaN;
    end
    sol = optimizeCbModel(LP,'min');
    if (sol.stat > 0)
        targetLowerBound(i) = sol.f;
    else
        targetLowerBound(i) = NaN;
    end
end

% Plot results
lineHandle=plot([biomassValues fliplr(biomassValues)],[targetUpperBound fliplr(targetLowerBound)],lineColor,'LineWidth',2);
axis tight;
hold on
pgon = polyshape([biomassValues fliplr(biomassValues)],[targetUpperBound fliplr(targetLowerBound)]);
plot(pgon,'FaceColor',lineColor,'FaceAlpha',0.1);
hold off

%ylabel([strrep(targetRxn,'_','-') ' (mmol/gDW h)']);
%xlabel('Growth rate (1/h)');

biomassValues = biomassValues';
targetValues = [targetLowerBound' targetUpperBound'];