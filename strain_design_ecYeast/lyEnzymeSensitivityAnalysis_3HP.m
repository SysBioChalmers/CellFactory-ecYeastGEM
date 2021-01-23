current = pwd;
%Load ecYeast model
load('ec3HP-asp.mat')
%Carbon sources and media types
Csources   = {'D-glucose'};
MediaTypes = {'Min'};
%Get enzymes indexes (mets)
enzIndxs = find(contains(model.mets,'prot_'));
%Avoid protein pool
enzIndxs = enzIndxs([1:963 965 966 967]);
%Map proteins to gene IDs
proteins     = strrep(model.mets(enzIndxs),'prot_','');
[~,geneIndx] = ismember(proteins,model.enzymes);
genes        = model.enzGenes(geneIndx);

objIndex = find(model.c);
cd (current)
for media = MediaTypes
    sensitivities = [];
    disp(media{1})
    cd (current)
       for source = Csources
        disp(source{1})
        uptakeRxn = strcat(source,' exchange (reversible)');
        [Model,~] = lychangeMedia_batch(model,uptakeRxn{1},media);
        sol = solveLP(Model);
        column = zeros(length(enzIndxs),1);
        if ~isempty(sol.f)
            disp(['Base growth :' num2str(-sol.f)])
            for i=1:length(enzIndxs)
                tempModel = Model;
            	FCC       = 0; 
                KcatIndexes = find(tempModel.S(enzIndxs(i),:));
                %Just get the indexes of the kinetic coefficients and avoid
                %the enzyme usage reaction indexes
                KcatIndexes = KcatIndexes(1:end-1);
                %Multiply all coefficients
                tempModel.S(enzIndxs(i),KcatIndexes) = model.S(enzIndxs(i),KcatIndexes)/1.001;
                newSol = solveLP(tempModel);
                if ~isempty(newSol.f)
                    %If second optimization is feasible, then fix the
                    %optimal objective value and minimize the total protein
                    %usage
                    tempModel.lb(objIndex) = -0.9999999*newSol.f;
                    tempModel.ub(objIndex) = -newSol.f;
                    tempModel.c            = zeros(length(tempModel.c),1);
                    tempModel.c(enzIndxs)  = -1;
                    newSol = solveLP(tempModel);
                    if ~isempty(newSol.f)
                        objValue = newSol.x(objIndex);
                        %Calculate obj control coefficient for the i-th
                        %enzyme
                        FCC      = 1000*(objValue-(-sol.f))/abs(sol.f);
                        if FCC<1E-6
                            FCC = 0;
                        end
                    end
                end
                column(i) = FCC;
            end
        end
        sensitivities = [sensitivities, column];
    end
    %Write the results for all the carbon sources in a given media type in
    %tab separated file
    MediaTable = table(sensitivities,'RowNames',genes);
    fileName   = ['KcatSensitivities_' media{1}, '.txt'];
    cd ../result_ecYeast/3HP
    writetable(MediaTable, fileName,'Delimiter','\t','QuoteStrings',true,'WriteRowNames',true);
end