function [Res] = multimodal_cell_match(cellmatch, meancellradius, clusteroffset,cludteroffsetid)
% This function determines the best matching paris of cells between in vivo
% ca2+ imaging and post hoc histology.
%
% Inputs:
%          cellmatch : metrics for matching in vivo recorded cells post hoc
%          in histology
%
%          meancellradius : the mean diameter of in vivo recorded cells
%
%          clusteroffset : a vector containing the distances of the
%          offsets for each cluster of in vivo recorded cells
%
%          cludteroffsetid : a vector of the same size as cells recorded
%          in in vivo Ca2+ imaging, assigning each cell to a cluster
%
% Outputs:  
%          cellmatchindex : indices for matched cells         
%                      
%          solveparameter : a matchingscore to quantify the ACCURRACY? of
%          the matched cells
%
%
% Function is written by Philip Anner (2020)

maxEucD = max(max(cellmatch.euclidean));
cellmatch.NormEuclidean = cellmatch.euclidean/maxEucD;
for i=1:size(cellmatch.pdf,1)
    prob = optimproblem('ObjectiveSense','max');
    
    % get matching parameters for each ca2+ cell
    pdf =  reshape(cellmatch.pdf(i,:),[],1);
    NormEuclidean= reshape(cellmatch.NormEuclidean(i,:),[],1);

    zdim = log( reshape(cellmatch.zdim(i,:),[],1));
    
    somajaccard =  reshape(cellmatch.somajaccard(i,:),[],1);
    allareajaccard=  reshape(cellmatch.allareajaccard(i,:),[],1);
    
    % define optimization parameters
    optpdf = optimvar('pdf',length(pdf),'Type','integer','LowerBound',0,'UpperBound',1);
    opteuclidean = optimvar('NormEuclidean', length(NormEuclidean) ,1,'Type','integer','LowerBound',0,'UpperBound',1);
    optzdim = optimvar('zdim', length(zdim) ,1,'Type','integer','LowerBound',0,'UpperBound',1);
    optsomajaccard =  optimvar('somajaccard', length(somajaccard) ,1,'Type','integer','LowerBound',0,'UpperBound',1);
    optallareajaccard =  optimvar('allareajaccard', length(allareajaccard) ,1,'Type','integer','LowerBound',0,'UpperBound',1);
    
    
    % define constraints
    % only find 1 match
    constmatchcellseuc = sum(opteuclidean)==1;
    prob.Constraints.opteuclidean = constmatchcellseuc;
    
    constmatchcellspdf = sum(optpdf)==1;
    prob.Constraints.optpdf = constmatchcellspdf;
    
    constmatchcellszdim = sum(optzdim)==1;
    prob.Constraints.optzdim = constmatchcellszdim;
    
    constmatchcellssomajaccard = sum(optsomajaccard)==1;
    prob.Constraints.optsomajaccard = constmatchcellssomajaccard;
    
    constmatchcellsallareajaccard = sum(optallareajaccard)==1;
    prob.Constraints.optallareajaccard = constmatchcellsallareajaccard;
    
    constEucD = NormEuclidean'*opteuclidean <=  nansum([meancellradius clusteroffset(cludteroffsetid(i))]);            
    prob.Constraints.maxEucD = constEucD;
      
    % Constraint euc dist > minimum value
    %make sure to select same cell
    constmatcheuclideanopt = optpdf == opteuclidean;
    prob.Constraints.consteucpdf = constmatcheuclideanopt;
    
    constmatchzdimopt = optzdim == opteuclidean;
    prob.Constraints.constzdimeuc = constmatchzdimopt;
    
    constmatchsomajaccardopt = optsomajaccard == opteuclidean;
    prob.Constraints.constsomajaccardeuc = constmatchsomajaccardopt;
    
    constmatchallareajaccardopt = optallareajaccard == opteuclidean;
    prob.Constraints.constallareajaccardeuc = constmatchallareajaccardopt;
    
    % optimization function
    prob.Objective = pdf'*optpdf + somajaccard'*optsomajaccard + allareajaccard'* optallareajaccard - NormEuclidean'*opteuclidean;
    
    %solve it
    [sol,p,c] = solve(prob);
    
    %find cell match
    if c ~= -2
        Res.cellmatchindex(i,1) = find(round(sol.NormEuclidean));
        Res.solveparameter(i,1)=p;
        Res.Normsolveparameter(i,1) = p/3;
    end
end

end