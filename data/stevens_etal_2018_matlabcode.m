% This is the primary script for running the evolutionary simulations
% To run this simulation, copy all scripts in a single directory, and
% create a directory named "data". 

%%%%%%%%%%%%%%%%%%%%%%%
% Clear and assign variables
%%%%%%%%%%%%%%%%%%%%%%%
clear;  % clear variables
clearvars;  % clear variables
numIndivs=90;   % number of individuals
numGenerations= 250;     % number of generations
memoryParam1=-1;    % lambda for memory function         (-1: no application of memory function/perfect recall)
memoryParam2=-1;    % psi for memory function	(decay rate)
numStrategies=9;    % number of strategies
numStructures=3;    % number of contact patterns (skew)
numMemories=7;      % number of forgetting rates
numMemoryTypes=2;   % number of memory error types 1: random, 2: forget
numRepetitions=1000;  % number of simulation repetitions
imperfectOwnActions=1;  % own action can be forgotten (if =1)
numChosen = 20;     % number of highest fitness agents used for truncation selection operator
deterministic = 0;  % flag for deterministic reproduction (0 = reproduce equally) or stochastic (1 = reproduce with uniform probability 1/p) truncation selection operator
bAnalysisOnly=0;    % flag for test analysis with limited condition, repetitions, and generations and no data recording (0=no; 1=yes)
bDistractor=0;      % flag for including additional distractor interactions (0=no; 1=yes)
numDistractionsPerIndiv=100;    % number of distractor interactions
fileDir='./data/';  % directory name for output: comment when running on Microsoft Windows
% fileDir='.\data\';  % directory name for output: uncomment when running on Microsoft Windows
tic;    % start timer
if bAnalysisOnly    % when running test analysis, assign parameters
    numGenerations=1;
    numRepetitions=1;
    numMemoryTypes=1;
    numMemories=1;
end

%%%%%%%%%%%%%%%%%%%%%%%
% Primary control code
%%%%%%%%%%%%%%%%%%%%%%%
for MemoryT=1:numMemoryTypes    % for each memory error type
    for MemoryI=1:numMemories   % for each forgetting rate
        if or((MemoryT==1), (MemoryI>1))   
            for structuresI=1:numStructures     % for each contact pattern
                pairStrengthSequence=zeros(numIndivs-1,1);
                switch structuresI  % Contact patterns
                    case 1  % High skew
                        pairStrengthSequence(1:10,1)=[33 23 15 10 7 5 3 2 1 1];                               % number of interactions with the 10 interaction partners (high skew)
                        structName='VerySkewed';
                        numContacts=10;
                    case 2  % Low skew
                        pairStrengthSequence(1:10,1)=[19 16 14 11 10 8 7 6 5 4];								% number of interactions with the10 interaction partners (low skew)
                        structName='LessSkewed';
                        numContacts=10;     
                    case 3  % No skew
                        pairStrengthSequence(1:10,1)=10;																		% 10 interactions with each of the 10 interaction partners (no skew)
                        structName='Equal';
                        numContacts=10;
                end
                switch MemoryI  % Probability of corrrect recall: p = λ * ((k + 1) ^ (− ψ)) 
											  %		 λ (lambda) starting point of the forgetting function 
											  % 		ψ (psi) decay rate of the forgetting function 
											  % 		k number of interactions between the target and the current interaction (k = 0 if the target interaction was the most recent interaction)
                    case 1						% no forgetting
                        Lambda=1;
                        Psi=0;
                          memName='L1P00';
                    case 2
                        Lambda=1;
                        Psi=0.10;
                          memName='L1P10';
                    case 3
                        Lambda=1;
                        Psi=0.25;
                        memName='L1P25';
                    case 4
                        Lambda=0.25;
                        Psi=0.50;
                        memName='L25P50';  
                    case 5  % Empirical forgetting rate
                        Lambda=0.8775247;
                        Psi=0.2289332;
                        memName='L87P22';      
                    case 6
                        Lambda=0.4;
                        Psi=0.25;
                        memName='L40P25';       
                    case 7
                        Lambda=0.6;
                        Psi=0.25;
                        memName='L60P25';       
                end
                switch MemoryT  % Memory error types
                    case 1  % Errors of commission
                        MemoryType='Random';
                        memTName='RND';
                    case 2  % Errors of omission
                        MemoryType='Forget';
                        memTName='FGT';
                end        
                typeProportions(1:numStrategies,1)=1/numStrategies;                   %proportion of agents playing each type of strategy in the starting population
                numTypes=size(typeProportions,1);													%number of different strategies still played in the population			 
                matrixSave=zeros(numRepetitions, numGenerations, numStrategies);				
                matrixOutcomeSave=zeros(numRepetitions, numGenerations, numStrategies);
                matrixInteractionSave=zeros(numRepetitions, numGenerations, numStrategies,numStrategies);                  
                matrixInteractionMeanSave=zeros(numRepetitions, numGenerations, numStrategies,numStrategies);
                matrixInteractionCountSave=zeros(numRepetitions, numGenerations, numStrategies,numStrategies);
                matrixCooperationRatioSave=zeros(numRepetitions, numGenerations);																% stores overall level of cooperation for all generations in all runs
                matrixCooperationRatioStratSave=zeros(numRepetitions, numGenerations,numStrategies);							% stores the ratio of strategies for all generations in all runs
                for ci = 1:numRepetitions                                                 % runs through all simulation repetitions
                    disp(strcat(int2str(ci),': ',num2str(toc)));													
                    TypeHist=zeros(numGenerations,numTypes);											% stores the proportion of agents playing each type of strategy for all generations this simulation run
                    OutcomeHist=zeros(numGenerations,numTypes);									% stores the score for each type of strategy for all generations this simulation run
                    InteractionHist=zeros(numGenerations,numTypes,numTypes);            % 
                    InteractionMeanHist=zeros(numGenerations,numTypes,numTypes);
                    CountHist=zeros(numGenerations,numTypes,numTypes);
                    CoopRatioHist=zeros(numGenerations,1);
                    CoopRatioHistStrat=zeros(numGenerations,numTypes);
                    playerTypes=  createPopulation(numIndivs,typeProportions);																			%contains the strategies (1-9) of all agents
                    pairConnections = createPairConnections( numIndivs, numContacts, pairStrengthSequence );               % creates a matrix specifying the number of interaction of each agent with each other agent  
                    for gi=1:numGenerations
                        disp(strcat('T',int2str(ci),'G',int2str(gi),': ',num2str(toc)));																				%outputs simulation number, generation number and time stamp
                        contactList=createContactListRandomOrderDistractor(pairConnections,bDistractor,numDistractionsPerIndiv);
                        if bAnalysisOnly																											% only for flag for test analysis with limited condition, repetitions, and generations and no data recording 
                            sumConnect=sum(pairStrengthSequence,1);
                            maxConnect=max(pairStrengthSequence,[],1);
                            [summary1, summary2]=analyzeContactList(contactList, numIndivs, numContacts,sumConnect,maxConnect);   % outputs matrix with interaction frequency of all agents with all other agents; and mean interaction number for each agent
                        else 
                            [indivOutcomes, typeMatrix, typeMeanMatrix, typeCountMatrix, cooperationRatio, cooperationRatioStrat]=playGames(contactList,playerTypes,Lambda,Psi,MemoryType,numStrategies,imperfectOwnActions);
                            a=cooperationRatioStrat(1,:);
                            if a<1
                                b=a;
                            end
                            % Selection operators (uncomment desired selection operator
                            [playerTypes, typeCount, typeOutcome]=generateNewPopulationRoulette(indivOutcomes, playerTypes,numTypes);
%                             [playerTypes, typeCount, typeOutcome]=generateNewPopulationSUS(indivOutcomes, playerTypes,numTypes);
%                             [playerTypes, typeCount, typeOutcome]=generateNewPopulationTruncation(indivOutcomes, playerTypes,numTypes,numChosen,deterministic);

                            TypeHist(gi,1:numTypes)=typeCount(1:numTypes,1)./numIndivs;
                            OutcomeHist(gi,1:numTypes)=typeOutcome(1:numTypes,1);
                            InteractionHist(gi,1:numTypes,1:numTypes)=typeMatrix(1:numTypes,1:numTypes);
                            InteractionMeanHist(gi,1:numTypes,1:numTypes)=typeMeanMatrix(1:numTypes,1:numTypes);
                            CountHist(gi,1:numTypes,1:numTypes)=typeCountMatrix(1:numTypes,1:numTypes);
                            CoopRatioHist(gi,1)=cooperationRatio;
                            CoopRatioHistStrat(gi,1:numTypes)=cooperationRatioStrat(1:numTypes);
                        end   
                    end
                    if ~bAnalysisOnly 
                        matrixSave(ci, 1:numGenerations, 1:numTypes)=TypeHist(1:numGenerations, 1:numTypes);
                        matrixOutcomeSave(ci, 1:numGenerations, 1:numTypes)=OutcomeHist(1:numGenerations, 1:numTypes);
                        matrixInteractionSave(ci,1:numGenerations,1:numTypes,1:numTypes)=InteractionHist(1:numGenerations, 1:numTypes, 1:numTypes);
                        matrixInteractionMeanSave(ci,1:numGenerations,1:numTypes,1:numTypes)=InteractionMeanHist(1:numGenerations, 1:numTypes, 1:numTypes);
                        matrixInteractionCountSave(ci,1:numGenerations,1:numTypes,1:numTypes)=CountHist(1:numGenerations, 1:numTypes, 1:numTypes);
                        matrixCooperationRatioSave(ci,1:numGenerations)=CoopRatioHist(1:numGenerations,1);
                        matrixCooperationRatioStratSave(ci,1:numGenerations,1:numTypes)=CoopRatioHistStrat(1:numGenerations,1:numTypes);
                    end
                end
                if ~bAnalysisOnly 
                    matrixInteractionSave=squeeze(nanmean(matrixInteractionSave,1));
                    matrixInteractionMeanSave=squeeze(nanmean(matrixInteractionMeanSave,1));
                    matrixInteractionCountSave=squeeze(nanmean(matrixInteractionCountSave,1));
                    partMatrix=zeros(numRepetitions, numGenerations, numTypes);
                    partMatrixMean=zeros(numStrategies,numGenerations);
                    partMatrixSurvival=zeros(numStrategies,numGenerations);
                    partMatrixDominance=zeros(numStrategies,numGenerations);
                    partOutcomeMean=squeeze(nanmean(matrixOutcomeSave,1));
                    meanTempRatioStratSave=squeeze(nanmean(matrixCooperationRatioStratSave,1));
                    forSummaryMatrix=zeros(numRepetitions,numTypes);
                    for ti=1:numTypes 
                        partMatrixInteraction=zeros(numGenerations,numStrategies);
                        partMatrixInteractionMean=zeros(numGenerations,numStrategies);
                        partMatrixCountMean=zeros(numGenerations,numStrategies);
                        fileName=strcat(fileDir,'SingleProportion','_M-',memName,'_E-',memTName,'_S-',structName,'_','T',int2str(ti),'.txt');
                        partMatrix=matrixSave(1:numRepetitions, 1:numGenerations, ti);
                        partMatrixMean(ti,1:numGenerations)=mean(partMatrix,1);
                        partMatrixSurvival(ti,1:numGenerations)=mean(partMatrix>0.00000000001,1);
                        partMatrixDominance(ti,1:numGenerations)=mean(partMatrix>0.9999999999,1);
                        dlmwrite(fileName, partMatrix, '\t');
                        forSummaryMatrix(1:numRepetitions,ti)=partMatrix(1:numRepetitions, numGenerations);
                        fileName=strcat(fileDir,'Interaction','_M-',memName,'_E-',memTName,'_S-',structName,'_','T',int2str(ti),'.txt');
                        partMatrixInteraction(1:numGenerations,1:numStrategies)=squeeze(matrixInteractionSave(1:numGenerations,ti,1:numStrategies));
                        dlmwrite(fileName, partMatrixInteraction, '\t');
                        fileName=strcat(fileDir,'InteractionMean','_M-',memName,'_E-',memTName,'_S-',structName,'_','T',int2str(ti),'.txt');
                        partMatrixInteractionMean(1:numGenerations,1:numStrategies)=squeeze(matrixInteractionMeanSave(1:numGenerations,ti,1:numStrategies));
                        dlmwrite(fileName, partMatrixInteractionMean, '\t');
                        fileName=strcat(fileDir,'InteractionCount','_M-',memName,'_E-',memTName,'_S-',structName,'_','T',int2str(ti),'.txt');
                        partMatrixCountMean(1:numGenerations,1:numStrategies)=squeeze(matrixInteractionCountSave(1:numGenerations,ti,1:numStrategies));
                        dlmwrite(fileName, partMatrixCountMean, '\t');  
                    end
                    tftSelection=[3; 4;  5;  8];
                    tftMatrix= forSummaryMatrix(:,tftSelection);
                    tftSumMatrix=sum(tftMatrix,2);  
                    fileName=strcat(fileDir,'EndStateProportions','_M-',memName,'_E-',memTName,'_S-',structName,'.txt');
                    dlmwrite(fileName, [forSummaryMatrix tftSumMatrix], '\t');  
                    shortSummaryMatrix=zeros(4,numTypes+1);
                    shortSummaryMatrix(1,1:numTypes)=mean(forSummaryMatrix,1);
                    shortSummaryMatrix(1,numTypes+1)=mean(tftSumMatrix,1);
                    shortSummaryMatrix(2,1:numTypes)=std(forSummaryMatrix,0,1);
                    shortSummaryMatrix(2,numTypes+1)=std(tftSumMatrix,0,1);
                    shortSummaryMatrix(3,1:numTypes)= meanTempRatioStratSave(numGenerations,1:numStrategies);
                    shortSummaryMatrix(3,numTypes+1)=    sum((shortSummaryMatrix(3,tftSelection).*shortSummaryMatrix(1,tftSelection)),2) / shortSummaryMatrix(1,numTypes+1);
                    shortSummaryMatrix(4,1)=mean(matrixCooperationRatioSave(:,numGenerations),1);
                    shortSummaryMatrix(4,2)=std(matrixCooperationRatioSave(:,numGenerations),0,1);
                    fileName=strcat(fileDir,'shortSummaryPaperXXXXX','_M-',memName,'_E-',memTName,'_S-',structName,'.txt');
                    dlmwrite(fileName, shortSummaryMatrix, '\t');
                    fileName=strcat(fileDir,'cooperationRatioStrat','_M-',memName,'_E-',memTName,'_S-',structName,'.txt');
                    partmatrixCooperationRatioStratSave(1:numGenerations,1:numStrategies)=meanTempRatioStratSave(1:numGenerations,1:numStrategies);
                    dlmwrite(fileName, partmatrixCooperationRatioStratSave, '\t');
                    fileName=strcat(fileDir,'cooperationRatio','_M-',memName,'_E-',memTName,'_S-',structName,'.txt');
                    dlmwrite(fileName, matrixCooperationRatioSave', '\t');
                    fileName=strcat(fileDir,'MeanProportion','_M-',memName,'_E-',memTName,'_S-',structName,'.txt');
                    dlmwrite(fileName, partMatrixMean', '\t');
                    fileName=strcat(fileDir,'MeanOutcome','_M-',memName,'_E-',memTName,'_S-',structName,'.txt');
                    dlmwrite(fileName, partOutcomeMean, '\t');
                    fileName=strcat(fileDir,'MeanSurvival','_M-',memName,'_E-',memTName,'_S-',structName,'.txt');
                    dlmwrite(fileName, partMatrixSurvival', '\t');
                    fileName=strcat(fileDir,'MeanDominance','_M-',memName,'_E-',memTName,'_S-',structName,'.txt');
                    dlmwrite(fileName, partMatrixDominance', '\t');
                    set(0,'DefaultAxesColorOrder',[1 0 0;0 1 0;0 0 1;0.5 0.5 0.5; 0 1 1; 1 1 0; 1 0 1],...
                      'DefaultAxesLineStyleOrder','-|--|:');
                    hFig = figure(1);% from figure
                    plot(partMatrixMean');
                    title('Type numbers per generation')
                    xlabel('Generation')
                    ylabel('Type Frequency')
                    legend('ALLC','ALLD','TFT','TF2T','GTFT','WSLS','RAND','CTFT','GRIM'); % Add a legend in the upper left:
                    Figname=strcat(fileDir,'fig_',memName,'_E-',memTName,'_',structName,'.png');
                    print(hFig, '-dpng', Figname);
                end
            end
        end
    end
end
disp(toc);  % display timer

%%%%%%%%%%%%%%%%%%%%%%%
% analyzeContactList
%%%%%%%%%%%%%%%%%%%%%%%
function  [summary1, summaryMean,summaryQuartiles] = analyzeContactList(contactList, numIndivs, numContacts,sumConnect,maxConnect)
    % outputs matrix with interaction frequency of all agents with all other agents; and mean interaction number for each agent
    % is called only when test analysis turned on (bAnalysisOnly=1: limited condition, repetitions, and generations and no data recording)
	%		input:
	% 	contactList: matrix of specifying the number of interaction of each agent with each other agent
	% 	numIndivs: number of agents in this generation
	% 	numContacts: number of interaction partners 
	% 	sumConnect: number of interactions of an agent
	% 	maxConnect: maximum number of interactions with another agent
	% output:
	%		summary1:  matrix with interaction frequency of all agents with all other agents 
	%     summaryMean: mean interaction number for each agent
	%     summaryQuartiles: outputs 25% and 75% quantiles for maximum number of interactions
	
    numDyads=size(contactList,1);						% 
    summary1=zeros(numIndivs,sumConnect);
    summaryQuartiles=zeros(numIndivs,3);
    countContacts=zeros(numIndivs,1);
    lastContact=zeros(numIndivs*1000+numIndivs,1)-1;
    contactCodes=contactList(:,1)*1000+contactList(:,2);
    for ci=1:numDyads
        if contactCodes(ci,1)>0
            if (lastContact(contactCodes(ci,1))<0)
               lastContact(contactCodes(ci,1))=ci; 
            else
                countContacts(contactList(ci,1))=countContacts(contactList(ci,1))+1;
                countContacts(contactList(ci,2))=countContacts(contactList(ci,2))+1;
                count1=0;
                count2=0;
                for di=lastContact(contactCodes(ci,1)):ci
                    if ((contactList(di,1)==contactList(ci,1))||(contactList(di,2)==contactList(ci,1)))
                        count1=count1+1;
                    end
                    if ((contactList(di,1)==contactList(ci,2))||(contactList(di,2)==contactList(ci,2)))
                        count2=count2+1;
                    end
                end
                 summary1( contactList(ci,1),countContacts(contactList(ci,1))  )=count1;
                 summary1( contactList(ci,2),countContacts(contactList(ci,2))  )= count2;
               lastContact(contactCodes(ci,1))=ci; 
            end
        end
    end
    maxContacts=max(countContacts);
    summaryMean=mean(summary1(:,1:maxContacts),2);
    summaryQuartiles(:,1:3)=quantile(summary1(:,1:maxContacts),3,2);
end

%%%%%%%%%%%%%%%%%%%%%%%
% createContactListRandomOrderDistractor
%%%%%%%%%%%%%%%%%%%%%%%
function [ connectionList ] = createContactListRandomOrderDistractor( pairConnections,bDistractor,numDistractionsPerIndiv )
    % reshuffles the number of interactions of each agent with all other agents and adds  distractor interactions
    % input:
	% 	pairConnections: matrix specifying the number of interaction of each agent with each other agent
	%     bDistractor: additional distractor interactions are included
	%     numDistractionsPerIndiv:  number of distractor interactions
	% output:
	%		connectionList: matrix of specifying the number of interaction of each agent with each other agent including distractor interactions (marked as "-1")
    
    numIndivs=size(pairConnections,1);
    sumConnections=sum(pairConnections,1);
    sumConnections=sum(sumConnections,2);
    connectionList=zeros(sumConnections, 2);
    assigned=0;
    for ca=1:numIndivs-1
         for cb=ca:numIndivs
              newConns=pairConnections(ca,cb);
              connectionListTemp(assigned+1:assigned+newConns,1)=ca;
              connectionListTemp(assigned+1:assigned+newConns,2)=cb;
              assigned=assigned+newConns;
         end
    end    
    if bDistractor         % if additional distractor interactions are included
        newList=zeros(size(connectionListTemp,1)+numIndivs*numDistractionsPerIndiv,1);
        newList(1:size(connectionListTemp,1),1:2)=connectionListTemp(1:size(connectionListTemp,1),1:2);
        newList(size(connectionListTemp,1)+1:size(connectionListTemp,1)+numIndivs*numDistractionsPerIndiv,1)=-1;        % add -1 into all distractor interaction slots
        for di=1:numIndivs
            newList( size(connectionListTemp,1)+(di-1)*numDistractionsPerIndiv:size(connectionListTemp,1)+di*numDistractionsPerIndiv,2)=di;
        end
        connectionListTemp=newList;
        newOrder=randperm(size(connectionListTemp,1));
        connectionList(1:size(connectionListTemp,1),1:2)=connectionListTemp(newOrder,1:2);
    else                        % if no additional distractor interactions are included
        newOrder=randperm(sumConnections);
        connectionList(1:sumConnections,1:2)=connectionListTemp(newOrder,1:2);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%
% createPairConnections
%%%%%%%%%%%%%%%%%%%%%%%
function [ pairConnections ] = createPairConnections( numIndivs, numContacts, pairStrengthSequence )
    %summary: create a pair connection strength matrix
    % will create a matrix of pair connections for
    %	numIndivs Individuals with numContacts contacts: the strength of each
    %	contact in sequence are contained in pairSterengthSequence. The
    %	resulting structure is deterministically created based on a
    %	round-robin tournament
	
     % input:
	% 	numIndivs: total number of agents
	%     numContacts: number of interaction partners
	%     pairStrengthSequence: number of interactions with each partner (high skew, low skew, now skew)
	% output:
	%		pairConnections: matrix containing the number of contacts of each agent with each other agent in one generation
	
    pairConnections=zeros(numIndivs, numIndivs);
    for pairIndex=1:numContacts
        covered=zeros(numIndivs,1);
        pairConnections(min(numIndivs, pairIndex),max(numIndivs, pairIndex))=pairStrengthSequence(pairIndex,1);
        covered(numIndivs)=1;
        covered(pairIndex)=1;
        goalNumber=mod(2*pairIndex, numIndivs-1);
        for testIndiv=1:numIndivs
            if covered(testIndiv,1)==0
                firstNumber=mod(testIndiv, numIndivs-1);
                secondNumber=goalNumber-firstNumber;
                while (secondNumber<1)||(covered(secondNumber,1)~=0)
                    secondNumber=secondNumber+numIndivs-1;
                end
                covered(testIndiv,1)=1;
                covered(secondNumber,1)=1;
                pairConnections(min(testIndiv, secondNumber),max(testIndiv, secondNumber))=pairStrengthSequence(pairIndex,1);
            end
        end
    end 
end

%%%%%%%%%%%%%%%%%%%%%%%
% createPopulation
%%%%%%%%%%%%%%%%%%%%%%%
function [ playerTypes ] = createPopulation( numIndivs,typeProportions )
    % summary: outputs line vector containing the strategy of every single agent
   %	input:
    %		numIndivs: total number of agents
    %		typeProportions: proportion of the population playing each strategy
    %	output:
	%		playerTypes: line vector containing the strategy of every agent (in randomized order)
    
    playerTypes=zeros(numIndivs,1);
    b=sum(typeProportions,2);
    b=sum(b,1);
    if round(b)==1
        %values are relative frequencies
        numTypes=size(typeProportions,1);
        currentPos=1;
        for pi= 1: numTypes
            for ci=1:round(typeProportions(pi,1)*numIndivs)
                playerTypes(currentPos,1)=pi;
                currentPos=currentPos+1;
            end
        end
    else
        % values are total frequencies
        numTypes=size(typeProportions,1);
        currentPos=1;
        for pi= 1: numTypes
            for ci=1:typeProportions(pi,1)
                playerTypes(currentPos,1)=pi;
                currentPos=currentPos+1;
            end
        end
    end
    newOrder=randperm(numIndivs);
    playerTypes(1:numIndivs)=playerTypes(newOrder);
end

%%%%%%%%%%%%%%%%%%%%%%%
% gameOutcome
%%%%%%%%%%%%%%%%%%%%%%%
function [ actionA,actionB ] = gameOutcome( typeA,typeB,historyA, historyB )
    % determines actions of both agents in a single interaction
   
   %input:
    %		typeA: strategy of agent A
    %		typeB: strategy of agent B
    %		historyA: past interactions of agent A with agent B
    %		historyB: past interactions of agent B with agent A
    %output:
	%		actionA: action of agent A this round (cooperate(1) or defect(2))
    %		actionB: action of agent Bthis round (cooperate(1) or defect(2))
    
    numPastGamesA=size(historyA,1);
    numPastGamesB=size(historyB,1);
    if numPastGamesA==0                         
        historyA=[];
    end
    if numPastGamesB==0
       historyB=[];  
    end
    actionA=reaction(historyA,typeA);         % next action of agent A in dependence of the past interaction history with the same partner and the agent's strategy
    actionB=reaction(historyB,typeB);		   % next action of agent B in dependence of the past interaction history with the same partner and the agent's strategy
end

%%%%%%%%%%%%%%%%%%%%%%%
% gamePayoff
%%%%%%%%%%%%%%%%%%%%%%%
function [ payoff ] = gamePayoff( actionA,actionB )
    % summary: outputs the payoff for agent A this round
	
   % input:
    %		actionA: action of agent A this round (cooperate(1) or defect(2))
	%		actionB: action of agent Bthis round (cooperate(1) or defect(2))
    % output:
	%		payoff: payoff for agent A:
	%			A cooperates, B cooperates: 3
	%			A cooperates, B defects: 0
	%			A defects, B cooperates: 5
	%			A defects, B defects: 1
		
    if actionA==1 && actionB==1
        payoff=3;
    elseif actionA==1 && actionB==2
        payoff=0;
    elseif actionA==2 &&  actionB==1
        payoff=5;
    elseif actionB==2 && actionA==2
        payoff=1;
    end  
end

%%%%%%%%%%%%%%%%%%%%%%%
% generateNewPopulationRoulette
%%%%%%%%%%%%%%%%%%%%%%%
function [ playerTypes, typeCount, typeOutcome] = generateNewPopulationRoulette( indivOutcomes, playerTypes ,numTypes)
    % Creates on new generation of agents based on the payoff in the current generation
    %   The number of agents playing a strategy in the new generation is chosen. 
	%   The robability of reproductive success is proportional to the total score achieved in a generation.
	%
    % input:
	%		indivOutcomes: total payoff for each agent
	%		playerTypes: contains the strategies (1-9) of all agents
	%		numTypes: number of different strategies still played in the population
    %		
    % output:
	%		playerTypes: strategies that survive to the new genration
	%		typeCount: number of agents playing each strategy in new generation
	%     typeOutcome: average payoff for an agent playing each of the surviving strategies in the last generation
    
	typeCount=zeros(numTypes,1); 
    typeOutcome=zeros(numTypes,1);
    typeCountold=zeros(numTypes,1); 
    numIndivs=size(playerTypes,1);
    totalFitness=sum(indivOutcomes,1);
    relativeFitness=indivOutcomes./totalFitness;
    fitnessSum=zeros(numTypes,1);
    for ci=1:numIndivs
        fitnessSum(playerTypes(ci,1))=fitnessSum(playerTypes(ci,1))+relativeFitness(ci,1);
        typeOutcome(playerTypes(ci,1),1)=typeOutcome(playerTypes(ci,1),1)+indivOutcomes(ci,1);
        typeCountold(playerTypes(ci,1),1)=typeCountold(playerTypes(ci,1),1)+1;
    end
    typeOutcome=typeOutcome./typeCountold;
    for ci=1:numIndivs
        ri=rand(1,1);
        typei=0;
        while ri>= 0
            typei=typei+1;
            ri=ri-fitnessSum(typei,1); 
        end
        playerTypes(ci,1)=typei;
        typeCount(typei,1)=typeCount(typei,1)+1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%
% generateNewPopulationSUS
%%%%%%%%%%%%%%%%%%%%%%%
function [ playerTypes, typeCount, typeOutcome] = generateNewPopulationSUS( indivOutcomes, playerTypes ,numTypes)
    % Creates on new generation of agents based on the payoff in the current generation
    %   The number of agents playing a strategy in the new generation is chosen. 
	%   The robability of reproductive success is proportional to the total score achieved in a generation.
    %
    % input:
	%		indivOutcomes: total payoff for each agent
	%		playerTypes: contains the strategies (1-9) of all agents
	%		numTypes: number of different strategies still played in the population
    %		
    % output:
	%		playerTypes: strategies that survive to the new genration
	%		typeCount: number of agents playing each strategy in new generation
	%     typeOutcome: average payoff for an agent playing each of the surviving strategies in the last generation

    typeCount=zeros(numTypes,1); 						%number of agents playing each strategy
    typeOutcome=zeros(numTypes,1);
    typeCountold=zeros(numTypes,1); 
    numIndivs=size(playerTypes,1);							%total number of agents
    totalFitness=sum(indivOutcomes,1);					%total fitness: sum of all agent's scores
    relativeFitness=indivOutcomes./totalFitness;		%proportional fitness: proportion of individual scores of total fitness
    fitnessSum=zeros(numIndivs,1);
    fitnessSum(1,1)=relativeFitness(1,1);   
    typeOutcome(playerTypes(1,1),1)=typeOutcome(playerTypes(1,1),1)+indivOutcomes(1,1);
    typeCountold(playerTypes(1,1),1)=typeCountold(playerTypes(1,1),1)+1;
    for ci=2:numIndivs
        fitnessSum(ci,1)=fitnessSum(ci-1,1)+relativeFitness(ci,1);
        typeOutcome(playerTypes(ci,1),1)=typeOutcome(playerTypes(ci,1),1)+indivOutcomes(ci,1);
        typeCountold(playerTypes(ci,1),1)=typeCountold(playerTypes(ci,1),1)+1;
    end
    typeOutcome=typeOutcome./typeCountold;
    constantDistance=1/numIndivs;
    searchValue=rand(1,1)*constantDistance;
    searchI=1;
    for ci=1:numIndivs
        while(fitnessSum(searchI)<searchValue)
            searchI=searchI+1;
        end
        searchValue=searchValue+constantDistance;
        playerTypes(ci,1)=playerTypes(searchI);
        typeCount(playerTypes(ci,1),1)=typeCount(playerTypes(ci,1),1)+1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%
% generateNewPopulationTruncation
%%%%%%%%%%%%%%%%%%%%%%%
function [ playerTypes, typeCount, typeOutcome] = generateNewPopulationTruncation( indivOutcomes, playerTypes ,numTypes, numChosen, deterministic)
    %   The number of agents playing a strategy in the new generation is chosen. 
	%   number of agents that can reproduce is truncated to the top (numChosen) agents  
	%   	If deterministic=0: reproduce equally  
   	%   	If deterministic=1: reproduce with uniform probability 1/p
	
    %input:
	%		indivOutcomes: total payoff for each agent
	%		playerTypes: contains the strategies (1-9) of all agents
	%		numTypes: number of different strategies still played in the population
    %		numChosen: number of Individuals from the Top that can reproduce 
    %  	deterministic: whether they reproduce evenly (not randomly chosen) or with uniform probabilities (1/p)
    %output:
	%		playerTypes: strategies that survive to the new genration
	%		typeCount: number of agents playing each strategy in new generation
	%     typeOutcome: average payoff for an agent playing each of the surviving strategies in the last generation
	
    typeCount=zeros(numTypes,1); 				%number of agents playing each strategy
    typeOutcome=zeros(numTypes,1);						
    typeCountold=zeros(numTypes,1); 						
    numIndivs=size(playerTypes,1);							%total number of agents
    totalFitness=sum(indivOutcomes,1);					%total fitness: sum of all agent's scores
    relativeFitness=indivOutcomes./totalFitness;	%proportional fitness: proportion of individual scores of total fitness
    fitnessSum=zeros(numTypes,1);
    for ci=1:numIndivs													
        fitnessSum(playerTypes(ci,1))=fitnessSum(playerTypes(ci,1))+relativeFitness(ci,1);
        typeOutcome(playerTypes(ci,1),1)=typeOutcome(playerTypes(ci,1),1)+indivOutcomes(ci,1);
        typeCountold(playerTypes(ci,1),1)=typeCountold(playerTypes(ci,1),1)+1; 
    end
    typeOutcome=typeOutcome./typeCountold;
    [sortOutcomes,IndexOutcomes]=sort(indivOutcomes,1);
    relevantIndivs(1:numChosen)=IndexOutcomes(numIndivs-numChosen+1:numIndivs,1);
    relevantTypes(1:numChosen)=playerTypes(relevantIndivs);
    if (deterministic==1)     
        fractionEven=floor(numIndivs/numChosen);
        fractionRemainder=numIndivs-fractionEven*numChosen;
        for ci=1:numChosen
            playerTypes((ci-1)*fractionEven+1:ci*fractionEven,1)=relevantTypes(1,ci);
            typeCount(relevantTypes(1,ci),1)=typeCount(relevantTypes(1,ci),1)+fractionEven;
        end
        if fractionRemainder>0
            playerTypes(numIndivs-fractionRemainder+1:numIndivs,1)...
              = relevantTypes(1,numChosen-fractionRemainder+1: numChosen);
            typeCount( relevantTypes(1,numChosen-fractionRemainder+1: numChosen),1)...
              = typeCount( relevantTypes(1,numChosen-fractionRemainder+1: numChosen),1)+1;
        end
    else      
        for ci=1:numIndivs 
            ri=randi(numChosen);
            playerTypes(ci,1)=relevantTypes(1,ri);
            typeCount( playerTypes(ci,1),1)=typeCount( playerTypes(ci,1),1)+1;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%
% Apply forgetting to memory or retrieve random chosen past action
%%%%%%%%%%%%%%%%%%%%%%%
function [ newHistoryA, newHistoryB ] = memoryTransform(currentTimeA, currentTimeB, oldHistory,memoryParam1,memoryParam2 ,memoryType,imperfectOwnActionsRandom)
    % creates new Past Interaction Memories for both agents in an interaction
    %  applies errors of omission (past interactions forgotten: recall of previous  memory instance) and/or errors of commission (Past interactions misremembered) (set by the 'MemoryType' handle)
    % input:
	%		currentTimeA: current total interaction number of agent A in this generation
	%		currentTimeB: current total interaction number of agent B in this generation
	%		oldHistory: past interaction history between agent A and agent B
	%		memoryParam1: lambda in  forgetting function
	%		memoryParam2: psi in forgetting function (decay rate)
	%		memoryType: Types of memory errors the agents commit
	%				 MemoryType='Random' : Errors of commission (Past interactions misremembered)
    % 			 MemoryType='Forget' : Errors of omission (past interactions forgotten: recall of previous  memory instance)
	%		imperfectOwnActionsRandom: agent imperfectly remembers own actions as well (if =1)
	% output:
	% 	newHistoryA: updated interaction history of agent A
	%		newHistoryB: updated interaction history of agent B
    lambda=memoryParam1;
    psi=memoryParam2;
    numInteractions= size(oldHistory,1);
    if numInteractions>0
        % empirical error rate
        % power function: p=1-92(1+n) ^-0.08 (p error: n: number of intervening
        % interractions
        RememberA(1:numInteractions,1)=lambda.*((currentTimeA-oldHistory(1:numInteractions,3)).^(-psi));      %rembered interaction history of agent A with agent B in agent A's "mind": forgetting function applied to interaction history
        RememberA=(rand(numInteractions,1)-RememberA(1:numInteractions,1))<0;												   %all remembered instances of agent A's interaction history with agent B marked								   
        RememberB(1:numInteractions,1)=lambda.*((currentTimeB-oldHistory(1:numInteractions,4)).^(-psi));		%rembered interaction history of agent B with agent A in agent B's "mind": forgetting function applied to interaction history
        RememberB=(rand(numInteractions,1)-RememberB(1:numInteractions,1))<0;													%all remembered instances of agent B's interaction history with agent A marked		
        if strcmp(memoryType,'Forget')											%Errors of omission (past interactions forgotten: recall of previous  memory instance)
            rmA= find( RememberA(:,1)>0);
            rmB= find( RememberB(:,1)>0);
            a=[1 2 3];
            b=[2 1 4];
            newHistoryA=oldHistory(rmA,a);
            newHistoryB=oldHistory(rmB,b);
        elseif strcmp(memoryType,'Random')										% Errors of commission (Past interactions misremembered)
            randomA=random('bino',1,0.5,numInteractions,1);
            randomB=random('bino',1,0.5,numInteractions,1);
            ForgA= (1-( RememberA(:,1)>0)).*randomA;
            ForgB= (1- ( RememberB(:,1)>0)).*randomB;
            if imperfectOwnActionsRandom								% agent imperfectly remembers own actions as well
                randomA1=random('bino',1,0.5,numInteractions,1);
                randomB1=random('bino',1,0.5,numInteractions,1); 
                ForgA1= (1-( RememberA(:,1)>0)).*randomA1;
                ForgB1= (1- ( RememberB(:,1)>0)).*randomB1;
            end
            a=[1 2 3];
            b=[2 1 4];
            newHistoryA=oldHistory(:,a);
            newHistoryB=oldHistory(:,b);
            newHistoryA(1:numInteractions,2 )= (newHistoryA(1:numInteractions,2)-1)+ForgA.*(1-2*(newHistoryA(1:numInteractions,2)-1)) +1;
            newHistoryB(1:numInteractions,2)= (newHistoryB(1:numInteractions,2)-1)+ForgB.*(1-2*(newHistoryB(1:numInteractions,2)-1)) +1;
            if imperfectOwnActionsRandom
               newHistoryA(1:numInteractions,1 )= (newHistoryA(1:numInteractions,1)-1)+ForgA1.*(1-2*(newHistoryA(1:numInteractions,1)-1)) +1;
               newHistoryB(1:numInteractions,1)= (newHistoryB(1:numInteractions,1)-1)+ForgB1.*(1-2*(newHistoryB(1:numInteractions,1)-1)) +1;
            end    
        end
    else
        newHistoryA=zeros(0,3);
        newHistoryB=zeros(0,3);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%
% playGames
%%%%%%%%%%%%%%%%%%%%%%%
function [ indivOutcomes, typeMatrix, typeMeanMatrix, typeCountMatrix,cooperationRatio,cooperationRatioStrat ] = playGames( connectionList,playerTypes,memoryParam1,memoryParam2,memoryType ,numStrategies,imperfectOwnActions)
    % Main function that simulates all agent interactions in a generation
    %   
    % input: 
	%		connectionList: interaction number of each agents with each other agent
	%		playerTypes: strategies of all agents
	%		memoryParam1: lambda for memory function
	%		memoryParam2: psi for memory function (decay rate)
	%		memoryType: Types of memory errors the agents commit
	%				 MemoryType='Random' : Errors of commission (Past interactions misremembered)
    % 			 MemoryType='Forget' : Errors of omission (past interactions forgotten: recall of previous  memory instance)
	%		numStrategies: number of strategies
	%		imperfectOwnActionsRandom: agent imperfectly remembers own actions as well (if =1)
	% output:
	%		indivOutcomes: total payoff for each agent
	%		typeMatrix: contains total payoff for matchups of each strategy with each other strategy
	%		typeMeanMatrix: contains average payoff per agent for matchups of each strategy with each other strategy
	%		typeCountMatrix: contains number of agents for matchups of each strategy with each other strategy
	%		cooperationRatio: overall ratio of cooperation within generation
	%		cooperationRatioStrat: ratio of cooperation (C) of each strategy
    numGames=size(connectionList,1);
    numIndivs=size(playerTypes,1);				% number of individuals
    indivOutcomes=zeros(numIndivs,1);
    typeMatrix=zeros(numStrategies);
    typeCountMatrix=zeros(numStrategies);
    outcomeList=zeros(numGames,2);
    playerGames=zeros(numIndivs,numIndivs,30);
    indivRounds=zeros(numIndivs,numIndivs,30);
    interactionCount=zeros(numIndivs,1);
    cooperationCount=0; 
    cooperationRatioStrat=zeros(numStrategies,1);					         % ratio of cooperation decision for each of the strategies played 
    for gameIndex=1:numGames
        playerA=min(connectionList(gameIndex,1),connectionList(gameIndex,2));
        playerB=max(connectionList(gameIndex,1),connectionList(gameIndex,2));
        if playerA<0 %distractorEvent
            interactionCount(playerB,1)=interactionCount(playerB,1) +1; 
        end     
        if playerA>0 % if incomplete gamelist, last numbers are 0,0 for players
            historyLength= playerGames(playerA,playerB,1);
            history=zeros(historyLength,4);
            interactionCount(playerA,1)=interactionCount(playerA,1) +1;
            interactionCount(playerB,1)=interactionCount(playerB,1) +1;
            history(1:historyLength,1:2)=outcomeList(playerGames(playerA,playerB,2:historyLength+1),1:2);
            history(1:historyLength,3)=indivRounds(playerA,playerB,1:historyLength); 
            history(1:historyLength,4)=indivRounds(playerB,playerA,1:historyLength); 
            if memoryParam1>-1          % if lambda for memory function > -1: apply memory function
                [historyA,historyB]=memoryTransform(interactionCount(playerA,1),interactionCount(playerB,1),history,memoryParam1,memoryParam2,memoryType,imperfectOwnActions);  
            else
                a=[1 2 3];
                b=[2 1 4];
                historyA=history(a);
                historyB=history(b);
            end
        [actionA,actionB]=gameOutcome(playerTypes(playerA,1),playerTypes(playerB,1),historyA, historyB);                    % determines actions of both agents in a single interaction
        outcomeList(gameIndex,1)=actionA;
        outcomeList(gameIndex,2)=actionB;
        if (actionA==1)
            cooperationCount=cooperationCount+1;
            cooperationRatioStrat(playerTypes(playerA,1),1)=cooperationRatioStrat(playerTypes(playerA,1),1)+1;
        end
        if (actionB==1)
            cooperationCount=cooperationCount+1;
            cooperationRatioStrat(playerTypes(playerB,1),1)=cooperationRatioStrat(playerTypes(playerB,1),1)+1;
        end
        playerGames(playerA,playerB,1)=playerGames(playerA,playerB,1)+1;
        playerGames(playerA,playerB,playerGames(playerA,playerB,1)+1)=gameIndex;
        indivRounds(playerA,playerB,playerGames(playerA,playerB,1))= interactionCount(playerA,1);
        indivRounds(playerB,playerA,playerGames(playerA,playerB,1))= interactionCount(playerB,1);
        payoffA=gamePayoff(actionA,actionB);
        payoffB=gamePayoff(actionB,actionA);
        indivOutcomes(playerA)=indivOutcomes(playerA)+payoffA ;
        indivOutcomes(playerB)=indivOutcomes(playerB)+ payoffB;
        typeCountMatrix(playerTypes(playerA,1),playerTypes(playerB,1))=typeCountMatrix(playerTypes(playerA,1),playerTypes(playerB,1))+1;
        typeCountMatrix(playerTypes(playerB,1),playerTypes(playerA,1))=typeCountMatrix(playerTypes(playerB,1),playerTypes(playerA,1))+1;
        typeMatrix(playerTypes(playerA,1),playerTypes(playerB,1))= typeMatrix(playerTypes(playerA,1),playerTypes(playerB,1))+payoffA;
        typeMatrix(playerTypes(playerB,1),playerTypes(playerA,1))= typeMatrix(playerTypes(playerB,1),playerTypes(playerA,1))+payoffB;
        end      
    end
    typeMeanMatrix=typeMatrix./typeCountMatrix;
    cooperationRatio=cooperationCount/(2*numGames);
    countInteractions=sum(typeCountMatrix,2);
    cooperationRatioStrat=cooperationRatioStrat./countInteractions;
end

%%%%%%%%%%%%%%%%%%%%%%%
% Agent's next reaction
%%%%%%%%%%%%%%%%%%%%%%%
function [ action ] = reaction( history,pType )
    % determines the next action of an agent in dependence of 
		% a) the past interaction history with the same partner
		% b) the agent's strategy
    %   input:
	%	         history: past interaction history with interaction partner
	%			 ptype: strategy of the agent
	%   output 
	%          action: agent's next action 
	%					1: cooperation 
	%					2: defection

    
    % enter memory-updater
    switch pType
        case 1 %Always Cooperate ALLC
            action=1;
        case 2 % Always Defect ALLD
            action=2;
        case 3 %TIT-for-TAT TFT
            numPreviousGames=size(history,1);    
            if numPreviousGames==0						%if first interaction with specific partner: cooperate
                action=1;
            else
                action=history(numPreviousGames,2);	%else copy partner's previous move
            end
        case 4 %Tit-For-Two-Tat TF2T
                numPreviousGames=size(history,1);
            if numPreviousGames<2						%if first or second interaction with specific partner: cooperate
                action=1;
            else
               if ( history(numPreviousGames,2)==1) || (history(numPreviousGames-1,2)==1)
                   action=1;											%if partner cooperated in any of the last two interactions: cooperate
               else
                   action=2;											%if partner defected in the last two interactions: defect
               end
            end
        case 5 %Generous Tit-For-Tat GTFT 
            numPreviousGames=size(history,1);	
            if numPreviousGames==0						%if first interaction with specific partner: cooperate
                action=1;
            else
                randR=rand(1,1);									
                if history(numPreviousGames,2)==1     %if partner cooperated last round: cooperate with 99% probability, defect with 0% probability
                    if randR<.99
                        action=1;
                    else
                        action=2;
                    end                    
                elseif randR<.33										%if partner defected last round: cooperate with 33% probability, defect with 66% probability									
                    action=1;
                else
                    action=2;
                end              
            end
        case 6 %Win-Stay-Lose-Shift (Pavlov) WSLS
            numPreviousGames=size(history,1);
            if numPreviousGames==0						%if first interaction with specific partner: cooperate
                action=1;
            elseif (history(numPreviousGames,1)==history(numPreviousGames,2))       %if the agent and the partner both cooperated or both defected last round: cooperate
                action=1;
            else																															%if the agent and the partner acted different from each other last round: cooperate
                action=2;						
            end
        case 7      %Random       RAND: cooperate with 50% probability, defect with  50% probability
 
            randR=rand(1,1);
            if randR<.50											
                    action=1;
            else
                    action=2;
            end     
        case 8 %Contrite Tit-For-Tat CTFT
            % cooperate-> cooperate, defect -> only cooperate
            % if you defected two rounds ago and partner did not
            numPreviousGames=size(history,1);
            if numPreviousGames==0							%if first interaction with specific partner: cooperate
                action=1;
            elseif (numPreviousGames>2)&&((history(numPreviousGames-1,2)==1)&&(history(numPreviousGames-1,1)==2))    %if two rounds prior both the partner cooperated and the agent defected: cooperate
                action=1;
            elseif (history(numPreviousGames,2)==1)								%if partner cooperated previous round: cooperate
                action=1;
            else																								%if partner defected last round and two rounds prior (the agent cooperated or the partner defected): defect
                action=2;
            end 
        case 9 %Grim Trigger GRIM
            % http://en.wikipedia.org/wiki/Grim_trigger
            numPreviousGames=size(history,1);
            if numPreviousGames==0							%if first interaction with specific partner: cooperate
                action=1;
            elseif (history(numPreviousGames,1)==2)||(history(numPreviousGames,2)==2)      %if partner defected last round or agent defected last round: defect
                action=2;
            else								%if agent and partner cooperated last round: cooperate
                action=1;
            end  
    end 
end


