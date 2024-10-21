
clear; clc; close all

% Set path and load the Cfg and Paths files
[projectPath,~,~]   = fileparts(mfilename('fullpath'));
projectPath         = extractBefore(projectPath,"3_analysis");
addpath (genpath (projectPath)); 

[Cfg, Paths, Log]       = getCfg(projectPath);
condNames               = Cfg.condNames;  
groups                  = Cfg.groups;


% Displaying version for the ITIs
% 1 = violin plots, 2 = median + std
ITIDisplay = 1;

%% Create the ITI file for the distance between median ITIs and possible meter periodicities

MetersName      = {'Duple';'Triple'};
possibleMeters  = {[400,800,1600];[600,1200]};

i = 1;

for iGroup = 1:2

    group   = groups{iGroup};
    load (fullfile(projectPath,['3_analysis/TappingAnalysis_',group,'.mat']));

    subjects        = fieldnames(ITI.(condNames{1}));
    idx2keep        = find(contains(subjects,'sub0'));
    subjects        = subjects(idx2keep);   

    for iSub = 1:length(subjects)
        for iCond = 1:length(condNames)

            for iTrial = 1:3 % length(fieldnames(ITI.(condNames{iCond}).(subjects{iSub})))-1

                try 
                    selectedITIs        = ITI.(condNames{iCond}).(subjects{iSub}).(['Trial',num2str(iTrial)]).ITI;
                catch % in one participant there is a missing tapping trial so I merge the two first trials to simulate a third one
                    selectedITIs        = [ITI.(condNames{iCond}).(subjects{iSub}).(['Trial',num2str(iTrial-2)]).ITI,...
                                           ITI.(condNames{iCond}).(subjects{iSub}).(['Trial',num2str(iTrial-1)]).ITI];
                end

                medianITI           = median(selectedITIs)*1000;                

                for iMeters = 1:2

                    % Table
                    Subject{i,1}        = [strrep(groups{iGroup},'icians',''),'_',...
                                           strrep(subjects{iSub},'sub','')]; 
                    Group{i,1}          = groups{iGroup};
                    Condition{i,1}      = condNames{iCond}; 
                    Trial{i,1}          = ['Trial',num2str(iTrial)];  
                    Meter{i,1}          = MetersName{iMeters};
                    FriedmanLine{i,1}   = [Subject{i,1},'_',Condition{i,1},'_',Trial{i,1}];

                    MedianITI{i,1}      = medianITI;
                    StdITI{i,1}         = std(selectedITIs)*1000;
                    Distance{i,1}       = min(abs(medianITI-possibleMeters{iMeters}));

                    i = i+1;
                end
            end
        end
    end
end

T = table(Subject, Group, Condition, Trial, Meter, FriedmanLine, MedianITI, StdITI, Distance);
T = sortrows(T,{'Group','Subject'}); 

% writetable(T,'StatisticalAnalysisTapping_DistanceMedianITI&PossibleMeters.xls')    

%% Create the ITI file for the distance between the 400 and 800 periodicities

periodicityNames      = {'400ms';'800ms'};
possiblePeriodicity   = {[400];[800]};

i = 1;

for iGroup = 1:2

    group   = groups{iGroup};
    load (fullfile(projectPath,['3_analysis/TappingAnalysis_',group,'.mat']));

    subjects        = fieldnames(ITI.(condNames{1}));
    idx2keep        = find(contains(subjects,'sub0'));
    subjects        = subjects(idx2keep);   

    for iSub = 1:length(subjects)
        for iCond = 1:length(condNames)

            for iTrial = 1:3 % length(fieldnames(ITI.(condNames{iCond}).(subjects{iSub})))-1

                try 
                    selectedITIs        = ITI.(condNames{iCond}).(subjects{iSub}).(['Trial',num2str(iTrial)]).ITI;
                catch % in one participant there is a missing tapping trial so I merge the two first trials to simulate a third one
                    selectedITIs        = [ITI.(condNames{iCond}).(subjects{iSub}).(['Trial',num2str(iTrial-2)]).ITI,...
                                           ITI.(condNames{iCond}).(subjects{iSub}).(['Trial',num2str(iTrial-1)]).ITI];
                end

                medianITI           = median(selectedITIs)*1000;                

                for iPeriodicities = 1:length(periodicityNames)

                    % Table
                    Subject{i,1}        = [strrep(groups{iGroup},'icians',''),'_',...
                                           strrep(subjects{iSub},'sub','')]; 
                    Group{i,1}          = groups{iGroup};
                    Condition{i,1}      = condNames{iCond}; 
                    Trial{i,1}          = ['Trial',num2str(iTrial)];  
                    Periodicity{i,1}    = periodicityNames{iPeriodicities};
                    FriedmanLine{i,1}   = [Subject{i,1},'_',Condition{i,1},'_',Trial{i,1}];

                    MedianITI{i,1}      = medianITI;
                    StdITI{i,1}         = std(selectedITIs)*1000;
                    Distance{i,1}       = min(abs(medianITI-possiblePeriodicity{iPeriodicities}));

                    i = i+1;
                end
            end
        end
    end
end

T = table(Subject, Group, Condition, Trial, Periodicity, FriedmanLine, MedianITI, StdITI, Distance);
T = sortrows(T,{'Group','Subject'}); 

% writetable(T,'StatisticalAnalysisTapping_DistanceMedianITI&DuplePeriodicities.xls')

%% Create the ITI file for the violin plots on R
i=1 ;
for iGroup = 1:2

    group   = groups{iGroup};
    load (fullfile(projectPath,['3_analysis/TappingAnalysis_',group,'.mat']));
    subjects        = fieldnames(ITI.(condNames{1}));
    idx2keep        = find(contains(subjects,'sub0'));
    subjects        = subjects(idx2keep); 

    for iSub = 1:length(subjects)
        for iCond = 1:length(condNames)


            subject         = subjects{iSub};
            condName        = condNames{iCond}; 
            trials          = fieldnames(ITI.(condName).(subject));

            for iTrial = 1:length(trials)-1
                for iITI = 1:length(ITI.(condName).(subject).(trials{iTrial}).ITI)

                    % Table
                    Subject{i,1}        = [strrep(groups{iGroup},'icians',''),'_',...
                                           strrep(subjects{iSub},'sub','')]; 
                    Group{i,1}          = groups{iGroup};
                    Condition{i,1}      = condNames{iCond}; 
                    Trial{i,1}          = ['Trial',num2str(iTrial)];  

                    data                = ITI.(condName).(subject).(trials{iTrial}).ITI(iITI);

                    if ~isequal(data,0) && data < 1.3
                        Iti(i,1)        = data;
                        i               = i+1;
                    end
                end
            end
        end
    end
end

T2 = table(Subject, Group, Condition, Trial, Iti);
T2 = sortrows(T2,{'Group','Subject'}); 

% writetable(T2,'NR2_ITI_violinPlots.csv') 


%% Figure preparation

figure('Position',[1 1 1500 700],'Color',[1 1 1]) 
       
plotPos = {[0.08 0.6 0.42 0.32],...
           [0.58 0.6 0.42 0.32],...
           [0.08 0.01 0.42 0.5],...
           [0.58 0.01 0.42 0.5]};    

fontsize    = 23;
linewidth   = 2;
Colors      = {[62, 100, 163]./256   , [34 55 89]./256       ,[62, 100, 163]./256  ; ... % bleu
               [245, 162, 97]./256   , [230, 108, 15]./256   ,[246, 175, 121]./256 ; ... % orange
               [42, 157, 143]./256   , [34, 129, 118]./256   ,[126, 221, 210]./256 ;...  % vert
               [42, 157, 143]./256   , [34, 129, 118]./256   ,[126, 221, 210]./256 };    % vert       

% musician and non-musician groups
annotation('textbox', [0.09,0.9,0.42,0.08], 'String','Musicians',...
           'fontsize',fontsize, 'EdgeColor','none', 'FontWeight','bold', 'HorizontalAlignment','center'); 
annotation('textbox', [0.58,0.9,0.42,0.08], 'String','Non-musicians',...
           'fontsize',fontsize, 'EdgeColor','none', 'FontWeight','bold', 'HorizontalAlignment','center'); 
       
% "A" and "B" annotations
annotation('textbox', [0.025,0.940,0.04,0.04], 'String','A',...
           'fontsize',fontsize+6, 'EdgeColor','none', 'FontWeight','bold', 'HorizontalAlignment','center'); % 0.009333333333333,0.939857142857143,0.04,0.04
       
annotation('textbox', [0.025,0.54,0.04,0.04], 'String','B',...
           'fontsize',fontsize+6, 'EdgeColor','none', 'FontWeight','bold', 'HorizontalAlignment','center');     % 0.014722222222222,0.530714285714286,0.04,0.04   
%% Load the data

    for iGroup = 1:length(groups)

        group   = groups{iGroup};
        data    = load (fullfile(projectPath,['3_analysis/TappingAnalysis_',group,'.mat']));

        Asynchronies.(group)  = data.Asynchrony;
        ITI.(group)         = data.ITI;
        Taps.(group)        = data.preprocessedTaps;

        for iCond = 1:length(condNames) 

            Asynchronies.(group).(condNames{iCond})= orderfields(Asynchronies.(group).(condNames{iCond}));
            ITI.(group).(condNames{iCond})         = orderfields(ITI.(group).(condNames{iCond}));
            Taps.(group).(condNames{iCond})        = orderfields(Taps.(group).(condNames{iCond}));
        end
    end
    
%% ITIs   

if ITIDisplay == 1
      
    maxNITIs = zeros(length(groups)*length(condNames),1);
    iMax     = 1;    

    for iGroup = 1:length(groups)
        
        axes('Position',plotPos{iGroup})
        iPos        = 1;
        group       = groups{iGroup};
        subjects    = fieldnames(Asynchronies.(group).(condNames{1}));
        idx2keep    = find(contains(subjects,'sub0'));
        subjects    = subjects(idx2keep);
      
        for iCond = [2,1,3,4]
            
            tempITIs = [];

            for iSubjects = 1:length(subjects)

                % extract all subjects' ITIs for each condition
                tempITIs = [tempITIs ITI.(group).(condNames{iCond}).(subjects{iSubjects}).AllTrials.ITI]; 
                tempITIs(find(tempITIs>1.7)) = [];
            end

            % Work around to have the right colors for each condition (I create another violin plot at -5 that will not be shown - comes from the fact that the violinplot function refuses to only plot one graph...)
            ITIs2plot = zeros(length(tempITIs),length(condNames))-5000;

            ITIs2plot(1:length(tempITIs),iPos) = tempITIs'.*1000;
            
            % get the blue color for ABAB2
            if iCond == 3;  iStim = 1;
            else;           iStim = iCond;
            end

            [~,~,~,~,~] = violin(ITIs2plot,'x',1:4,'facecolor',Colors{iStim,1},... 
                                      'edgecolor',Colors{iStim,2},...
                                      'facealpha',0.4,...
                                      'mc',[],'medc',[],'bw',10); hold on 
                                  
          maxNITIs(iMax) = length(tempITIs);
          iMax = iMax +1;   
          iPos = iPos+1; 
        end
        
        % Equal violin plot spread across both groups (doesn't work...)
%         fakeITIs = zeros(max(maxNITIs),4)-5000;
%         [~,~,~,~,~] = violin(fakeITIs,'x',1:4,'facecolor',[1 1 1],... 
%                                   'edgecolor',[1 1 1],...
%                                   'mc',[],'medc',[],'bw',10);        
                                                      
        box off
        set(gca,'TickDir','out','fontsize',fontsize,'LineWidth',linewidth,...
                'ylim',[0 1300],'ytick',[0:400:1200],...
                'xticklabel',{[]}, 'xlim',[0.5 4.5])                              

        for iLines = 0:200:1600
            if iLines == 400 || iLines == 800
                line([0 10],[iLines iLines],'LineStyle','--','Color',[0.5 0.5 0.5],'LineWidth',1.75)
            else
                line([0 10],[iLines iLines],'LineStyle',':','Color',[0.6 0.6 0.6],'LineWidth',1.75)
            end
        end
        
        h = gca;
        h.XAxis.Visible = 'off';    

        if iGroup == 1
            ylabel('ITIs (in ms)','fontsize',fontsize+2, 'FontWeight', 'bold')  
        else
            h.YAxis.Visible = 'off';
        end    
    end    
end



if ITIDisplay == 2
    for iGroup = 1:length(groups)
        
        axes('Position',plotPos{iGroup})
        
        group           = groups{iGroup};
        subjects        = fieldnames(Asynchronies.(group).(condNames{1}));
        idx2keep        = find(contains(subjects,'sub0'));
        subjects        = subjects(idx2keep);
        
        % plot position parameters
        iPos            = 1;
        range           = 0.7;
        intervals       = range/length(subjects);
            
        for iCond = [2,1,3,4]
            
             % get the blue color for ABAB2
            if iCond == 3;  iStim = 1;
            else;           iStim = iCond;
            end           
            
            for iSubjects = 1:length(subjects)
                
                % get median and std (but should also try with the 5-95 percentiles or standard error of the mean)
                % The SD quantifies scatter — how much the values vary from one another. The SEM quantifies how precisely you know the true mean of the population.
                ITIs = ITI.(group).(condNames{iCond}).(subjects{iSubjects}).AllTrials.ITI .* 1000;
                ITIs(find(ITIs>1700)) = [];
                
                medITI      = median(ITIs);
                stdITI      = std(ITIs);
                prctile5    = interp1(linspace(1/numel(ITIs),1,numel(ITIs)), sort(ITIs), 0.05);
                prctile95   = interp1(linspace(1/numel(ITIs),1,numel(ITIs)), sort(ITIs), 0.95);
                                            
                xPos = (iPos-range/2) + (iSubjects-1) * intervals;
                scatter(xPos, medITI,'MarkerEdgeColor',Colors{iStim,2},'MarkerFaceColor',Colors{iStim,1},'MarkerFaceAlpha',0.4,'SizeData',50); 
                hold on
                
                %line([xPos xPos], [medITI-stdITI,  medITI+stdITI],'Color',Colors{iStim,1},'LineWidth',0.75)
                line([xPos xPos], [prctile5 prctile95],'Color',Colors{iStim,1},'LineWidth',0.75)
               
            end
            iPos = iPos + 1;
        end
               
        box off
        set(gca,'TickDir','out','fontsize',fontsize,'LineWidth',linewidth,...
                'ylim',[0 1600],'ytick',[0:400:1600],...
                'xticklabel',{[]}, 'xlim',[0.5 4.5],'XTick',[1:4])
            
        h = gca;
        h.XAxis.Visible = 'off';    

        if iGroup == 1
            ylabel('ITIs (in ms)','fontsize',fontsize+2, 'FontWeight', 'bold', 'Position', [0.07, 800,-1])  
        else
            h.YAxis.Visible = 'off';
        end    
            
        for iLines = 0:200:1600
            if iLines == 400 || iLines == 800
                line([0 10],[iLines iLines],'LineStyle','--','Color',[0.5 0.5 0.5],'LineWidth',1.75)
            else
                line([0 10],[iLines iLines],'LineStyle',':','Color',[0.6 0.6 0.6],'LineWidth',1.75)
            end
        end     
    end
end
          
        
%% Asynchronies
    
for iGroup = 1:length(groups)
         
        axes('Position',plotPos{2+iGroup})
    
        group           = groups{iGroup};
        subjects        = fieldnames(Asynchronies.(group).(condNames{1}));
        idx2keep        = find(contains(subjects,'sub0'));
        subjects        = subjects(idx2keep);        
    
        % store data 
        for iSubjects = 1:length(subjects)    
            for iCond = 1:length(condNames)
                
                Asynch(iSubjects,iCond) = mean(Asynchronies.(group).(condNames{iCond}).(subjects{iSubjects}).AllTrials.R);
            end
        end    

        iPos = 1;
        
        
        for iCond = [2,1,3,4] 
            x               = repmat(iPos,length(subjects),1); 
            jitterAmount    = 0.08;
            jitterValuesX   = 2*(rand(size(x))-0.5)*jitterAmount;
            
            if iCond == 3
                iStim = 1;
            else 
                iStim = iCond;
            end
            
            scatter(x+jitterValuesX ,Asynch(:,iCond),... %'Marker',marker,...
                            'MarkerEdgeColor',Colors{iStim,2},...
                            'MarkerFaceColor',Colors{iStim,1},...
                            'MarkerFaceAlpha',0.4,...
                            'SizeData',75)

            hold on
            iPos = iPos + 1;
        end        
 
       
        % boxplot (ordering made in the boxplot function)
        iBox = 1;
        for iCond = 1:length(condNames)
            boxAsynch(iBox:iBox+length(subjects)-1) = Asynch(:,iCond);
            boxCond(iBox:iBox+length(subjects)-1)   = repmat(iCond,length(subjects),1);
            
            iBox = iBox+length(subjects);
        end
       
        clear Asynch 
        
        h = boxplot(boxAsynch,boxCond,'Colors','k','Widths',0.3,'Positions',[2,1,3,4],'symbol', '','Labels','','Whisker',0);    %'Positions',xPos{iGroup} 
        set(h,{'linew'},{linewidth-0.25})
        set(gca, 'xlim',[0.5 4.5])

    set(gca,'TickDir','out','LineWidth',linewidth,'fontsize',fontsize,'YTick',[0 1],'XTick',[1:4],'YLim',[0,1])
    box off
    
    % multi-row tilted labels
    row1 = {'Medium Pattern' 'Long Pattern' 'Long Pattern' 'No Pattern'};
    row2 = {'Repetition' 'Repetition #1' 'Repetition #2' 'Repetition'};    
    labelArray = [row1; row2];
    labelArray = strjust(pad(labelArray),'center'); 
    tickLabels = strtrim(sprintf('%s\\newline%s\n', labelArray{:}));
    ax = gca(); 
    ax.XTick = 1:4;
    ax.TickLabelInterpreter = 'tex';
    ax.XTickLabel = tickLabels; 
    ax.FontWeight = 'bold';
    xtickangle(45)
    
    if iGroup == 1
        ylabel('Tapping Stability','fontsize',fontsize+2, 'FontWeight', 'bold', 'Position', [0.07, 0.5,-1])  
    else
        h = gca;
        h.YAxis.Visible = 'off';    
    end  
            
    clear boxAsynch boxCond   
end
  
cd(fullfile(projectPath,'4_figures'))
figName = 'Fig3_Tapping Analysis';
print([figName,'.jpg'],'-djpeg','-r600') 
exportgraphics(gcf,[figName,'.eps'],...   
    'ContentType','vector',...
    'BackgroundColor','none')    