% load image registration parameters

% settings for Image registration
Options=struct('Similarity','cc','Registration','Both','Penalty',0.0005,'MaxRef',2,'Grid',[],'Spacing',[8 8],'MaskMoving',[],'MaskStatic',[],'Verbose',0,'Points1',[],'Points2',[],'PStrength',[],'Interpolation','Linear','Scaling',[1 1]);

%% Process tail vein injection recording (RECDSA)
% Detect blood vessels in vivo
InVivo_DetectBloodVessels

%% Load in vivo Ca2+ imaging sessions 
%% Session 1 - 2411 
[REC2411] = preprocessRecording('.\Data\In vivo Ca2+\2411\',0,1);

% Define landmarks for image registration
if exist('Sampledata','var')
    [REC2411.Registration.movingPoints,REC2411.Registration.fixedPoints] = cpselect(REC2411.minrec,RECDSA.minrec,Sampledata.REC2411.Registration.movingPoints,Sampledata.REC2411.Registration.fixedPoints, 'Wait', true);
else
    [REC2411.Registration.movingPoints,REC2411.Registration.fixedPoints] = cpselect(REC2411.minrec,RECDSA.minrec,'Wait', true);
end 
REC2411.Registration.movingPoints=fliplr(REC2411.Registration.movingPoints);
REC2411.Registration.fixedPoints=fliplr(REC2411.Registration.fixedPoints);
REC2411.Registration.PStrength = ones(size(REC2411.Registration.fixedPoints,1),1);
REC2411.Registration.PStrength = REC2411.Registration.PStrength .* 0.95;
Options.Points1=REC2411.Registration.movingPoints;
Options.Points2=REC2411.Registration.fixedPoints;
Options.PStrength=REC2411.Registration.PStrength;

[REC2411.Ireg,REC2411.O_trans,REC2411.Spacing,~,~,~] = image_registration(REC2411.minrec,RECDSA.minrec,Options);

% Assess registration result
% implay(cat(3,REC2411.Ireg,RECDSA.minrec))

% Align all images of  the recording
REC2411_Reg = registerrecording(REC2411);
% The original struct is not needed anymore
clear REC2411
%% Session 2 - 2811 
[REC2811] = preprocessRecording('.\Data\In vivo Ca2+\2811\',0,1);

% Define landmarks for image registration
if exist('Sampledata','var')
    [REC2811.Registration.movingPoints,REC2811.Registration.fixedPoints] = cpselect(REC2811.minrec,RECDSA.minrec,Sampledata.REC2811.Registration.movingPoints,Sampledata.REC2811.Registration.fixedPoints, 'Wait', true);
else
    [REC2811.Registration.movingPoints,REC2811.Registration.fixedPoints] = cpselect(REC2811.minrec,RECDSA.minrec, 'Wait', true);
end
REC2811.Registration.movingPoints=fliplr(REC2811.Registration.movingPoints);
REC2811.Registration.fixedPoints=fliplr(REC2811.Registration.fixedPoints);
REC2811.Registration.PStrength = ones(size(REC2811.Registration.fixedPoints,1),1);
REC2811.Registration.PStrength = REC2811.Registration.PStrength .* 0.95;
Options.Points1=REC2811.Registration.movingPoints;
Options.Points2=REC2811.Registration.fixedPoints;
Options.PStrength=REC2811.Registration.PStrength;
[REC2811.Ireg,REC2811.O_trans,REC2811.Spacing,~,~,~] = image_registration(REC2811.minrec,RECDSA.minrec,Options);

% Assess registration result
%implay(cat(3,REC2811.Ireg,RECDSA.minrec))

% Align all images of  the recording
REC2811_Reg = registerrecording(REC2811);
% The original struct is not needed anymore
clear REC2811

%% Session 3 - 1412
[REC1412] = preprocessRecording('.\Data\In vivo Ca2+\1412\',0,1);

% Define landmarks for image registration
if exist('Sampledata','var')
    [REC1412.Registration.movingPoints,REC1412.Registration.fixedPoints] = cpselect(REC1412.minrec,RECDSA.minrec,Sampledata.REC1412.Registration.movingPoints,Sampledata.REC1412.Registration.fixedPoints, 'Wait', true);
else
    [REC1412.Registration.movingPoints,REC1412.Registration.fixedPoints] = cpselect(REC1412.minrec,RECDSA.minrec, 'Wait', true);
end
REC1412.Registration.movingPoints=fliplr(REC1412.Registration.movingPoints);
REC1412.Registration.fixedPoints=fliplr(REC1412.Registration.fixedPoints);
REC1412.Registration.PStrength = ones(size(REC1412.Registration.fixedPoints,1),1);
REC1412.Registration.PStrength = REC1412.Registration.PStrength .* 0.95;
Options.Points1=REC1412.Registration.movingPoints;
Options.Points2=REC1412.Registration.fixedPoints;
Options.PStrength=REC1412.Registration.PStrength;
[REC1412.Ireg,REC1412.O_trans,REC1412.Spacing,~,~,~] = image_registration(REC1412.minrec,RECDSA.minrec,Options);

% Assess registration result
%implay(cat(3,REC1412.Ireg,RECDSA.minrec))

% Align all images of  the recording
REC1412_Reg = registerrecording(REC1412);
% The original struct is not needed anymore
clear REC1412
%% merge sessions
combinedcells = trackcells(0.6, REC2411_Reg, REC2811_Reg, REC1412_Reg);

%% Manually verify tracked cells across recordings
% show aligned cells across recordings
%allregisteredcells = cat(3,REC2411_Reg.allIC,REC2811_Reg.allIC,REC1412_Reg.allIC);
%implay(allregisteredcells)

% Generate a graph to interactively evaluate tracked cells across
% recordings
G = construct_graph(combinedcells, REC2411_Reg, REC2811_Reg,REC1412_Reg);

% Plot the interactive GUI for assessing automatic cell tracking across
% experiments
G = plotinteractivegraph(G);

% Label the same cells across recordings, that have not been identified 
% automatically, if the sample data is loaded
if exist('Sampledata','var')
    G = addedge(G,{'Session1-IC trace 28','Session3-IC trace 12','Session1-IC trace 6','Session2-IC trace 6','Session3-IC trace 42','Session2-IC trace 51'},{'Session3-IC trace 12','Session1-IC trace 28','Session2-IC trace 6','Session1-IC trace 6','Session2-IC trace 51','Session3-IC trace 42'}, 0);
end

% Calculate transitive closure to consistently match cells across all
% recordings
H = transclosure(G);
% remove self loops
H = rmedge(H, 1:numnodes(H), 1:numnodes(H));

% Visually validate tracked cells across recordings
H = plotinteractivegraph(H);

%% Calculate the datastructures of tracked cells across all recordings
Allcells = mergecells(H);

% We don't need the structs anymore. Delete them to save memory
clear G H combinedcells Options

% his movie helps to detect if all instances of a cell across recordings
% was correctly identified
AllmatchedcellsMovie = cumsum(Allcells.image,3);
%implay(AllmatchedcellsMovie)