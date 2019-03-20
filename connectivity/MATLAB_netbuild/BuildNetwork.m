close all;
clear all;
clc;

% This script will generate connectivity based on pre defined cell
% positions and connectivity rules.

% Set random seed
seed = 123421;
rng(seed);

% Load positions
% dirname = './position.dat';
% fid = fopen(dirname,'r');
% pos = textscan(fid,'%f%f%f%f%f','HeaderLines',1);
% pos = cell2mat(pos); % convert from cell to a matrix
% pos = sortrows(pos,1); % sort by cell ID

%%% TODO:
%%% 1) Reciprocal connectivity
%%% 2) Assign weights

%%%%%%%%%%%%%%%%%%%%%%%%% Network parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_cells = 10286; 
num_pyrA = 4115; 
num_pyrC = 4114; 
num_axo = 411; 
num_pv = 1646;

if (num_pyrA + num_pyrC + num_axo + num_pv ~= num_cells)
    'WARNING: cell numbers do not add up'
    break
end

% geometry
length = 1400;
width = 1400;
height = 200;

% Conduction velocity
cond_vel = 500; %um/ms

% Connectivity
PN2FSI = 0.12;
FSI2PN = 0.34;
PN2PN = 0;
FSI2FSI = 0.26;
PN2AXOrecip_perc = 0.4; %reciprocal

% Weights (mean)
PN2AXOmean = 0.00235; PN2AXOstd = 0.1*PN2AXOmean; PN2AXODIST='normal';
PN2PVmean = 0.002; PN2PVstd = 0.1*PN2PVmean; PN2PVDIST='normal';
PN2PNmean = 0; PN2PNstd = 0.1*PN2PNmean; PN2PNDIST='normal';

PV2PVmean = 0.004; PV2PVstd = 0.1*PV2PVmean; PV2PVDIST='normal';
PV2AXOmean = 0.008; PV2AXOstd = 0.1*PV2AXOmean; PV2AXODIST='normal';
PV2PNmean = 0.008; PV2PNstd = 0.1*PV2PNmean; PV2PNDIST='normal';

AXO2PVmean = 0.003; AXO2PVstd = 0.1*AXO2PVmean; AXO2PVDIST='normal';
AXO2AXOmean = 0.008; AXO2AXOstd = 0.1*AXO2AXOmean; AXO2AXODIST='normal';
AXO2PNmean = 0.003; AXO2PNstd = 0.2*AXO2PNmean; AXO2PNDIST='lognormal';


%%%%%%%%%%%%%%%%%%%%%%% Initialize matrices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% matrices for outputs (connections, weights, and delays)
connout = NaN(num_cells,1000);
weightout = NaN(num_cells,1000);
delayout = NaN(num_cells,1000);

PN2AXOrecip_conns = NaN(num_pyrA+num_pyrC,1000);

% indices of cells
pyrA = int32(floor([0,num_pyrA-1]));
pyrC = int32(floor([num_pyrA,num_pyrA+num_pyrC-1]));
axo = int32(floor([num_pyrA+num_pyrC, num_pyrA+num_pyrC+num_axo-1]));
bask = int32(floor([num_pyrA+num_pyrC+num_axo, num_cells-1]));

% Generate distributions for weights
if PN2AXOmean~=0 & PN2AXODIST=='normal'
    PN2AXOwgtpdf = PN2AXOstd.*randn(10000,1) + PN2AXOmean;
    % Get rid of samples outside a physiological range
    PN2AXOwgtpdf(PN2AXOwgtpdf<0 | PN2AXOwgtpdf>5*PN2AXOmean) = [];
end
if PN2PVmean~=0 & PN2PVDIST=='normal'
    PN2PVwgtpdf = PN2PVstd.*randn(10000,1) + PN2PVmean;
    PN2PVwgtpdf(PN2PVwgtpdf<0 | PN2PVwgtpdf>5*PN2PVmean) = [];
end
if PN2PNmean~=0 & PN2PNDIST=='normal'
    PN2PNwgtpdf = PN2PNstd.*randn(10000,1) + PN2PNmean;
    PN2PNwgtpdf(PN2PNwgtpdf<0 | PN2PNwgtpdf>5*PN2PNmean) = [];
end

if PV2AXOmean~=0 & PV2AXODIST=='normal'
    PV2AXOwgtpdf = PV2AXOstd.*randn(10000,1) + PV2AXOmean;
    % Get rid of samples outside a physiological range
    PV2AXOwgtpdf(PV2AXOwgtpdf<0 | PV2AXOwgtpdf>5*PV2AXOmean) = [];
end
if PV2PNmean~=0 & PV2PNDIST=='normal'
    PV2PNwgtpdf = PV2PNstd.*randn(10000,1) + PV2PNmean;
    PV2PNwgtpdf(PV2PNwgtpdf<0 | PV2PNwgtpdf>5*PV2PNmean) = [];
end
if PV2PVmean~=0 & PV2PVDIST=='normal'
    PV2PVwgtpdf = PV2PVstd.*randn(10000,1) + PV2PVmean;
    PV2PVwgtpdf(PV2PVwgtpdf<0 | PV2PVwgtpdf>5*PV2PVmean) = [];
end

if AXO2PNmean~=0 & AXO2PNDIST=='lognormal'
    mu = log((AXO2PNmean^2)/sqrt(AXO2PNstd^2+AXO2PNmean^2));
    sigma = sqrt(log(AXO2PNstd^2/(AXO2PNmean^2)+1));
    AXO2PNwgtpdf = lognrnd(mu,sigma,1,10000);
    % Get rid of samples outside a physiological range
    AXO2PNwgtpdf(AXO2PNwgtpdf<0 | AXO2PNwgtpdf>5*AXO2PNmean) = [];
end
if AXO2PVmean~=0 & AXO2PVDIST=='normal'
    AXO2PVwgtpdf = AXO2PVstd.*randn(10000,1) + AXO2PVmean;
    AXO2PVwgtpdf(AXO2PVwgtpdf<0 | AXO2PVwgtpdf>5*AXO2PVmean) = [];
end
if AXO2AXOmean~=0 & AXO2AXODIST=='normal'
    AXO2AXOwgtpdf = AXO2AXOstd.*randn(10000,1) + AXO2AXOmean;
    AXO2AXOwgtpdf(AXO2AXOwgtpdf<0 | AXO2AXOwgtpdf>5*AXO2AXOmean) = [];
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STEP 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Generate connectivity based on percentages in literature and store
%%%% them in CONNOUT. This is a pre-centric matrix (the rows are 
%%%% presynaptic cells). However, the network model needs the rows to be 
%%%% POSTsynaptic cells. This will be done in step 2.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate positions
pos = generate_positions(num_pyrA+num_pyrC,num_axo,num_pv,length,width,height,seed);

cand_pos = zeros(2,3);
for i=0:num_cells-1
    
    radius = 300; %um
    cell_pos = pos(i+1,2:4); %projecting cell's position
    cand_pos(1,:) = cell_pos + radius; % space around projecting cell
    cand_pos(2,:) = cell_pos - radius;
    
    x2 = pos(:,2); y2 = pos(:,3); z2 = pos(:,4);
    x1 = cell_pos(1); y1 = cell_pos(2); z1 = cell_pos(3);
    
    % possible connections
    possible_conns_ID = pos(sqrt((x2-x1).^2+(y2-y1).^2+(z2-z1).^2)<=radius);
    
    % cannot connect to itself
    possible_conns_ID(possible_conns_ID==i)=[];
    
    %scatter3(pos(:,2),pos(:,3),pos(:,4),'k');hold on;
    %scatter3(pos(possible_conns_ID+1,2),pos(possible_conns_ID+1,3),pos(possible_conns_ID+1,4),'g.');hold on;
    %scatter3(cell_pos(1),cell_pos(2),cell_pos(3),'r.')
    %xlim([0 1400]);ylim([0 1400]);zlim([0 1400]);
    
    % First, get all possible connections based on distance (above)
    % Then break those IDs into types
    possible_INT_ID = possible_conns_ID(possible_conns_ID>=axo(1) & possible_conns_ID<=num_cells-1);
    possible_PN_ID = possible_conns_ID(possible_conns_ID>=pyrA(1) & possible_conns_ID<=pyrC(2));
    
    % If presynaptic cell is pyramidal, do the following
    if i>=pyrA(1) && i<=pyrC(2) 
        % Find connections
        num_PN = ceil(PN2PN*size(possible_PN_ID,1));
        num_INTS = ceil(PN2FSI*size(possible_INT_ID,1));
        outgoingPNconns = sortrows(datasample(possible_PN_ID,num_PN,'Replace',false));
        outgoingINTconns = sortrows(datasample(possible_INT_ID,num_INTS,'Replace',false));
        outgoingAXOconns = outgoingINTconns(outgoingINTconns>axo(1) & outgoingINTconns<axo(2));
        outgoingPVconns = outgoingINTconns(outgoingINTconns>bask(1) & outgoingINTconns<bask(2));
        
        outgoingconns = cat(1,outgoingPNconns,outgoingINTconns); %concatenate
        connout(i+1,1) = num_INTS+num_PN; %number of connections in row
        connout(i+1,2) = i; %GID of presynaptic cell
        connout(i+1,3:size(outgoingconns,1)+2) = outgoingconns';
        
        % Some need to be reciprocal, store those here
%         PN2AXOrecip_num = ceil(PN2AXOrecip_perc*size(outgoingAXOconns,1));
%         if PN2AXOrecip_num>0
%             PN2AXOrecip_conns(i+1,1:PN2AXOrecip_num) = sortrows(datasample(outgoingAXOconns,PN2AXOrecip_num,'Replace',false))';
%         end
        
        % Now get axonal delay
        pos_outgoingconns = pos(ismember(pos(:,1),outgoingconns),:);
        dist = (pos_outgoingconns(:,2:4)-repmat(cell_pos,size(pos_outgoingconns,1),1)).^2;
        dist = sqrt(sum(dist,2)); % distance in microns
        delays = dist/cond_vel; % divide by conduction velocity
        delays(delays<0.1) = 0.1; % set minimum delay (should be > dt in NEURON)
        
        delayout(i+1,1) = num_INTS; %number of connections in row
        delayout(i+1,2) = i; %GID of presynaptic cell
        delayout(i+1,3:size(outgoingconns,1)+2) = delays';
        
        % Finally, assign weights
%         weights = zeros(size(outgoingconns,1),2); %initialize
%         weights(:,1) = outgoingconns;
%         
%         outgo_AXO = (weights(:,1)>=axo(1) & weights(:,1)<=axo(2));
%         N = size(weights(outgo_AXO,1),1);
%         
%         weights(outgo_AXO,2) = bigsample(1:N);
    end
    
    % If presynaptic cell is FSI, do this
    if i>=axo(1) && i<=bask(2) 
        %Find connections
        num_PN = ceil(FSI2PN*size(possible_PN_ID,1));
        num_INTS =  ceil(FSI2PN*size(possible_INT_ID,1));
        outgoingPNconns = sortrows(datasample(possible_PN_ID,num_PN,'Replace',false));
        outgoingINTconns = sortrows(datasample(possible_INT_ID,num_INTS,'Replace',false));
        outgoingAXOconns = outgoingINTconns(outgoingINTconns>axo(1) & outgoingINTconns<axo(2));
        outgoingPVconns = outgoingINTconns(outgoingINTconns>bask(1) & outgoingINTconns<bask(2));
        
        outgoingconns = cat(1,outgoingPNconns,outgoingINTconns); %concatenate
        connout(i+1,1) = num_INTS+num_PN; %number of connections in row
        connout(i+1,2) = i; %GID of presynaptic cell
        connout(i+1,3:size(outgoingconns,1)+2) = outgoingconns';
        
        % Now get axonal delay
        pos_outgoingconns = pos(ismember(pos(:,1),outgoingconns),:);
        dist = (pos_outgoingconns(:,2:4)-repmat(cell_pos,size(pos_outgoingconns,1),1)).^2;
        dist = sqrt(sum(dist,2)); % distance in microns
        delays = dist/cond_vel; % divide by conduction velocity
        delays(delays<0.1) = 0.1; % set minimum delay (should be > dt in NEURON)
        
        delayout(i+1,1) = num_INTS; %number of connections in row
        delayout(i+1,2) = i; %GID of presynaptic cell
        delayout(i+1,3:size(outgoingconns,1)+2) = delays';
    end
    
end

% create adjacency matrix and visualize connections of particular cell for
% verification
init = 9968; % presynaptic cell of interest (gid)
adj_effpyr = zeros(num_cells,num_cells);
for i=3:3+connout(init+1,1)-1
    k = connout(init+1,i)+1; % gid needs 1 added
    adj_effpyr(init+1,k) = 1;
end

gplot(adj_effpyr,pos(:,2:4),'b');
xlim([0 1400]);ylim([0 1400]);zlim([0 1400]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STEP 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Take the pre-centric matrix and convert it to a post-centric one.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% afferentconn = NaN(num_cells,1000);
% for i=0:num_cells-1
%     % Get all connections to postsynaptic cell (skip first two columns)
%     afferentconn_inds = find(connout(:,3:end)==i);
%     
%     % These indices are numbers from
%     % 1 to (size(connout,1)*size(connout,2)), numbered down the columns
%     afferentconn(i+1,1) = length(afferentconn_inds);
%     afferentconn(i+1,2) = i;
%     afferentconn(i+1,3:length(afferentconn_inds)+2) = mod(afferentconn_inds,size(connout,1))';
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STEP 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Write the output files.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Connectivity output file
% fileID = fopen('connectivity.dat','w');
% fprintf(fileID,'%d\n',num_cells); %first line should be number of cells
% 
% for i=0:num_cells-1
%     fprintf(fileID,'%d\t',afferentconn(i+1,~isnan(afferentconn(i+1,:))));
%     fprintf(fileID,'\n');
% end
% 
% fclose(fileID);

% Position output file
% fileID = fopen('positions.dat','w');
% fprintf(fileID,'%d\n',num_cells); %first line should be number of cells
% 
% for i=0:num_cells-1
%     fprintf(fileID,'%d\t',pos(i+1,:));
%     fprintf(fileID,'\n');
% end
% 
% fclose(fileID);
%% Trust but verify (connectivity)

% Pyramidal -> FSI




% Break out positions by cell type
% pyrA_pos = pos(find(pos(:,1)>=pyrA(1) & pos(:,1)<=pyrA(2)),1:4);
% pyrC_pos = pos(find(pos(:,1)>=pyrC(1) & pos(:,1)<=pyrC(2)),1:4);
% chn_pos = pos(find(pos(:,1)>=axo(1) & pos(:,1)<=axo(2)),1:4);
% bask_pos = pos(find(pos(:,1)>=bask(1) & pos(:,1)<=bask(2)),1:4);
% 
% figure
% hold all
% scatter3(pyrA_pos(:,2),pyrA_pos(:,3),pyrA_pos(:,4),'m.');
% scatter3(pyrC_pos(:,2),pyrC_pos(:,3),pyrC_pos(:,4),'r.');
% scatter3(chn_pos(:,2),chn_pos(:,3),chn_pos(:,4),'b.');
% scatter3(bask_pos(:,2),bask_pos(:,3),bask_pos(:,4),'g.');