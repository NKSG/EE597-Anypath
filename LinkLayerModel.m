%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%											
%   "Copyright (c) 2004, 2005 The University of Southern California"				
%   All rights reserved.								
%											
%   Permission to use, copy, modify, and distribute this script and its		
%   documentation for any purpose, without fee, and without written agreement is	
%   hereby granted, provided that the above copyright notice, the following		
%   two paragraphs and the names in the credits appear in all copies of this software.		
%											
%   NO REPRESENTATIONS ARE MADE ABOUT THE SUITABILITY OF THE SCRIPT FOR ANY		
%   PURPOSE. IT IS PROVIDED "AS IS" WITHOUT EXPRESS OR IMPLIED WARRANTY. NO 
%   LIABILITY IS ASSUMED BY THE DEVELOPERS.
%											
%											
%   Authors:		Marco Zuniga, Rahul Urgaonkar
%   Director:       Prof. Bhaskar Krishnamachari
%   Autonomous Networks Research Group, University of Southern California
%   http://ceng.usc.edu/~anrg
%
%   Contact: marcozun@usc.edu
%   Previous Version: 1.0, 2004/07/02
%   Current Version 1.1
%   Date last modified: 2005/12/20							
%											
%   Anything following a "%" is treated as a comment.					
%											
%											
%   Description:									
%	LINKLAYERMODEL  generates an instance of the connectivity graph of a wireless 					
%                   sensor network  
%       All input parameters (except one, explained below) should be entered in inputFile.m.
%       Example
%       >> [ topology, prrM] = linklayermodel;
%           topology    provides the x and y coordinates in a Nx2 matrix where
%                       N is the number of nodes, first column provides X
%                       coordinates and second the Y coordinates.
%           prrM        Is the link quality matrix, quality between 0 and 1
%
%       Different topologies (Grid, Uniform, Random, File )can be generated in inputFile.m
%           if an specific topology is desired, it has to be given as an
%           input argument (Nx2 matrix). Example
%       >> [ topology, prrM] = linklayermodel(desiredTopology);
%											
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ topology, prrM] = LinkLayerModel(varargin)

numInputs = length(varargin);
%if (numInputs > 2)
%	error('Error: only input argument can be topology');
%end

% load input file
[nameInFile] = inputFile;
load(nameInFile);
NUMBER_OF_NODES=varargin{2};
TOPOLOGY=varargin{1};
PATH_LOSS_EXPONENT=varargin{3};
SHADOWING_STANDARD_DEVIATION=varargin{4};
MODULATION=varargin{5};
ENCODING=varargin{6};

if ( (TOPOLOGY == 4) && (numInputs == 0) )
	error('Error: enter topology as an argument');
elseif ( (TOPOLOGY == 4) && (numInputs == 1) )
    [inVar] = deal(varargin{1});
    topology = inVar;
end

% some security checks
if(PATH_LOSS_EXPONENT < 0)
	error('Error: value of PATH_LOSS_EXPONENT must be positive');
end
if(SHADOWING_STANDARD_DEVIATION < 0)
	error('Error: value of SHADOWING_STANDARD_DEVIATION must be positive');
end
if(PL_D0 < 0)
	error('Error: value of PL_D0 must be positive');
end
if(D0 < 0)
	error('Error: value of D0 must be positive');
end
if(PREAMBLE_LENGTH < 0)
	error('Error: value of PREAMBLE_LENGTH must be positive');
end
if(FRAME_LENGTH < 0)
	error('Error: value of FRAME_LENGTH must be positive');
end
if(NUMBER_OF_NODES < 0)
	error('Error: value of NUMBER_OF_NODES must be positive');
end
if(TERRAIN_DIMENSIONS_X*TERRAIN_DIMENSIONS_Y < 0)
	error('Error: value of dimensions must be positive');
end
if( (COVM(1,2) ~= COVM(2,1)) || ( max(COVM(1,1), COVM(2,2)) < COVM(1,2) ) )
  	error('Error: COV must be positive-definite');
end
    
area = TERRAIN_DIMENSIONS_X*TERRAIN_DIMENSIONS_Y;

%%%%%%%%%%%%%%%%%%%%%
% Create Topology
%%%%%%%%%%%%%%%%%%%%%

if      (TOPOLOGY == 1) % GRID
    
    if (GRID_UNIT < D0)
        error('value of GRID_UNIT must be greater than D0');
    end
    
    if (sqrt(NUMBER_OF_NODES) ~= round(sqrt(NUMBER_OF_NODES)))
        error('Number of nodes should be a perfect square');
    end
    
    for i=1:NUMBER_OF_NODES
        % X coordinate
        topology(i,1) = GRID_UNIT*rem(i-1, sqrt(NUMBER_OF_NODES));
        % Y coordinate
        topology(i,2) = GRID_UNIT*floor((i-1)/(sqrt(NUMBER_OF_NODES)));
    end
    
elseif  (TOPOLOGY == 2) % UNIFORM
    
    if (( TERRAIN_DIMENSIONS_X < 0) | (TERRAIN_DIMENSIONS_Y < 0))
        error('Terrain dimensions must be positive');
    end
    
    cellLength = sqrt(area/NUMBER_OF_NODES);
    nodesX = ceil(TERRAIN_DIMENSIONS_X/cellLength);
    cellLength = TERRAIN_DIMENSIONS_X/nodesX;
    % 1.4 (below) is an arbitrary number chosen to decrease the probability
    % that nodes get closer than D0 to one another.
    if (cellLength < D0*1.4)
        error('UNIFORM topology: density too high, increase area');
    end
    for i=1:NUMBER_OF_NODES
        topology(i,1) = ((rem(i-1, nodesX))*cellLength) + rand(1)*cellLength;
        topology(i,2) = ((floor((i-1)/nodesX))*cellLength) + rand(1)*cellLength;
        % rest of for loop checks that no internode distance is smaller than d0
        wrongPlacement = 1;
        while (wrongPlacement == 1)
            for j = 1:i
                xdist = topology(i,1) - topology(j,1);
                ydist = topology(i,2) - topology(j,2);
                dist = sqrt(xdist^2 + ydist^2);
                if ((dist < D0) & (i ~= j))  
                    topology(i,1) = ((rem(i-1, nodesX))*cellLength) + rand(1)*cellLength;
                    topology(i,2) = ((floor((i-1)/nodesX))*cellLength) + rand(1)*cellLength;
                    wrongPlacement = 1;
                    break;
                end
            end       
            if (j==i)
                wrongPlacement = 0;
            end         
        end
    end

elseif  (TOPOLOGY == 3) % RANDOM
    
    if (( TERRAIN_DIMENSIONS_X < 0) | (TERRAIN_DIMENSIONS_Y < 0))
        error('Terrain dimensions must be positive');
    end
    cellLength = sqrt(area/NUMBER_OF_NODES);
    % 1.4 (below) is an arbitrary number chosen to decrease the probability
    % that nodes get closer than D0 to one another.
    if (cellLength < D0*1.4)
        error('RANDOM topology: density too high, increase area');
    end
    for i=1:NUMBER_OF_NODES
        topology(i,1) = TERRAIN_DIMENSIONS_X*rand(1);
        topology(i,2) = TERRAIN_DIMENSIONS_Y*rand(1);
        % rest of for loop checks that no internode distance is smaller than d0
        wrongPlacement = 1;
        while (wrongPlacement == 1)
            for j=1:i
                xdist = topology(i,1) - topology(j,1);
                ydist = topology(i,2) - topology(j,2);
                dist = sqrt(xdist^2 + ydist^2);
               	if ((dist < D0) & (i ~= j))
               		topology(i,1) = TERRAIN_DIMENSIONS_X*rand(1);
               		topology(i,2) = TERRAIN_DIMENSIONS_Y*rand(1);
               		wrongPlacement = 1;
               		break;
                end
            end 
            if (j==i)
         	    wrongPlacement = 0;
            end
        end
    end
    
elseif(TOPOLOGY==4)

% The input matrix topology has dimensions NUMBER_OF_NODESx2 
% where
%       topology(i,1) denotes the x coordinates of node i
%       topology(i,2) denotes the y coordinate of node i

    [p q] = size(topology);
    if (p ~= NUMBER_OF_NODES)
        error('Number of nodes in file topology.m does not agree with NUMBER_OF_NODES');
    end
    if (q ~= 2)
        error('Wrong format of file topology.m');
    end
    
    for i=1:NUMBER_OF_NODES
        for j=1:i
            xdist = topology(i,1) - topology(j,1);
            ydist = topology(i,2) - topology(j,2);
            dist = sqrt(xdist^2 + ydist^2);
            if ( (dist < D0) & (i ~= j) )
                error('Error: topology.m contains internode distances less than D0');
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine Output Power and Noise Floor
%       both follow a gaussian distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:NUMBER_OF_NODES
    if  (ASYMMETRY == 1)
        % Cholesky Decomposition is used to generate multivariate random
        % variables:
        %   covariance matrix COVM = T' x T
        %   T is a 2x2 upper triangular 
        %       P_T   = P_T + T(1,1)* rn1
        %       P_N   = P_N + T(1,2) * rn1 + T(2,2) * rn2
        %  where rn1 and rn2 are normal(0,1) random variables.
        
        
        T11 = sqrt(COVM(1,1));
        T12 = COVM(1,2)/sqrt(COVM(1,1));
        T21 = 0;
        T22 = sqrt( (COVM(1,1)*COVM(2,2)-COVM(1,2)^2) / COVM(1,1) );
       
        rn1 = randn;
        rn2 = randn;
       
        noiseFloor(i)  = NOISE_FLOOR + T11 * rn1;
        outputPower(i) = OUTPUT_POWER + T12 * rn1 + T22 * rn2;
    else
        
        noiseFloor(i)  = NOISE_FLOOR;
        outputPower(i) = OUTPUT_POWER;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Obtain RSSI
%       use topology information and
%       channel parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:NUMBER_OF_NODES
    for j=i+1:NUMBER_OF_NODES
        xdist = topology(i,1) - topology(j,1);
        ydist = topology(i,2) - topology(j,2);
        dist = sqrt(xdist^2 + ydist^2);
        % mean rssi decay dependent on distance 
        pathLoss = - PL_D0 - (10*PATH_LOSS_EXPONENT*(log(dist/D0)/log(10))) + randn * SHADOWING_STANDARD_DEVIATION;
        % assymetric links are given by running two different
		% R.V.s for each unidirectional link.
		%   NOTE: this approach is not accurate, assymetry is due mainly to
		%   to hardware imperfections and not for assymetric paths
        rssi(i,j) = outputPower(i) + pathLoss;
        rssi(j,i) = outputPower(j) + pathLoss;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Obtain Prob. of bit Error
%       use rssi and modulation chosen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:NUMBER_OF_NODES
    for j=1:NUMBER_OF_NODES
        if (i==j)
            pe(i,j) = 0;
        else
            snr = ( 10^((rssi(i,j) - noiseFloor(j))/10) ) / .64;  % division by .64 converts from Eb/No to RSSI
                                                                  % this is specific for each radio (read paper: Data-rate(R) / Bandwidth-noise(B)) 
            if (MODULATION == 1)    % NCASK
                pe(i,j) = 0.5*( exp(-0.5*snr) + Q( sqrt(snr) ) );
            elseif(MODULATION == 2) % ASK
                pe(i,j) = Q( sqrt(snr/2) );
            elseif(MODULATION == 3) % NCFSK
                pe(i,j) = 0.5*exp(-0.5*snr);
            elseif(MODULATION == 4) % FSK
                pe(i,j) = Q( sqrt(snr) );
            elseif(MODULATION == 5) % BPSK
                pe(i,j) = Q( sqrt(2*snr) );
            elseif(MODULATION == 6) % DPSK
                pe(i,j) = 0.5*exp(-snr);
            else
                error('MODULATION is not correct');
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Obtain PRR
%   use prob. of error and encoding scheme
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:NUMBER_OF_NODES
    for j = 1:NUMBER_OF_NODES
        if (i == j)
            prrM(i,j) = 1;
        else
            preseq = (1-pe(i,j))^(8*PREAMBLE_LENGTH);
            if (ENCODING == 1)      % NRZ
                prrM(i,j) = preseq*((1-pe(i,j))^(8*(FRAME_LENGTH-PREAMBLE_LENGTH)));
            elseif (ENCODING == 2)  % 4B5B
                prrM(i,j) = preseq*((1-pe(i,j))^(8*1.25*(FRAME_LENGTH-PREAMBLE_LENGTH)));
            elseif (ENCODING == 3)  % MANCHESTER
                prrM(i,j) = preseq*((1-pe(i,j))^(8*2*(FRAME_LENGTH-PREAMBLE_LENGTH)));
                %fprintf('x=%d,y=%d,prrm=%f\n',i,j,prrM(i,j));
            elseif (ENCODING == 4)  % SECDED
                prrM(i,j) = ((preseq*((1-pe(i,j))^8)) + (8*pe(i,j)*((1-pe(i,j))^7)))^((FRAME_LENGTH-PREAMBLE_LENGTH)*3);
            else
                error('ENCODING is not correct');
            end
        end
    end
end
