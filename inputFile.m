%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%											
%   "Copyright (c) 2004 The University of Southern California"				
%   All rights reserved.								
%											
%   Permission to use, copy, modify, and distribute this software and its		
%   documentation for any purpose, without fee, and without written agreement is	
%   hereby granted, provided that the above copyright notice, the following		
%   two paragraphs and the author appear in all copies of this software.		
%											
%   NO REPRESENTATIONS ARE MADE ABOUT THE SUITABILITY OF THE SOFTWARE FOR ANY		
%   PURPOSE. IT IS PROVIDED "AS IS" WITHOUT EXPRESS OR IMPLIED WARRANTY.		
%											
%   Neither the software developers, the Autonomous Network Research Group		
%   (ANRG), or USC, shall be liable for any damages suffered from using this		
%   software.										
%											
%   Author:		Marco Zuniga 
%   Director: Prof. Bhaskar Krishnamachari
%   Autonomous Networks Research Group, University of Southern California
%   http://ceng.usc.edu/~anrg/
%   Contact: marcozun@usc.edu
%
%   Date last modified:	2004/06/23 marcozun						
%											
%   Anything following a "%" is treated as a comment.					
%											
%											
%   Description:									
%	This file contains the required parameters to generate an			
%	instance of the Link Layer Model						
%											
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [nameFile] = inputFile(topologyFile); 

%%%%%%%%%%%%%%%%%%%%%%
%		     
% Channel Parameters 
%		     
%%%%%%%%%%%%%%%%%%%%%%
%
% The channel model is based on the Log-Normal Shadowing Path Loss Model
%  PATH_LOSS_EXPONENT		    is an adimensional constant
%  SHADOWING_STANDARD_DEVIATION	is in dB
%  PL_D0			            is in dB, is the close-in reference pathloss
%  D0				            is in m, is the close-in reference distance
% all values should be positive

PATH_LOSS_EXPONENT = 4.7;
SHADOWING_STANDARD_DEVIATION = 3.2;
PL_D0 = 55.0;
D0 = 1.0;

%%%%%%%%%%%%%%%%%%%%
%		   
% Radio Parameters 
%		   
%%%%%%%%%%%%%%%%%%%%
%
% The radio model is based on probability-of-error expressions for
% channels considering only white gaussian noise, no dynamic effects
% are included, i.e. node mobility or highly-dynamic environments. 
%

% Modulation options
% ------------------
% NCASK	1 (Non Coherent Amplitude Shift Keying)
% ASK	2 (Amplitude Shift Keying)	
% NCFSK	3 (Non Coherent Frequency Shift Keying)
% FSK	4 (Frequency Shift Keying)
% BPSK	5 (Binary Phase Shift Keying)
% DPSK	6 (Differential Phase Shift Keying)

MODULATION = 3;

% Encoding options
% ----------------
% NRZ		    1 (No Return to Zero)
% 4B5B		    2 (4-Bit Data Symbol, 5-Bit Code)
% MANCHESTER    3 (Manchester)
% SECDED	    4 (Single Error Detection Double Error Correction)	

ENCODING = 3;

% Radio Output Power in dBm

OUTPUT_POWER = -7.0;

% Noise Floor in dBm

NOISE_FLOOR = -105.0;

% Receiver-Transmitter Correlation, Asymmetric links
%   First row of COV should have noise floor (rx) stats
%   Second row should have output power (tx) stats
%   COV = [Srx Srx-tx; Stx-rx Stx]
% The values presented below are for MICA2 radios.

COVM = [3.7 -3.3; -3.3 6.0];

% Asymmetry = 0 does not consider COV, i.e. symmetric links
% Asymmetry = 1 considers asymmetric links

ASYMMETRY = 1;

% Preamble and Frame Size in bytes

PREAMBLE_LENGTH = 2;
FRAME_LENGTH = 50;

%%%%%%%%%%%%%%%%%%%%%%%
%		      
% Topology Parameters 
%		      
%%%%%%%%%%%%%%%%%%%%%%%
%
% Terrain dimensions and nodes' position are in meters
%

% Number of nodes

NUMBER_OF_NODES = 400;
% Physical terrain (meters)
% the density (NUMBER_OF_NODES / area) can not be higher than
% 0.5 nodes / m^2.

TERRAIN_DIMENSIONS_X = 25.0;
TERRAIN_DIMENSIONS_Y = 25.0;

% Topology Options
% ----------------
%
% GRID    1
% UNIFORM 2
% RANDOM  3
% FILE	  4
%
% Choose the desired topology and place the corresponding number in TOPOLOGY (below)
% GRID: Node placement starts at (0, 0). ONLY for this option the GRID_UNIT variable
%  is required (meters). The number of nodes has to be square of an integer.
% UNIFORM: Based on the number of nodes, the physical terrain is divided 
%  into a number of cells. Within each cell, a node is placed randomly.
% RANDOM: Nodes are placed randomly within the physical terrain.
% FILE: Position of nodes is read as a matrix argument (read instructions for
% LinkLayerModel function)
%  The format of the matrix should be:
%
%  topology = [ 
%  Xcoordinate Ycoordinate
%  ];
%

TOPOLOGY = 2;
GRID_UNIT = 2.0;

%%%%%%%%%%%%%%%
%	      
% Output Data 
%	      
%%%%%%%%%%%%%%%
%
% Two formats are avaible a regular one, with parenthesis and colons. And another
% which provides a MATLAB format, for the MATLAB option set the MATLAB_FORMAT 
% variable to 1
% 

save nameFile;
nameFile = 'nameFile';