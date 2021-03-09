%Linear Algebraic Nodal Analysis (LANA) Algorithm: Example 2
%
%	This implementation is essentially the Linear Algebraic Nodal Analysis
%	algorithm described in:
%  
%	Jeffrey A. Anderson, Electrify the linear-systems problem for our
%	students, PRIMUS, XX (2020), pp. nnn-mmm.
%
% For more information, please see the LANA Support Website at
%
% http://www.appliedlinearalgebra.com/blog/for-teachers/linear-algebra-laboratory-exercises/electrify-linear-systems
%
% This script is designed to demonstrate how to use the LANA algorithm 
% to model relatively complicated circuits using linear algebra. In this
% script file, we run through each of the 14 steps of the algorithm and
% store all relevant matrices, vectors, and data based on each step in our
% algorithm. For a breakdown of each step in the LANA algorithm, please see
% the website listed above and look for the "Step-by-step guide to the LANA
% algorithm (.pdf) document."
%
% Outputs:
%	output1 - Description
%	output2 - Description
%
% Example: Run the following line of code in MATLAB's Command Window: 
%   LANA_Algorithm_Example_2
%
% Other m-files required: none
% MAT-files required: none
%
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2

% Author: Jeff Anderson, PhD
% Work address: 
%     Foothill College
%     12345 El Monte Road
%     Los Altos Hills, CA 94022
% email: andersonjeff@fhda.edu
% Website: http://www.appliedlinearalgebra.com/
% YouTube Channel Name: appliedlinearalgebra.com
% YouTube URL: https://www.youtube.com/channel/UCPGNxQ2p2eNrmaIsU4zuf2Q
% Created on: Monday 04/20/2020
% Last revision: Monday 04/20/2020
% Version: 1.0

%------------- BEGIN CODE --------------

%STEP 1: IDENTIFY AND LABEL THE ENTIRE SET OF NODES IN CIRCUIT

%%STEP 2: MODEL THE CIRCUIT AS A DIRECTED GRAPH
  %
  %STEP 2A: TRACK THE DIMENSIONS OF KEY FEATURES IN CIRCUIT
  %
  %Identify column dimension of incidence matrix (Nodes)
  ng = 7;   %Total number of nodes in circuit (including ground node)

  %Identify row dimension of incidence matrix (Edges/Elements)
  mr = 8;           %Number of resistors in circuit
  mv = 2;           %Number of dc voltage sources in circuit
  mi = 2;           %Number of dc current sources in circuit
  m = mr + mv + mi; %Total number of elements (edges) in circuit

%  
%%STEP 3: CREATE THE ENTIRE INCIDENCE MATRIX 
%

  %Create the resistor subblock of the entire incidence matrix 
  Arg = [ 1     0     0      0     0     0     -1;
         -1     1     0      0     0     0      0;
          0     1     0      0     0     0     -1;
          0     1    -1      0     0     0      0;
          0     0     1      0     0    -1      0;
          0     0     0      0     0     1     -1;
          0     0     1     -1     0     0      0;
          0     0     0      0     1    -1      0];

  %Create the voltage source subblock of the entire incidence matrix 
  Avg = [ 0     1     0      0     0     0     -1;
          0     0     0      1    -1     0      0];

  %Create the current source subblock of the entire incidence matrix 
  Aig = [ 0    -1     0      1     0     0      0;
          0     0     0      0    -1     0      1];

  %Create the entire incidence matrix using pre-defined subblocks
  Ag = [Arg; Avg; Aig];


%  
%STEP 5: STATE THE ENTIRE SET OF CIRCUIT EQUATIONS 
%
  %  
  %STEP 5A: SET UP RIGHT-HAND SIDE VECTORS FOR KVLS IN POTENTIAL FORM
  %
 
  %Vv = [4.97 ; 4.97 ];  %Measured values provided by voltage sources
  Vv = [5; 5]          %Known ideal values provided by voltage sources
  
  Ii = [2.5; 2.5];      %Known ideal values provided by voltage sources

  %  
  %STEP 5B: SET UP MATRICES FOR BRANCH CONSTITUENT RELATIONS
  %
  
  %Resistor matrix: diagonal matrix encoding resistor values
  % In our electronics learning laboratory kit, we build circuits using
  % exclusively 1k resistors. Moreover, we design all equipement to work in
  % units of volts (V), miliamps (mA), and kilohms (k) leading to unit eq
  %
  %       1V = 1mA * 1k
  %
  % Thus, all resistance values are written in terms of kilohms so that a
  % value of r1 = 1 is in kilohms. For users who switch to resistors with
  % other values, please note this design feature and change the entries of
  % the matrix below appropriately. 
  R = diag([1, 1, 1, 1, 1, 1, 1, 1]);

  %Conductance matrix: inverse of resistor matrix
  G = inv(R);

%
%STEP 7: GROUND THE CIRCUIT
%
  %
  %STEP 7A: CHOOSE A GROUND NODE FOR THE CIRCUIT
  %
  D0 = eye(ng, ng);
  D0 = D0(:,[1, 2, 3, 4, 5, 6]);

  %
  %STEP 7B: GROUND THE CIRCUIT EQUATIONS USING DEFLATION
  %
  Ar0 = Arg*D0;
  Av0 = Avg*D0;
  Ai0 = Aig*D0;
  A0 = Ag*D0;
  
%
%STEP 10: ELIMINATE THE NODE DEPENDENCIES THAT ARISE DUE TO VOLTAGE SOURCES
%

  %Create a particular solution to GLSP
  pg = [0; Vv(1,1); 0; 0; -Vv(2,1); 0; 0];

  %Create the deflation matrix for the entire set of voltage source KVLs
  Dvg  =    [  1,    0,    0,    0;
               0,    0,    0,    0; 
               0,    1,    0,    0;
               0,    0,    1,    0;
               0,    0,    1,    0;
               0,    0,    0,    1;
               0,    0,    0,    0];

  D = D0'*Dvg;
  Ar = Ar0*D;
  Ai = Ai0*D;
  x0  = D0'*pg;
  b  = -Ar0*x0;


K = Ar'*G*Ar;
f = Ar'*G*b - Ai'*Ii;

u = K\f;
u0 = x0 + D*u

% Kt = [    Ar0'*G*Ar0,    Av0';
%             Av0,        zeros(2,2)];
%         
% ft = [- Ai0'*Ii; Vv];
% Kt\ft
