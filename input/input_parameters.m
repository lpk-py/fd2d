%==========================================================================
% path where seismic sources are located
%==========================================================================

source_path='../../input/sources/forward/';

%==========================================================================
% set basic simulation parameters
%==========================================================================

Lx=200.0;   % model extension in x-direction [m]
Lz=250.0;   % model extension in y-direction [m]

nx=200;     % grid points in x-direction
nz=250;     % grid points in z-direction

dt=0.0001;  % time step [s]
nt=650;     % number of iterations

order=4;    % finite-difference order (2 or 4)

%==========================================================================
% model type
%==========================================================================

model_type=1;

% 1=homogeneous 
% 2=homogeneous with localised density perturbation
% 3=layered medium
% 4=layered with localised density perturbation
% 5=vertical gradient medium
% 6=vertical gradient medium with localised density perturbation

% 'initial'= read initial model for waveform inversion (mu_initial, rho_initial)

%==========================================================================
% source-time function
%==========================================================================

f_min=10.0;     % minimum frequency [Hz]
f_max=5000.0;    % maximum frequency [Hz]

%==========================================================================
% simulation mode
%==========================================================================

simulation_mode='forward_correlation';

% 'forward'                 regular forward simulation
% 'forward_correlation'     forward simulation where Fourier transform is computed on-the-fly (preparation to compute correlation function)
% 'correlation'             compute correlation functions as output, requires Fourier-transformed Green function to be present

%==========================================================================
% receiver positions
%==========================================================================

%- a large number of receivers in an open rectangular configuration
%rec_x=[50.0 50.0 50.0 50.0 50.0 50.0 50.0 50.0 50.0  50.0  50.0  50.0  60.0 70.0 80.0 90.0 100.0 110.0 120.0 130.0 60.0  70.0  80.0  90.0  100.0 110.0 120.0 130.0];
%rec_z=[70.0 80.0 90.0 100.0 110.0 120.0 130.0 140.0 150.0 160.0 170.0 180.0 70.0 70.0 70.0 70.0 70.0  70.0  70.0  70.0  180.0 180.0 180.0 180.0 180.0 180.0 180.0 180.0];

%- just one receiver
%rec_x=[30.0];
%rec_z=[125.0];

%- a large number of receivers in a closed rectangular configuration
rec_x=[50.0  50.0  50.0  50.0  50.0   50.0    70.0  90.0 110.0 130.0   70.0  90.0 110.0 130.0  150.0 150.0 150.0 150.0 150.0  150.0];
rec_z=[70.0  90.0 110.0 130.0 150.0  170.0    70.0  70.0  70.0  70.0  170.0 170.0 170.0 170.0   70.0  90.0 110.0 130.0 150.0  170.0];

%==========================================================================
% absorbing boundaries
%==========================================================================

width=30.0;     % width of the boundary layer in km

absorb_left=1;  % absorb waves on the left boundary
absorb_right=1; % absorb waves on the right boundary
absorb_top=1;   % absorb waves on the top boundary
absorb_bottom=1;% absorb waves on the bottom boundary

%==========================================================================
% make wavepropagation movie
%==========================================================================

make_movie='yes';                   % 'yes' or 'no'
movie_file='../../output/wavemovie.mp4';   % output file name 