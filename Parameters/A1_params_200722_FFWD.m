%%
% Parameter function required to run A1_freq_sweep_fn.m
% Assumes a small amount of DSI is inherited from Layer 4 by setting 
% param.ffwd_factor = 0.8;
%
% First argument is sweep speed (oct/ms)
% Second argument is one of the following:
%   'Default' (e.g., param = A1_params_Kato_200211(10*10^-3,'Default');
%   'SOM Block'
%   'PV Block'
%   'nonISN'
%   'SOM Projections'
%   'Threshold'
%
% See A1_run_all_sims.m for examples of each
%
% Parameters are optimized for speed and size of saved data, but may result
% in numerical artifacts.
%   For great precision uncomment the following lines below (adjustments
%   are also needed to the tolerance of ode45)
%       param.dx = 0.01;
%       param.tspan = [0:0.1:param.stim_stop+param.extra_time];
%
%
% This code was writen by Gregory Handy (2020)
% Please email ghandy@uchicago.edu with any questions
%%
function [ param ] = A1_params_200722_FFWD( varargin )

if nargin == 1
    error('Must have at least two inputs: Sweep speed and parameter set');
end

%% Stimulus parameters
if isnumeric(varargin{1})
    sweep_speed_up = varargin{1};   % units: oct/ms
    
    % Positive is a upward sweep, negative is a downward sweep
    param.sweep_speed_up = sweep_speed_up;
    param.sweep_speed_down = -sweep_speed_up;
else
    error('First input must be the numerical sweep speed (oct/ms)');
end
    
% starting frequencies for the up and downward directions
param.f_0_up = 4;
param.f_0_down = 64;

param.stim_delay = 500;
param.stim_stop = log2(param.f_0_down/param.f_0_up)*1/sweep_speed_up+param.stim_delay;
param.extra_time = 1000;
            
% peak amplitude of the response
param.b_amp = [10, 7];
param.ffwd_factor = 0.8;
% stimulus width
param.sigma = [0.5, 0.5];

%% Spatial parameters

% param.dx = 0.01; % use this for greater numerical precision
param.dx = 0.1;
param.x = [2:param.dx:6];
param.Nx = length(param.x);


%% Neuronal parameters
param.Npop = 3;

param.threshold = -1;
param.inhibitory_curr = [0, 0, 0];

% connectivity matrix
param.W = [1.2 -1.5 -1.5;
           1.2 -1.2 -1.2;
           2.5  0  0.0];

% projection widths
% Warning: numerical errors if lambda is taken to be too small
% (better approximation would be to take it to 0)
param.lambda = [0.5 0.5 2;
                0.5 0.5 2;
                2   0   0];
            
%% Normalization factor for connection strength 
% Normalize by the analytical solution, which is the integral of 
% exp(-(x-y)^2/lambda^2) 
% (i.e., the erf function with an extra scaling factor)
param.scaling_factor = zeros(3,3,length(param.x));
for i = 1:3
    for j = 1:3
        if param.lambda(i,j) ~=0
            param.scaling_factor(i,j,:) = 0.5*param.lambda(i,j)*sqrt(pi)*...
                (erf((-min(param.x)+param.x)/param.lambda(i,j))-...
                erf((-max(param.x)+param.x)/param.lambda(i,j)));
        end
    end
end

%% time constants for rate equations and EPSP/IPSP filters
param.tau_A = [10, 10, 10]; % time constant in ms
param.tau_S = [0, 0, 100]; % time constant in ms

%% ODE45 time parameters (in ms)

% param.tspan = [0:0.1:param.stim_stop+param.extra_time]; % use this for greater numerical precision
param.tspan = [0:1:param.stim_stop+param.extra_time];

%% Store the steady state

param.delta_r0 = zeros(param.Nx*param.Npop*2,1);

%%
param.color_scheme = [0.000 0.447 0.741
    0.850 0.325 0.098
    0.929 0.694 0.125];

%% Consider different datasets
if strcmp(varargin{2},'Default') ~= 1
    if nargin == 2
        error('For non-default parameter sets, user must also indicate additional parameter value');
    elseif strcmp(varargin{2},'SOM Block')
        
        if nargin ~= 4
            error(['SOM block requires two inputs: the amount to weaken the connection,' ...
                newline 'and the strength of inhibitory current']);
        else
            SOM_block_factor = varargin{3};
            SOM_inhibitory_curr = varargin{4};
            
            if isnumeric(SOM_block_factor) ~=1 || SOM_block_factor > 1 || SOM_block_factor < 0
                error('Third input should be a number between 0 and 1');
            end
            
            if isnumeric(SOM_inhibitory_curr) ~=1 || SOM_inhibitory_curr > 0 
                error('Fourth input should be negative');
            end
            
            param.inhibitory_curr(3) =  SOM_inhibitory_curr; % -0.15;
            param.W(3,1) = param.W(3,1)*SOM_block_factor;
            param.W(1,3) = param.W(1,3)*SOM_block_factor;
            
            % Quickly find the steady state using the non-spatial model
            [~,delta_r_quick] = ode15s(@(t,r)A1_quick_SS_ODE(t,r,param),[0 1000],[0 0 0]);
            param.delta_r0 = [delta_r_quick(end,1)*ones(param.Nx,1);
                delta_r_quick(end,1)*ones(param.Nx,1);
                delta_r_quick(end,2)*ones(param.Nx,1);
                delta_r_quick(end,2)*ones(param.Nx,1);
                delta_r_quick(end,3)*ones(param.Nx,1);
                delta_r_quick(end,3)*ones(param.Nx,1)];
        end
        
    elseif strcmp(varargin{2},'PV Block')
        
        if nargin ~= 4
            error(['PV block requires two inputs: the amount to weaken the connection,' ...
                newline 'and the strength of inhibitory current']);
        else
            
            PV_block_factor = varargin{3};
            PV_inhibitory_curr = varargin{4};
            
            if isnumeric(PV_block_factor) ~=1 || PV_block_factor > 1 || PV_block_factor < 0
                error('Third input should be a number between 0 and 1');
            end
            
            if isnumeric(PV_inhibitory_curr) ~=1 || PV_inhibitory_curr > 0 
                error('Fourth input should be negative');
            end
            
            param.inhibitory_curr(2) =  PV_inhibitory_curr; %-0.15;
            param.b_amp(2) = param.b_amp(2)*PV_block_factor;
            param.W(2,1) = param.W(2,1)*PV_block_factor;
            param.W(2,2) = param.W(2,2)*PV_block_factor;
            param.W(2,3) = param.W(2,3)*PV_block_factor;
            param.W(1,2) = param.W(1,2)*PV_block_factor;
            
            % Quickly find the steady state using the non-spatial model
            [~,delta_r_quick] = ode15s(@(t,r)A1_quick_SS_ODE(t,r,param),[0 1000],[0 0 0]);
            param.delta_r0 = [delta_r_quick(end,1)*ones(param.Nx,1);
                delta_r_quick(end,1)*ones(param.Nx,1);
                delta_r_quick(end,2)*ones(param.Nx,1);
                delta_r_quick(end,2)*ones(param.Nx,1);
                delta_r_quick(end,3)*ones(param.Nx,1);
                delta_r_quick(end,3)*ones(param.Nx,1)];
        end
        
    elseif strcmp(varargin{2},'nonISN')
        
        nonISN_factor = varargin{3};
        
        if isnumeric(nonISN_factor) ~=1 || nonISN_factor > 1 || nonISN_factor < 0
            error('Third input should be a number between 0 and 1');
        end
        
        param.W = param.W*nonISN_factor;
        param.b_amp = param.b_amp*nonISN_factor;
        param.threshold = -0.1;
        
    elseif strcmp(varargin{2},'SOM Projections')
        
        SOM_projection = varargin{3};
        if isnumeric(SOM_projection) ~=1 || SOM_projection > 1 || SOM_projection < 0
            error('Third input should be a number between 0 and 1');
        end
        
        param.lambda(1,3) = SOM_projection;
        param.lambda(2,3) = SOM_projection;
        param.lambda(3,1) = SOM_projection;
        
    elseif strcmp(varargin{2},'Threshold')
        
        new_threshold = varargin{3};
        if isnumeric(new_threshold) ~=1 || new_threshold>0
            error('Third input should be a number less than 0');
        end
        
        param.threshold = new_threshold;
    elseif strcmp(varargin{2},'Short Sweep')
        
        f_0_up = varargin{3};
        f_0_down = varargin{4};
        
        if isnumeric(f_0_up) ~=1 || f_0_up > 64 || f_0_up < 4
            error('Third input should be a number between 4 and 64');
        end
        
        if isnumeric(f_0_down) ~=1 || f_0_down > 64 || f_0_down < f_0_up
            error('Fourth input should less than 64 and greater than the second input');
        end
        
        param.f_0_up = f_0_up;
        param.f_0_down = f_0_down;
        
        param.stim_stop = log2(param.f_0_down/param.f_0_up)*1/sweep_speed_up+param.stim_delay;
    else
        error('Parameter set not reconginized');
    end
end

end