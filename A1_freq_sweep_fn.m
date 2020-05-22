%%
% Wilson and Cowan, spatial rate model of the A1 tonotopic map using
% a threshold-linear approximation for the firing rate response function
%
% Runs both an upward and downward sweep
%
% Three populations: E, PV, and SOM (but can be generalized to more)
%
% Main options: Default, SOM_Block, PV_Block, nonISN, and SOM_Projections
% See A1_run_all_sims.m for examples of how to call these options
%
% Also, one can change the 'Default' parameter set by adjusting values in 
% A1_params_Kato_200211.m found in the Parameters folder
%
% This code was writen by Gregory Handy (2020)
% Please email gregoryhandy@pitt.edu with any questions
%%
function [ param, delta_r_up, delta_r_down, DSI_values] = A1_freq_sweep_fn(varargin)

restoredefaultpath;
folder = fileparts(which('A1_freq_sweep_fn.m')); 
addpath(genpath(folder));
rmpath(folder)

if nargin == 1
    error('Must have at least two inputs: Sweep speed and parameter set');
end

%% Store the wave speed
if isnumeric(varargin{1})
    sweep_speed_up = varargin{1};   % units: oct/ms
else
    error('First input must be the numerical sweep speed (oct/ms)');
end

%% Load the correct parameter set
param_set = varargin{2};
if strcmp(param_set,'Default') ~= 1
    if nargin == 2
        error('For non-default parameter sets, user must also indicate additional parameter value');
    elseif strcmp(param_set,'SOM Block')
        
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
            param = A1_params_Kato_200211(sweep_speed_up,param_set,SOM_block_factor,SOM_inhibitory_curr);
        end
    elseif strcmp(param_set,'PV Block')
        
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
            
            param = A1_params_Kato_200211(sweep_speed_up,param_set,PV_block_factor,PV_inhibitory_curr);
            
        end
        
    elseif strcmp(param_set,'nonISN')
        
        if nargin ~=3
            error('nonISN requires just one input: the amount to weaken the connectivity');
        else
            
            nonISN_factor = varargin{3};
            
            if isnumeric(nonISN_factor) ~=1 || nonISN_factor > 1 || nonISN_factor < 0
                error('Third input should be a number between 0 and 1');
            end
            
            param = A1_params_Kato_200211(sweep_speed_up,param_set,nonISN_factor);
        end
        
    elseif strcmp(param_set,'SOM Projections')
        
        if nargin ~=3
            error(['SOM Projections requires just one input: the spread of SOM incoming'...
                newline 'and outgoing connections']);
        else
            SOM_projection = varargin{3};
            if isnumeric(SOM_projection) ~=1 || SOM_projection > 1 || SOM_projection < 0
                error('Third input should be a number between 0 and 1');
            end
            
            param = A1_params_Kato_200211(sweep_speed_up,param_set,SOM_projection);
        end
        
    elseif strcmp(param_set,'Threshold')
        
        if nargin ~=3
            error(['Threshold requires just one input: the new threshold']);
        else
            new_threshold = varargin{3};
            if isnumeric(new_threshold) ~=1 || new_threshold>0
                error('Third input should be a number less than 0');
            end
            
            param = A1_params_Kato_200211(sweep_speed_up,param_set,new_threshold);
        end
        
    else
        error('Parameter set not reconginized');
    end
else
    param = A1_params_Kato_200211(sweep_speed_up,param_set);
end

%% Run the upward and downward sweeps
% Preallocation
delta_r = zeros(length(param.tspan), length(param.x)*3*2,2);

% Starting frequency depends on upward or downward direction
f_0 = [param.f_0_up, param.f_0_down];
sweep_speeds = [param.sweep_speed_up, param.sweep_speed_down];

% Tolerance level (can omit for computation speed, at the cost of
% numerically noisy DSI calculation
%options2 = odeset('RelTol', 1e-5, 'AbsTol', 1e-6);

% Uncomment these lines and change the for to a parfor to run in parallel
% delete(gcp('nocreate'))
% parpool(2);
for cc = 1:2
    [~,delta_r(:,:,cc)] = ode45(@(t,r)A1_spatial_sweep_6_ODE(t,r,param.x,...
        param.dx, param.Nx, param.Npop, param.b_amp, param.W, ...
        param.tau_A, param.sigma, param.lambda, param.scaling_factor, ...
        sweep_speeds(cc),f_0(cc), param.stim_stop,param.stim_delay,...
        param.tau_S,param.inhibitory_curr,param.threshold),...
        param.tspan,param.delta_r0); %#ok<*PFBNS>
end

%%
delta_r_up = delta_r(:,:,1);
delta_r_down = delta_r(:,:,2);

[DSI_values,~,~] = DSI_calc(param.tspan,param.x,delta_r_up,delta_r_down);

end

