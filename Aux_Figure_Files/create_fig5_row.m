%%
% Create the rows that appear in Figure 5
%
% This code was writen by Gregory Handy (2020)
% Please email gregoryhandy@pitt.edu with any questions
%%
function [] = create_fig5_row(num_param_sets, DSI_default, DSI_param_set_1, DSI_param_set_2,...
    delta_r_up_param_set_2,delta_r_down_param_set_2, param, total_rows, ...
    row_num, fig_handle)

%%  Pick a BF to plot
x_interested = find(2.^param.x > 5,1);

%% Load the incoming excitatory input for the E population at BF specified

[~, ~, recurrent_exc_up, recurrent_exc_down] =...
    recurrent_input(param, x_interested, delta_r_up_param_set_2,...
    delta_r_down_param_set_2);

%% Load upward and downward feeforward input at specific BF

[~, ~, ffwd_total_up, ffwd_total_down] =...
    ffwd_input(param, x_interested);

%% Calculate the ratio of total/ffwd excitatory input

total_up = recurrent_exc_up+ffwd_total_up;
total_down = recurrent_exc_down+ffwd_total_down;

% normalize by the ffwd input
x_temp = [1, 2];
data_up = [ffwd_total_up,total_up]/ffwd_total_up;
data_down = [ffwd_total_down,total_down]/ffwd_total_up;

%% Plot the data
figure(fig_handle);
subplot(total_rows,3,1+3*(row_num-1))
plot(param.x,DSI_default,'linewidth',1.5,'color',param.color_scheme(1,:))
hold on
if num_param_sets == 3
    plot(param.x,DSI_param_set_1,'linewidth',1.5,'color',param.color_scheme(2,:))
    plot(param.x,DSI_param_set_2,'linewidth',1.5,'color',param.color_scheme(3,:))
else
    plot(param.x,DSI_param_set_2,'linewidth',1.5,'color',param.color_scheme(2,:))
end

ylim([-1 1])
xlim([2.27 5.73])
set(gca,'fontsize',16)
xticks([2 3 4 5 6])
xticklabels({'4','8','16','32','64'})
ylabel('DSI')
xlabel('BF')

subplot(total_rows,3,2+3*(row_num-1))
t_start = 450;
t_end = 850;
plot(param.tspan,delta_r_up_param_set_2(:,x_interested)-param.threshold,'linewidth',1.5,'color',param.color_scheme(num_param_sets,:))
hold on
plot(param.tspan,delta_r_down_param_set_2(:,x_interested)-param.threshold,'--','linewidth',1.5,'color',param.color_scheme(num_param_sets,:))
set(gca,'fontsize',16)
xlabel('Time (ms)')
xlim([t_start t_end])
ylim([-0.25 max(max(delta_r_up_param_set_2(:,x_interested), delta_r_down_param_set_2(:,x_interested))-param.threshold)])
ylabel('r (Hz)')

subplot(total_rows,3,3+3*(row_num-1))
hold off
plot(x_temp,data_up,'*-','markersize',16,'linewidth',1.5,'color',param.color_scheme(1,:))
hold on
plot(x_temp,data_down,'.-','markersize',16,'linewidth',1.5,'color',param.color_scheme(2,:))
set(gca,'fontsize',16)
xlim([0.8 2.2])
xticks([1 2])
ylabel('EPSC Ratio')
xticklabels({'Ffwd','Total'})
ylim([0.99 max([data_down, data_up])*1.01])


end
