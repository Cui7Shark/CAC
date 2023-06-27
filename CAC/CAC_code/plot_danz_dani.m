clc;close all;clear;
load('data_di_ro.mat');
load('data_di.mat');
%%
%%时间
T_ro = t_ro * 1e9;
T_di = t_di * 1e9;
%角度
A30_di = 61;
A90_di = 181;
A60_di = 121;
A120_di = 241;
A240_di = 481;
A300_di = 601;

A30_ro = 60;
A60_ro = 120;
A90_ro = 180;
A120_ro = 240;
A240_ro = 480;
A300_ro = 600;
%时间
T_ro_17ns = 112;
T_ro_50ns = 301;
T_ro_100ns = 369;

T_di_17ns = 119;
T_di_50ns = 306;
T_di_100ns = 360;
%%
%%Di-旋转
[u_ro_pro_V17,~] = smoothdata(u_ro(T_ro_17ns,k1_ro(:,2)) - u_ro(T_ro_17ns,k1_ro(:,1)) , 'sgolay', 16);
[u_ro_pro_V50, ~] = smoothdata(u_ro(T_ro_50ns,k1_ro(:,2)) - u_ro(T_ro_50ns,k1_ro(:,1)) , 'sgolay', 12);
[u_ro_pro_V100, ~] = smoothdata(u_ro(T_ro_100ns,k1_ro(:,2)) - u_ro(T_ro_100ns,k1_ro(:,1)) , 'sgolay', 12);

[u_ro_pro_N17,~] = smoothdata(u_ro(T_ro_17ns,a_ro+1:a_ro+a1_ro), 'gaussian', 18);
[u_ro_pro_N50, ~] = smoothdata(u_ro(T_ro_50ns,a_ro+1:a_ro+a1_ro), 'gaussian', 12);
[u_ro_pro_N100, ~] = smoothdata(u_ro(T_ro_100ns,a_ro+1:a_ro+a1_ro), 'gaussian', 12);
%%%%%%%%%%%%%%%%%%%
%%%%%  Di  %%%%%%%%
%%%%%%%%%%%%%%%%%%%
%%17ns TMP di 数据平滑处理
tmp_60_di = u_di(T_di_17ns,k1_di(A60_di,2)) - u_di(T_di_17ns,k1_di(A60_di,1));
tmp_120_di = u_di(T_di_17ns,k1_di(A120_di,2)) - u_di(T_di_17ns,k1_di(A120_di,1));
tmp_240_di = u_di(T_di_17ns,k1_di(A240_di,2)) - u_di(T_di_17ns,k1_di(A240_di,1));
tmp_300_di = u_di(T_di_17ns,k1_di(A300_di,2)) - u_di(T_di_17ns,k1_di(A300_di,1));

[u_di_pro_V17,~] = smoothdata(u_di(T_di_17ns,k1_di(:,2)) - u_di(T_di_17ns,k1_di(:,1)) , 'sgolay', 16);

u_di_pro_V17(1,A60_di) = tmp_60_di;
u_di_pro_V17(1,A120_di) = tmp_120_di;
u_di_pro_V17(1,A240_di) = tmp_240_di;
u_di_pro_V17(1,A300_di) = tmp_300_di;
%%17ns N
[u_di_pro_N17,~] = smoothdata(u_di(T_di_17ns,a_di+1:a_di+a1_di), 'gaussian', 18);
%%50ns TMP
tmp_60_di = u_di(T_di_50ns,k1_di(A60_di,2)) - u_di(T_di_50ns,k1_di(A60_di,1));
tmp_120_di = u_di(T_di_50ns,k1_di(A120_di,2)) - u_di(T_di_50ns,k1_di(A120_di,1));
tmp_240_di = u_di(T_di_50ns,k1_di(A240_di,2)) - u_di(T_di_50ns,k1_di(A240_di,1));
tmp_300_di = u_di(T_di_50ns,k1_di(A300_di,2)) - u_di(T_di_50ns,k1_di(A300_di,1));

[u_di_pro_V50, ~] = smoothdata(u_di(T_di_50ns,k1_di(:,2)) - u_di(T_di_50ns,k1_di(:,1)) , 'sgolay', 12);

u_di_pro_V50(1,A60_di) = tmp_60_di;
u_di_pro_V50(1,A120_di) = tmp_120_di;
u_di_pro_V50(1,A240_di) = tmp_240_di;
u_di_pro_V50(1,A300_di) = tmp_300_di;
%%50ns N
tmp_N50_60_di = u_di(T_di_50ns,a_di+A60_di);
tmp_N50_120_di = u_di(T_di_50ns,a_di+A120_di);
tmp_N50_240_di = u_di(T_di_50ns,a_di+A240_di);
tmp_N50_300_di = u_di(T_di_50ns,a_di+A300_di);

[u_di_pro_N50, ~] = smoothdata(u_di(T_di_50ns,a_di+1:a_di+a1_di), 'gaussian', 12);

u_di_pro_N50(1,A60_di) = tmp_N50_60_di;
u_di_pro_N50(1,A120_di) = tmp_N50_120_di;
u_di_pro_N50(1,A240_di) = tmp_N50_240_di;
u_di_pro_N50(1,A300_di) = tmp_N50_300_di;

u_di_pro_N50(1,A60_di-1) = tmp_N50_60_di;
u_di_pro_N50(1,A120_di-1) = tmp_N50_120_di;
u_di_pro_N50(1,A240_di-1) = tmp_N50_240_di;
u_di_pro_N50(1,A300_di-1) = tmp_N50_300_di;

u_di_pro_N50(1,A60_di+1) = tmp_N50_60_di;
u_di_pro_N50(1,A120_di+1) = tmp_N50_120_di;
u_di_pro_N50(1,A240_di+1) = tmp_N50_240_di;
u_di_pro_N50(1,A300_di+1) = tmp_N50_300_di;
%%100ns
tmp_60_di = u_di(T_di_100ns,k1_di(A60_di,2)) - u_di(T_di_100ns,k1_di(A60_di,1));
tmp_120_di = u_di(T_di_100ns,k1_di(A120_di,2)) - u_di(T_di_100ns,k1_di(A120_di,1));
tmp_240_di = u_di(T_di_100ns,k1_di(A240_di,2)) - u_di(T_di_100ns,k1_di(A240_di,1));
tmp_300_di = u_di(T_di_100ns,k1_di(A300_di,2)) - u_di(T_di_100ns,k1_di(A300_di,1));

[u_di_pro_V100, ~] = smoothdata(u_di(T_di_100ns,k1_di(:,2)) - u_di(T_di_100ns,k1_di(:,1)) , 'sgolay', 12);

u_di_pro_V100(1,A60_di) = tmp_60_di;
u_di_pro_V100(1,A120_di) = tmp_120_di;
u_di_pro_V100(1,A240_di) = tmp_240_di;
u_di_pro_V100(1,A300_di) = tmp_300_di;
%%N
tmp_N100_60_di = u_di(T_di_100ns,a_di+A60_di);
tmp_N100_120_di = u_di(T_di_100ns,a_di+A120_di);
tmp_N100_240_di = u_di(T_di_100ns,a_di+A240_di);
tmp_N100_300_di = u_di(T_di_100ns,a_di+A300_di);

[u_di_pro_N100, ~] = smoothdata(u_di(T_di_100ns,a_di+1:a_di+a1_di), 'gaussian', 12);

u_di_pro_N100(1,A60_di) = tmp_N100_60_di;
u_di_pro_N100(1,A120_di) = tmp_N100_120_di;
u_di_pro_N100(1,A240_di) = tmp_N100_240_di;
u_di_pro_N100(1,A300_di) = tmp_N100_300_di;

u_di_pro_N100(1,A60_di-1) = tmp_N100_60_di;
u_di_pro_N100(1,A120_di-1) = tmp_N100_120_di;
u_di_pro_N100(1,A240_di-1) = tmp_N100_240_di;
u_di_pro_N100(1,A300_di-1) = tmp_N100_300_di;

u_di_pro_N100(1,A60_di+1) = tmp_N100_60_di;
u_di_pro_N100(1,A120_di+1) = tmp_N100_120_di;
u_di_pro_N100(1,A240_di+1) = tmp_N100_240_di;
u_di_pro_N100(1,A300_di+1) = tmp_N100_300_di;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)   
hold on
box on;

plot(T_ro, u_ro(:, k1_ro(A30_ro,2)) - u_ro(:,k1_ro(A30_ro,1)), 	'Color', '#A2142F',  'LineWidth', 2); % ro-30 
plot(T_ro, u_ro(:, k1_ro(A60_ro,2)) - u_ro(:,k1_ro(A60_ro,1)), 	'Color', '#EDB120',  'LineWidth', 2); % ro-60 
plot(T_ro, u_ro(:, k1_ro(A90_ro,2)) - u_ro(:,k1_ro(A90_ro,1)),  'Color', '#0072BD',  'LineWidth', 2); % ro-90

plot(T_di, u_di(:, k1_di(A30_di,2)) - u_di(:,k1_di(A30_di,1)), 	'--', 'Color', '#A2142F',  'LineWidth', 2); % di-30
plot(T_di, u_di(:, k1_di(A60_di,2)) - u_di(:,k1_di(A60_di,1)), 	'--','Color', '#EDB120',  'LineWidth', 2); % di-60
plot(T_di, u_di(:, k1_di(A90_di,2)) - u_di(:,k1_di(A90_di,1)),  '--','Color', '#0072BD',  'LineWidth', 2); % di-90

xlabel('Time (ns)', 'FontWeight','bold');
ylabel('TMP (V)', 'FontWeight','bold');
legend({'$\mathbf {A-30^{\circ}}$', '$\mathbf {A-60^{\circ}}$', '$\mathbf {A-90^{\circ}}$', ...
    '$\mathbf {B-30^{\circ}}$', '$\mathbf {B-60^{\circ}}$', '$\mathbf {B-90^{\circ}}$'}, 'Interpreter','latex', 'FontSize', 15);
axis([0 200 0 1.6]) % 设置坐标轴的范围  
set(gca, 'XTick',(0:50:200)) % 设置x坐标轴的刻度 
set(gca, 'YTick',(0:0.4:1.6))  % 设置y坐标轴的刻度
set(gca, 'FontSize',20, 'FontWeight', 'bold','LineWidth', 2) % 设置坐标轴字体是 8
set(gcf,'unit','centimeters','position',[5 5 17 17*(3/4)]);
hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
hold on
box on;
plot(T_ro, u_ro(:, a_ro+A30_ro), 'Color', '#A2142F',  'LineWidth', 2); 
plot(T_ro, u_ro(:, a_ro+A60_ro), 'Color', '#EDB120',  'LineWidth', 2); 
plot(T_ro, u_ro(:, a_ro+A90_ro), 'Color', '#0072BD',  'LineWidth', 2); 

plot(T_di, u_di(:, a_di+A30_di), '--', 'Color', '#A2142F',  'LineWidth', 2);
plot(T_di, u_di(:, a_di+A60_di), '--', 'Color', '#EDB120',  'LineWidth', 2);
plot(T_di, u_di(:, a_di+A90_di), '--','Color', '#0072BD', 'LineWidth', 2);
xlabel('Time (ns)', 'FontWeight','bold');
ylabel('N (m^{-2})', 'FontWeight','bold');
set(gca,'YScale','log');
axis([0 200 10^9 10^16]) % 设置坐标轴的范围  
set(gca, 'XTick',(0:50:200)) % 设置x坐标轴的刻度
set(gca, 'FontSize',20, 'FontWeight', 'bold','LineWidth', 2) % 设置坐标轴字体是 10
set(gcf,'unit','centimeters','position',[5 5 17 17*(3/4)]);
legend({'$\mathbf {A-30^{\circ}}$', '$\mathbf {A-60^{\circ}}$', '$\mathbf {A-90^{\circ}}$', ...
    '$\mathbf {B-30^{\circ}}$', '$\mathbf {B-60^{\circ}}$', '$\mathbf {B-90^{\circ}}$'}, 'Interpreter','latex', 'FontSize', 15);
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
figure(3) %%电压-弧度 17ns
hold on;
box on;
plot(k1_ro(:,4), u_ro_pro_V17, 'Color',[223/255 122/255 094/255], 'LineWidth', 2); % ro-17ns-V
plot(k1_di(:,4), u_di_pro_V17, '-.', 'Color',[223/255 122/255 094/255], 'LineWidth', 2); % di-17ns-V

xlabel('$\mathbf {Angle (\theta)}$',  'Interpreter','latex', 'FontSize', 20,'FontWeight','bold');
ylabel('TMP (V)', 'FontWeight','bold');

legend({'$\mathbf {A-17ns}$', '$\mathbf {B-17ns}$'}, 'Interpreter','latex', 'FontSize', 15);
axis([0 360 -1.8 1.8]) % 设置坐标轴的范围  
set(gca,'XTick', [0 60 120 180 240 300 360]);
set(gca, 'FontSize',20, 'FontWeight', 'bold','LineWidth', 2) % 设置坐标轴字体是 10
set(gcf,'unit','centimeters','position',[5 5 17 17*(3/4)]);

tag1_x = 0:360;
tag1_y = ones(length(tag1_x),1)';
plot(tag1_x, tag1_y, '--k','LineWidth', 2,'HandleVisibility','off');
tag2_x = 0:360;
tag2_y = -ones(length(tag2_x),1)';
plot(tag2_x, tag2_y, '--k','LineWidth', 2,'HandleVisibility','off');
hold off;
%% 50ns
figure(4) 
hold on;
box on;
plot(k1_ro(:,4), u_ro_pro_V50, 'Color',	[060/255 064/255 091/255], 'LineWidth', 2);   %ro-50ns-V
plot(k1_di(:,4), u_di_pro_V50, '-.', 'Color', [060/255 064/255 091/255], 'LineWidth', 2);   %di-50ns-V

xlabel('$\mathbf {Angle (\theta)}$',  'Interpreter','latex', 'FontSize', 20,'FontWeight','bold');
ylabel('TMP (V)', 'FontWeight','bold');
legend({'$\mathbf {A-50ns}$', '$\mathbf {B-50ns}$'}, 'Interpreter','latex', 'FontSize', 15);
axis([0 360 -1.8 1.8]) % 设置坐标轴的范围  
set(gca,'XTick', [0 60 120 180 240 300 360]);
set(gca, 'FontSize',20, 'FontWeight', 'bold','LineWidth', 2) % 设置坐标轴字体是 10
set(gcf,'unit','centimeters','position',[5 5 17 17*(3/4)]);

tag1_x = 0:360;
tag1_y = ones(length(tag1_x),1)';
plot(tag1_x, tag1_y, '--k','LineWidth', 2,'HandleVisibility','off');
tag2_x = 0:360;
tag2_y = -ones(length(tag2_x),1)';
plot(tag2_x, tag2_y, '--k','LineWidth', 2,'HandleVisibility','off');
hold off;
%% 100 ns
figure(5) 
hold on;
box on;
plot(k1_ro(:,4), u_ro_pro_V100, 'Color', [130/255 178/255 154/255], 'LineWidth', 2);   %ro-100ns-V
plot(k1_di(:,4), u_di_pro_V100, '-.', 'Color', [130/255 178/255 154/255], 'LineWidth', 2);   %di-100ns-V

xlabel('$\mathbf {Angle (\theta)}$',  'Interpreter','latex', 'FontSize', 20,'FontWeight','bold');
ylabel('TMP (V)', 'FontWeight','bold');
legend({'$\mathbf {A-100ns}$', '$\mathbf {B-100ns}$'}, 'Interpreter','latex', 'FontSize', 15);
axis([0 360 -1.8 1.8]) % 设置坐标轴的范围  
set(gca,'XTick', [0 60 120 180 240 300 360]);
set(gca, 'FontSize',20, 'FontWeight', 'bold','LineWidth', 2) % 设置坐标轴字体是 10
set(gcf,'unit','centimeters','position',[5 5 17 17*(3/4)]);

tag1_x = 0:360;
tag1_y = ones(length(tag1_x),1)';
plot(tag1_x, tag1_y, '--k','LineWidth', 2,'HandleVisibility','off');
tag2_x = 0:360;
tag2_y = -ones(length(tag2_x),1)';
plot(tag2_x, tag2_y, '--k','LineWidth', 2,'HandleVisibility','off');
hold off;

%%

figure(6) %孔-弧度 17ns 50ns
hold on
box on;
plot(k1_ro(:,4),u_ro_pro_N17, 'Color',  [223/255 122/255 094/255], 'LineWidth', 2);   

plot(k1_di(:,4),u_di_pro_N17, '-.', 'Color', [223/255 122/255 094/255], 'LineWidth', 2);   

%%==========================================================================%%
tag1_x = 0:360;
tag1_y = zeros(length(tag1_x),1)';
tag1_y(1:end) = 1e14;
plot(tag1_x, tag1_y, '--k','LineWidth', 2,'HandleVisibility','off');
tag2_x = 0:360;
tag2_y = zeros(length(tag2_x),1)';
tag2_y(1:end) = 1e14;
plot(tag2_x, tag2_y, '--k','LineWidth', 2,'HandleVisibility','off');
%%===========================================================================%%
xlabel('$\mathbf {Angle (\theta)}$',  'Interpreter','latex', 'FontSize', 20,'FontWeight','bold');
ylabel('N (m^{-2})', 'FontWeight','bold');
shading interp 
axis([0 360 10^9 10^16]) % 设置坐标轴的范围  
legend({'$\mathbf {A-17ns}$', '$\mathbf {B-17ns}$'}, 'Interpreter','latex', 'FontSize', 15);
set(gca,'XTick', [0 60 120 180 240 300 360]);
set(gca, 'FontSize',20, 'FontWeight', 'bold','LineWidth', 2) % 设置坐标轴字体是 10
set(gcf,'unit','centimeters','position',[5 5 17 17*(3/4)]);
set(gca,'YScale','log');

%% 50 ns
figure(7) %孔-弧度 17ns 50ns
hold on
box on;

plot(k1_ro(:,4),u_ro_pro_N50, 'Color',	[060/255 064/255 091/255], 'LineWidth', 2);  

plot(k1_di(:,4),u_di_pro_N50, '-.', 'Color', [060/255 064/255 091/255], 'LineWidth', 2); 

%%==========================================================================%%
tag1_x = 0:360;
tag1_y = zeros(length(tag1_x),1)';
tag1_y(1:end) = 1e14;
plot(tag1_x, tag1_y, '--k','LineWidth', 2,'HandleVisibility','off');
tag2_x = 0:360;
tag2_y = zeros(length(tag2_x),1)';
tag2_y(1:end) = 1e14;
plot(tag2_x, tag2_y, '--k','LineWidth', 2,'HandleVisibility','off');
%%===========================================================================%%
xlabel('$\mathbf {Angle (\theta)}$',  'Interpreter','latex', 'FontSize', 20,'FontWeight','bold');
ylabel('N (m^{-2})', 'FontWeight','bold');
shading interp 
axis([0 360 10^9 10^16]) % 设置坐标轴的范围  
legend({'$\mathbf {A-50ns}$', '$\mathbf {B-50ns}$'}, 'Interpreter','latex', 'FontSize', 15);
set(gca,'XTick', [0 60 120 180 240 300 360]);
set(gca, 'FontSize',20, 'FontWeight', 'bold','LineWidth', 2) % 设置坐标轴字体是 10
set(gcf,'unit','centimeters','position',[5 5 17 17*(3/4)]);
set(gca,'YScale','log');

%% 100 ns
figure(8) %孔-弧度 17ns 50ns
hold on
box on;

plot(k1_ro(:,4),u_ro_pro_N100, 'Color',	[130/255 178/255 154/255], 'LineWidth', 2);  

plot(k1_di(:,4),u_di_pro_N100, '-.', 'Color', [130/255 178/255 154/255]', 'LineWidth', 2); 
%%==========================================================================%%
tag1_x = 0:360;
tag1_y = zeros(length(tag1_x),1)';
tag1_y(1:end) = 1e14;
plot(tag1_x, tag1_y, '--k','LineWidth', 2,'HandleVisibility','off');
tag2_x = 0:360;
tag2_y = zeros(length(tag2_x),1)';
tag2_y(1:end) = 1e14;
plot(tag2_x, tag2_y, '--k','LineWidth', 2,'HandleVisibility','off');
%%===========================================================================%%
xlabel('$\mathbf {Angle (\theta)}$',  'Interpreter','latex', 'FontSize', 20,'FontWeight','bold');
ylabel('N (m^{-2})', 'FontWeight','bold');
shading interp 
axis([0 360 10^9 10^16]) % 设置坐标轴的范围  
legend({'$\mathbf {A-100ns}$', '$\mathbf {B-100ns}$'}, 'Interpreter','latex', 'FontSize', 15);
set(gca,'XTick', [0 60 120 180 240 300 360]);
set(gca, 'FontSize',20, 'FontWeight', 'bold','LineWidth', 2) % 设置坐标轴字体是 10
set(gcf,'unit','centimeters','position',[5 5 17 17*(3/4)]);
set(gca,'YScale','log');



