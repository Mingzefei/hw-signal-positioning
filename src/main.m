close all;
clear;
%% 开启并行
my_par = parpool; 

%% 数据写入

% 定义存储采集数据的矩阵
n_signals = zeros(8,8,4000); % 正常状态下数据
d_signals = zeros(8,8,4000); % 损伤状态下数据

% 循环读取文档中的采集数据
for i = 1:7
    % 读取无损伤采集数据
    filename = sprintf("..\\data\\raw\\a%dtone5p150k", i);
    data = load(filename);
    for j = i+1:8
        n_signals(i,j,:) = data(j,:);
    end
    % 读取损伤采集数据
    filename = sprintf("..\\data\\raw\\da%dtone5p150k", i);
    data = load(filename);
    for j = i+1:8
        d_signals(i,j,:) = data(j,:);
    end
end

% 绘制 1->2 信号（有、无损伤）
n_data_12 = reshape(n_signals(1,2,:),1,4000);
d_data_12 = reshape(d_signals(1,2,:),1,4000);
f1 = figure(1);
t = 1:length(n_data_12);
plot(t,d_data_12,'r','LineWidth',1.5)
plot(t,n_data_12,'LineWidth',1.5)
set(gca, 'linewidth',0.8)
saveas(f1,'..\\reports\\figures\\由1到2的原始信号.png')
close(f1);
n_data_12 = reshape(n_signals(1,2,1:150),1,150);
d_data_12 = reshape(d_signals(1,2,1:150),1,150);
f1 = figure(1);
t = 1:length(n_data_12);
plot(t,d_data_12,'r','LineWidth',1.5)
plot(t,n_data_12,'LineWidth',1.5)
set(gca, 'linewidth',0.8)
saveas(f1,'..\\reports\\figures\\由1到2的原始信号（前150）.png')
close(f1);

%% 数据预处理

for i = 1:7
    for j = i+1:8
        % 去平均
        n_signals(i,j,:) = n_signals(i,j,:) - mean(n_signals(i,j,:));
        d_signals(i,j,:) = d_signals(i,j,:) - mean(d_signals(i,j,:));
        % 归一化
        n_signals(i,j,:) = n_signals(i,j,:) / max(abs(n_signals(i,j,:)));
        d_signals(i,j,:) = d_signals(i,j,:) / max(abs(d_signals(i,j,:)));
    end
end

% 绘制 1->2 信号（有、无损伤）
n_data_12 = reshape(n_signals(1,2,:),1,4000);
d_data_12 = reshape(d_signals(1,2,:),1,4000);
f2 = figure(2);
t = 1:length(n_data_12);
plot(t,d_data_12,'r','LineWidth',1.5)
plot(t,n_data_12,'LineWidth',1.5)
saveas(f2,'..\\reports\\figures\\由1到2的预处理信号.png')
close(f2);
n_data_12 = reshape(n_signals(1,2,1:150),1,150);
d_data_12 = reshape(d_signals(1,2,1:150),1,150);
f2 = figure(2);
t = 1:length(n_data_12);
plot(t,d_data_12,'r','LineWidth',1.5)
plot(t,n_data_12,'LineWidth',1.5)
saveas(f2,'..\\reports\\figures\\由1到2的预处理信号（前150）.png')
close(f2);

% 绘制 5->6 信号（有、无损伤）
n_data_56 = reshape(n_signals(5,6,1:150),1,150);
d_data_56 = reshape(d_signals(5,6,1:150),1,150);
f2 = figure(2);
t = 1:length(n_data_56);
plot(t,d_data_56,'r','LineWidth',1.5)
plot(t,n_data_56,'LineWidth',1.5)
saveas(f2,'..\\reports\\figures\\由5到6的预处理信号（前150）.png')
close(f2);

%% 算反映损伤信息的特征参数

% 选取前N个数据
N = 150;
n_signals = n_signals(:,:,1:N);
d_signals = d_signals(:,:,1:N);
% 定义存储相关系数的矩阵
rho = zeros(8,8);
DI = zeros(8,8);
for i = 1:7
    for j = i+1:8
        r_temp = corrcoef(n_signals(i,j,:),d_signals(i,j,:));
        rho(i,j) = r_temp(1,2);
        DI(i,j) = 1 - rho(i,j);
    end
end
% 定义存储概率的矩阵
length = 130*2; % 待求正方形区域边长，单位：mm
resolution = 2; % 分辨率为2，单位：mm
p_xy = zeros(length/resolution,length/resolution);
pos_x = zeros(length/resolution,length/resolution);
pos_y = zeros(length/resolution,length/resolution);
for i = 1:length/resolution
    parfor j = 1:length/resolution
        x = -length/2 + resolution*i;
        y = -length/2 + resolution*j;
        p_xy(i,j) = P(DI,x,y,1,1.05);
        pos_x(i,j) = x;
        pos_y(i,j) = y;
    end
end

%% 生成图像

f3 = figure(3);
contourf(pos_x,pos_y,p_xy)
xlabel('x');  % 设置 x 轴标签
ylabel('y');  % 设置 y 轴标签
colorbar;
saveas(f3,'..\\reports\\figures\\概率图-beta-105-alpha-100.png')
close(f3);

%% 敏感性分析

% beta 
N = 150;
alpha = 1;
beta_list = [1.01, 1.02, 1.03, 1.05, 1.07, 1.10, 1.15,1.25,1.35,1.55,1.75,2.10];
H_beta_list = zeros(size(beta_list));
for i = 1:numel(beta_list)
    beta = beta_list(i);
    [p_xy,pos_x,pos_y] = Analyse(n_signals,d_signals,N,alpha,beta);
    % 归一化
    p_xy_normalized = p_xy / sum(p_xy(:));
    % 计算信息熵
    H_xy = -sum(p_xy_normalized(:) .* log2(p_xy_normalized(:)), 'omitnan');
    H_beta_list(i) = H_xy;
end
f_beta = figure();
plot(beta_list, H_beta_list, 'o-','LineWidth',1.5)
xlabel('$\beta$','Interpreter', 'latex')
ylabel('$H(\beta)$','Interpreter', 'latex')
xlim([0.9,2.4]);
ylim([12.8,14.5]);
set(gca, 'linewidth',0.8)
saveas(f_beta,'..\\reports\\figures\\信息熵-beta.png')
close(f_beta);

% alpha 
N = 150;
beta = 1.02;
alpha_list = [1, 2, 3, 4, 5, 7, 9, 12, 15,20];
H_alpha_list = zeros(size(alpha_list));
for i = 1:numel(alpha_list)
    alpha = alpha_list(i);
    [p_xy,pos_x,pos_y] = Analyse(n_signals,d_signals,N,alpha,beta);
    % 归一化
    p_xy_normalized = p_xy / sum(p_xy(:));
    % 计算信息熵
    H_xy = -sum(p_xy_normalized(:) .* log2(p_xy_normalized(:)), 'omitnan');
    H_alpha_list(i) = H_xy;
end
f_alpha = figure();
plot(alpha_list, H_alpha_list, 'o-', 'LineWidth',1.5)
xlabel('$\alpha$','Interpreter', 'latex')
ylabel('$H(\alpha)$','Interpreter', 'latex')
xlim([0,22]);
ylim([0,15]);
set(gca, 'linewidth',0.8)
saveas(f_alpha,'..\\reports\\figures\\信息熵-alpha.png')
close(f_alpha);

%% 结束并行

delete(my_par)

%% 定义损伤概率 P
function p = P(DI,x,y,alpha,beta)
    theta = [6 5 4 3 2 1 0 -1 -2]*2*pi/8; % 8 个传感器的角坐标
    d = 240; % 传感器分布圆的直径，单位：mm
    s_x = d/2*cos(theta);
    s_y = d/2*sin(theta);
    p = 0;
    for i = 1:7
        for j = i+1:8
           pos_i = [s_x(i) s_y(i)];
           pos_j = [s_x(j) s_y(j)];
           pos_xy = [x y];
           RD = (pdist2(pos_i,pos_xy,'euclidean')+pdist2(pos_j,pos_xy,'euclidean'))...
               /pdist2(pos_i,pos_j,'euclidean');
           R = min(RD, beta);
           p = p + DI(i,j)*(beta-R)/(beta-1);
        end
    end
    p = power(p, alpha);
end

%% 定义敏感性分析函数 Analyse
function [p_xy,pos_x,pos_y] = Analyse(n_signals,d_signals,N,alpha,beta)
    n_signals = n_signals(:,:,1:N);
    d_signals = d_signals(:,:,1:N);
    % 定义存储相关系数的矩阵
    rho = zeros(8,8);
    DI = zeros(8,8);
    for i = 1:7
        for j = i+1:8
            r_temp = corrcoef(n_signals(i,j,:),d_signals(i,j,:));
            rho(i,j) = r_temp(1,2);
            DI(i,j) = 1 - rho(i,j);
        end
    end
    % 定义存储概率的矩阵
    length = 130*2; % 待求正方形区域边长，单位：mm
    resolution = 2; % 分辨率为2，单位：mm
    p_xy = zeros(length/resolution,length/resolution);
    pos_x = zeros(length/resolution,length/resolution);
    pos_y = zeros(length/resolution,length/resolution);
    for i = 1:length/resolution
        parfor j = 1:length/resolution
            x = -length/2 + resolution*i;
            y = -length/2 + resolution*j;
            p_xy(i,j) = P(DI,x,y,alpha,beta);
            pos_x(i,j) = x;
            pos_y(i,j) = y;
        end
    end
    % 绘图
    f = figure();
    contourf(pos_x,pos_y,p_xy,14)
    xlabel('x');  % 设置 x 轴标签
    ylabel('y');  % 设置 y 轴标签
    title(sprintf('$\\beta=$%.2f, $\\alpha=$%.2f', beta, alpha),'Interpreter', 'latex');
    colorbar;
    saveas(f,sprintf('..\\reports\\figures\\概率图-beta-%d-alpha-%d.png', beta*100, alpha*100))
    close(f);
end

