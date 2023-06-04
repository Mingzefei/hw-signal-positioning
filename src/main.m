%% 数据写入

% 定义存储采集数据的矩阵
n_signals = zeros(8,8,4000); % 正常状态下数据
d_signals = zeros(8,8,4000); % 损伤状态下数据

% 循环读取文档中的采集数据
for i = 1:7
    % 读取无损伤采集数据
    filename = sprintf("..\\data\\raw\\a%dtone5p150k", i);
    data = load(filename);
    for j = i:8
        % 写入元组
        n_signals(i,j,:) = data(j,:);
    end
    % 读取损伤采集数据
    filename = sprintf("..\\data\\raw\\da%dtone5p150k", i);
    data = load(filename);
    for j = i:8
        % 写入元组
        d_signals(i,j,:) = data(j,:);
    end
end

% 绘制 1->2 信号（有、无损伤）
n_data_12 = reshape(n_signals(1,2,:),1,4000);
d_data_12 = reshape(d_signals(1,2,:),1,4000);
figure;
t = 1:length(n_data_12);
plot(t,d_data_12,'r')
plot(t,n_data_12,'b')
saveas(gcf,'..\\reports\\figures\\由1到2的原始信号.png')

%% 数据预处理

for i = 1:7
    for j = i:8
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
figure;
t = 1:length(n_data_12);
plot(t,d_data_12,'r')
plot(t,n_data_12,'b')
saveas(gcf,'..\\reports\\figures\\由1到2的预处理信号.png')

%% 算反映损伤信息的特征参数

% 选取前N个数据
N = 150;
n_signals = n_signals(:,:,1:N);
d_signals = d_signals(:,:,1:N);
% 定义存储相关系数的矩阵
rho = zeros(8,8);
DI = zeros(8,8);
for i = 1:7
    for j = i:8
        r_temp = corrcoef(n_signals(i,j,:),d_signals(i,j,:));
        rho(i,j) = r_temp(1,2);
        DI(i,j) = 1-abs(rho(i,j));
    end
end
% 定义存储概率的矩阵
p_xy = zeros(320,320);
for x = -160:160
    for y = -160:160
        p_xy(i,j) = P(DI,x,y,1,1.02);
    end
end
%% 生成图像


%% 定义损伤概率 P
function p = P(DI,x,y,alpha,beta)
    theta = [6 5 4 3 2 1 0 -1 -2]*2*pi/8; % 8 个传感器的角坐标
    d = 240; % 传感器分布圆的直径，单位：mm
    s_x = d/2*cos(theta);
    s_y = d/2*sin(theta);
    p = 0;
    for i = 1:7
        for j = i:8
           pos_i = [s_x(i) s_y(i)];
           pos_j = [s_x(j) s_y(j)];
           pos_xy = [x y];
           RD = (pdist2(pos_i,pos_xy,'euclidean')+pdist2(pos_j,pos_xy,'euclidean'))/pdist2(pos_i,pos_j,'euclidean');
           R = min(RD, beta);
           p = p + DI(i,j)*(beta-R)/(beta-1);
        end
    end
    p = power(p, alpha);
end