function createfigure(X1, Y1)
%CREATEFIGURE(X1, Y1)
%  X1:  x 数据的向量
%  Y1:  y 数据的向量

%  由 MATLAB 于 05-Jun-2023 16:21:10 自动生成

% 创建 figure
figure1 = figure;

% 创建 axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% 创建 plot
plot(X1,Y1,'DisplayName','传感器1到传感器2','Parent',axes1,'Color',[0 0 1]);

% 取消以下行的注释以保留坐标区的 X 范围
% xlim(axes1,[0 4000]);
% 取消以下行的注释以保留坐标区的 Y 范围
% ylim(axes1,[-2.4431947284708 1.0568052715292]);
box(axes1,'on');
hold(axes1,'off');
% 创建 legend
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.636749199593712 0.837448559670783 0.248968357807013 0.0624142618001403],...
    'Interpreter','latex',...
    'FontSize',12,...
    'EdgeColor',[1 1 1]);

% 创建 axes
axes2 = axes('Parent',figure1,...
    'Position',[0.431911966987621 0.156378600823049 0.383768913342504 0.279835390946499]);
hold(axes2,'on');

% 创建 plot
plot(X1,Y1,'DisplayName','前150个数据点','Color',[0 0 1]);

% 取消以下行的注释以保留坐标区的 X 范围
% xlim(axes2,[0 151.595889801654]);
% 取消以下行的注释以保留坐标区的 Y 范围
% ylim(axes2,[-1.43338867712711 1]);
box(axes2,'on');
hold(axes2,'off');
% 创建 legend
legend2 = legend(axes2,'show');
set(legend2,...
    'Position',[0.437291127286565 0.372246560870038 0.21733149492593 0.0514403279303523],...
    'FontSize',12,...
    'EdgeColor',[1 1 1]);

% 创建 rectangle
annotation(figure1,'rectangle',...
    [0.132049518569464 0.345679012345679 0.0261348005502063 0.504115226337449],...
    'Color',[1 0 0],...
    'LineWidth',1.5);

% 创建 arrow
annotation(figure1,'arrow',[0.171939477303989 0.392022008253095],...
    [0.390946502057613 0.308641975308642],'Color',[1 0 0]);

