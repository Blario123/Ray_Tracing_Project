function plot_variance(X1, Y1, X2, Y2)

%
PLOT_VARIANCE(X1, Y1, X2, Y2
)
%  X1:
vector of
x data
%  Y1:
vector of
y data
%  X2:
vector of
x data
%  Y2:
vector of
y data

%  Auto-
generated by
MATLAB on
08-Jan-2022 20:31:49

%
Create figure
figure('WindowState','fullscreen');

%
Create axes
axes1 = axes;
hold(axes1,
'on');

%
Create loglog
loglog(X1, Y1,
'SeriesIndex',1,...
'DisplayName','Variance of the estimator',...
'MarkerSize',8,...
'LineWidth',1,...
'Color',[0 0 0]);

%
Create loglog
loglog(X2, Y2,
'DisplayName','Line with gradient of -1','LineWidth',1,...
'Color',[1 0 0]);

%
Create ylabel
ylabel('Variance','Interpreter','latex');

%
Create xlabel
xlabel('N','Interpreter','latex');

%
Create title
title('Variance vs N','Interpreter','latex');

%
Uncomment the
following line
to preserve
the X
-
limits of
the axes
%
xlim(axes1,
[10 1000]);
box(axes1,
'on');
hold(axes1,
'off');
%
Set the
remaining axes

properties
set(axes1,

'FontSize',20,'LabelFontSizeMultiplier',1.5,'XMinorTick','on',...
'XScale','log','YMinorTick','on','YScale','log');
%
Create legend
legend1 = legend(axes1, 'show');
set(legend1, ...
'Position',[0.618055555555553 0.758017492711368 0.268518518518518 0.148688046647227],...
'FontSize',20,...
'FontName','Helvetica Neue');
end
