function [] = slider_plot_tempresid(MaxTemp,colIdx,residMatrix,locationsFormatted)
% Plot different plots according to slider location.

DatesBeingForecasted = datenum(cat(1,MaxTemp(:,colIdx).DateBeingForecasted));
NumDaysAhead = colIdx - 1;
datemin = min(DatesBeingForecasted);
datemax = max(DatesBeingForecasted);
residmin = ceil(min(cat(1,residMatrix(:))));
residmax = ceil(max(cat(1,residMatrix(:))));

figure('units','normalized',...
    'position',[0.05 0.05 0.6 0.8],...
    'menubar','none',...
    'name','slider_plot',...
    'numbertitle','off',...
    'resize','on');
axes('unit','normalized',...
    'position',[0.1 0.2 0.8 0.7]);
LN = plot(DatesBeingForecasted,residMatrix(:,1));
title(locationsFormatted{1});
uicontrol('style','slide',...
    'unit','normalized',...
    'position',[0.1 0.05 0.8 0.05],...
    'min',1,'max',size(residMatrix,2),'val',1,...
    'sliderstep',[1/(size(residMatrix,2)-1) 10/(size(residMatrix,2)-1)],...
    'callback',{@sl_call,residMatrix,locationsFormatted,LN});
xlim([datemin,datemax])
datetick('x','mmm-yyyy')
xlabel(['Date Being Forecasted (',sprintf('%d days ahead)',NumDaysAhead)])
ylim([residmin,residmax])
ylabel('Residual')

function [] = sl_call(varargin)
% Callback for the slider.
[h,S,locations,LN] = varargin{[1,3,4,5]};  % calling handle and data structure.
set(LN,'ydata',S(:,round(get(h,'value'))))
title(locations{round(get(h,'value'))})