function goodplot(papersize, margin, fontsize,ADJUST_FIG_ONLY)
% function which produces a nice-looking plot
% and sets up the page for nice printing
% fonttype = 'Arial';
fonttype = 'Times New Roman';
% fonttype = 'Helvetica';
% fonttype = 'CMU Serif Roman';
if nargin == 0
    papersize = [6 2.5];
    margin = [0.01 0.01];
    fontsize = 16;
    ADJUST_FIG_ONLY = false;
elseif nargin == 1
    margin = [0.01 0.01];
    fontsize = 16;
    ADJUST_FIG_ONLY = false;
elseif nargin == 2
    fontsize = 16;
    ADJUST_FIG_ONLY = false;
elseif nargin == 3
    ADJUST_FIG_ONLY = false;
end
if  ~ADJUST_FIG_ONLY
    ax = gca;
    ax.FontSize = fontsize-3;
    set(get(gca,'xlabel'),'FontName',fonttype,'FontSize', fontsize);
    set(get(gca,'ylabel'),'FontName',fonttype,'FontSize', fontsize);
    set(get(gca,'title'),'FontName',fonttype,'FontSize', fontsize-1);
    H=findobj(gca,'Type','text');
    set(H,'FontName',fonttype,'FontSize', fontsize-2);
    
    % box off; axis square;
    % set(gca,'LineWidth',2);
    % set(gca,'DefaultLineLineWidth',2)
%     set(gca,'FontName',fonttype);
%     set(gca,'FontSize',fonstsize);
    % set(gca,'FontWeight','Bold');
    % set(gcf,'color','w');
    set(gca,'Box','on');
end
% hLegend = findobj(gca, 'Type', 'Legend');
legend boxoff  % remove the border of the legend.

% use plot(xx,'HandleVisibility','off') to hide the legend
h = findobj(gcf,'Tag','legend');
set(h,'FontSize',fontsize-2);
for i=1:length(h)
    if ~isempty(h(1).String) && strcmp(h(i).String{1},'data1')  % indicates automatically generated legend
        set(h(i),'Visible','off');
    else 
        set(h(i),'FontSize',fontsize-1) % 15 is default
        h(i).ItemTokenSize = [13.4400    8.0640]*1.3; %h(i).ItemTokenSize*0.8;
    end
end

set(gcf,'PaperUnits','inches');
set(gcf,'PaperSize', papersize);
set(gcf,'PaperPosition',[margin(1)-0.22 margin(2) papersize(1)+0.75-margin(1) papersize(2)+0.1-margin(2)]); % left bottom width height
set(gcf,'PaperPositionMode','Manual');
% for print 
%print -painters -dpdf -r150 ModelSetFR.pdf
end