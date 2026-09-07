%%%set environment%%%
clear all;
%clf;
set(0,'DefaultFigureWindowStyle','normal');
set(0,'defaulttextinterpreter','latex');

%Plot Format
set(0,'DefaultTextFontsize',18, ...
    'DefaultTextFontname','Arial', ...
    'DefaultTextFontWeight','normal', ...
    'DefaultTextFontname','Arial', ...
    'DefaultAxesFontsize',18, ...
    'DefaultAxesFontname','Arial', ...
    'DefaultLineLineWidth', 1.5)

set(0,'defaulttextinterpreter','latex')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
datatype = 'casestudy';

rootdir = './fig_fmax';
Figdir = './fig/colorcontour';
FileFormat = 'pdf';

FaultLength = 64.8456;

SaveFigure = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R0 = 1192.1598;
N_L = R0;

%%

D = load(sprintf('./sensordataset_%s/plotfmax_%s.mat', datatype, datatype));

% derive amplification
cx = D.plotdata.cx * N_L;
cy = D.plotdata.cy * N_L;
%famp = D.plotdata.fmax ./ C.plotdata.fmax;

%extract artifacts
%%famp(famp>30) = 15;
%FMAX IS ALREADY IN REAL DIMENSION
%
% d1 = log10(D.plotdata.fmax);
% mn = min(d1(:));
% rng = max(d1(:))-mn;
% d = 1+63*(d1-mn)/rng; % Self scale data
% L = [1 2 5 10 20 50 70];
% %L = [1 5 10 50 100];
% % Choose appropriate
% % or somehow auto generate colorbar labels
% ll = 1+63*(log10(L)-mn)/rng;


%%

fig = figure(1);
fig.Units = 'normalized';
fig.Position = [0 1 0.6 0.6];
clf(fig,'reset'); cla(fig,'reset'); hold on;
axis equal;

gridsize = 0.1 * R0;
domainY = 8.8075 * R0;

[xq,yq] = meshgrid(1*N_L:gridsize:25*N_L, -domainY:gridsize:domainY);
vq = griddata(cx,cy,log10(D.plotdata.fmax),xq,yq, 'v4');

contourresolution = 100;
[~, hp] = contourf(xq/1e3,yq/1e3,vq,contourresolution,'LineWidth',0.1);
set(hp,'LineColor','none');


%%
contourresolution_l= log10([10 10]);

[hph, hpl] = contourf(xq/1e3,yq/1e3,vq,contourresolution_l,...
   'LineColor', 'k', 'LineWidth',0.5);
set(hpl,'Fill','off');
%%
% hp = pcolor(xq,yq,vq);
% set(hp,'EdgeAlpha',0);
% hold on;
% hs = scatter(cx/1e3,cy/1e3,20,'b' , 'x');
% hs.MarkerEdgeAlpha = 0.4;

%%
% for i = 1:D.plotdata.Cnet_count
%     x1 = D.plotdata.Cnet_table(i).x1/1e3;
%     y1 = D.plotdata.Cnet_table(i).y1/1e3;
%     x2 = D.plotdata.Cnet_table(i).x2/1e3;
%     y2 = D.plotdata.Cnet_table(i).y2/1e3;
%
%     for j = 1:length(x1)
%         %plot cracks
%         hf = plot([x1(j) x2(j)],[y1(j) y2(j)], '-', 'Color', 'k', 'LineWidth',1);
%         hf.Color(4)=0.4;
%     end
% end

%axis
XLimit = [0 25] * N_L/1e3;
YLimit = [-domainY domainY]/1e3;

ax1 = gca;
xlabel(ax1,'$x$ (km)');
ylabel(ax1,'$y$ (km)');
ax1.Layer = 'top';
ax1.XLim = XLimit;
ax1.YLim = YLimit;

%importColormapFromParaview('parajet.json', 0, 1, 0)
importColormapFromParaview('pararainbow.json', 0, 0, 0);

cmin = log10(0.1);
cmax = log10(100);
caxis([cmin cmax]);
L = [0.1 1 2 5 10 20 50 100]; %(m/s)
Ltick = log10(L);
h = colorbar('southoutside','YTick',Ltick,'YTickLabel',L);
ylabel(h, 'Maximum cut-off frequency (Hz)');

%%
%cbpos = get(h,'position');
%set(h, 'Position', [1.125*cbpos(1) cbpos(2)+0.08 cbpos(3) 0.85*cbpos(4)]);
%set(ax1,'outerposition',[0.1,0.1,0.85,0.85]);

pos = get(gca,'pos');
set(gca,'pos',[pos(1) pos(2) 0.9*pos(3) 0.9*pos(4)]);

%next plot axis above for Nc
ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');

xlabel(ax2,'$x/R_0$');
ylabel(ax2,'$y/R_0$');
linkaxes([ax1,ax2])
axis equal
ax2.XLim = XLimit;
ax2.YLim = YLimit;

%xaxis cohesive zone resolution
nondimpoints_x = [0 5 10 15 20 25] * (N_L/1e3) ; %km
nondimvalues = round(nondimpoints_x ./ (N_L/1e3),1);
ax2.XTick = nondimpoints_x;
ax2.XTickLabel = num2cell(nondimvalues);

%yaxis cohesive zone resolution
nondimpoints_y = linspace(-8,8,5)*R0/1e3; %km
nondimvalues_y = round(nondimpoints_y ./ (N_L/1e3),1);
ax2.YTick = nondimpoints_y;
ax2.YTickLabel = num2cell(nondimvalues_y);


%%
if SaveFigure
    fodir = [Figdir,'/',FileFormat];
    mkdir(fodir);
    foname = sprintf('fmax_presentation%s.%s',datatype, FileFormat);
    set(gcf, 'Color', 'none');
    %set(ax1, 'Visible', 'off');
    %set(ax2, 'Visible', 'off');
    export_fig([fodir,'/',foname],'-nocrop','-r0');
end
