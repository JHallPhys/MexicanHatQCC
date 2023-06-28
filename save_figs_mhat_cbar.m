clear all

% Custom heat maps
% inferno=inferno();
% viridis=viridis();
% return

input_2=input('Do you want to save the figs [y/n]:  ','s');

nrow=2;
ncol=2;
string_d=strcat('_',num2str(nrow),'x',num2str(ncol));
% Get the figure dimensions for page
text_width=13.7; % paper width
[dim_out,width_out,height_out]=dfig(nrow,ncol,text_width);
width=dim_out;
height=dim_out;
str_ex='.png';

if input_2=='y'
% Save the current directory
parent_d = cd;    
cd '/Users/joe/Desktop/Mexican_Hat_Project_Kate/Notes_Mhat/figs'
fnum=length(findobj('type','figure'));
for f=1:fnum

disp(strcat('figure ',num2str(f)))
input_2=input('[y/n]:  ','s');

if input_2=='y'
% Save fig 1 
fname1 = input('type file name fig :','s'); 
% fname1=strcat(fname1,string_d);
figure(f)
fze=5;
set(gca,'FontSize',fze)
box on
% xlabel('q')
% ylabel('p')
% xticks([0 0.5 1])
% yticks([-0.5 0 0.5])
% xticks([-1 0 1])
% yticks([-1 0 1])
% xticks([-1.5 0 1.5])
% yticks([-1.5 0 1.5])
title(' ')
colorbar

set(gca,'FontSize',fze);
set(gca,'Units','normalized');
yhere = get(gca,'ylabel'); % handle to the label object
ypos = get(yhere,'position'); % get the current position property of the label
paxis_inner=get(gca,'Position');
paxis_outer=get(gca,'OuterPosition');

eps=0.02; % small parameter to create a little bit of a boundary for y
lshift=paxis_inner(1)+ypos(1); % find the diference between the ylabel and edge
lshift=lshift-eps;
set(gca,'Position',[paxis_inner(1) paxis_inner(2) paxis_inner(3) paxis_inner(4)])
paxis_inner=get(gca,'Position');  % Get the new information
% =========================================================================
% Fit the colorbar to the the right side of figure    
% =========================================================================
cbar_xpos=paxis_inner(3)+paxis_inner(1); % x coordinate of upper limit of x-axis
cbar_ypos=paxis_inner(2); % y coordinate of upper limit of x-axis 
%==========================================================================
% Height of the colorbar
%==========================================================================
cbar_height=paxis_inner(4); % same size as the y-axis
yhere = get(gca,'ylabel'); % handle to the label object
ypos = get(yhere,'position'); % get the current position property of the label
lshift=paxis_inner(1)+ypos(1);
cbar_width=0.02;
%==========================================================================
% Place the colourbar (asssumes don't have colour
%==========================================================================
c = colorbar('eastoutside');
set(c,'Position',[cbar_xpos cbar_ypos cbar_width cbar_height])
L=cellfun(@(x)sprintf('%.3f',x),num2cell(get(c,'xtick')),'Un',0);
set(c,'xticklabel',L)
% set(c,'FontSize',4) % This is how we keep font size the same
set(c,'Units','centimeters');
cbar_width=get(c,'Position');
cbar_width=cbar_width(3);
set(c,'Units','normalized');
% ==========================================================================
% ==========================================================================
set(gcf,'units','centimeters','position',[0,0,width+0.15*width,height])
% return
ax = gca;
% exportgraphics(ax,strcat(fname1,str_ex),'Resolution',600) % High Quality
exportgraphics(ax,strcat(fname1,str_ex),'Resolution',300) % Lower quality



saveas(gcf,fname1)
% fig = gcf;
% fig.PaperUnits = 'centimeters';
% fig.PaperPosition = [0  0 width height];
% print(fname1,'-dpng','-r600')

end

end
cd(parent_d)
end


close all