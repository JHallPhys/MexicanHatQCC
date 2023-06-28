clear all

% Add cumstom string based on size of fig

nrow=1;
ncol=1;
fnt_scale=0.1;
string_d=strcat('_',num2str(nrow),'x',num2str(ncol));

% Get the figure dimensions for page
text_width=13.7; % paper width 
[dim_out,width_out,height_out]=dfig(nrow,ncol,text_width);
width=dim_out
height=dim_out

% return
% Save the figs
% input_3=input('Is it going abobve fig with cbar?  ','s');
input_2=input('Do you want to save the figs [y/n]:  ','s');
if input_2=='y'
    
% Save the current directory
parent_d = cd;    
cd '/Users/joe/Desktop/Mexican_Hat_Project_Kate/Notes_Mhat/figs'


fnum=length(findobj('type','figure'));
for f=1:fnum

disp(strcat('figure ',num2str(f)))
input_2=input('[y/n]:  ','s');

if input_2=='y'
% Save fig
fname1 = input('type file name fig :','s'); 
% fname1=strcat(fname1,string_d);
figure(f)
fze=5;
% fze=8;
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
% xticks([0 200 400 600 800 1000])
% yticks([0 0.2 0.4 0.6 0.8 1])
saveas(gcf,fname1)
fig = gcf;
fig.PaperUnits = 'centimeters';
% width=4.5;
% height=4.5;
fig.PaperPosition = [0 0 width height];
% fig.PaperPosition = [0 0 width+fnt_scale*width height]; % reformat
% if input_3=='y'
% set(gcf,'units','centimeters','position',[0,0,width+0.15*width,height])
% end
print(fname1,'-dpng','-r300')

end

end
cd(parent_d)
end

close all
