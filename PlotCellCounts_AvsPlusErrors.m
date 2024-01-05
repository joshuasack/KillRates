% Version 07/08/2021
% plot cell counts over time
load colorblind_colormap/colorblind_colormap
% type "colornames" to see the names of each of the 12 colors.
% colornames =
%   12×1 cell array
% 
%     {'blue'      }
%     {'red'       }
%     {'yellow'    }
%     {'darkgray'  }
%     {'black'     }
%     {'orange'    }
%     {'magenta'   }
%     {'teal'      }
%     {'darkblue'  }
%     {'darkgreen' }
%     {'cyan'      }
%     {'darkorchid'}

% Set the path where the data is
path = 'CTL28000/dataCounts/';
RunName ='aveCounts.dat';



% specify the time points
TimePoints = 0:10:300;
% generate a matrix to store the cell counts: one row for each time point,
% 5 columns, where the first one is the time step.
% The file 'aveCounts' has six columns with two rows for each time point.
% The first row are the average counts, the second the standard errors.
% See the file ``aveStats" for the names of each column:
% loMHC, hiMHC ,ctlprerec, ctlpostrec, nklprerec, nklpostrec
CellCounts = zeros(length(TimePoints),5); 
StdErrors = zeros(length(TimePoints),5);
Data = load([path,RunName]);
for i=1:length(TimePoints)
    TimePoint = TimePoints(i);
    CellCounts(i,1)=TimePoint;
    CellCounts(i,2:5)= Data((2*(i-1)+1),1:4 );
    StdErrors(i,1)=TimePoint;
    StdErrors(i,2:5) = Data((2*i),1:4);
end
save([path,'/CellCounts.mat'],'CellCounts','StdErrors')

%% Make a new figure: not logarithmic, no pre-recognition
orng = '#FF8800';  % orange color
magnta = colorblind(7,:); % magenta
drkorchid = colorblind(12,:); % dark orchid
drkblue = colorblind(6,:); % dark blue

f = figure
%semilogy(CellCounts(:,1),CellCounts(:,2),'linewidth',2);
%hold on
plot(CellCounts(:,1),CellCounts(:,2),'linewidth',3,'Color',drkblue);
hold on
plot(CellCounts(:,1),CellCounts(:,3),'linewidth',3,'Color',drkorchid);
plot(CellCounts(:,1),CellCounts(:,5),'linewidth',3,'LineStyle','-- ','Color',magnta);
% add errors 
upper_LoMHC = CellCounts(:,2)+StdErrors(:,2);
lower_LoMHC = max((CellCounts(:,2)-StdErrors(:,2)),eps); %avoid negative values
fill_time = [CellCounts(:,1);flipud(CellCounts(:,1))];
fill_border_LoMHC = [upper_LoMHC;flipud(lower_LoMHC)];
fill(fill_time,fill_border_LoMHC,drkblue,'FaceAlpha',0.1, ...
    'EdgeColor',drkblue)
upper_HiMHC = CellCounts(:,3)+StdErrors(:,3);
lower_HiMHC = max((CellCounts(:,3)-StdErrors(:,3)),eps); %avoid negative values
fill_time = [CellCounts(:,1);flipud(CellCounts(:,1))];
fill_border_HiMHC = [upper_HiMHC;flipud(lower_HiMHC)];
fill(fill_time,fill_border_HiMHC,drkorchid,'FaceAlpha',0.1, ...
    'EdgeColor',drkorchid)
upper_post = CellCounts(:,5)+StdErrors(:,5);
lower_post = max((CellCounts(:,5)-StdErrors(:,5)),eps);
fill_border_post = [upper_post;flipud(lower_post)];
fill(fill_time,fill_border_post,magnta,'FaceAlpha',0.05, ...
    'EdgeColor',magnta)
axis([0, max(CellCounts(:,1)), 1, 80])
xlabel('Time','interpreter','latex','fontsize',20)
ylabel('Cell Counts','interpreter','latex','fontsize',20)
legend('Low MHC tumor','High MHC tumor',...
    'CTL post-recognition',...
    'Location','southeast', ...
    'fontsize',18, 'interpreter','latex')
title('Average Cell Counts: CTL only','interpreter','latex','FontSize',24)
%set(gcf, 'PaperPositionMode', 'auto','PaperOrientation','landscape');
%print(gcf, ['CellCounts',RunName,'.pdf'], '-dpdf', '-fillpage');
exportgraphics(f,[path,'AverageCellCounts_wstd','.pdf'],'ContentType','vector')
%close(f)

%% WE CAN SEE THE CHANGES BETTER USING A LOG SCALE FOR THE CELL COUNTS
orng = '#FF8800';  % orange color

f = figure
semilogy(CellCounts(:,1),CellCounts(:,2),'linewidth',2);
hold on
semilogy(CellCounts(:,1),CellCounts(:,3),'linewidth',2);
semilogy(CellCounts(:,1),CellCounts(:,4),'linewidth',2);
semilogy(CellCounts(:,1),CellCounts(:,5),'linewidth',2);
% add errors to pre- and post- recognition
upper_pre = CellCounts(:,4)+StdErrors(:,4);
lower_pre = max((CellCounts(:,4)-StdErrors(:,4)),eps); %avoid negative values
fill_time = [CellCounts(:,1);flipud(CellCounts(:,1))];
fill_border_pre = [upper_pre;flipud(lower_pre)];
fill(fill_time,fill_border_pre,'y','FaceAlpha',0.1, ...
    'EdgeColor',orng)
upper_post = CellCounts(:,5)+StdErrors(:,5);
lower_post = max((CellCounts(:,5)-StdErrors(:,5)),eps);
fill_border_post = [upper_post;flipud(lower_post)];
fill(fill_time,fill_border_post,'m','FaceAlpha',0.05, ...
    'EdgeColor','m')
axis([0, max(CellCounts(:,1)), 1, 1e5])
xlabel('Time','interpreter','latex','fontsize',20)
ylabel('Cell Counts','interpreter','latex','fontsize',20)
legend('Low MHC tumor','High MHC tumor','CTL pre-recognition',...
    'CTL post-recognition',...
    'Location','southeast', ...
    'fontsize',12, 'interpreter','latex')
title('Average Cell Counts: CTL only','interpreter','latex','FontSize',24)
%set(gcf, 'PaperPositionMode', 'auto','PaperOrientation','landscape');
%print(gcf, ['CellCounts',RunName,'.pdf'], '-dpdf', '-fillpage');
exportgraphics(f,[path,'AverageCellCounts_wstd','.pdf'],'ContentType','vector')
%close(f)

% Plot Hi MHC - Low MHC
% f2 = figure
% plot(CellCounts(:,1),(CellCounts(:,3)-CellCounts(:,2)),'linewidth',2);
% hold on
% plot([0,300],[0,0],'k')
% xlabel('Time','interpreter','latex','fontsize',20)
% ylabel('Cell Counts','interpreter','latex','fontsize',20)
% title('High MHC minus Low MHC Tumor Cells: CTL only','interpreter','latex','FontSize',20)
% exportgraphics(f2,['HiMinusLowMHC',RunName,'.pdf'],'ContentType','vector')

