function MakeHexGridColorRGB_Fig10(X,Y,M1,M2,M3,M4,M5,M6,M7,M8,TimePoint,ColorMat)
% Draw an array representing the cell types in a hexagonal grid
% The vectors X and Y give the centers of each grid element
% ColorMat is a 3 by 3 array with RGB colors of 
% Tumor, CTL and NK in each row, respectively
%% compare cell types in each grid element
% M1 tumor lowMHC
% M2 tumor hiMHC
% M3 CTL pre recognition of a tumor
% M4 NK pre recognition of a tumor
% M5 CTL using fasL to kill tumor
% M6 CTL using perforin to kill tumor
% M7 NK using fasL to kill tumor
% M8 NK using perforin to kill tumor
Cmax = max(M3+M5+M6);
Nmax = max(M4+M7+M8);
Emax = max(Cmax,Nmax);
T = min(.5*M1+M2,1);
%-%-% SCALE BY MAX NUMBER IN MATRIX 
% if Cmax == 0
%     CTL = (M3+M5+M6);  % this should be zero in this case
% else
%     CTL = (M3+M5+M6)/Cmax;
% end
% if Nmax == 0
%     NK = M4+M7+M8;
% else
%     NK = (M4+M7+M8)/Nmax;
% end
%-%-% SCALE BY MAX EFFECTOR CELL NUMBER
if Emax == 0
    NK = 0;
    CTL = 0;
else
    NK = log(M4+M7+M8)/log(Emax); % use log to get more color at low end
    CTL = log(M3+M5+M6)/log(Emax);
end;
%-%-% DO NOT SCALE
%CTL = max((M3+M5+M6)/15,1);
%NK = max((M4+M7+M8)/15,1);
%-%-% Use relative size of effector cells
% Emax = max(max(M3+M5+M6+M4+M7+M8));
% if Emax ==0
%     CTL = (M3+M5+M6); 
%     NK = M4+M7+M8;
% else
%      CTL = (M3+M5+M6)/Emax;
%      NK = (M4+M7+M8)/Emax;
% end
%T = 1 - min(.5*M1+M2,1); % tumor color lightener
%CTL = 1 - min((M3+M5+M6)/Cmax,1);  % CTL  color lightener
%NK = 1 - min((M4+M7+M8)/Nmax,1);   % NK  color lightener
Tumor_col = ColorMat(1,:);
CTL_col = ColorMat(2,:);
NK_col = ColorMat(3,:);
% Now the color matrices should all be between 0 and 1, with 1 the lightest

%% color will use RGB values based on values given for T, CTL and NK
f = figure('visible','off'); % create a new figure but do not display it
s = (X(2)-X(1))/2;
% vertex offsets, going counterclockwise from the top of the hexagon
xoff = [0,-s,-s,0,s,s];
yoff = [3/2*s,s/2,-s/2,-3/2*s,-s/2,s/2];
hold on
for i=1:length(X)
    if mod(i,1000)==1
        i
    end
    % generate vertices
    xv = X(i)+xoff;
    yv = Y(i)+yoff;
    % Tcol = Tumor_col + T(i)*(1-Tumor_col);
    % CTLcol = CTL_col + CTL(i)*(1-CTL_col);
    % NKcol = NK_col + NK(i)*(1-NK_col);
    % USE CELL TYPE WITH MAX NUMBER:
    [Cell,Ind] = max([T(i),CTL(i),NK(i)]);
    col = ColorMat(Ind,:) + (1-Cell)*(1-ColorMat(Ind,:));
    % USE WEIGHTED SUM OF EFFECTOR CELL COLORS
    %SumCell = CTL(i)+NK(i);
    %if SumCell ==0
    %   EFF_col = [1,1,1];  % background is white
    %else
    %    TotCol = 1/SumCell*(T(i)*Tumor_col + CTL(i)*CTL_col + NK(i)*NK_col);
    %    col = TotCol + (1-min(SumCell,1))*(1-TotCol);
    %end;
    fill(xv,yv,col,'linestyle','none');
    %shading('flat')   
end
hold off
axis('square')
axis('off')
title(['Time Step ',num2str(TimePoint)])
print(gcf, ['TimeStep',num2str(TimePoint),'.pdf'], '-dpdf', '-fillpage');
close(f)
    