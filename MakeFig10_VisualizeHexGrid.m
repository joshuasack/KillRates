% Version 01/13/2021
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
NK_col = colorblind(6,:);  % orange
CTL_col = colorblind(9,:); % semi-darkblue
Tum_col = colorblind(2,:);  % red
ColorMat = [Tum_col;CTL_col;NK_col];
%Bckgrnd = colorblind(3,:); % yellow background 
%Bckgrnd = [1,1,1]; % white background
%% Set path to load data
% Set the path where the data is
path = 'InVivo_grow/data';
%% loop through time points

% load the data at the specified time point
TimePoints = 0:50:200;
for k=1:length(TimePoints)
    TimePoint = TimePoints(k);
    % For in vivo runs
    Data = load([path,'/inVivoTime',num2str(TimePoint),'Pos.dat']);
    % For in vitro runs
    %Data = load([path,'/inVitroPosExp1Time',num2str(TimePoint),'.dat']);
    % extract 8 matrices: M1 through M8 with various types of cell counts
    % also returns n, the size of the matrices (which are nxn)
    ExtractCounts9
    % get the  xy-coordinates for the centers of the hexagons
    [X,Y] = gethexind(n);
    MakeHexGridColorRGB_Fig10(X,Y,M1,M2,M3,M4,M5,M6,M7,M8,TimePoint,ColorMat);

end
