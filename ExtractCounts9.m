% extract cell counts for visualization
% Assume that Data is a "position" matrix generated from the code
% e.g.: load ctlPos.dat
%       Data = ctlPos
%
% Modified April 12, 2017
% The matrix Data contains 18 columns - see biofunctions file
% 1 - lowMHC tumor
% 2 - hiMHC tumor
% 3:4, 7:8 - CTL not yet recognized a tumor
% 11:12, 15:16 - NK not yet recognized a tumor
% 5:6 - CTL, fasL
% 13:14 - NK, fasL
% 9:10 - CTL, perforin
% 17:18 - NK, perforin

% M1 tumor lowMHC
% M2 tumor hiMHC
% M3 CTL pre recognition of a tumor
% M4 NK pre recognition of a tumor
% M5 CTL using fasL to kill tumor
% M6 CTL using perforin to kill tumor
% M7 NK using fasL to kill tumor
% M8 NK using perforin to kill tumor

M1 = Data(:,1);
M2 = Data(:,2);
M3 = sum(Data(:,[3:4,7:8]),2);
M4 = sum(Data(:,[11:12,15:16]),2);
M5 = sum(Data(:,5:6),2);
M6 = sum(Data(:,9:10),2);
M7 = sum(Data(:,13:14),2);
M8 = sum(Data(:,17:18),2);
n = sqrt(length(Data(:,1))); % size of the computational grid

