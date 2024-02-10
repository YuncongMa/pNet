function Data=fHCP_To_32K(Data0,Index_MW,varargin)
% By Yuncong Ma, Jan. 20, 2020
% Data=fHCP_To_32K(Data,Index_MW)
% Map Data from HCP to FS_32K 
% Data=fHCP_To_32K(Data,Index_MW,{'Initialization',0})
% Initialization could be 0 or NaN for values in middle wall

Options.Initialization=0;
Options=fOption('fHCP_To_32K',Options,varargin);
if isempty(Options)
    return;
end

if Options.Initialization==0
    Data=zeros([length(Index_MW),size(Data0,2)]);
elseif isnan(Options.Initialization)
    Data=nan([length(Index_MW),size(Data0,2)]);
else
    error('Error in fHCP_To_32K: Optional settings Initialization could only be 0 or NaN.\n');
end
Data(Index_MW==0,:)=Data0;

end