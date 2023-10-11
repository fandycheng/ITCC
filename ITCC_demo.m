%% Demo of ITCC detection and tracking alogrithm

%==================================================================================================
%/       Author: Franklin Cheng (fandycheng@ust.hk)
%/  Last Update: October 10, 2023
%==================================================================================================
    
Your_Source_Codes_Path = '/home/tfchengac/';   %/ Adapt this path to where you have placed the source codes
addpath(genpath([Your_Source_Codes_Path, '/ITCC_Source_Codes_github']));

%/ Read OLR data (ERA5) 
OLR_name      = 'OLR'; 
masterfolder  = '/disk/r059/tfchengac/Trimonsoon_ITCC/prcssd_data_4plotting/';  %/ Adapt this path to where you have placed the data
OLR_filepath  = [masterfolder, sprintf('%s_pentad_clim_1979-2020.mat', OLR_name)];
load(OLR_filepath, 'S');
flds          = fieldnames(S);
dataset       = [];
for f = 1:length(flds)
    dataset.(OLR_name).(flds{f}) = S.(flds{f});
end
clear S

%/ Settings
OLR_thres            = 220;                     %/ Maximum OLR threshold (in W m^-2)
track_mode           = 'area_ellip';            %/ Track the ITCC area by fitting an ellipse on it
lon_IOWP             = [40, 210];               %/ ITCC detection range (in lon)
lat_IOWP             = [-25, 25];               %/ ITCC detection range (in lat)
seek_N               = 20;                      %/ Find N objects that satisfy the OLR thres
largest_n            = 20;                      %/ Select the n largest patches from the N objects
min_Area_thres       = 6e5;                     %/ Minimum area threshold of the fitted ellipse (in km^2)
min_Majlen_geo_thres = 1500;                    %/ Minimum major-axis length of the fitted ellipse (in km)
remove_highland      = 1;                       %/ [0 or 1] Turn it on to avoid detection ITCC over highland (>3000 m)
plot_ITCC            = 0;                       %/ [0 or 1] Outline ITCC on the map requires M_Map Matlab package. Set 0 by default.
pinpoint_grids       = 0;                       %/ [0 or 1] Pinpoint the grids of the ITCC/LCCs on the map
derive_attrs         = 1;                       %/ [0 or 1] Output 'LCC_attrs' 
t_loop               = 1:73;                    %/ 73 calendar pentads (climatology)
LCC_attrs            = cell(length(t_loop), 6); %/ Attributes: centroid, bndry_data, BW_matrix, length scale, Ellip Area, Total Grid Area
CISO_attrs           = cell(length(t_loop), 6); %/ Attributes: centroid, bndry_data, BW_matrix, length scale, Ellip Area, Total Grid Area
ITCC_attrs           = [];

lon    = dataset.(OLR_name).lon; 
lat    = dataset.(OLR_name).lat;
    
%/ Run ITCC_detection.m (to generate 'LCC_attrs', a list of attributes of Large-scale convective cells (LCCs))
for t = t_loop
    fprintf('*** t = %d ... ***\n', t)
    OLR_2D = dataset.(OLR_name).pentad_clim(:,:,t); %/ 2D OLR data
    
    LCC_attrs = ITCC_detection('OLR_2D',  OLR_2D,  'lon',  lon, 'lat',  lat,...
                               'domain_lon_range', lon_IOWP,  'domain_lat_range',  lat_IOWP,...
                               'tot_timestamp', length(t_loop), 'timestamp', t, 'track_mode', track_mode,...
                               'derive_attrs', derive_attrs, 'LCC_attrs', LCC_attrs, 'ITCC_attrs', ITCC_attrs,...
                               'OLR_thres', OLR_thres, 'seek_N', seek_N, 'largest_n', largest_n, 'min_Area_thres', min_Area_thres, 'min_Majlen_geo_thres', min_Majlen_geo_thres,...
                               'plot_ITCC', plot_ITCC, 'remove_highland', remove_highland, 'pinpoint_grids', pinpoint_grids);
end

%/ Run ITCC_tracking.m (with the full list of 'LCC_attrs')
ITCC_attrs = ITCC_tracking('lon', lon, 'lat', lat, 'LCC_attrs', LCC_attrs);

Centr      = ITCC_attrs{:, 1};      %/ Centroid of the ellipse(s)
bndry_data = ITCC_attrs{:, 2};      %/ Boundary of the ellipse(s)
output_BW  = ITCC_attrs{:, 3};      %/ Other attr: The detected LCC grids (logical matrix).
MajLen_geo = ITCC_attrs{:, 4};      %/ Other attr: Elliptical length scale (km)
Ellip_Area = ITCC_attrs{:, 5};      %/ Other attr: Elliptical area (km^2)
Grid_Area  = ITCC_attrs{:, 6};      %/ Other attr: Total grid area (km^2)

%% Below is just for double-checking, please ignore it.
% %/  - Load ITCC_attrs
% masterfolder    = '/disk/r059/tfchengac/Trimonsoon_ITCC/';
% data_folder     = strcat(masterfolder, 'prcssd_data_4plotting/');
% ITCC_attrs_filename = [data_folder,...
%     sprintf('ITCC_attrs_%s_pentad_clim_P1-73_area_ellip_IOWP%s_40to210E_-25to25N_10Largest.mat', 'OLR', '_v7')];
% S = load(ITCC_attrs_filename, 'ITCC_attrs');
% ITCC_attrs_ori = S.ITCC_attrs;
% 
% isequal([ITCC_attrs{:,1}], [ITCC_attrs_ori{:,1}])
% isequal([ITCC_attrs{:,4:5}], [ITCC_attrs_ori{:,4:5}])
% 
% bndry = [ITCC_attrs{:,2}];
% bndry_ori = [ITCC_attrs_ori{:,2}];
% isequal([bndry{:}], [bndry_ori{:}])
% a = [bndry{:}] - [bndry_ori{:}];
% max(abs(a), [], 'all')
% 
% isequal([ITCC_attrs{:,3}], [ITCC_attrs_ori{:,3}])