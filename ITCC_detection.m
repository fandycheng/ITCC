%%
function LCC_attrs = ITCC_detection(varargin)
  
    pnames = {    'OLR_2D',          'lon',            'lat',     'domain_lon_range',     'domain_lat_range',...
           'tot_timestamp',    'timestamp',     'track_mode',         'derive_attrs',            'LCC_attrs',       'ITCC_attrs', 'show_ITCC_only', 'Dlon', 'Dlat',...
               'OLR_thres',       'seek_N',      'largest_n',       'min_Area_thres', 'min_Majlen_geo_thres',    'min_MSE_thres',...
                  'MSE_2D',      'MSE_lon',        'MSE_lat',      'remove_highland',       'pinpoint_grids',...
               'plot_ITCC',    'recursive',   'track_colmap',           'markersize',               'linewi',           'linest', };
    
    dflts  = {          [],             [],               [],                     [],                     [],...
                        [],             [],               [],                      0,                     [],                 [],                 0,        [],      [],...
                       220,             20,               [],                     [],                     [],                 [],...
                        [],             [],               [],                      0,                      0,...
                         1,              0,               [],                      6,                      3,                '-' };
    
    [               OLR_2D,            lon,              lat,       domain_lon_range,       domain_lat_range,...
             tot_timestamp,      timestamp,       track_mode,           derive_attrs,              LCC_attrs,         ITCC_attrs,  show_ITCC_only,    Dlon,  Dlat,...
                 OLR_thres,         seek_N,        largest_n,         min_Area_thres,   min_Majlen_geo_thres,      min_MSE_thres,...
                    MSE_2D,        MSE_lon,          MSE_lat,        remove_highland,         pinpoint_grids,...
                 plot_ITCC,      recursive,     track_colmap,             markersize,                 linewi,             linest, ] = ...
                                            internal.stats.parseArgs(pnames, dflts, varargin{:});
    %%
    %==================================================================================================
    %/       Author: Franklin Cheng (fandycheng@ust.hk)
    %/  Last Update: October 10, 2023
    %/
    %/  DESCRIPTION:
    %/          'recursive':     1]: then it loops to draw all the bndry or dots till the given timestamp. 
    %/                           0]: draw current bndry or dot
    %/      
    %/          'OLR_thres':     By default, we filter OLR < 220 W m**-2. 
    %                            It indicates tropical convection (Webster et al. 1998; Arkin and Meisner 1987)
    %==================================================================================================
    
    if isempty(track_mode)   warning('[ITCC_detection]: No valid input of track_mode. Skip tracking minOLR.');    return;   end

    if recursive
        t_list        = 1:timestamp;
        track_colmap  = hsv(tot_timestamp);                  %/ Always use hsv as it has a closed-loop coloring.
    else
        t_list        = timestamp;
        if isempty(track_colmap)
            track_colmap  = repmat([0 0 0], largest_n-1, 1); %/ All in grey.
        end
    end
    
    if derive_attrs
        if unique(diff(lat)) < 0     %/ <- mind this!
            lat    = flip(lat, 1);
            OLR_2D = flip(OLR_2D, 2);
        end

        [lon_2D, lat_2D] = meshgrid(lon, lat);
        lon_2D = lon_2D';  lat_2D = lat_2D'; 

        filtered_data_2D = OLR_2D;  %/ copy.

        %/ [Step 1] Filter the domain if given.
        if ~isempty(domain_lon_range) && ~isempty(domain_lat_range)
            filtered_data_2D = OLR_2D; 

            cond_lon_2D = (lon_2D >= domain_lon_range(1) & lon_2D <= domain_lon_range(2));
            cond_lat_2D = (lat_2D >= domain_lat_range(1) & lat_2D <= domain_lat_range(2));

            filtered_data_2D(~cond_lon_2D | ~cond_lat_2D) = nan;  %/ fulfilling either logical array will set nan. == mask data
        else
            domain_lon_range = [min(lon), max(lon)];
            domain_lat_range = [min(lat), max(lat)];
        end

        %/ [Step 2] Remove highland grids - by default
        if remove_highland
            %/ 2D interpolation from high-res topo to 1x1 topo that fits our data.
            m_proj('Miller Cylindrical','longitudes',[0 360], 'latitudes', [-90 90]);
            [topo, topo_lon_2D, topo_lat_2D] = m_tbase([-179.5, 179.5, -90, 90]); %/ 1x1 topo  m_elev.    
            topo = topo';               %/ make sure in (lon,lat) dim

            topo_lon = topo_lon_2D(1,:)';
            topo_lat = topo_lat_2D(:,1);
            if diff(topo_lat(1:2)) < 0
                topo_lat = flip(topo_lat, 1);
                topo = flip(topo, 2);
            end

            ind_nve = find(topo_lon < 0);
            ind_pve = find(topo_lon >= 0);

            %/ Shift grids from (-180, 180) to (0, 360]
            topo_lon = topo_lon([ind_pve;ind_nve]);
            topo_lon(topo_lon < 0) = topo_lon(topo_lon < 0) + 360;
            topo = topo([ind_pve;ind_nve],:);
            [topo_lon_2D, topo_lat_2D] = meshgrid(topo_lon, topo_lat);
            topo_interp = interp2(topo_lon_2D, topo_lat_2D, topo', lon_2D', lat_2D', 'linear')';  %/ make sure 1st dim is lat!
            
            %/ Mask out plateau (>3000 m)
            cond_rm_topo = (topo_interp > 3000);
            filtered_data_2D(cond_rm_topo) = nan;
        end

        %/ [Step 3] Fitler OLR < OLR_thres W m^-2
        filtered_data_2D(filtered_data_2D >= OLR_thres) = nan;
        
        %/ Check if any grid cells left after filtering
        if isempty(find(~isnan(filtered_data_2D)))
           warning('!!! [ITCC_detection] No grid cells to track after filtering. Skip tracking. !!!')  
           %/ NOTE: no need to set empty LCC_attrs since they are the inputs.
           return;
        end
        
        %/ [Step 4] Find the largest object(s) (by grids or by fitting an ellipse)
        if isequal(track_mode, 'area')
            error('%s track_mode has been abandoned!', track_mode);
%             %/ Turn into binary data 
%             logical_2D = filtered_data_2D;
%             logical_2D(isnan(logical_2D)) = 0;
%             logical_2D = logical(logical_2D);
% 
%             bndry_data = get_bndry_from_logi('logical_2D', logical_2D, 'bndry_lon', lon, 'bndry_lat', lat, 'outputcell', 1,...
%                                              'map_lon_lower', domain_lon_range(1), 'map_lon_upper', domain_lon_range(2), 'map_lat_lower', domain_lat_range(1), 'map_lat_upper', domain_lat_range(2),...
%                                              'glb_data_mode', 0, 'hole_mode', 'noholes', 'draw_rings', 0);
%             bndry_len        = cellfun(@(x) length(x), bndry_data);
%             [~, ind_largest] = max(bndry_len);
% 
%             bndry_data_list{timestamp} =  bndry_data{ind_largest};

        elseif isequal(track_mode, 'area_ellip')
            field_data = filtered_data_2D;
            [Centr, ~, ~, ~, ~, bndry_data, Ellip_Area, Grid_Area, MajLen_geo, output_BW] = fit_ellipse_BW_geo('field_data', field_data,  'lon', lon, 'lat', lat,...
                                                                                   'seek_N', seek_N, 'largest_n', largest_n, 'draw_ellipse', 0,  'create_fig', 0);
            
            %/ [Step 5] Remove "small-size" LCCs (if required).
            %/          if both thresholds are not empty, remove the ones that satisfy none of the criteria.
            if ~isempty(min_Area_thres) && ~isempty(min_Majlen_geo_thres) 
                ind = find(Ellip_Area < min_Area_thres & MajLen_geo < min_Majlen_geo_thres);
                Centr(ind,:)      = [];
                bndry_data(ind,:) = [];
                Ellip_Area(ind,:) = [];
                Grid_Area(ind,:)  = [];
                MajLen_geo(ind,:) = [];
                output_BW(:,:,ind)= [];
                
            elseif ~isempty(min_Area_thres)
                ind = find(Ellip_Area < min_Area_thres);
                Centr(ind,:)      = [];
                bndry_data(ind,:) = [];
                Ellip_Area(ind,:) = [];
                Grid_Area(ind,:)  = [];
                MajLen_geo(ind,:) = [];
                output_BW(:,:,ind)= [];
                
            elseif ~isempty(min_Majlen_geo_thres)
                ind = find(MajLen_geo < min_Majlen_geo_thres);
                Centr(ind,:)      = [];
                bndry_data(ind,:) = [];
                Ellip_Area(ind,:) = [];
                Grid_Area(ind,:)  = [];
                MajLen_geo(ind,:) = [];
                output_BW(:,:,ind)= [];
            end
            
            %/ [Step 6, deprecated] Remove LCCs with regional mean <MSE> lower than min_MSE_thres 
            if ~isempty(min_MSE_thres)
                ind = [];
                for i = 1:length(bndry_data)
                    bndry_data_bc = bndry_data{i};
                    MSE_2D_masked = masked_by_bndry('X', MSE_2D, 'lon', MSE_lon, 'lat', MSE_lat, 'bndry_data', bndry_data_bc);
                    
                    A = calc_grid_area('lon', MSE_lon, 'lat', MSE_lat);
                    cond = isnan(MSE_2D_masked);
                    A(cond) = nan;
                    
                    aw_reg_mean_MSE = nansum(MSE_2D_masked.*A, 'all')./nansum(A, 'all');
                    if aw_reg_mean_MSE < min_MSE_thres
                        ind = [ind; i];
                        fprintf('LCC with an area-weight regional mean <MSE> = %g (< %g) is detected, and has been removed!\n', aw_reg_mean_MSE, min_MSE_thres)
                    end
                end
                Centr(ind,:)      = [];
                bndry_data(ind,:) = [];
                Ellip_Area(ind,:) = [];
                Grid_Area(ind,:)  = [];
                MajLen_geo(ind,:) = [];
                output_BW(:,:,ind)= [];
            end
            
            %/ Appending to the input list of points and boundaries
            LCC_attrs{timestamp, 1} = Centr;              %/ Centroid of the ellipse(s)
            LCC_attrs{timestamp, 2} = bndry_data;         %/ Boundary of the ellipse(s)
            LCC_attrs{timestamp, 3} = output_BW;          %/ Other attr: The detected LCC grids (logical matrix).
            LCC_attrs{timestamp, 4} = MajLen_geo;         %/ Other attr: Elliptical length scale (km)
            LCC_attrs{timestamp, 5} = Ellip_Area;         %/ Other attr: Elliptical area (km^2)
            LCC_attrs{timestamp, 6} = Grid_Area;          %/ Other attr: Total grid area (km^2)
        end
    else
        fprintf('!!! Detected that derive_attrs == 0. No LCC_attrs will be derived. !!!\n');
    end
    
    %/======= Plotting =======%
    if plot_ITCC 
        %/ Tracked area 
        if isequal(track_mode, 'area') || isequal(track_mode, 'area_ellip')
        bndry_alpha = 1;
        %/ If recursive = 1, then it loops to draw all the cells till the given timestamp
        for t = t_list
            %/ Check if we want to show all LCCs or only the ITCC.
            if show_ITCC_only == 0
                n_list = 1:length(LCC_attrs{t,2});
                %/ If n_largest > 1, there could be > 1 objected fitted at each timestamp.
                for n = n_list
                    if recursive == 0
                        color = track_colmap(n,:);
                    else
                        color = track_colmap(t,:);
                    end

                    h1 = m_line(LCC_attrs{t,2}{n,1}(:,1), LCC_attrs{t,2}{n,1}(:,2), 'marker','none', 'color', color,...
                           'linewi', linewi, 'linest', linest, 'markersize', markersize, 'markerfacecolor', 'none');
                    h1.Color(4) = bndry_alpha;
                    hold on;

                    if isequal(track_mode, 'area_ellip')
                        %/ Draw centroids
                        m_line(LCC_attrs{t,1}(n, 1), LCC_attrs{t,1}(n, 2), 'marker','o', 'color', 'none', 'linewi', linewi,...
                           'linest','none','markersize', markersize*4, 'markerfacecolor', color);
                    end

                    if recursive == 0 && pinpoint_grids
                        output_BW = LCC_attrs{timestamp, 3}; 
                        cond = logical(output_BW(:,:,n));

                        if unique(diff(lat)) > 0     %/ NOTE: Since the output_BW is generated by descending lat.
                            lat     = flip(lat, 1);
                        end
                        [lon_2D, lat_2D] = meshgrid(lon, lat);
                        lon_2D = lon_2D';  lat_2D = lat_2D'; 

                        BW_pt_lon = lon_2D(cond);
                        BW_pt_lat = lat_2D(cond);

                        m_line(BW_pt_lon, BW_pt_lat, 'marker','.', 'color', color, 'linewi', linewi,...
                               'linest','none','markersize', markersize*2, 'markerfacecolor', color);
                    end
                end
            end
            
            if ~isempty(ITCC_attrs) && recursive == 0
                fprintf('*** Non-empty ITCC_attrs is found. The ITCC will be outlined in red. ***\n');
                
                color = [255 0 0]./255;  
                h1 = m_line(ITCC_attrs{t,2}{:}(:,1), ITCC_attrs{t,2}{:}(:,2), 'marker','none', 'color', color,...
                       'linewi', linewi, 'linest', linest, 'markersize', markersize, 'markerfacecolor', 'none');
                h1.Color(4) = bndry_alpha;
                hold on;
                
                if isequal(track_mode, 'area_ellip')
                    %/ Draw centroids
                    m_line(ITCC_attrs{t,1}(1), ITCC_attrs{t,1}(2), 'marker','o', 'color', 'none', 'linewi', linewi,...
                       'linest','none','markersize', markersize*4, 'markerfacecolor', color);
                end

                if recursive == 0 && pinpoint_grids
                    output_BW = ITCC_attrs{timestamp, 3}; 
                    cond = logical(output_BW);
                    if unique(diff(lat)) > 0     %/ NOTE: Since the output_BW is generated by descending lat.
                        lat     = flip(lat, 1);
                    end
                    [lon_2D, lat_2D] = meshgrid(lon, lat);
                    lon_2D = lon_2D';  lat_2D = lat_2D'; 
                    BW_pt_lon = lon_2D(cond);
                    BW_pt_lat = lat_2D(cond);
                    m_line(BW_pt_lon, BW_pt_lat, 'marker','.', 'color', color, 'linewi', linewi,...
                           'linest','none','markersize', markersize*2, 'markerfacecolor', color);
                end
                
                %/ Draw the box region around the ITCC's centroid, where we did the cross-section following it.
                if show_ITCC_only && ~isempty(Dlon) && ~isempty(Dlat)
                    ITCC_centr = ITCC_attrs{t,1};
                    bndry_cross_section = [ITCC_centr(1)-Dlon, ITCC_centr(2)-Dlat;
                                           ITCC_centr(1)-Dlon, ITCC_centr(2)+Dlat;
                                           ITCC_centr(1)+Dlon, ITCC_centr(2)+Dlat;
                                           ITCC_centr(1)+Dlon, ITCC_centr(2)-Dlat;
                                           ITCC_centr(1)-Dlon, ITCC_centr(2)-Dlat;];
                    m_line(bndry_cross_section(:,1), bndry_cross_section(:,2), 'marker','none', 'color', 'y', 'linewi', linewi,...
                           'linest','--','markersize', markersize*2, 'markerfacecolor', 'none');
                end
            end
        end
        drawnow; pause(0.05);
        end
    
        %/ Tracked line
        if isequal(track_mode, 'line')
            if derive_attrs
                [ind_lon_minOLR, ind_lat_minOLR] = find(filtered_data_2D == min(reshape(filtered_data_2D, [], 1)));

                %/ Restore to the true lon, lat
                lon_minOLR = lon(ind_lon_minOLR);  
                lat_minOLR = lat(ind_lat_minOLR);

                %/ Appending to the input list of points
                LCC_attrs{timestamp,1} = [lon_minOLR, lat_minOLR];
            end

            %/ If recursive = 1, then it loops to draw a color-gradient line till the given timestamp
            for t = t_list
                %/ Draw nodes
                m_line(LCC_attrs{t,1}(:, 1), LCC_attrs{t,1}(:, 2), 'marker','o', 'color', 'none', 'linewi', linewi,...
                       'linest','none','markersize', markersize*4, 'markerfacecolor', track_colmap(t,:));

                %/ Draw tracking segments from the beginning when nn > 1
                if t ~= 1
                    listOfpts_bc = cat(1, LCC_attrs{t-1:t,1});

                    m_line(listOfpts_bc(:,1), listOfpts_bc(:,2), 'marker','none', 'color', track_colmap(t,:), 'linewi', linewi,...
                           'linest', linest,'markersize', markersize, 'markerfacecolor', 'none');
                end

                if t == tot_timestamp
                    %/ Close the loop
                    listOfpts_bc = cat(1, LCC_attrs{[t, 1],1});

                    m_line(listOfpts_bc(:, 1), listOfpts_bc(:, 2), 'marker','none', 'color', track_colmap(1,:), 'linewi', linewi,...
                           'linest', linest,'markersize', markersize, 'markerfacecolor', 'none');
                end
                hold on;
            end
            drawnow; pause(.5) %/ Important! Sometimes it has not draw the ellipse before the next fig shows up.
        end
    end
end