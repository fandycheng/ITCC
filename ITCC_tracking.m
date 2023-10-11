%%
function ITCC_attrs = ITCC_tracking(varargin)

    pnames = {'lon', 'lat', 'LCC_attrs'};
    dflts  = cell(length(pnames),1);
    
    [          lon,   lat,   LCC_attrs] = ...
                         internal.stats.parseArgs(pnames, dflts, varargin{:});
    %%
    %==================================================================================================
    %/       Author: Franklin Cheng (fandycheng@ust.hk)
    %/  Last Update: October 10, 2023
    %/
    %/ DESCRIPTION: The tracking of ITCC will begin from P18 in the
    %/              Indochina-Sumatra region ('IDC-Sumatra-new') by first
    %/              seeking the largest one, and then track it using
    %/              area-overlapping (AOL) method.
    %/ INPUT:
    %/          'LCC_attrs': The attributes of the fitted LCCs. (from ITCC_detection.m)
    %/ OUTPUT:
    %/         'ITCC_attrs': The attributes of the ITCC tracked.
    %==================================================================================================
    
    LCC_centr_list = LCC_attrs(:,1);
    LCC_bndry_list = LCC_attrs(:,2);
    ind_ITCC       = nan(length(LCC_bndry_list), 1);
    A              = calc_grid_area('lon',lon,'lat',lat);

    %/ The full cycle of ITCC begins in P18, at least based on climatology.
    start_pentad = 18;  
    ITCC_cycle = [start_pentad:length(LCC_bndry_list), 1:start_pentad-1];

    %/ Read the IDC-Sumatra (IS) box regions (springtime fast transition region)
    dom_name = 'IDC-Sumatra-new'; zm_or_mm = 1;
    [dom_box_lon, dom_box_lat, ~, ~] = hovmoller_box_reg('dom_name', dom_name, 'zm_or_mm', zm_or_mm); 
    dom_bndry = [dom_box_lon(1) dom_box_lat(1); dom_box_lon(1) dom_box_lat(2); dom_box_lon(2) dom_box_lat(2);...
                 dom_box_lon(2) dom_box_lat(1); dom_box_lon(1) dom_box_lat(1);];
    dom_IS = ones(length(lon),length(lat));
    dom_IS = masked_by_bndry('X', dom_IS, 'lon', lon, 'lat', lat, 'bndry_data', dom_bndry, 'output_logical', 1);

    %/ The AOL algorithm
    prev_t = nan;
    for i = 1:length(ITCC_cycle)
        curr_t = ITCC_cycle(i);
        
        %/ Current LCCs (t)
        curr_bndry_list = LCC_bndry_list{curr_t};
        n_LCCs          = length(curr_bndry_list);
        curr_LCCs       = nan(length(lon),length(lat),n_LCCs);
        for j = 1:n_LCCs
            X = ones(length(lon),length(lat));
            curr_LCCs(:,:,j) = masked_by_bndry('X', X, 'lon', lon, 'lat', lat, 'bndry_data', curr_bndry_list{j});
        end
        
        %/ First, start with the IS region (springtime fast transition region)
        if curr_t == ITCC_cycle(1)  flag_slct_largest = 1;  end

        if flag_slct_largest
            fprintf('*** curr_t = %2d [Searching ITCC in IDC-Sumatra] ***\n', curr_t)
            %/ If no LCC at curr_t = ITCC_cycle(i), will skip to the next time step.
            if isempty(LCC_bndry_list{curr_t}) 
                prev_t = curr_t; %/ update the prev time.
                continue;
            end
            
            %/ Find the largest ellipse that encompasses the IS region.
            [B, I]            = max(squeeze(nansum(curr_LCCs.*dom_IS, [1,2])));
            ind_ITCC(curr_t)  = I;
            
            %/ If the ellipse's centroid is outside within the IS region, 
            %/ continue to search for the largest one within the region.
            [in,~] = inpolygon(LCC_centr_list{curr_t}(I,1), LCC_centr_list{curr_t}(I,2), dom_bndry(:,1), dom_bndry(:,2));  
%             [in,~] = inpoly2([LCC_centr_list{curr_t}(I,:)], dom_bndry);  %/ inpoly2 is 600xx faster, but it's error-prone.
            if in 
                flag_slct_largest = 0; %/ turn it off (we've now successfully located the ITCC over the IS region)
            end
%             [LCC_centr_list{curr_t}(I,:)]
            
        %/ Then, perform ITCC tracking based on AOL.
        else
            fprintf('*** curr_t = %2d [Tracking ITCC] ***\n', curr_t)
            %/ If no LCC is found at the given time, 
            %/ we continue by choosing the largest one at the next timestep.
            if isempty(LCC_bndry_list{curr_t}) 
                flag_slct_largest = 1;
                prev_t = curr_t; %/ update the prev time.
                continue;
            end

            %/ Previous ITCC (t-1)
            X               = ones(length(lon),length(lat));
            prev_bndry_list = LCC_bndry_list{prev_t};
            prev_ITCC_bndry = prev_bndry_list{ind_ITCC(prev_t)};
            prev_ITCC       = masked_by_bndry('X', X, 'lon', lon, 'lat', lat, 'bndry_data', prev_ITCC_bndry);
            
            %/ Find the intersected grids, then sum their areas up.
            AOL = prev_ITCC.*curr_LCCs.*A;  %/ make use of the property: 1*nan = nan
            AOL = squeeze(nansum(AOL,[1,2]));

            %/ If no overlapped area at all, we assign the ITCC as the largest one.
            [~,I] = max(AOL);
            ind_ITCC(curr_t) = I;
        end
        prev_t = curr_t; %/ update prev_t.
    end

    %/ Store as ITCC_attrs (for the later convenience).
    ITCC_attrs = cell(size(LCC_attrs));
    for k = 1:size(LCC_attrs,2)
        A = cell(length(ind_ITCC),1);
        for i = 1:length(ind_ITCC)
            if k == 3 %/ corresponds to the BW_matrix (lon,lat,n_DCC)
                A{i} = squeeze(LCC_attrs{i,k}(:,:,ind_ITCC(i)));
            else
                A{i} = LCC_attrs{i,k}(ind_ITCC(i),:);
            end
        end
        ITCC_attrs(:,k) = A;
    end
    
    fprintf('!!! ITCC Tracking has completed successfully !!!\n')
end