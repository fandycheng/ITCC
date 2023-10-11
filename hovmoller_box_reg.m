%%
function [dom_box_lon, dom_box_lat, border_lines, time_lines] = hovmoller_box_reg(varargin)

    pnames = {'project_name', 'dom_name', 'zm_or_mm', 'select_field'};
    dflts  = {        'ITCC',         [],         [],       'pentad'};
    
    [           project_name,   dom_name,   zm_or_mm,  select_field ] = ...
                         internal.stats.parseArgs(pnames, dflts, varargin{:});
    %%
    %==================================================================================================
    %/       Author: Franklin Cheng (fandycheng@ust.hk)
    %/  Last Update: October 10, 2023
    %/
    %/         NOTE: This function is tailor-made for consistently plotting
    %/               homvoller diagrams for ITCC study.
    %==================================================================================================
    
    dom_box_lon = []; dom_box_lat =[]; border_lines = []; time_lines = [];
    
    %/ use which_study to distinguish the domain definitions. Avoid inconsistency.
    if isequal(project_name, 'trimonsoon')
        date_mmdd_OneYr = date_array_gen('years', 1979, 'st_month', 1, 'st_day', 1, 'ed_month', 12, 'ed_day', 31, 'output_date_format', 'yyyymmdd');
        date_mmdd_OneYr = mod(date_mmdd_OneYr, 1e4);  %/ 365 days (ignore leap day) -> consistent with pentad field processing.
        
        if isequal(dom_name, 'EAM')
           if zm_or_mm == 1
               dom_box_lon  = [110 135];
               dom_box_lat  = [15   45];
               
           elseif zm_or_mm == 2
               dom_box_lon  = [110 140];
               dom_box_lat  = [20   40];
           else
               dom_box_lon  = [110 140];
               dom_box_lat  = [20   45];  
           end
           %/ mark the clim. onset dates
           clim_onset_dates = [212, 402, 519, 616, 722, 907, 1022, 1210];
           time_lines       = findismember_loop(date_mmdd_OneYr, clim_onset_dates);
           
        elseif isequal(dom_name, 'IM')
           %/ same lon lat range for zm or mm.
           dom_box_lon  = [60   90];
           dom_box_lat  = [9    30];  
           
           %/ mark the clim. onset dates
           clim_onset_dates = [317, 516, 608, 906, 1010, 1226];
           time_lines       = findismember_loop(date_mmdd_OneYr, clim_onset_dates);
            
       elseif isequal(dom_name, 'IM_5N-30N')
           %/ same lon lat range for zm or mm.
%            dom_box_lon  = [60   90];
           dom_box_lon  = [70   90];
           dom_box_lat  = [5    30];
           
           %/ mark the clim. onset dates
           clim_onset_dates = [317, 516, 608, 906, 1010, 1226];
           time_lines       = findismember_loop(date_mmdd_OneYr, clim_onset_dates);
           
        elseif isequal(dom_name, 'AuM')
           %/ same lon lat range for zm or mm.
           dom_box_lon  = [110 130];
           dom_box_lat  = [-15  -5];
           
           %/ mark the clim. onset dates
           clim_onset_dates = [310, 330, 1112, 1213];
           time_lines       = findismember_loop(date_mmdd_OneYr, clim_onset_dates);
        
        else
            error('No preset domain for %s in the set of %s!', dom_name, project_name)
        end
    end
    
    if isequal(project_name, 'ITCC')
        
        if isequal(dom_name, 'MC_10Sto5N')
           dom_box_lon  = [30 190];
           dom_box_lat  = [-10 5];
           border_lines = []; %/AS, IN, BoB

        elseif isequal(dom_name, 'SA')
           if zm_or_mm == 1
               dom_box_lon  = [ 52 87];
               dom_box_lat  = [-15 35];
               border_lines = [-10 10 28];  %/ Tropics, South Asia

           elseif zm_or_mm == 2
               dom_box_lon  = [40 105];
               dom_box_lat  = [10  23];
               border_lines = [52; 73; 87; 98]; %/AS, IN, BoB
           end

        elseif isequal(dom_name, 'MJO') && zm_or_mm == 2
           dom_box_lon  = [ 40  105];
           dom_box_lat  = [-10  10];
           border_lines = [];

        elseif isequal(dom_name, 'MJO-north') && zm_or_mm == 2
           dom_box_lon  = [ 40  105];
           dom_box_lat  = [  0  10];
           border_lines = [];

        elseif isequal(dom_name, 'MJO-south') && zm_or_mm == 2
           dom_box_lon  = [ 40  105];
           dom_box_lat  = [-10    0];
           border_lines = [];

        elseif isequal(dom_name, 'WP-EA')
           if zm_or_mm == 1
               dom_box_lon  = [118  138];
               dom_box_lat  = [-20   40];
               border_lines = [6; 31];       %/ MC, Phillipine Sea, Japan

           elseif zm_or_mm == 2
               dom_box_lon  = [70  180];
               dom_box_lat  = [10   20];     %/ Follow Fig. 10b in Wang and Xu 1997.
    %            dom_box_lat  = [15   20];     %/ Follow Fig. 10b in Wang and Xu 1997.
    %            dom_box_lat  = [12   22];     %/ Follow Fig. 10b in Wang and Xu 1997.

               border_lines = [];       
           end

        elseif isequal(dom_name, 'SA-IDC-MC') && zm_or_mm == 1
           dom_box_lon  = [70  150];
           dom_box_lat  = [-25  40];    
           border_lines = [-10, 6, 27];       %/ southmost and northmost of MC, Southern flank of Himalayas 

        elseif isequal(dom_name, 'SA-IDC') && zm_or_mm == 2
           dom_box_lon  = [70  150];
           dom_box_lat  = [6    27];    
           border_lines = [95];      %/ westmost of MC or IDC

%       elseif isequal(dom_name, 'IDC-Sumatra') && zm_or_mm == 1   %<--- springtime fast transition region.
%            dom_box_lon  = [90   108];
%            dom_box_lat  = [-10   27];    
%            border_lines = [8];              %/ a rough lat division of IDC and MC.
%            
%            if contains(select_field, 'pentad')
%                time_lines   = [24, 33, 51, 61]; %/ two time periods for transition in Sumatra-Indochina
%                
%            elseif contains(select_field, 'daily')
%                time_lines   = [116, 165, 251, 305]; %/ two time periods for transition in Sumatra-Indochina
%            end
           
        elseif isequal(dom_name, 'IDC-Sumatra-new') && zm_or_mm == 1   %<--- springtime fast transition region.
           dom_box_lon  = [89   110];
           dom_box_lat  = [-9    25];   
           border_lines = [8];              %/ a rough lat division of IDC and MC.
           
           if contains(select_field, 'pentad')
               time_lines   = [26, 32, 51, 60]; %/ two time periods for transition in Sumatra-Indochina
               
           elseif contains(select_field, 'daily')
               time_lines   = [116, 165, 251, 305]; %/ two time periods for transition in Sumatra-Indochina
           end
           
        elseif isequal(dom_name, 'IDC-Sumatra-new-v2') && zm_or_mm == 1   %<--- springtime fast transition region.
           dom_box_lon  = [90   110]; %<- narrower to avoid the Indian CISO wave.
           dom_box_lat  = [-4    20];   
           border_lines = [8];              %/ a rough lat division of IDC and MC.
           
           if contains(select_field, 'pentad')
               time_lines   = [26, 32, 51, 60]; %/ two time periods for transition in Sumatra-Indochina
               
           elseif contains(select_field, 'daily')
               time_lines   = [116, 165, 251, 305]; %/ two time periods for transition in Sumatra-Indochina
           end
       
       elseif isequal(dom_name, 'IDC-Sumatra-new-v3') && zm_or_mm == 1   %<--- springtime fast transition region.
           dom_box_lon  = [90   110]; %<- narrower to avoid the Indian CISO wave.
           dom_box_lat  = [-6    18];   
           border_lines = [6];              %/ a rough lat division of IDC and MC.
           
           if contains(select_field, 'pentad')
               time_lines   = [26, 32, 51, 60]; %/ two time periods for transition in Sumatra-Indochina
               
           elseif contains(select_field, 'daily')
               time_lines   = [116, 165, 251, 305]; %/ two time periods for transition in Sumatra-Indochina
           end
           
        elseif isequal(dom_name, 'MC') && zm_or_mm == 2
           dom_box_lon  = [90  160];
           dom_box_lat  = [-8   2];    
           border_lines = [120];     

        elseif isequal(dom_name, 'ISM')
           dom_box_lon  = [62  92];
           dom_box_lat  = [6   27];

        elseif isequal(dom_name, 'IDC')
           dom_box_lon  = [92  110];
           dom_box_lat  = [6   27];

          if isequal(select_field(1:6), 'pentad')
               time_lines   = [24, 33, 51, 61]; %/ two time periods for transition in Sumatra-Indochina
               
          elseif isequal(select_field(1:5), 'daily')
               time_lines   = [116, 165, 251, 305]; %/ two time periods for transition in Sumatra-Indochina
          end
           
        elseif isequal(dom_name, 'Meiyu')
           dom_box_lon  = [110  140];
           dom_box_lat  = [27   35];

        elseif isequal(dom_name, 'MC-west')
           dom_box_lon  = [95  125];
           dom_box_lat  = [-10    6];

        elseif isequal(dom_name, 'MC-east')
           dom_box_lon  = [125  152];
           dom_box_lat  = [-10    0];
        else
            error('code not set!');
        end
    end
    
    %/ check if anything is output
    if isempty(dom_box_lon) && isempty(dom_box_lat)
         error('Nothing is output! Check the input ''project_name''!');
    end
end