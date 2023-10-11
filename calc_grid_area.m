%%
function grid_area = calc_grid_area(varargin)

    pnames = {'lon', 'lat'};
    dflts = repmat([], 1, length(pnames));
    [lon, lat] = internal.stats.parseArgs(pnames, dflts, varargin{:});
    %%
    %==================================================================================================
    %/       Author: Franklin Cheng (fandycheng@ust.hk)
    %/  Last Update: October 10, 2023
    %/
    %/       Output: area (m^2)
    %==================================================================================================
    
    r_earth = 6.371e6;
    if isvector(lon) && isvector(lat)
        if size(lat, 2) == 1         lat = lat';         end   %/ column to row array

        dx = abs(diff(lon(1:2)));
        dy = abs(diff(lat(1:2)));

        %/ Since the meridional length of the grid (almost) does not vary with latitude
        grid_y = r_earth*dy/180*pi;
        grid_y_array = repmat(grid_y, length(lon), 1); %/ column vector, round() is to make the number integer.

        R = r_earth*cosd(lat);
        grid_x_array = R.*dx/180*pi;               %/ row vector
        grid_area    = grid_y_array*grid_x_array;
    else
        %/ A bit complicated to compute the area of uneven grid cells
        dx      = diff(lon,[],1);
        dx_half = reshape([dx'/2;dx'/2], size(lon,2), [])';
        dx_half = [dx_half(1,:); dx_half; dx_half(end,:)];
        dx_adj  = reshape(sum(reshape(dx_half,2,[]),1), size(lon,1), []);
        R       = r_earth*cosd(lat);
        grid_x  = R.*dx_adj/180*pi;               %/ row vector

        dy      = diff(lat,[],2);
        dy_half = reshape([dy/2; dy/2], size(lat,1), []);
        dy_half = [dy_half(:,1), dy_half, dy_half(:,end)];
        dy_adj  = reshape(sum(reshape(dy_half',2,[]),1), size(lat,1), []);

        grid_y    = r_earth*dy_adj/180*pi;
        grid_area = grid_x.*grid_y;
    end

end
