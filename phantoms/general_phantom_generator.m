function phantom = general_phantom_generator(pixel_size, base_dims, layers, vessels)
%GENERAL_PHANTOM_GENERATOR Create 3-D layered phantom with vessels.
%   phantom = GENERAL_PHANTOM_GENERATOR(pixel_size, base_dims, layers, vessels)
%   builds a volume where the x-axis is along the vessel direction.
%
%   pixel_size - pixel size in millimetres
%   base_dims  - [nx ny nz] base dimensions before scaling by 0.1/pixel_size
%   layers     - struct array with fields:
%                   thickness : thickness of the layer in base units. Use [] for
%                               the remaining volume.
%                   value     : voxel value to assign to this layer
%   vessels    - struct array with fields:
%                   center_y  : centre coordinate in y (base units)
%                   center_z  : centre coordinate in z (base units)
%                   radius_y  : radius along y (base units)
%                   radius_z  : radius along z (base units, optional, default =
%                               radius_y)
%                   value     : voxel value of the vessel
%
%   The function returns a 3-D matrix "phantom" of size
%   round(base_dims*(0.1/pixel_size)).

fac = 0.1 / pixel_size;
 dims = round(base_dims .* fac);
 nx = dims(1); ny = dims(2); nz = dims(3);
 phantom = ones(nx, ny, nz);

% Create layered structure
current = 1;
for i = 1:length(layers)
    if isempty(layers(i).thickness)
        thick = nz - current + 1; % rest of the volume
    else
        thick = round(layers(i).thickness * fac);
    end
    last = min(current + thick - 1, nz);
    phantom(:, :, current:last) = layers(i).value;
    current = last + 1;
    if current > nz
        break;
    end
end

% Add vessels
for i = 1:length(vessels)
    cy = round(vessels(i).center_y * fac);
    cz = round(vessels(i).center_z * fac);
    ry = round(vessels(i).radius_y * fac);
    if isfield(vessels(i),'radius_z') && ~isempty(vessels(i).radius_z)
        rz = round(vessels(i).radius_z * fac);
    else
        rz = ry;
    end
    val = vessels(i).value;
    for x = 1:nx
        for z = 1:nz
            for y = 1:ny
                if ((y - cy)^2)/(ry^2) + ((z - cz)^2)/(rz^2) <= 1
                    phantom(x,y,z) = val;
                end
            end
        end
    end
end
end
