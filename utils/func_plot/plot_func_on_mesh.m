function [] = plot_func_on_mesh(S,f)
if nargin < 2
    trimesh(S.surface.TRIV, S.surface.X, S.surface.Y, S.surface.Z, ...
        'FaceColor','interp', 'EdgeColor', 'none','FaceAlpha',0.5); axis equal; axis off;
else
    trimesh(S.surface.TRIV, S.surface.X, S.surface.Y, S.surface.Z,f, ...
        'FaceColor','interp', 'EdgeColor', 'none','FaceAlpha',0.89); axis equal; axis off;
end
if isfield(S,'name')
    title(S.name,'Interpreter','none')
end
end