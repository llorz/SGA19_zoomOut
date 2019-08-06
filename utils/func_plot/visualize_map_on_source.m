function [] = visualize_map_on_source(S1, S2, T12)
[g1,g2,g3] = set_mesh_color(S2);
    f1 = g1(T12);
    f2 = g2(T12);
    f3 = g3(T12);
    trimesh(S2.surface.TRIV, S2.surface.X, S2.surface.Y, S2.surface.Z, ...
    'FaceVertexCData', [g1, g2, g3],...
    'FaceColor', 'interp', 'EdgeColor', 'none'); axis equal; axis off; 
end


function [g1,g2,g3] = set_mesh_color(S)
g1 = normalize_function(0,1,S.surface.X);
g2 = normalize_function(0,1,S.surface.Y);
g3 = normalize_function(0,1,S.surface.Z);
g1 = reshape(g1,[],1);
g2 = reshape(g2,[],1);
g3 = reshape(g3,[],1);
end

function fnew = normalize_function(min_new,max_new,f)
    fnew = f - min(f);
    fnew = (max_new-min_new)*fnew/max(fnew) + min_new;
end