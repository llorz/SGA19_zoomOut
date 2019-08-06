clc; close all; clear;
addpath(genpath('utils/'));
addpath('func_main/');
mesh_dir = 'data/';

s1_name = 'tr_reg_017';
%% read the mesh and compute the LB basis
S1 = MESH.MESH_IO.read_shape([mesh_dir, s1_name]);
S1 = MESH.compute_LaplacianBasis(S1, 100);
%% set the initial map
B1 = S1.evecs(:,1:4);
figure(1);
for i = 1:4
    subplot(1,4,i); plot_func_on_mesh(S1, B1(:,i)); title([s1_name, ': LB',num2str(i)]); view([0,90])
end
% we can see that the 3rd and 4th LB function show the left-to-right
% symmetry, therefore, we can flip the sign of the 3rd and 4th LB
C_ini = diag([1,1,-1,-1]);
T_ini = fMAP.fMap2pMap(B1,B1,C_ini);

figure(2);
subplot(1,2,1); visualize_map_on_source(S1, S1, T_ini); title('Source'); view([0,90])
subplot(1,2,2); visualize_map_on_target(S1, S1, T_ini); title('The initial symmetric map'); view([0,90])

%you can also try the code provided here https://github.com/llorz/SGA18_orientation_BCICP_code
% to compute a symmetric point-wise using the orientation-reversing term
% then instead of using the provided BCICP but the zoomOut to refine it :D
%% apply zoomOut
para.k_init = 20;
para.k_step = 5;
para.k_final = 100;
tic
T11 = zoomOut_refine(S1.evecs, S1.evecs, T_ini, para);
t = toc;
fprintf('ZoomOut runtime: %.2f sec\n',t)

% fast version 
para.num_samples = 500;
tic
T11_fast = zoomOut_refine_fast(S1, S1, T_ini, para,1);
t = toc;
fprintf('ZoomOut with sampling runtime: %.2f sec\n',t)

%% visualize the maps
figure(3);
subplot(1,4,1); visualize_map_on_source(S1, S1, T_ini); title('Source'); view([0,90]);
subplot(1,4,2); visualize_map_on_target(S1, S1, T_ini); title('Initialization'); view([0,90])
subplot(1,4,3); visualize_map_on_target(S1, S1, T11); title('zoomOut'); view([0,90])
subplot(1,4,4); visualize_map_on_target(S1, S1, T11_fast); title('zoomOut (fast)'); view([0,90])