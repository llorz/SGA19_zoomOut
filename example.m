clc; close all; clear;
addpath(genpath('utils/'))
%% some parameters
% fmap size: k2-by-k1
k1 = 100;
k2 = 100;
% params to compute WKS descriptors
% #WKS = length(1:numSkip:numTimes)
numTimes = 100;
numSkip = 100;
% relative weights of different terms to optimize a functional map
para.a = 2e-1;
para.b = 1e-2;
para.c = 8e-4;  % weight for the Laplacian mask term!
para.alpha = 1e-1;
% params to pre-process the meshes
meshOptions = {'IfComputeGeoDist',false,'IfComputeLB',true,'IfComputeNormals',true,'numEigs',100};
% params to visualize the maps
plotOptions = {'IfShowCoverage',false,'OverlayAxis','y','cameraPos',[0,90]};
%% read the mesh
mesh_dir = 'data/';
s1_name = 'tr_reg_012';
s2_name = 'tr_reg_015';
S1 = MESH.MESH_IO.read_shape([mesh_dir, s1_name]);
S2 = MESH.MESH_IO.read_shape([mesh_dir, s2_name]);
%% preprocess the meshes
S1 = MESH.preprocess(S1,meshOptions{:});
S2 = MESH.preprocess(S2,meshOptions{:});

%% compute the WKS descriptors
B1 = S1.evecs(:,1:k1); B2 = S2.evecs(:,1:k2);
Ev1 = S1.evals(1:k1); Ev2 = S2.evals(1:k2); 

fct1_all = fMAP.waveKernelSignature(B1, Ev1, S1.A, numTimes);
fct2_all = fMAP.waveKernelSignature(B2, Ev2, S2.A, numTimes);
fct1 = fct1_all(:,1:numSkip:end);
fct2 = fct2_all(:,1:numSkip:end);

%% optimize the functional map using the standard or the complex resolvent Laplacian term
[C12, M_old] = compute_fMap_complRes(S1,S2,B1,B2,Ev1,Ev2,fct1,fct2,para, 'standard');
[C12_slant, M_slant] = compute_fMap_complRes(S1,S2,B1,B2,Ev1,Ev2,fct1,fct2,para, 'slant');
[C12_new, M_new] = compute_fMap_complRes(S1,S2,B1,B2,Ev1,Ev2,fct1,fct2,para, 'complRes');
T21 = fMAP.fMap2pMap(B1,B2,C12);
T21_slant = fMAP.fMap2pMap(B1,B2,C12_slant);
T21_new = fMAP.fMap2pMap(B1,B2,C12_new);

%%
% visualize the computed maps
figure(1);
subplot(1,3,1);
MESH.PLOT.visualize_map_colors(S2,S1,T21,plotOptions{:}); title('standard Mask');
subplot(1,3,2);
MESH.PLOT.visualize_map_colors(S2,S1,T21_slant,plotOptions{:}); title('slanted Mask');
subplot(1,3,3);
MESH.PLOT.visualize_map_colors(S2,S1,T21_new,plotOptions{:}); title('complex resolvent Mask');

% visualize the mask 
figure(2);
subplot(1,3,1); imagesc(M_old); axis square;  title('Standard Laplacian Mask');
subplot(1,3,2); imagesc(M_slant); axis square;  title('Standard Laplacian Mask');
subplot(1,3,3); imagesc(M_new); axis square;  title('Complex Resolvent Laplacian Mask');