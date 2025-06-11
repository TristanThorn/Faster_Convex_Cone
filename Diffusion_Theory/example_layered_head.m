clear all
addpath(genpath('/home/zz111/Spectral_unmixing/NIRFAST-9.1'));
%%mesh loading
% we can create a mesh from a segmented mha file:
%image2mesh_gui
% and then follow the guide in powerpoint 'Head Model Tutorial.pptx'
% OR use one made for you already...NOTE we will use a spectral mesh
mesh = load_mesh('digital_vessel_mimicCC_0.2_nirfast_mesh');

% for this example, we will look at only 2 wavelengths
mesh.wv = mesh.wv(2:end);

% We have 2 source located at
mesh.source.coord
% We have 6 detectors, located at
mesh.meas.coord

%% parameter assignment
% assign tissue properties: Note we are assigning HbO, Hb, Water and
% scattering...see powerpoint 'Head Model Tutorial.pptx' for values
% skin (region 1)
values.sa = 2;
values.sp = 0.5;
values.HbO =  0.045;
values.deoxyHb = 0.015;
values.Water = 0.5;
mesh = set_mesh(mesh,1,values);
% Bone (region 2)
values.sa = 1.4;
values.sp = 1.4;
values.HbO =  0.0392;
values.deoxyHb = 0.0098;
values.Water = 0.15;
mesh = set_mesh(mesh,2,values);
% CSF (region 3)
values.sa = 0.5;
values.sp = 0.5;
values.HbO =  0.0009;
values.deoxyHb = 0.0001;
values.Water = 0.99;
mesh = set_mesh(mesh,3,values);
% Brain (regions 4 & 5)
values.sa = 0.76;
values.sp = 0.76;
values.HbO =  0.054;
values.deoxyHb = 0.022;
values.Water = 0.78;
mesh = set_mesh(mesh,4,values);
mesh = set_mesh(mesh,5,values);
clear values;

% let us look at the cross section of HbO
% To do this, we will interpolate values onto a uniform grid, at a cross
% section
[x,y,z]=meshgrid(min(mesh.nodes(:,1)):1:max(mesh.nodes(:,1)),100,min(mesh.nodes(:,3)):1:max(mesh.nodes(:,3)));
F = TriScatteredInterp(mesh.nodes(:,1),mesh.nodes(:,2),mesh.nodes(:,3),mesh.conc(:,1));
val = F(x,y,z);
figure; surf(squeeze(x),squeeze(y),squeeze(z),squeeze(val));
view(0,0); shading interp; colorbar; axis tight

%% how to generate modelled data
% now let us generate boundary data
% 1. CW mode
tic
data_cw = femdata(mesh,0);
toc

% data.paa(:,1) shows intensity for wavelength 1
% data.paa(:,3) shows intensity for wavelength 2
% data.paa(:,5) shows intensity for wavelength 3
% note even columns of data.paa are zero as we have no frequency
% information in CW

% 2. FD mode at 100 MHz
data_FD = femdata(mesh,100);

% data.paa(:,1) shows intensity for wavelength 1
% data.paa(:,3) shows intensity for wavelength 2
% data.paa(:,5) shows intensity for wavelength 3
% data.paa(:,2) shows phase for wavelength 1
% data.paa(:,4) shows phase for wavelength 2
% data.paa(:,6) shows phase for wavelength 3

%% simulate data for a hypoxia model in Brain
% now let us generate some forward data for some change in
% oxygenation for brain, ranging from initial 71% down to 41% and up
% we will for simplicity assume CW data
delta_SO2 = [0.71:-0.1:0.41 0.41 0.41 0.51:0.1:0.71]; % note SO2 is in fractions

% now let us update properties of brain and calculate simulated data
for i = 1 : size(delta_SO2,2)
    HbT = 0.076;
    values.HbO =  HbT.*delta_SO2(i);
    values.deoxyHb = HbT-values.HbO;
    mesh_anom = set_mesh(mesh,4,values);
    mesh_anom = set_mesh(mesh_anom,5,values);
    data_anom{i} = femdata(mesh_anom,0);
end

%% Parameter Recover--Jacobian!
% now let us recover SO2 change from initial state using a Jacobian (also known as weight
% or sensitivity matrix) based on rest state, for each tissue type
% Please see NIRFAST paper for structure of this Jacobian matrix:
J = jacobian(mesh,0);

% let us look at the Jacobian for 3 sources at a wavelength for Hb
% To do this, we will interpolate values onto a uniform grid, at a cross
% section
nn = size(mesh.nodes,1);
[x,y,z]=meshgrid(min(mesh.nodes(:,1)):1:max(mesh.nodes(:,1)),100,min(mesh.nodes(:,3)):1:max(mesh.nodes(:,3)));
F = TriScatteredInterp(mesh.nodes(:,1),mesh.nodes(:,2),mesh.nodes(:,3),J(4,1:nn)'./mesh.support);
val = F(x,y,z);
figure; surf(squeeze(x),squeeze(y),squeeze(z),squeeze(val));
view(0,0); shading interp; colorbar; axis tight; title('1st source and detector');
F = TriScatteredInterp(mesh.nodes(:,1),mesh.nodes(:,2),mesh.nodes(:,3),J(5,1:nn)'./mesh.support);
val = F(x,y,z);
figure; surf(squeeze(x),squeeze(y),squeeze(z),squeeze(val));
view(0,0); shading interp; colorbar; axis tight; title('2nd source and detector');
F = TriScatteredInterp(mesh.nodes(:,1),mesh.nodes(:,2),mesh.nodes(:,3),J(6,1:nn)'./mesh.support);
val = F(x,y,z);
figure; surf(squeeze(x),squeeze(y),squeeze(z),squeeze(val));
view(0,0); shading interp; colorbar; axis tight; title('3rd source and detector');

% This Jacobian relates each detector measurement to a small change to each
% node of the head model, therefore to get a Jacobian for each tissue type

%% Parameter Reconstruction
% we will need to sum appropriate nodal values -- See NIRFAST paper for
% structure of Jacobian, but in brief
% J = [wavelength_1 and Hb0, wavelength_1 for Hb, wavelength_1 for water
%      wavelength_2 and Hb0, wavelength_2 for Hb, wavelength_2 for water
%      wavelength_3 and Hb0, wavelength_3 for Hb, wavelength_3 for water];
nn = size(mesh.nodes,1); % get number nodes

% find nodes in each region
ind1 = find(mesh.region==1); % Skin (region 1)
ind2 = find(mesh.region==2); % Bone (region 2)
ind3 = find(mesh.region==3); % CSF (region 3)
ind4_5 = find(mesh.region>=4); % Brain (regions 4 and 5)

% sum appropriate nodes to get tissue specific value
J_layers = [sum(J(:,ind1),2) sum(J(:,ind1+nn),2) sum(J(:,ind1+nn+nn),2) ...
            sum(J(:,ind2),2) sum(J(:,ind2+nn),2) sum(J(:,ind2+nn+nn),2) ...
            sum(J(:,ind3),2) sum(J(:,ind3+nn),2) sum(J(:,ind3+nn+nn),2) ...
            sum(J(:,ind4_5),2) sum(J(:,ind4_5+nn),2) sum(J(:,ind4_5+nn+nn),2)];
        
% We can also assume a homogenous head and recover parameters for the whole
% head!
J_homog = [sum(J(:,ind1),2) sum(J(:,ind1+nn),2) sum(J(:,ind1+nn+nn),2)] + ...
            [sum(J(:,ind2),2) sum(J(:,ind2+nn),2) sum(J(:,ind2+nn+nn),2)] + ...
            [sum(J(:,ind3),2) sum(J(:,ind3+nn),2) sum(J(:,ind3+nn+nn),2)] + ...
            [sum(J(:,ind4_5),2) sum(J(:,ind4_5+nn),2) sum(J(:,ind4_5+nn+nn),2)];

% Or we can ignore everything BUT the brain
J_brain = [sum(J(:,ind4_5),2) sum(J(:,ind4_5+nn),2) sum(J(:,ind4_5+nn+nn),2)];
      
% let us calculate an inverse to the homogenous Jacobian and recover the
% bulk SO2 -- 
lambda = 1e-5;
Hess = J_homog*J_homog';
reg = eye(size(Hess)).*max(diag(Hess)).*lambda;
invJ_homog = J_homog'*inv(Hess+reg);
Hess = J_layers*J_layers';
reg = eye(size(Hess)).*max(diag(Hess)).*lambda;
invJ_layers = J_layers'*inv(Hess+reg);
Hess = J_brain*J_brain';
reg = eye(size(Hess)).*max(diag(Hess)).*lambda;
invJ_brain = J_brain'*inv(Hess+reg);

% now let us recover parameters based on the modelled data
% Note: We use log of data!
data_ref = [data_cw.paa(:,1);data_cw.paa(:,3)];
for i = 1 : size(delta_SO2,2)
    data_change(:,i) = [data_anom{i}.paa(:,1);data_anom{i}.paa(:,3)];
    recon_layers(:,i) = invJ_layers*log(data_change(:,i)./data_ref);
    recon_homog(:,i) = invJ_homog*log(data_change(:,i)./data_ref);    
    recon_brain(:,i) = invJ_brain*log(data_change(:,i)./data_ref);    
    % and add baseline to the calculated change
    recon_layers(:,i) = recon_layers(:,i)+...
        [mesh.conc(ind1(1),:) mesh.conc(ind2(1),:) mesh.conc(ind3(1),:) mesh.conc(ind4_5(1),:) ]';
    recon_homog(:,i) = recon_homog(:,i)+mesh.conc(ind4_5(1),:)';
    recon_brain(:,i) = recon_brain(:,i)+mesh.conc(ind4_5(1),:)';    
end

%% plot recovered SO2
% let us extract recovered SO2 parameters 
% note: these recovered parameters will be relative change in SO2 from
% baseline
recon_skin_SO2 = recon_layers(1,:)./sum(recon_layers(1:2,:),1);
recon_bone_SO2 = recon_layers(4,:)./sum(recon_layers(4:5,:),1);
recon_csf_SO2 = recon_layers(7,:)./sum(recon_layers(7:8,:),1);
recon_brain_SO2 = recon_layers(10,:)./sum(recon_layers(10:11,:),1);
recon_homog_SO2 = recon_homog(1,:)./sum(recon_homog(1:2,:),1)
recon_brainonly_SO2 = recon_brain(1,:)./sum(recon_brain(1:2,:),1)

% Plot results
figure; hold on
plot(delta_SO2,'k');
plot(recon_homog_SO2,'ro');
plot(recon_brain_SO2,'rx');
plot(recon_brainonly_SO2,'r+');
xlabel('measurement number');
ylabel('Calculated SO_2 (%)');
legend('True','Homogeneous model','Heterogeneous model (Brain)',...
    'Heterogeneous model (Brain Only)');