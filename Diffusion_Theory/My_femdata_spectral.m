function data = My_femdata_spectral(mesh,frequency,save_num,save_fluence)

% data = femdata_spectral(mesh,frequency,wv)
%
% Calculates data (phase and amplitude) for a given
% spectral mesh at a given frequency (MHz).
% outputs phase and amplitude in structure data
% and mesh information in mesh
%
% mesh is the input mesh (variable or filename)
% frequency is the modulation frequency (MHz)
% wv is optional wavelength array
if nargin < 4
    save = false;
end

if nargin < 3
    save_num = [];
end

parallel = parallel_init();

% error checking
if frequency < 0
    errordlg('Frequency must be nonnegative','NIRFAST Error');
    error('Frequency must be nonnegative');
end

% If not a workspace variable, load mesh
if ischar(mesh)== 1
  mesh = load_mesh(mesh);
end


wv_array = sort(mesh.wv);
linki = 1:length(wv_array);

nwv = length(wv_array);

%**************************************************************
% Initiate log file
% fid_log = fopen([mesh.name,'.log'],'w');
% fprintf(fid_log,'**************************\n\nForward Model Log\n\n*************************\n\n');
% fprintf(fid_log,'Forward Mesh  = %s\n',mesh.name);
% fprintf(fid_log,'Wavelength Array  =  %dnm\n',wv_array);

%****************************************************************
% Run femdata for each wavelength and save data as .paa

data.paa = []; data.wv=[];
  disp(sprintf('Calculating data for: %g nm',(wv_array(1))))

  % Extract coefficient
  for region = 1:8
      ind = find(mesh.region==region);
      mua_spec = mesh.mua_spec{region};
      musr_spec = mesh.musr_spec{region};
      kappa_spec = mesh.kappa_spec{region};
      mesh.mua(ind) = mua_spec(1);
      mesh.mus(ind) = musr_spec(1);
      mesh.kappa(ind) = kappa_spec(1);
  end
  
  %**************************************************************** 
  % if sources are not fixed, move sources depending on mus
  if mesh.source.fixed == 0
    mus_eff = mesh.mus;
    [mesh]=move_source(mesh,mus_eff,3);
    clear mus_eff
  end
  
  mesh_temp = mesh;
  % build link file for first wv from first wv of function call.
  mesh_temp.link = [mesh.link(:,1:2) mesh.link(:,linki(1)+2)];
  
  [data_single_wv,junk] = femdata_stnd(mesh_temp,frequency,wv_array(1),save_num, save_fluence);
  data.paa = [data.paa, data_single_wv.paa];
  data.wv = [data.wv wv_array(1)];
  data.link = mesh_temp.link;
  clear data_single_wv junk
  
  % PARALLEL

parfor i = 2 : nwv
  data_single_wv=[];
  mus_eff=[];
  disp(sprintf('Calculating data for: %g nm',(wv_array(i))))


  %****************************************************************
  % calculate absorption and scattering coefficients from concetrations and
  % scattering parameters a and b
  mesh_temp=mesh;
  mesh_temp.link = [mesh.link(:,1:2) mesh.link(:,linki(i)+2)];

  for region = 1:8
      ind = find(mesh.region==region);
      mua_spec = mesh.mua_spec{region};
      musr_spec = mesh.musr_spec{region};
      kappa_spec = mesh.kappa_spec{region};
      mesh_temp.mua(ind) = mua_spec(i);
      mesh_temp.mus(ind) = musr_spec(i);
      mesh_temp.kappa(ind) = kappa_spec(i);
  end

  % if sources are not fixed, move sources depending on mus
  if mesh_temp.source.fixed == 0
    mus_eff = mesh_temp.mus;
    [mesh_temp]=move_source(mesh_temp,mus_eff,3);
   % clear mus_eff
  end

  [data_single_wv,mesh_temp] = femdata_stnd(mesh_temp,frequency,wv_array(i),save_num,save_fluence);
  data_amp(:,i) = data_single_wv.paa(:,1);
  data_phase(:,i)= data_single_wv.paa(:,2);
  data_wv(i) = wv_array(i);
  data_link(:,i-1) = mesh.link(:,linki(i)+2);
  % clear data_single_wv
end

for i=2:nwv
    data.paa=[data.paa,data_amp(:,i), data_phase(:,i)];
    fluence = sum(data_amp,1);
%     plot(fluence)
    data.wv(i)=data_wv(i);
    data.link(:,i+2) = data_link(:,i-1);
end
    
if isfield(mesh,'R')
    mesh=rmfield(mesh,'R');
end