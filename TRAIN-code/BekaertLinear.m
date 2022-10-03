function ph_tropo_linear=BekaertLinear(hgt, lon, lat, phuw, height, width, save_path,TRAINMATLABPATH)
% [ph_tropo_linear] = aps_linear(save_path)
% Scipt to compute the tropospheric delay map from a linear relation
% between phase and topography, optional a non-deforming polygon can be
% specified by setting crop_flag to 'y' in the parms_aps list.
% The computed tropospheric delay is in the same units as the
% inputed interferogram phases.
%
% All required inputs will be read from the aps_parm list. This includes:
% non_defo_flag 'y' or 'n' to use a non-deformaing region.
%               Polygon of the non-deforming area. By default 'n' the whole
%               interferogram is used. Change to 'y' by using setparm_aps.
%               Note that this variable is a matrix with in its columns the
%               longitude and latitude of the non-deforming area.
% hgt_matfile   Path to the interferograms
%               Interferogram phases, given as a column vector or matrix
%               with in its columens the different interferograms.
%               Stamps structure will automatically be recognised. Use
%               setparm_aps so change the data path pointing to the .mat
%               file.
% hgt_matfile   Path to the heights file.
%               Colum vector of the topography.
%               Stamps structure will automatically be recognised. Use
%               setparm_aps so change the data path pointing to the .mat
%               file.
% ll_matfile    Path to the longitude and latitude file
%               Matrix with in its columns the longitude and latitude.
%               Stamps structure will automatically be recognised. Use
%               setparm_aps so change the data path pointing to the .mat
%               file. In case the poly argument is specified, both need to
%               have the same units.
%
% OUTPUTS:
% ph_tropo_linear   The topography correlated delays either estimated from a
%                   linear relationship over the whole interferogram of using
%                   a non-deforming region. By default the output of this
%                   function is stored in the 'tca2.mat' or 'tca_sb2.mat'
%                   for StaMPS SM and SB option. In case no StaMPs
%                   structure is used the data is saved in 'tca2.mat'.
%
%     Copyright (C) 2015  Bekaert David - University of Leeds
%     Email: eedpsb@leeds.ac.uk or davidbekaert.com
%
%     This program is free software; you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation; either version 2 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License along
%     with this program; if not, write to the Free Software Foundation, Inc.,
%     51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
%
% Bekaert David - Leeds University 2013
% Modifcations:
% DB:   04/2013     Include SB functionality for stamps
% DB:   04/2013     Uset setparm_aps and getparm_aps to get the processign
%                   parameters
% DB:   06/2015     Suppress the output, as its saved anyway.
% DB:   01/2016     Include non-Stamps support
% DB:   02/2016     Remove inlfuence of NaN's
% DB:   12/2016     Revert to aps_save to allow saving of large variables


addpath(TRAINMATLABPATH);

test_fig = 0;   % debug flag which plots the scatter cloud and the
                % estimated line
n_fig_line = 7; % number of ifgs per row for the plots
fontsize =10;   % figure fontsize

phuw = reshape(phuw,height*width,1);
hgt = reshape(hgt,height*width,1);
lon = reshape(lon,height*width,1);
lat = reshape(lat,height*width,1);
lonlat = [lon lat];
lonlat = reshape(lonlat,height*width,2);

if nargin <1 || isempty(save_path)
    save_path = [pwd '\'];
end

stamps_processed = 'n';

if strcmp(getparm_aps('non_defo_flag',1),'n')
    non_defo_flag = 0;
else
    non_defo_flag = 1;
end

if strcmp(stamps_processed,'y')
   load psver
else
    psver = 2;
end

% loading the data
%phuw = load(phuw_matfile);
%lonlat = load(ll_matfile);
%hgt = load(hgt_matfile);
if strcmp(stamps_processed,'y')
   load psver
else
    psver = 2;
end
%phuw = phuw.phuw_matfile(:,1);%phuw.ph_uw;
%lonlat = lonlat.ll_matfile;%lonlat.lonlat;
%hgt = hgt.hgt(:,1);
%% Loading of the data
% file names of the output data

apsname = [save_path filesep 'tca' num2str(psver) '.mat'];
apssbname = [save_path filesep 'tca_sb' num2str(psver) '.mat'];


% the number of interferograms
n_dates = size(phuw,2);


%% use a non-deforming area
if non_defo_flag==1
    non_defo = load('non_defo.mat');
    poly = non_defo.poly;

    % search those points within the non-deforming polygon
    ix_points = [1:size(hgt,1)]';
    ix_temp = inpolygon(lonlat(:,1),lonlat(:,2),poly(:,1),poly(:,2));
    ix_points=ix_points(ix_temp);
    ixnon_points=[1:size(hgt,1)]';
    ixnon_points=ixnon_points(ix_points);


    clear ix_temp
else
    % use all points
    ix_points = [1:size(hgt,1)]';
    ixnon_points=[];
end


%% correct for DEM error
DEM_corr = getparm_aps('powerlaw_DEM_corr',1);
% geting the number of interferograms and points from the phase matrix
n_interferograms= size(phuw,2);
n_points = size(phuw,1);
if strcmp(DEM_corr,'y') && n_interferograms>5;
   % these are erros that scale with perpendicular baseline
   % loading the perpendicualr baseline information
   bperp_matfile = getparm_aps('bperp_matfile',1);
   bperp = load(bperp_matfile);
    if strcmp(stamps_processed,'y')
      bperp = bperp.bperp;
   end
   % checking the size of bperp
   if size(bperp,2) > 1
      bperp = bperp';
      if size(bperp,2) >1
          error('myApp:argChk', ['bperp is not a vector,... \nAbort,... \n'])
      end
   end
   % estimating the correlated errors
   DEM_corr_e = lscov(bperp,phuw')';

   % removing DEM correlated errors
   phuw = phuw-repmat(bperp',n_points,1).*repmat(DEM_corr_e,1,n_interferograms);
   clear A
else
    if strcmp(DEM_corr,'y') && n_interferograms<=5
        fprintf('Not enough interferograms to make a reliable estimate for the DEM error \n')
        DEM_corr = 'n';
    end
    DEM_corr_e = zeros([n_points 1]);
end





%% Compute the linear relation between phase and topography for each interferogram
% and compute the tropospheric delay for the full interferogram from it.
% initialisation of the tropospheric delay matrix
ph_tropo_linear = zeros([size(hgt,1) n_dates]);
hgt_range = [min(hgt) max(hgt)]';
if hgt_range(2)>10
  % height are in m
  hgt = hgt/1000;       % km
  hgt_range = hgt_range/1000;
end

xy = (llh2local(lonlat',mean(lonlat)))';
xy = [[1:size(xy,1)]' xy];
ps.xy = xy;
ps.n_ifg = size(phuw,2);
ps.n_ps = size(phuw,1);


ix_points_or = ix_points;
ixnon_points_or = ixnon_points;
for k=1:n_dates
    ixnon_points = ixnon_points_or;
    ix_points = ix_points_or;


    % removing NaN's from the phase
    ixnon_points = unique([ixnon_points ; find(isnan(phuw(:,k))==1)]);

    ix_points = [1:n_points]';
    ix_points(ixnon_points)=[];


    % Setting up the design matrix
    A = [hgt(ix_points) ones([length(ix_points) 1])];


    % compute the linear relation
    coeff = lscov(A,phuw(ix_points,k));
    % computation of the delay
    ph_tropo_linear(:,k) = [hgt ones(size(hgt))]*coeff;

    % set those pixels not used back to NaN
    ph_tropo_linear(isnan(phuw(:,k)),k)=NaN;


end

save 'ph_tropo_linear.mat' ph_tropo_linear;
%return ph_tropo_linear;
end
