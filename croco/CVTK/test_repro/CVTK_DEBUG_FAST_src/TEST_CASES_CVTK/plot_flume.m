%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Make plot from the results of the FLUME test case
% 
%  Further Information:  
%  http://www.crocoagrif.org/
%  
%  This file is part of CROCOTOOLS
%
%  CROCOTOOLS is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation; either version 2 of the License,
%  or (at your option) any later version.
%
%  CROCOTOOLS is distributed in the hope that it will be useful, but
%  WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software
%  Foundation, Inc., 59 Temple Place, Suite 330, Boston,
%  MA  02111-1307  USA
%
%  Ref: Penven, P., L. Debreu, P. Marchesiello and J.C. McWilliams,
%       Application of the ROMS embedding procedure for the Central 
%      California Upwelling System,  Ocean Modelling, 2006.
%
%  Patrick Marchesiello, IRD 2017
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
%================== User defined parameters ===========================
%
% --- model params ---
%
fname = 'flume_his.nc'; % croco file name

morph_fac = 24;         % morphological factor (sediment.in)
morph_cpl =  1;         % feedback to currents

makepdf    = 0;         % make pdf file
add_plot  =  0;         % addition plots
%
%======================================================================
%
% ---------------------------------------------------------------------
% --- get grid from numerical model ---
% ---------------------------------------------------------------------
%
yindex = 3; % Mm=3 with NS periodic conditions
%
nc=netcdf(fname,'r');
tindex  =length(nc{'scrum_time'}(:)); % reads last record
tindex12=min(tindex,7);               % 12 h
tindex24=min(tindex,13);              % 24 h
tindex0 =min(tindex,2);               %  2 h
%
% horizontal grid
hr=squeeze(nc{'h'}(yindex,:));
xr=squeeze(nc{'x_rho'}(yindex,:));
hu=0.5*(hr(1:end-1)+hr(2:end));
xu=0.5*(xr(1:end-1)+xr(2:end));
L=length(hr);
%
% new bathy from last record
bed0=squeeze(nc{'bed_thick'}(1,1,yindex,:));
bed =squeeze(nc{'bed_thick'}(tindex,1,yindex,:));
bed=bed-bed0;
hnew=hr-bed;
if morph_cpl,
 h=hnew;
else
 h=hr;
end
%
% vertical grid
N=length(nc('s_rho'));
theta_s=nc.theta_s(:); 
theta_b=nc.theta_b(:); 
hc=nc.hc(:); 
zeta=squeeze(nc{'zeta'}(tindex,yindex,:));
zr=squeeze(zlevs(h,zeta,theta_s,theta_b,hc,N,'r',2));
zw=squeeze(zlevs(h,zeta,theta_s,theta_b,hc,N,'w',2));
zru=0.5*(zr(:,1:end-1)+zr(:,2:end));
zwu=0.5*(zw(:,1:end-1)+zw(:,2:end));
%
xr2d=repmat(xr,[N 1]);
xu2d=repmat(xu,[N 1]);
xw2d=repmat(xr,[N+1 1]);
D   =zw(N+1,:)-zw(1,:);
D2d =repmat(D,[N 1]);
Du  =zwu(N+1,:)-zwu(1,:);
Du2d=repmat(Du,[N 1]);

% ---------------------------------------------------------------------
% --- read/compute numerical model fields (index 1) ---
% --------------------------------------------------------------------
time=morph_fac*nc{'scrum_time'}(tindex)/86400;

% ... Minimum depth in wetting-drying scheme
Dcrit=nc{'Dcrit'}(:);

% ... zonal velocity ...                         ---> xu,zu
u=squeeze(nc{'u'}(tindex,:,yindex,:));

% ... vertical velocity ...                      ---> xr,zw
w=squeeze(nc{'w'}(tindex,:,yindex,:));

% ... sediment concentration ...                 ---> xr,zr
t=squeeze(nc{'sand_1'}(tindex,:,yindex,:));

% ... total viscosity
Akv=squeeze(nc{'AKv'}(tindex,:,yindex,:));

% ... viscosity due to wave breaking ...
Akb=squeeze(nc{'Akb'}(tindex,:,yindex,:));

% ... new bathy at 12h and 24h ...  
bed0=squeeze(nc{'bed_thick'}(1,1,yindex,:));
bed =squeeze(nc{'bed_thick'}(tindex12,1,yindex,:));
bed=bed-bed0;
hnew12=bed-hr;
bed =squeeze(nc{'bed_thick'}(tindex24,1,yindex,:));
bed=bed-bed0;
hnew24=bed-hr;

% ... wave setup ...  
sup=squeeze(nc{'zeta'}(tindex0,yindex,:)); % init time
sup(hr<0)=sup(hr<0)+hr(hr<0)-Dcrit;

% ... u undertow ...
ubot =squeeze(nc{'u'}(tindex0,3,yindex,:)); % init time

% ... hrms ...  
hrms =squeeze(nc{'hrm'}(tindex0,yindex,:)); % init time

close(nc)

%============================================================
% --- plot ---
%=============================================================
%
zeta(D<Dcrit)=NaN;
u(Du2d<Dcrit)=NaN;
sup(D<=max(0.1,Dcrit))=NaN;
sup(hr<0)=NaN;
%
%-----------------------------------
% Eulerian velocities u,v,w
%-----------------------------------
%
figure('Position',[100 150 500 750])

subplot(2,1,1)
cmin=-0.1; cmax=0.1; nbcol=40;
cint=(cmax-cmin)/nbcol;
map=colormap(jet(nbcol));
map(nbcol/2  ,:)=[1 1 1];
map(nbcol/2+1,:)=[1 1 1];
colormap(map);
contourf(xu2d,zru,u,[cmin:cint:cmax]); hold on
shading flat; colorbar;
if morph_cpl,
 plot(xr,-hnew,'color','k','LineWidth',3);
else
 plot(xr,-hr,'k',xr,-hnew,'k--','LineWidth',3);
end
plot(xr,zeta,'color','g','LineWidth',3);
grid on
axis([-Inf Inf -0.6 0.05])
caxis([cmin cmax])
thour=floor(time*24);
title(['FLUME: U at Time ',num2str(thour),' hour'])
hold off

subplot(2,1,2)
cmin=-3.e-3; cmax=3.e-3; nbcol=40;
cint=(cmax-cmin)/nbcol;
map=colormap(jet(nbcol));
map(nbcol/2  ,:)=[1 1 1];
map(nbcol/2+1,:)=[1 1 1];
colormap(map);
contourf(xr2d,zr,w,[cmin:cint:cmax]); hold on
shading flat; colorbar;
if morph_cpl,
 plot(xr,-hnew,'color','k','LineWidth',3);
else
 plot(xr,-hr,'k',xr,-hnew,'k--','LineWidth',3);
end
plot(xr,zeta,'color','g','LineWidth',3);
grid on
axis([-Inf Inf -0.6 0.05])
caxis([cmin cmax])
thour=floor(time*24);
title(['FLUME: W at Time ',num2str(thour),' hour'])
hold off

if makepdf
 export_fig -transparent flume_u.pdf
end
%
%-------------------------------------------------
% Turbulent and wave-induced Viscosity/diffusivity
%-------------------------------------------------
%
if add_plot % ---
figure('Position',[100 150 500 750])
subplot(2,1,1)
cmin=0; cmax=1.e-3; nbcol=40;
cint=(cmax-cmin)/nbcol;
map=colormap(jet(nbcol));
colormap(map);
contourf(xw2d,zw,Akv,[cmin:cint:cmax]); hold on
shading flat; colorbar;
plot(xr,-h,'color','k','LineWidth',3);
grid on
axis([-Inf Inf -0.6 0.05])
caxis([cmin cmax])
thour=floor(time*24);
title(['FLUME: Akv at Time ',num2str(thour),' hour'])
hold off

subplot(2,1,2)
cmin=0; cmax=1.e-3; nbcol=40;
cint=(cmax-cmin)/nbcol;
map=colormap(jet(nbcol));
colormap(map);
contourf(xw2d,zw,Akb,[cmin:cint:cmax]); hold on
shading flat; colorbar;
plot(xr,-h,'color','k','LineWidth',3);
grid on
axis([-Inf Inf -0.6 0.05])
caxis([cmin cmax])
thour=floor(time*24);
title(['FLUME: Akb at Time ',num2str(thour),' hour'])
hold off

if makepdf
 export_fig -transparent flume_Kv.pdf
end
end % add_plot ---
%
%======================================================================
%            VALIDATION USING FLUME EXPERIMENT 
%               BY ROELVINK & STIVE 1989
%======================================================================
%
% ---------------------------------------------------------
% --- flume data from Roelvink & Stive 1989 ---
% ---------------------------------------------------------
%
xd=[0.:0.5:32];

% --- bathy at 0h, 12h & 24h on xd grid ---
hd_0 = [ ...
-0.5998  -0.6008  -0.6019  -0.6028  -0.6037  -0.6045  -0.6051  -0.6043 ...
-0.6036  -0.6029  -0.6021  -0.6013  -0.6006  -0.5861  -0.5656  -0.5467 ...
-0.5318  -0.5189  -0.5072  -0.4962  -0.4860  -0.4776  -0.4740  -0.4633 ...
-0.4501  -0.4371  -0.4250  -0.4116  -0.3982  -0.3878  -0.3773  -0.3668 ...
-0.3563  -0.3440  -0.3303  -0.3176  -0.3060  -0.2941  -0.2816  -0.2683 ...
-0.2516  -0.2365  -0.2230  -0.2127  -0.2048  -0.1957  -0.1832  -0.1705 ...
-0.1565  -0.1426  -0.1306  -0.1186  -0.1067  -0.0871  -0.0670  -0.0583 ...
-0.0496  -0.0378  -0.0258  -0.0137  -0.0017  +0.0114  +0.0260  +0.0367 ...
+0.0451 ];
hd_12 = [ ...
-0.6003  -0.5984  -0.5973  -0.5973  -0.5973  -0.5974  -0.5974  -0.5975 ...
-0.5975  -0.5976  -0.5976  -0.5977  -0.5977  -0.5759  -0.5592  -0.5433 ...
-0.5275  -0.5104  -0.4969  -0.4857  -0.4746  -0.4672  -0.4584  -0.4393 ...
-0.4211  -0.4025  -0.3807  -0.3627  -0.3467  -0.3309  -0.3150  -0.3076 ...
-0.3077  -0.3067  -0.3019  -0.2962  -0.2890  -0.2796  -0.2702  -0.2621 ...
-0.2540  -0.2444  -0.2358  -0.2313  -0.2256  -0.2117  -0.2017  -0.1962 ...
-0.1874  -0.1765  -0.1668  -0.1588  -0.1518  -0.1417  -0.1306  -0.1236 ...
-0.1163  -0.1013  -0.0906  -0.0469  -0.0204  +0.0044  +0.0241  +0.0370 ...
+0.0472 ];
hd_24 = [ ...
-0.5972  -0.5997  -0.5999  -0.5999  -0.6000  -0.6004  -0.6016  -0.6026 ...
-0.6023  -0.6005  -0.5981  -0.5961  -0.5928  -0.5721  -0.5469  -0.5328 ...
-0.5158  -0.4963  -0.4844  -0.4726  -0.4589  -0.4315  -0.4008  -0.3664 ...
-0.3335  -0.3341  -0.3357  -0.3370  -0.3383  -0.3382  -0.3308  -0.3213 ...
-0.3152  -0.3098  -0.3035  -0.2958  -0.2880  -0.2811  -0.2745  -0.2691 ...
-0.2646  -0.2597  -0.2525  -0.2435  -0.2347  -0.2257  -0.2170  -0.2096 ...
-0.2030  -0.1957  -0.1868  -0.1774  -0.1674  -0.1575  -0.1490  -0.1412 ...
-0.1320  -0.1210  -0.1122  -0.1086  -0.1042  -0.0389  +0.0200  +0.0375 ...
+0.0484 ];
hd_12(1:15)=NaN;
hd_24(1:15)=NaN;

% --- sea level (wave setup) on xd grid ---
sup_d = 1.e-3*[ ...
 0.4571   0.4498   0.4089   0.3681   0.3278   0.2897   0.2516   0.2160 ...
 0.2142   0.2124   0.2106   0.2089   0.2071   0.2053   0.2035   0.2017 ...
 0.1590  -0.0242  -0.2074  -0.3769  -0.3786  -0.3804  -0.3822  -0.3840 ...
-0.3858  -0.3876  -0.3894  -0.3911  -0.3929  -0.3947  -0.2894  -0.0807 ...
 0.1914   0.5881   0.9681   1.3378   1.7074   2.0942   2.4840   2.8738 ...
 3.2190   3.5608   3.9026   4.3259   4.7549   5.1842   5.6246   6.0649 ...
 6.5115   6.9871   7.4627   8.1402   8.8747  10.5380  12.4662  14.1726 ...
15.3697  16.4699  17.4170  18.8660  18.8642  18.8624  18.8606  18.8589 ...
    NaN ];

% --- Hrms ---
xd_hrms = [ ...
 3.02   4.91   7.93   9.86  12.77  15.00  17.38  19.61 ...
20.14  22.06  22.55  24.67  25.05  27.01 ];
hrms_d = [ ...
0.124   0.122   0.121   0.120   0.123   0.121   0.112   0.100 ...
0.094   0.087   0.081   0.068   0.065   0.047 ];

% --- Undertow ---
xd_ubot = [ ...
 9.85  14.76  17.00  17.00  19.68  22.15  24.58 24.58  27.04 31.00 ];
ubot_d = 1.e-2*[ ...
-0.48  -2.93  -4.60  -5.51  -6.64  -5.92  -5.92 -5.39 -3.85 0.06   ];
%
%----------------------------------------------------------
%  Plot Hrms + undertow + bed evolution
%----------------------------------------------------------
%
figure('Position',[800 150 700 700])
%
% hrms ...
%
subplot(3,1,1)
xmin=0; xmax=Inf; zmin=-0.05; zmax=0.15;
plot(xr,hrms,'b',xd_hrms,hrms_d,'b*','LineWidth',2); hold on
errorbar(xd_hrms,hrms_d,0.01*ones(size(hrms_d)),'o');
legend('Model','Flume','Location','SouthWest')
grid on
axis([xmin xmax zmin zmax])
thour=floor(time*24);
ylabel('Depth [m]','Fontsize',15)
title('Hrms','Fontsize',15)
set(gca,'Fontsize',15)
hold off
%
% Undertow ...
%
subplot(3,1,2)
xmin=0; xmax=Inf; zmin=-Inf; zmax=Inf;
plot(xu,100*ubot,'b',xd_ubot,100*ubot_d,'b*','LineWidth',2); hold on;
errorbar(xd_ubot,100*ubot_d,ones(size(ubot_d)),'o');
legend('Model','Flume', ...
      'Location','SouthWest')
grid on
axis([xmin xmax zmin zmax])
thour=floor(time*24);
ylabel('U [cm/s]','Fontsize',15)
title('Undertow','Fontsize',15)
set(gca,'Fontsize',15)
hold off
%
% bed evolution 12h + 24h ...
%
subplot(3,1,3)
xmin=0; xmax=Inf; zmin=-hr(1); zmax=max(zeta);
plot(xr,-hr,'k',xr,hnew12,'b',xd,hd_12,'b--', ...
                xr,hnew24,'m',xd,hd_24,'m--', ...
                xr,zeta,'g','LineWidth',2);
legend('0h','Model 12h','Flume 12h', ...
            'Model 24h','Flume 24h', ...
            'Location','NorthWest')
grid on
axis([xmin xmax zmin zmax])
thour=floor(time*24);
xlabel('Cross-shore distance [m]','Fontsize',15)
ylabel('Depth [m]','Fontsize',15)
title('Bed Evolution - 0h/12h/24h','Fontsize',15)
set(gca,'Fontsize',15)

if makepdf
 export_fig -transparent flume_valid.pdf
end
%
%-------------------------------------------
%  Bed evolution 12h + 24h
%-------------------------------------------
%
if add_plot % ---
figure('Position',[300 100 700 400])
xmin=0; xmax=Inf; zmin=-hr(1); zmax=max(zeta);
plot(xr,-hr,'k',xr,hnew12,'b',xd,hd_12,'b--', ...
                xr,hnew24,'m',xd,hd_24,'m--', ...
                xr,zeta,'g','LineWidth',2);
legend('0h','Model 12h','Flume 12h', ...
            'Model 24h','Flume 24h', ...
            'Location','NorthWest')
grid on
axis([xmin xmax zmin zmax])
thour=floor(time*24);
xlabel('Cross-shore distance [m]','Fontsize',15)
ylabel('Depth [m]','Fontsize',15)
title('Bed Evolution - 0h/12h/24h','Fontsize',15)
set(gca,'Fontsize',15)
hold off

if makepdf
 export_fig -transparent flume_valid_bed.pdf
end
end % add_plot
%
%-------------------------------------------
%  Wave setup
%-------------------------------------------
%
if add_plot % ---
figure('Position',[500 100 700 400])
xmin=5; xmax=Inf; zmin=-0.06; zmax=0.05;
plot(xr,sup,'b',xd,sup_d,'b--',xr,-hr/10,'k','LineWidth',2);
legend('Setup Model','Setup Exp','Bathy/10')
grid on
axis([xmin xmax zmin zmax])
thour=floor(time*24);
title(['FLUME: Setup at Time ',num2str(thour),' hour'])
hold off

if makepdf
 export_fig -transparent flume_valid_sup.pdf
end
end % add_plot ---
