%%
% Definitions


close all
clear all


elecmark='bo';
marksize=5;
az=0; % Azimuth for the plot 
el=90; % elongation for thew plot
fns=14; % font size
cmin=0;cmax=0; % plot range can be set here, 0 automatic choosed
logme=1; % linear plots -> 0 logarithmic else
saveplot=1; % if save the plot adjust and press space

fp=fopen('tmp.meshname','r');
gridelem=fscanf(fp,'%s',1);
fclose(fp);

fp=fopen('tmp.elecname','r');
elecfile=fscanf(fp,'%s',1);
fclose(fp);

fp=fopen('tmp.lastmod','r');
modfile=fscanf(fp,'%s',1);
fclose(fp);

checkme = exist ('tmp.fenster','file');
if checkme~=0
    fp=fopen('tmp.fenster','r');
    fenster=fscanf(fp,'%s',1);
    fclose(fp);
end
checkme = exist ('tmp.range','file');
if checkme~=0
    fp=fopen('tmp.range','r');
    cmin=fscanf(fp,'%f',1);
    cmax=fscanf(fp,'%f',1);
    fclose(fp);
end
%%
% Knoten und Elemente einlesen 
% Filenamen bestimmen

% [gridelem,wdir] = uigetfile('*','Grid/element file (from CRTomo)');
% ex     = cd (wdir);
fp     = fopen(gridelem,'r'); % Standard file handler

% Knoten und Elementtypen bestimmen
sanz   = fscanf (fp,'%d',1); % Anzahl der Knoten
typanz = fscanf (fp,'%d',1); % Anzahl der Elementtypen
nA     = fscanf (fp,'%d',1); % Bandrbreite der FE Matrix (wird nicht gebraucht)

mtyp(1:typanz)=0;nelanz(1:typanz)=0;selanz(1:typanz)=0; % Felder allozieren
snr(1:sanz)=0;   sx(1:sanz)=0;      sy(1:sanz)=0;

% Elementgrößen bestimmen
nel=0;
for i=1:typanz
    mtyp(i)    = fscanf (fp,'%d',1); % Elementtyp 
    if (mtyp(i)==3)
        ntri=3; %Dreieckselement
    end
    if (mtyp(i)==8)
        ntri=4; % Viereckselement
    end
    nelanz(i)  = fscanf (fp,'%d',1); % Anzahl der jeweiligen Elemente
    selanz(i)  = fscanf (fp,'%d',1); % Anzahl der Elementknoten/Element
    nel        = nel+nelanz(i); % Gesamtzahl der Elemente
end
sprintf('Reading in %d Nodes and %d Elements',sanz,nel);
% Knoten einlesen
for i=1:sanz
    snr(i) = fscanf (fp,'%d',1); % Knotennummer..
    sx(i)  = fscanf (fp,'%f',1); % x-Koordinate
    sy(i)  = fscanf (fp,'%f',1); % y-Koordinate
end
sxp=zeros(sanz,1);
syp=zeros(sanz,1);
snrp=zeros(sanz,1);
sprintf('rearranging node numbers \n');

[snrs,perm]=sort(snr);
for i=1:sanz
    sxp(i)=sx(snr(i));
    syp(i)=sy(snr(i));
    snrp(i)=snr(snr(i));
end
% Element Felder definieren
msel=max(max(selanz));
nrel=zeros(nel,msel);
fntri=find(selanz==ntri);
nelem=nelanz(fntri); % Anzahl der Hauptelemente
TRI=zeros(nelem,ntri);n=0;
vertices=zeros(sanz,2);
for i=1:typanz
    for j=1:nelanz(i)
        n=n+1;
        for k=1:selanz(i)
            a = fscanf (fp,'%d',1);
            nrel(n,k) = a;
            if (selanz(i)==ntri)
                TRI(j,k)=a;
            end
        end
    end
end
fclose (fp);
vertices(:,1)=sxp;
vertices(:,2)=syp;

%%
% Einlesen der Modellwerte (Knotenbasiert)

% [modfile,pdir] = uigetfile('*','Model file (from CRTomo)');
% cd (pdir);
fp=fopen(modfile,'r');
line=fgetl(fp);
nm=sscanf(line,'%d',1);
if (nm~=nelem)
    sprintf('There seems something wrong since the Element numbers %d \n',nelem);
    sprintf('do not match the number of Model cells %d !!!\n',nm); 
end
rho(1:nm)=0;phase(1:nm)=0;
for i=1:nm
    rho(i)=fscanf(fp,'%f',1);
    phase(i)=fscanf(fp,'%f',1);
end
fclose(fp);
%%
% Einlesen der Elektroden
fp=fopen(elecfile,'r');
line=fgetl(fp);
nelec=sscanf(line,'%d',1);
elecpos=zeros(nelec,2);
for i=1:nelec
    ie=fscanf(fp,'%d',1);
    elecpos(i,1)=sxp(ie);
    elecpos(i,2)=syp(ie);
end
fclose(fp);
%%
% open figure with name
name=sprintf('CRTomo model');
if checkme ~= 0
    name = sprintf('%s %s',name,fenster);
end
fig=figure('Name',name,'Numbertitle','off');

%%
% Plotten der Potentiale mit fill(X,Y,Z,C)
clf;
% plot
%trisurf(TRI,Vx,Vy,volt1,'edgecolor','none');
title(name,'fontsize',fns);

romin=min(rho); % absolutes minimum
romax=max(rho); % absolutes maximum
if (romin==romax)
    rho(round(nm/4))=1;(romin+romax)/2;
    rho(round(nm/2))=100;(romin+romax)/2;
end
if logme~=0
    rolog=log10(rho); % log trafo
else
    rolog=rho;
end

rolmin=min(rolog); % absolutes minimum
rolmax=max(rolog); % absolutes maximum

if (cmin==0)
    cmin=rolmin;
end
if (cmax==0)
    cmax=rolmax;
end
sprintf('Plot range:: %f\t%f\n',cmin,cmax);
patch('Faces',TRI,'Vertices',vertices,'CData',rolog','FaceColor','flat','Edgecolor','none')
caxis([cmin cmax])

%set(gca,'DataAspectRatio',[1 1 1],'fontsize',fns,'TickDir','out')
set(gca,'fontsize',fns,'TickDir','out')
xlabel('x [m]','fontsize',fns)
ylabel('z [m]','fontsize',fns)
axis tight
%if logme~=0
%    % Colortable neu definieren
%    ni=(rolmax-rolmin)/nd; %increment
%    
%    S=char(ones(nd+1,5));
%    for i=1:nd+1
%        if (logme~=0)
%            S(i,:)=sprintf('%5.f',10^(rolmin+(i-1)*ni));
%        else
%            S(i,:)=sprintf('%5.f',rolmin+(i-1)*ni);
%        end
%    end 
%    ct=cellstr(S);
%end

%colormap(jet(2*nd+1));
h=colorbar('vert');

set(h,'fontsize',fns)
set(h,'XaxisLocation','top')
if logme~=0
%    set(h,'YlimMode','manual');
%    set(h,'Ylim',[10^cmin 10^cmax]);
%    set(h,'Ytick','log')
%    set(h,'Yscale','log')
%    set(h,'YTickLabelMode','manual')
%    set(h,'YTickLabel',ct)
  set(get(h,'xlabel'),'String','log_{10}(\rho) [\Omega m]','fontsize',fns)
else
  set(get(h,'xlabel'),'String','\rho [\Omega m]','fontsize',fns)
end

view(az,el);
hold on
for i=1:nelec
    plot(elecpos(i,1),elecpos(i,2),elecmark,...
        'MarkerEdgeColor','k','MarkerFaceColor','k',...
        'MarkerSize',marksize);
end
if (saveplot==1)
    %pause
    [p,fls,app,m]=fileparts(modfile);
    fleps=strcat(fls,app,'.eps');
    set(gcf,'PaperPositionMode','auto');
    print('-depsc2','-r400',fleps)
%    print('-dpdf','-r400',flpdf) %%%geht nicht
    %print('-depsc2','-r600',file_save)
end
close all
clear all
