%%
% Definitions

elecmark='bo';
marksize=5;
az=0; % Azimuth for the plot 
el=90; % elongation for thew plot
fns=16; % font size
saveplot=1; % if save the plot adjust and press space

fp=fopen('tmp.meshname','r');
gridelem=fscanf(fp,'%s',1);
fclose(fp);

fp=fopen('tmp.elecname','r');
elecfile=fscanf(fp,'%s',1);
fclose(fp);

checkme = exist ('inv.lastmod','file');
if checkme~=0
  fp=fopen('inv.lastmod','r');
modfile=fscanf(fp,'%s',1);
fclose(fp);
 else
   
 end

fenster='';
checkme = exist ('tmp.fenster','file');
if checkme~=0
    fp=fopen('tmp.fenster','r');
    fenster=fgetl(fp);
    fclose(fp);
end
fenstert='';
checkme = exist ('tmp.fenstert','file');
if checkme~=0
    fp=fopen('tmp.fenstert','r');
    fenstert=fgetl(fp);
    fclose(fp);
end
cmin=0;camx=0;
checkme = exist ('tmp.crange','file');
if checkme~=0
    fp=fopen('tmp.crange','r');
    cmin=fscanf(fp,'%f',1);
    cmax=fscanf(fp,'%f',1);
    fclose(fp);
end
xmin=0;xmax=0;
checkme = exist ('tmp.xrange','file');
if checkme~=0
    fp=fopen('tmp.xrange','r');
    xmin=fscanf(fp,'%f',1);
    xmax=fscanf(fp,'%f',1);
    fclose(fp);
end
ymin=0;ymax=0;
checkme = exist ('tmp.yrange','file');
if checkme~=0
    fp=fopen('tmp.yrange','r');
    ymin=fscanf(fp,'%f',1);
    ymax=fscanf(fp,'%f',1);
    fclose(fp);
end
scrsz=[];
cbarn='';
checkme = exist ('tmp.cbarn','file');
if checkme~=0
    fp=fopen('tmp.cbarn','r');
    cbarn=fgetl(fp);
    scrsz=fscanf(fp,'%d',1);
    fclose(fp);
end
%%
% Knoten und Elemente einlesen 
% Filenamen bestimmen
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

[p,fls,appi,m]=fileparts(modfile);
fp=fopen(modfile,'r');
line=fgetl(fp);
nm=sscanf(line,'%d',1);
if (nm~=nelem)
    sprintf('There seems something wrong since the Element numbers %d \n',nelem);
    sprintf('do not match the number of Model cells %d !!!\n',nm); 
end
pha=strmatch(appi,'.pha');
mag=strmatch(appi,'.mag');
modl=strmatch(appi,'.modl');

if length(pha) ~= 0 | length(mag) ~= 0
    rho(1:nm)=0;
    for i=1:nm
        a=fscanf(fp,'%f',2);
        rho(i)=fscanf(fp,'%f',1);
    end
else
    rho(1:nm)=0;
    for i=1:nm
        rho(i)=fscanf(fp,'%f',1);
        a=fscanf(fp,'%f',1);
    end
end
fclose(fp);
%%
if length(cbarn) == 0
    if length(pha) ~= 0
        cbarn='Phase -[mRad]';
    elseif length(mag) ~= 0
        cbarn='log_{10}(\rho) [\Omega m]';
    elseif length(modl) ~= 0
        cbarn='\rho} [\Omega m]';  
    else
        cbarn='Parameter [Units]';
    end
end
fp=fopen('tmp.cbarn','w');
fprintf(fp,'%s \n',cbarn);
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
if length(fenster) ~= 0
    name = fenster;
end
if length(fenstert) ~= 0
    name = sprintf('%s\n%s',name,fenstert);
end
fp=fopen('tmp.fenster','w');
fprintf(fp,'%s \n',name);
fclose(fp);

%scrsz = get(0,'ScreenSize');
%fig=figure('Name',name,'Position',[1 scrsz(4) scrsz(3) scrsz(4)],'Numbertitle','off');
%if scrsz ~= 0
%else
fig=figure('Name',name,'Numbertitle','off');
%end
%%
% Plotten
clf;
title(name,'fontsize',fns);

romin=min(rho); % absolutes minimum
romax=max(rho); % absolutes maximum

if cmin~=cmax
    climits=[cmin cmax];
else
    climits=[romin romax];
end

patch('Faces',TRI,'Vertices',vertices,'CData',rho','FaceColor','flat','Edgecolor','none')
%patch('Faces',TRI,'Vertices',vertices,'CData',rho','FaceColor','flat')
caxis(climits)

set(gca,'fontsize',fns,'TickDir','out')
xlabel('x [m]','fontsize',fns)
ylabel('z [m]','fontsize',fns)
axis tight
axis image
if xmin~=xmax
    xlimits=[xmin xmax];
else
    xlimits=xlim;
end
xlim(xlimits);

if ymin~=ymax
    ylimits=[ymin ymax];
else
    ylimits=ylim;
end
ylim(ylimits);

h=colorbar('vert');

set(h,'fontsize',fns)
set(h,'XaxisLocation','top')
%cbarn=sprintf('$$\\mathsf{%s}$$',cbarn);
cbarn2=sprintf('\n%s\n',cbarn);
%set(get(h,'xlabel'),'interpreter','latex','String',cbarn2,'fontsize',fns)
set(get(h,'xlabel'),'String',cbarn2,'fontsize',fns)
view(az,el);
hold on
for i=1:nelec
    plot(elecpos(i,1),elecpos(i,2),elecmark,...
        'MarkerEdgeColor','k','MarkerFaceColor','k',...
        'MarkerSize',marksize);
end

%set(fig,'PaperPositionMode','auto');
print('-depsc2','-r400',strcat(fls,appi,'.eps'));
%print('-dpdf','-r400',strcat(fls,appi,'.pdf'));
print('-dpng','-r400',strcat(fls,appi,'.png'));
close (fig);

%write out some config files..
[p,fls,app,m]=fileparts(modfile);

fp=fopen('tmp.cbarn2','w');
fprintf(fp,'%s \n',cbarn2);
fclose(fp);

fp=fopen('tmp.crange','w');
fprintf(fp,'%f\t%f\n',climits);
fclose(fp);

fp=fopen('tmp.xrange','w');
fprintf(fp,'%f\t%f\n',xlimits);
fclose(fp);

fp=fopen('tmp.yrange','w');
fprintf(fp,'%f\t%f\n',ylimits);
fclose(fp);
