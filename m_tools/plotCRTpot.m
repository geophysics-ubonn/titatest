%%
% Definitions
az=0; % Azimuth for the plot 
el=90; % elongation for thew plot
fns=14; % font size
cmin=0;cmax=0; % plot range can be set here, 0 automatic choosed
logme=0; % linear plots -> 0 logarithmic else
saveplot=1;

%%
% Knoten und Elemente einlesen 

% Filenamen bestimmen
[gridelem,wdir] = uigetfile('*','Grid/element file (from CRTomo)');
ex     = cd (wdir);
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

% Knoten einlesen
for i=1:sanz
    snr(i) = fscanf (fp,'%d',1); % Knotennummer..
    sx(i)  = fscanf (fp,'%f',1); % x-Koordinate
    sy(i)  = fscanf (fp,'%f',1); % y-Koordinate
end

% Element Felder definieren
msel=max(max(selanz));
nrel=zeros(nel,msel);
nelem=nelanz(find(selanz==ntri)); % Anzahl der Hauptelemente
TRI=zeros(nelem,ntri);n=0;
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

%%
% Einlesen der Potentialwerte (Knotenbasiert)
[modfile,pdir] = uigetfile('*.dat','Potential file (from CRTomo)');
cd (pdir);
fp=fopen(modfile,'r');
Volt(1:sanz)=0;Vx(1:sanz)=0;Vy(1:sanz)=0;Phase(1:sanz)=0;
for i=1:sanz
    Vx(i)=fscanf(fp,'%f',1);
    Vy(i)=fscanf(fp,'%f',1);
    Volt(i)=fscanf(fp,'%f',1);
    Phase(i)=fscanf(fp,'%f',1);
end
%%
% Open Figure with name
name=sprintf('Calculated potential');
figure('Name',name,'Numbertitle','off');

%%
% Plotten der Potentiale
clf;
% Nullen raussschmeissen für die Trafo
[indzero]=Volt~=0;
wo_n_v=Volt(indzero);
[indzero]=Volt==0;
w_n_v=Volt(indzero);
% Log10 Tranformation
wo_n_v_amin=min(abs(wo_n_v)); % absolutes minimum
wo_n_v_amax=max(abs(wo_n_v)); % absolutes maximum
wo_n_v_n=wo_n_v/wo_n_v_amin; % potential normalized to minimum
wo_n_v_l=sign(wo_n_v_n).*log10(abs(wo_n_v_n)); % log trafo
w_n_v=ones(length(w_n_v))*min(wo_n_v_l); % die nullen aufs minimum projezieren

if (logme~=0)
    volt1=[wo_n_v_l;w_n_v]; % transformierte voltages zum plotten
else
    volt1=Volt;
end

% plot
trisurf(TRI,Vx,Vy,volt1)%,'edgecolor','none');
%trisurf(TRI,Vx,Vy,volt1);
title(name,'fontsize',fns);

%set(gca,'DataAspectRatio',[1 1 1],'fontsize',fns,'TickDir','out')
set(gca,'fontsize',fns,'TickDir','out')
xlabel('x [m]','fontsize',fns)
ylabel('z [m]','fontsize',fns)

% Colortable neu definieren
if logme~=0
    nd=4;
    lomin=log10(abs(wo_n_v_amin));
    lomax=log10(abs(wo_n_v_amax));
else
    lomax=max(abs(Volt));
    lomin=-lomax;
end

if (cmax==0) 
    cmax=lomax;
end
if logme~=0
    ni=(lomax-lomin)/nd; %increment
    S=char(ones(2*nd+1,8));
    for i=1:nd+1
        S(i,:)=sprintf('%8.2g',10^(lomax-(i-1)*ni));
    end
    for i=1:nd
        S(nd+i+1,:)=sprintf('%8.2g',-10^(lomin+i*ni));
    end
    
    ct=cellstr(S)
end
sprintf('Plot range:: %f\t%f\n',-cmax,cmax)

caxis([-cmax cmax]) ; 

%colormap(jet(2*nd+1));
h=colorbar('vert');
if (logme~=0)
    set(h,'YtickMode','manual')
    set(h,'YTickLabel',ct)
end
set(h,'fontsize',fns)
set(h,'XaxisLocation','top')
col_lab='[V]';
set(get(h,'xlabel'),'String',col_lab,'fontsize',fns)
view(az,el);

if (saveplot==1)
    pause
    fls='my_cr_pot';
    fleps=strcat(fls,'.eps');
    set(gcf,'PaperPositionMode','auto');
    print('-depsc2','-r400',fleps)
%    print('-dpdf','-r400',flpdf) %%%geht nicht
    %print('-depsc2','-r600',file_save)
end