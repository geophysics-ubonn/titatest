%%
clear all;
close all;
% Definitions
az=0; % Azimuth for the plot 
el=90; % elongation for thew plot
fns=14; % font size
cmin=0;cmax=0; % plot range can be set here, 0 automatic choosed
logme=0; % linear plots -> 0 logarithmic else
saveplot=1; % if save the plot adjust and press space

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
sprintf('Reading in %d Nodes and %d Elements',sanz,nel);
% Knoten einlesen
for i=1:sanz
    snr(i) = fscanf (fp,'%d',1); % Knotennummer..
    sx(i)  = fscanf (fp,'%f',1); % x-Koordinate
    sy(i)  = fscanf (fp,'%f',1); % y-Koordinate
end
sxp(1:sanz)=sx;
syp(1:sanz)=sy;
snrp=snr;
if (snr(1)~=1)

    sprintf('rearranging node numbers \n');

    [snrs,perm]=sort(snr);
    for i=1:sanz
        sxp(i)=sx(snr(i));
        syp(i)=sy(snr(i));
        snrp(i)=snr(snr(i));
    end
end
% Element Felder definieren
msel=max(max(selanz));
nrel=zeros(nel,msel);
fntri=find(selanz==ntri);
nelem=nelanz(fntri); % Anzahl der Hauptelemente
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
vertices(:,1)=sxp;
vertices(:,2)=syp;

%%
% Einlesen der Potentialwerte (Knotenbasiert)
[modfile,pdir] = uigetfile('*','Model file (from CRTomo)');
cd (pdir);
fp=fopen(modfile,'r');
line=fgetl(fp);
nm=sscanf(line,'%d',1);
if (nm~=nelem)
    sprintf('There seems something wrong since the Element numbers %d \n',nelem);
    sprintf('do not match the number of Model cells %d !!!\n',nm); 
end
rho(1:nm)=0;phase(1:nm)=0;midx(1:nm)=0;midy(1:nm)=0;
for i=1:nm
    rho(i)=fscanf(fp,'%f',1);
    phase(i)=fscanf(fp,'%f',1);
    spx=0;spy=0;
    for k=1:ntri
        spx=spx+sxp(nrel(i,k));
        spy=spy+syp(nrel(i,k));
    end
    midx(i)=spx/ntri;
    midy(i)=spy/ntri;
end
%%
rho(:)=100.;
mx=0.36511864859458693;
my=0.36701755939812991;
auge1m=[0.23 0.55]; % mittelpunkt erstes "auge"
auge2m=[0.47 0.55];
mundx=[0.2 0.5];
mundy=[0.17 0.33];

auge1m=[0.23-mx 0.55-my]; % mittelpunkt erstes "auge"
auge2m=[0.47-mx 0.55-my];
mundx=[0.2 0.5]-mx;
mundy=[0.17 0.33]-my;

ra=0.05;

count=0;
for i=1:nm
    
    txa1=midx(i)-auge1m(1);
    tya1=midy(i)-auge1m(2);

    txa2=midx(i)-auge2m(1);
    tya2=midy(i)-auge2m(2);
    
    rma1=sqrt(txa1^2+tya1^2);
    rma2=sqrt(txa2^2+tya2^2);
    lx=( mundx(1) <= midx(i) );
    lx=( lx && ( midx(i)<=mundx(2) ) );
    ly=( mundy(1) <= midy(i) );
    ly=( ly && ( midy(i)<=mundy(2) ) );
    lt=( lx && ly );
    if ( ( rma1 <= ra ) || ( rma2 <= ra ) || lt)
        sprintf('Element %d mittelpunkt %f %f hat abstand %f %f'...
            ,i,midx(i),midy(i),rma1,rma2)
        count=count+1;
        rho(i)=1.0;
    end
end
sprintf('Counted %d\n',count)
%%
% open figure with name
name=sprintf('CRTomo model')
figure('Name',name,'Numbertitle','off');

%%
% Plotten der Potentiale mit fill(X,Y,Z,C)
clf;
% plot
%trisurf(TRI,Vx,Vy,volt1,'edgecolor','none');
title(name,'fontsize',fns);

romin=min(rho); % absolutes minimum
romax=max(rho); % absolutes maximum

if (romin==romax) 
    rho(round(nm/4))=romin;(romin+romax)/2;
    rho(round(nm/2))=romax;(romin+romax)/2;
end
if (logme~=0)
    nd=8; %Anzahl der Ticks
    rolog=log10(rho); % log trafo
else
    nd=10; %Anzahl der Ticks
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
sprintf('Plot range:: %f\t%f\n',cmin,cmax)
patch('Faces',TRI,'Vertices',vertices,'CData',rolog','FaceColor','flat')%,'Edgecolor','none')
caxis([cmin cmax])

%set(gca,'DataAspectRatio',[1 1 1],'fontsize',fns,'TickDir','out')
set(gca,'fontsize',fns,'TickDir','out')
xlabel('x [m]','fontsize',fns)
ylabel('z [m]','fontsize',fns)

if logme~=0
    % Colortable neu definieren
    ni=(rolmax-rolmin)/nd; %increment
    
    S=char(ones(nd+1,5));
    for i=1:nd+1
        if (logme~=0)
            S(i,:)=sprintf('%5.f',10^(rolmin+(i-1)*ni));
        else
            S(i,:)=sprintf('%5.f',rolmin+(i-1)*ni);
        end
    end 
    ct=cellstr(S);
end

%colormap(jet(2*nd+1));
h=colorbar('vert');
if (logme~=0)
    set(h,'YtickMode','manual')
    set(h,'YTickLabel',ct)
end
set(h,'fontsize',fns)
set(h,'XaxisLocation','top')
set(get(h,'xlabel'),'String','\rho [\Omega m]','fontsize',fns)
view(az,el);

fp=fopen('gen_cr.model','w');
fprintf(fp,'%d\n',nm)
for i=1:nm
    fprintf(fp,'%f\t%f\n',rho(i),0.0);
end
fclose(fp);
fp=fopen('gen_cr.mag','w');
fprintf(fp,'%d\n',nm)
for i=1:nm
    fprintf(fp,'%f\t%f\t%f\n',midx(i),midy(i),log10(rho(i)));
end
fclose(fp);
saveplot=1
if (saveplot==1)
    pause
    fls='my_cr_mod';
    fleps=strcat(fls,'.eps');
    set(gcf,'PaperPositionMode','auto');
    print('-depsc2','-r400',fleps)
%    print('-dpdf','-r400',flpdf) %%%geht nicht
    %print('-depsc2','-r600',file_save)
end
