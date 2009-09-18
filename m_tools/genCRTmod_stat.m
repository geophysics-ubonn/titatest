%%

close all;
% Definitions
elecmark='bo';
marksize=5;
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
disp(sprintf('Reading in %d Nodes and %d Elements',sanz,nel));
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

    disp(sprintf('rearranging node numbers \n'));

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
% Einlesen der Elektroden
[elecfile,pdir] = uigetfile('*','Electrode file (from CRTomo)');
cd (pdir);
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
% Perparing model section..
nm=nelem;
rho(1:nm)=100;
midx(1:nm)=0;midy(1:nm)=0;
nx=0;
ny=0;
for i=1:nm
    spx=0;spy=0;
    for k=1:ntri
        spx=spx+sxp(nrel(i,k));
        spy=spy+syp(nrel(i,k));
    end
    midx(i)=spx/ntri;
    midy(i)=spy/ntri;
    if (midx(i)==midx(1))
        ny=ny+1;
    end
    if (midy(i)==midy(1))
        nx=nx+1;
    end
end
sprintf('nx:%d\n ny:%d\n nx*ny:%d nm:%d\n',nx,ny,nx*ny,nm)
%%
% open figure with nam%%
% Now statistical model
logme=1;
saveplot=1;

[exim,myepscr]=system('which my_epscrop');
myepscr=strcat(myepscr,' ./');

scrsz = get(0,'ScreenSize');
myfig=figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)/2],'Numbertitle','off');

hs=100;
deca=3;

if (isunix), dirtrail='/'; else; dirtrail='\';end
decadir=sprintf('%s%d_decades%s',pdir,deca,dirtrail);
exi=exist(decadir,'dir');
if (exi~=0)
    rmdir(decadir,'s');
end
mkdir(decadir);
cd(decadir);

lhs=log10(hs);
if (logme~=0)
    cmin = lhs-deca/2;
    cmax = lhs+deca/2;
else
    cmin=10^(lhs-deca/2);
    cmax=10^(lhs+deca/2);
end
    
    
xmi=min(midx);
xma=max(midx);
ymi=min(abs(midy));
yma=max(abs(midy));
dxi=(xma-xmi)/(nx-1);
dyi=(yma-ymi)/(ny-1);
x=[0:1:(nx-1)].*dxi+xmi;
y=[0:1:(ny-1)].*dyi+ymi;

mean_r=0;
Ix=[1 1 5 5];
Iy=[1 0.1 5 0.5];

V=visim_init(x,y);
V.rseed=2;
V.Va.ang1=90;   % Rotation angle of dip(clockwise from north)
V.gmean=mean_r;
V.gvar=1;

for iscale=1:length(Ix)
    mydir=sprintf('%sIx_%.1f_Iy_%.1f%s',decadir,Ix(iscale),Iy(iscale),dirtrail);
    mystr=sprintf('Ix_%.1f_Iy_%.1f',Ix(iscale),Iy(iscale));
    if (exist(mydir,'dir')~=0)
        disp(sprintf('remove %s\n',mydir));
        rmdir(mydir,'s')
    end
    mkdir(mydir);
    V.Va.a_hmax=Ix(iscale); % maximum correlation length
    V.Va.a_hmin=Iy(iscale); % minimum correlation length
    
    for its=1:4
        clf
        cd(mydir);
        V.Va.it=its;
        V=visim(V);    % run visim;
        % plot
        fp=fopen('visim.dat','r');
        visd=fscanf(fp,'%f',nm);
        fclose(fp);
        % plot
        vmin=min(visd);
        vmax=max(visd);
        
        d=deca/(vmax-vmin);
        
        rho=10.^(lhs+visd*d);
        
        [mytit,ctype]=visim_format_variogram(V);
        
        pldir=sprintf('%s%s%s',mydir,ctype,dirtrail);
        mkdir(pldir);
        cd(pldir);

        romin=min(rho); % absolutes minimum
        romax=max(rho); % absolutes maximum
        
        if (romin==romax)
            rho(round(nm/4))=romin+(romin+romax)/2;
            rho(round(nm/2))=romax-(romin+romax)/2;
        end
        if (logme~=0)
            rolog= log10(rho); % log trafo
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
        
        disp(sprintf('Plot range:: %f\t%f\n',cmin,cmax));
        patch('Faces',TRI,'Vertices',vertices,'CData',rolog,'FaceColor','flat','Edgecolor','none')
        caxis([cmin cmax])
        
        %set(gca,'DataAspectRatio',[1 1 1],'fontsize',fns,'TickDir','out')
        set(gca,'fontsize',fns,'TickDir','out')
        xlabel('x [m]','fontsize',fns)
        ylabel('z [m]','fontsize',fns)
        axis tight
        
        %if logme~=0
        % Colortable neu definieren
        %     ni=(rolmax-rolmin)/nd; %increment
        %     rostep = (romax-romin)/nd;
        %     S=char(ones(nd,5));
        %     for i=1:nd
        %         if (logme~=0)
        %             S(i,:)=sprintf('%5.f',10^(rolmin+(i-1)*ni));
        %         else
        %             S(i,:)=sprintf('%5.f',rolmin+(i-1)*ni);
        %         end
        %     end
        %     ct=cellstr(S);
        % end
        
        mytit=sprintf('%s / \\rho_{min}=%.f, \\rho_{max}=%.f, \\rho_{mean}=%.f',mytit,romin,romax,mean(rho));
        title(mytit,'fontsize',fns) ;
        
        %colormap(jet(2*nd+1));
        h=colorbar('vert');
        if (logme~=0)
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
        set(h,'fontsize',fns)
        set(h,'XaxisLocation','top')
        view(az,el);
        hold on
        for i=1:nelec
            plot(elecpos(i,1),elecpos(i,2),elecmark,...
                'MarkerEdgeColor','k','MarkerFaceColor','k',...
                'MarkerSize',marksize);
        end
        
        name=sprintf('%s_%s.model',mystr,ctype);
        fp=fopen(name,'w');
        fprintf(fp,'%d\n',nm);
        for i=1:nm
            fprintf(fp,'%f\t%f\n',rho(i),0.0);
        end
        fclose(fp);
        name=sprintf('%s_%s.mag',mystr,ctype);
        fp=fopen(name,'w');
        fprintf(fp,'%d\n',nm);
        for i=1:nm
            fprintf(fp,'%f\t%f\t%f\n',midx(i),midy(i),log10(rho(i)));
        end
        fclose(fp);
        if (saveplot==1)
            fls=sprintf('%s_%s',mystr,ctype);
            fleps=strcat(fls,'.eps');
            set(myfig,'PaperPositionMode','auto');
            print('-depsc2','-r400',fleps)
            %    print('-dpdf','-r400',flpdf) %%%geht nicht
            %print('-depsc2','-r600',file_save)
            if (exim==0), unix(myepscr); end
        end
    end
end
disp(sprintf('done'));