function [] = migrate_scan(mode,writeconf)

%mode=drawIntra
%mode=drawExtra
%mode=drawZi
%mode=drawPM %PM for pupil mask
%mode=drawCi  %Ci are the coefficients of polynomial representing the off-axis distortion

d=5;
g=2; %gap=2
t=g*d-1;
nz = 22;
gs = 401; %grid sampling
pn = 10;
flow=1.07;
fhigh =1.3;

ifigure=0;
%% draw intra focal images 
if (strcmp(mode,'drawIntra')==1)
    ifigure=ifigure+1;figure(ifigure); clf;
    for i=1:g:t
        for j=1:g:t
            k=(d-(i+1)/2)*d+(j+1)/2;
            filename=['wcs/conf/migrate/Image_' num2str(i) '_' num2str(j) '_intra.txt'];
            im=load(filename);
            im=flipud(im);
            subplot(d,d,k);
            imagesc(im);axis xy;axis square;axis off;
            fitswrite(im,['~/temp/image_' sprintf('%02d',((t-i+1)*d+j)) '.fits']);
        end
    end
end

%% draw extra focal images 
if (strcmp(mode,'drawExtra')==1)
    ifigure=ifigure+1;figure(ifigure); clf;
    for i=1:g:t
        for j=1:g:t
            k=(d-(i+1)/2)*d+(j+1)/2;
            im=load(['wcs/conf/migrate/Image_' num2str(i) '_' num2str(j) '_extra.txt']);
            im=flipud(im);
            subplot(d,d,k);
            imagesc(im);axis xy;axis square;axis off;
        end
    end
end

%% draw the Zernike coefficients along the 45 degree
if (strcmp(mode,'drawZi')==1)
    ifigure=ifigure+1;figure(ifigure);clf;
    zi=zeros(t,nz);
    for i=1:t
        zi(i,:)=load(['wcs/conf/migrate/z_' num2str(i) '_' num2str(i) '.txt']);
    end
    for i=1:nz
        dd=ceil(sqrt(nz));
        subplot(dd,dd,i);
        plot(1:t,zi(:,i),'.-r');
    end
end

%% draw ca,ra,cb,rb along 45 degree line
if (strcmp(mode,'drawPM')==1)
    ca=zeros(t,1);
    ra=zeros(t,1);
    cb=zeros(t,1);
    rb=zeros(t,1);
    vigR=zeros(t,1);
    ifigure=ifigure+1;figure(ifigure);clf;
    for i=1:t
        pm=zeros(gs,gs);
%         a=load(['wcs/conf/migrate/Grid_' num2str(i) '_' num2str(i) '_in.txt']);
        a=load(['~/wavefront/TestImages/LSST_chip_scan/Grid_' num2str(i) '_' num2str(i) '_in.txt']);
        pm(int32( (((1+a(:,3))*(gs-1)/2)-1)*gs+(1+a(:,4))*(gs-1)/2 ))=1;
        vigR(i)=sum(sum(pm));
        
        %calculate ca,ra
        lineout=diag(pm);
        x0=(find(lineout((gs-1)/2+1:end), 1, 'last' )+0.5)/((gs-1)/2);
        x1=0.1;
        lineout=pm(:,x1*(gs-1)/2+(gs-1)/2+1);
        y1=(find(lineout((gs-1)/2+1:end),1,'last')+0.5)/((gs-1)/2);
        ca(i)=0.7071*(x1^2+y1^2-2*x0^2)/(x1+y1-2*x0);
        ra(i)=1.414*x0-ca(i);
        
        %calculate cb,rb
        lineout=diag(pm);
        x0=(find(lineout(1:(gs-1)/2+1), 1, 'last' )+0.5)/((gs-1)/2)-1;
        x1=-0.15;
        lineout=pm(:,x1*(gs-1)/2+(gs-1)/2+1);
        y1=(find(lineout(1:(gs-1)/2+1),1,'last')+0.5)/((gs-1)/2)-1;
        cb(i)=0.7071*(x1^2+y1^2-2*x0^2)/(x1+y1-2*x0);
        rb(i)=cb(i)-1.414*x0;
        
        %draw and check pupil shapes
        r=0.61;
        rr=1;
        th=0:360;th=th/180*pi;
        xx=rr*cos(th);
        yy=rr*sin(th);
        x=r*cos(th);
        y=r*sin(th);
        
        th=53:127;
        th=th/180*pi;
        th=th-45./180*pi;
        x1=ra(i)*cos(th)+ca(i)*0.7071;
        y1=ra(i)*sin(th)+ca(i)*0.7071;
        
        th=18:162;th=(th+180)/180*pi;
        th=th-45./180*pi;
        x2=rb(i)*cos(th)+cb(i)*0.7071;
        y2=rb(i)*sin(th)+cb(i)*0.7071;
        
        dd=ceil(sqrt(t));    subplot(dd,dd,i);
        %     subplot(t,t,(t-i)*t+i);
        myaxis=-1:2/(gs-1):1;
        imagesc(myaxis,myaxis,pm);axis xy;axis square;
        hold on; scatter(x,y,'g.');scatter(xx,yy,'g.');scatter(x1,y1,'g.');scatter(x2,y2,'g.');
        hold off;xlim([-1 1]);ylim([-1 1]);axis off;
        
    end
    ifigure=ifigure+1;figure(ifigure);
    subplot(2,2,1);plot(ca,'.-r');title('Center of Upper Left Circle (ca)');
    subplot(2,2,2);plot(ra,'.-r');title('Radius of Upper Left Circle (ra)');
    subplot(2,2,3);plot(cb,'.-r');title('Center of Lower Left Circle (cb)');
    subplot(2,2,4);plot(rb,'.-r');title('Radius of Lower Left Circle (rb)');
    

    if writeconf
        pmparam=zeros(t,5);
        pmparam(:,1)=linspace(flow,fhigh,t);
        pmparam(:,2)=ca;
        pmparam(:,3)=ra;
        pmparam(:,4)=cb;
        pmparam(:,5)=rb;
        save('wcs/conf/mask_migrate.txt','pmparam','-ascii');
    end
    
    %for normalization of vignetting ratio (vignetting function)
    pm=zeros(gs,gs);
    [xp yp]=meshgrid(-1:1/((gs-1)/2):1);    
    rp=sqrt(xp.^2+yp.^2);
    pm(rp>=r& rp<=1)=1;
    vigR=vigR/sum(sum(pm));
    
    ifigure=ifigure+1;figure(ifigure);
    plot(linspace(flow,fhigh,t),vigR,'-ro','LineWidth',2,'MarkerSize',8);
    xlim([flow fhigh]);
    ylim([0.5 1]);
    xlabel('field angle');
    ylabel('pupil vignetting ratio (by area)');
    if writeconf
        pmparam=zeros(t,2);
        pmparam(:,1)=linspace(flow,fhigh,t);
        pmparam(:,2)=vigR;
        save('wcs/conf/vig_func.txt','pmparam','-ascii');
    end
end

%% the polynomial coefficients:  c1-66  
if (strcmp(mode,'drawCi')==1)
    pn=(pn+1)*(pn+2)/2;
    cn=zeros(t,pn);
    nx=7;
    ny=8;
    % intra x'
    ifigure=ifigure+1;figure(ifigure);clf;
    for i=1:t    
        a1=load(['wcs/conf/migrate/Grid_' num2str(i) '_' num2str(i) '_in.txt']);
        
        c0=zeros(pn,1);
        c0(2)=0.5;
        x=a1(:,3)';
        y=a1(:,4)';
        xy = [x; y];
        z = a1(:,nx)'; %x prime
        
        cn(i,:)=lsqcurvefit(@poly10_2D,c0,xy,z);%,'Display','off');
        zfit=poly10_2D(cn(i,:),xy);
        dd=ceil(sqrt(t));    subplot(dd,dd,i);
        scatter3(a1(:,3),a1(:,4),100*(a1(:,nx)-zfit'),30,100*(a1(:,nx)-zfit'),'filled');xlabel('x');ylabel('y');zlabel('diff (in pixel)');
    end
    cnw=zeros(t,pn+2);
    fx =linspace(flow,fhigh,t);
    cnw(:,1:2)=[fx' fx'];
    cnw(:,3:pn+2)=cn;
    save('wcs/conf/cxin.txt','cnw','-ascii');
    ifigure=ifigure+1;figure(ifigure);clf;    
    for i=1:pn
        subplot(11,ceil(pn/11),i);
        plot(cn(:,i),'-ro','LineWidth',2,'MarkerEdgeColor','k', 'MarkerFaceColor','g', 'MarkerSize',8);
    end
    
    % intra y'
    ifigure=ifigure+1;figure(ifigure);clf;
    for i=1:t
        a1=load(['wcs/conf/migrate/Grid_' num2str(i) '_' num2str(i) '_in.txt']);
        
        c0=zeros(pn,1);
        c0(2)=0.5;
        x=a1(:,3)';
        y=a1(:,4)';
        xy = [x; y];
        z = a1(:,ny)'; %y prime
        
        cn(i,:)=lsqcurvefit(@poly10_2D,c0,xy,z);%,'Display','off');
        zfit=poly10_2D(cn(i,:),xy);
        dd=ceil(sqrt(t));    subplot(dd,dd,i);
        scatter3(a1(:,3),a1(:,4),100*(a1(:,ny)-zfit'),30,100*(a1(:,ny)-zfit'),'filled');xlabel('x');ylabel('y');zlabel('diff (in pixel)');
    end
    cnw=zeros(t,pn+2);
    fx =linspace(flow,fhigh,t);
    cnw(:,1:2)=[fx' fx'];
    cnw(:,3:pn+2)=cn;
    save('wcs/conf/cyin.txt','cnw','-ascii');
    ifigure=ifigure+1;figure(ifigure);clf;
    for i=1:pn
        subplot(11,ceil(pn/11),i);
        plot(cn(:,i),'-ro','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','g', 'MarkerSize',8);
    end
   % extra x'
   ifigure=ifigure+1;figure(ifigure);clf;
    for i=1:t    
        a1=load(['wcs/conf/migrate/Grid_' num2str(i) '_' num2str(i) '_ex.txt']);
        
        c0=zeros(pn,1);
        c0(2)=-0.5;
        x=a1(:,3)';
        y=a1(:,4)';
        xy = [x; y];
        z = a1(:,nx)'; %x prime
        
        cn(i,:)=lsqcurvefit(@poly10_2D,c0,xy,z);%,'Display','off');
        zfit=poly10_2D(cn(i,:),xy);
        dd=ceil(sqrt(t));    subplot(dd,dd,i);
        scatter3(a1(:,3),a1(:,4),100*(a1(:,nx)-zfit'),30,100*(a1(:,nx)-zfit'),'filled');xlabel('x');ylabel('y');zlabel('diff (in pixel)');
    end
    cnw=zeros(t,pn+2);
    fx =linspace(flow,fhigh,t);
    cnw(:,1:2)=[fx' fx'];
    cnw(:,3:pn+2)=cn;
    save('wcs/conf/cxex.txt','cnw','-ascii');
    ifigure=ifigure+1;figure(ifigure);clf;
    for i=1:pn
        subplot(11,ceil(pn/11),i);
        plot(cn(:,i),'-ro','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',8);
    end
    
    % extra y'
    ifigure=ifigure+1;figure(ifigure);clf;
    for i=1:t
        a1=load(['wcs/conf/migrate/Grid_' num2str(i) '_' num2str(i) '_ex.txt']);
        
        c0=zeros(pn,1);
        c0(2)=-0.5;
        x=a1(:,3)';
        y=a1(:,4)';
        xy = [x; y];
        z = a1(:,ny)'; %y prime
        
        cn(i,:)=lsqcurvefit(@poly10_2D,c0,xy,z);%,'Display','off');
        zfit=poly10_2D(cn(i,:),xy);
        dd=ceil(sqrt(t));    subplot(dd,dd,i);
        scatter3(a1(:,3),a1(:,4),100*(a1(:,ny)-zfit'),30,100*(a1(:,ny)-zfit'),'filled');xlabel('x');ylabel('y');zlabel('diff (in pixel)');
    end
    cnw=zeros(t,pn+2);
    fx =linspace(flow,fhigh,t);
    cnw(:,1:2)=[fx' fx'];
    cnw(:,3:pn+2)=cn;
    save('wcs/conf/cyex.txt','cnw','-ascii');
    ifigure=ifigure+1;figure(ifigure);clf;
    for i=1:pn
        subplot(11,ceil(pn/11),i);
        plot(cn(:,i),'-ro','LineWidth',2, 'MarkerEdgeColor','k', 'MarkerFaceColor','g', 'MarkerSize',8);
    end
        
    
end
end
