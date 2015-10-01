function [] = cornerCoordiTrans_chipCenter(mode)

% This routine is for modeling of the LSST off-axis distortion

%mode =
% drawEdge
% drawCorner
% fitEdgeIntra
% fitEdgeExtra
% fitCornerIntraPoly8
% fitCornerExtraPoly8
% fitCornerIntraPoly10
% fitCornerExtraPoly10
% fitEdgeIntraZ22
% fitEdgeExtraZ22
% fitEdgeIntraZA22
% fitEdgeExtraZA22

if (~isempty(strfind(mode,'Edge') ))
    a=load('wcs/conf/LSSTA_OffAxis_XY_Normal.txt');
%     a=load('raytracingData/LSSTA_OffAxis_XY_Normal_vignetting.txt');
elseif (~isempty(strfind(mode,'Corner')))
    a1=load('wcs/conf/comp_1.0mm_intra.txt');
    a2=load('wcs/conf/comp_1.0mm_extra.txt');
else
    fprintf('Wrong mode, please check!');
    return;
end

nx=7;
ny=8;
if (strncmpi(mode,'draw',4)==1)
 
    r=0.613;
    rr=1;
    t=0:360;t=t/180*pi;
    xx=rr*cos(t);
    yy=rr*sin(t);
    x=r*cos(t);
    y=r*sin(t);
    
    r1=1.23;
    t=53:0.2:127;
    offcenter=0.30;
    
    t=t/180*pi;
    x1=r1*cos(t);
    y1=r1*sin(t)-offcenter;
    if (~isempty(strfind(mode,'drawCorner')))
        t=t-45./180*pi;
        x1=r1*cos(t)-offcenter*0.7071;
        y1=r1*sin(t)-offcenter*0.7071;
    end
    
    r2=0.57;
    offcenter=0.11;
    
    t=18:162;t=(t+180)/180*pi;
    x2=r2*cos(t);
    y2=r2*sin(t)-offcenter;
    if (~isempty(strfind(mode,'drawCorner')))
        t=t-45./180*pi;
        x2=r2*cos(t)-offcenter*0.7071;
        y2=r2*sin(t)-offcenter*0.7071;
    end
    
    figure(5);clf;
    scatter(a1(:,3),a1(:,4),'.b');axis square;
    hold on; scatter(x,y,'r.');scatter(xx,yy,'r.');scatter(x1,y1,'r.');scatter(x2,y2,'r.');
    hold off;xlim([-1 1]);ylim([-1 1]);
    
%     return;
    
    r=0.613;
    rr=1;
    t=0:360;
    xx=rr*cos(t);
    yy=rr*sin(t);
    x=r*cos(t);
    y=r*sin(t);
    
    figure(1);clf;
    subplot(1,2,1);
    scatter(a1(:,3),a1(:,4),'.b');axis square;
    hold on; scatter(x,y,'r.');scatter(xx,yy,'r.');
    hold off;
    
    subplot(1,2,2);
    scatter(a1(:,nx),a1(:,ny),'.b');axis square;
    
    figure(2);clf;
    subplot(1,2,1);
   if (~isempty(strfind(mode,'drawEdge')))
       b=load('~/largeData/ZemaxImages/LSST_N/z4_0.00_intra.txt');
   elseif (~isempty(strfind(mode,'drawCorner')))
       b=load('~/largeData/ZemaxImages/LSST_NE/z4_0.00_intra.txt');
   end
    b=flipud(b);
    imagesc(b);axis square;axis xy;
    hold on;scatter(a1(:,nx)*100+64,a1(:,ny)*100+64);axis square;
    hold off;
    
    subplot(1,2,2);
    if (~isempty(strfind(mode,'drawEdge')))
        c=load('~/largeData/ZemaxImages/LSST_N/z4_0.00_extra.txt');
    elseif (~isempty(strfind(mode,'drawCorner')))
        c=load('~/largeData/ZemaxImages/LSST_NE/z4_0.00_extra.txt');
    end
    c=flipud(c);
    imagesc(c);axis square;axis xy;
    hold on;scatter(a2(:,nx)*100+64,a2(:,ny)*100+64);axis square;
    hold off;

    figure(3);clf;
    scatter3(a1(:,3),a1(:,4),a1(:,nx),'.r');xlabel('x');ylabel('y');
    
    figure(4);clf;
    scatter3(a1(:,3),a1(:,4),a1(:,ny),'.r');xlabel('x');ylabel('y');

    figure(5);clf;
    scatter3(a1(:,3),a1(:,4),a2(:,nx),'.r');xlabel('x');ylabel('y');
    
    figure(6);clf;
    scatter3(a2(:,3),a2(:,4),a2(:,ny),'.r');xlabel('x');ylabel('y');
   
    figure(7);clf;
    subplot(1,2,1); scatter3(a1(:,3),a1(:,4),a1(:,11)*100,30,a1(:,11)*100);xlabel('x');ylabel('y');zlabel('deviation of x^\prime from ideal focal point (in pixel)');
    subplot(1,2,2); scatter3(a1(:,3),a1(:,4),a1(:,12)*100,30,a1(:,12)*100);xlabel('x');ylabel('y');zlabel('deviation of y^\prime from ideal focal point (in pixel)');
%% fitting
elseif (strncmpi(mode,'fit',3)==1 && ~isempty(strfind(mode,'Intra')))

    if (~isempty(strfind(mode,'Poly8')))
        nTerms=(8+1)*(8+2)/2;
        c0=zeros(nTerms,1);
        c0(2)=0.5;
    elseif (~isempty(strfind(mode,'Poly10')))
        nTerms=(10+1)*(10+2)/2;
        c0=zeros(nTerms,1);
        c0(2)=0.5;
    elseif (~isempty(strfind(mode,'Z22'))  )
        nTerms=22;
    elseif (  ~isempty(strfind(mode,'ZA22')) )
        nTerms=22;
        e=0.61;
    end

    x=a1(:,3)';
    y=a1(:,4)';
    xy = [x; y];
    z = a1(:,nx)'; %x prime
    if (~isempty(strfind(mode,'Poly8')))
        c=lsqcurvefit(@poly8_2D,c0,xy,z);%,'Display','off');
        zfit=poly8_2D(c,xy);
    elseif (~isempty(strfind(mode,'Poly10')))
        c=lsqcurvefit(@poly10_2D,c0,xy,z);%,'Display','off');
        zfit=poly10_2D(c,xy);
    elseif (  ~isempty(strfind(mode,'ZA22')) )
        c=ZernikeAnnularFit(z,x,y,nTerms,e);
        zfit=ZernikeAnnularEval(c,x,y,e);
    end
    display(c);
    
    figure(3);clf;
    scatter3(a1(:,3),a1(:,4),a1(:,nx),'.r');xlabel('x');ylabel('y');
    hold on;
    scatter3(a1(:,3),a1(:,4),zfit,'+b');
    hold off;
    scatter3(a1(:,3),a1(:,4),100*(a1(:,nx)-zfit'),30,100*(a1(:,nx)-zfit'),'filled');xlabel('x');ylabel('y');zlabel('diff (in pixel)');
    for i=1:nTerms
        fprintf('%e ',c(i));
    end
    fprintf('\n');
    %

%     c0=[0 0 0.5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    if (~isempty(strfind(mode,'Poly8')))
        nTerms=(8+1)*(8+2)/2;
        c0=zeros(nTerms,1);
        c0(2)=0.5;
    elseif (~isempty(strfind(mode,'Poly10')))
        nTerms=(10+1)*(10+2)/2;
        c0=zeros(nTerms,1);
        c0(2)=0.5;
    elseif (~isempty(strfind(mode,'Z22')) ||  ~isempty(strfind(mode,'ZA22')) )
        nTerms=22;
        c0=zeros(nTerms,1);
    end
    
    xy = [a1(:,3)'; a1(:,4)'];
    z = a1(:,ny)'; %y prime
    if (~isempty(strfind(mode,'Poly8')))
        c=lsqcurvefit(@poly8_2D,c0,xy,z);%,'Display','off');
        zfit=poly8_2D(c,xy);
    elseif (~isempty(strfind(mode,'Poly10')))
        c=lsqcurvefit(@poly10_2D,c0,xy,z);%,'Display','off');
        zfit=poly10_2D(c,xy);  
    end
    display(c);

    figure(4);clf;
    scatter3(a1(:,3),a1(:,4),a1(:,ny),'.r');xlabel('x');ylabel('y');
    hold on;
    scatter3(a1(:,3),a1(:,4),zfit,'+b');
    hold off;
    scatter3(a1(:,3),a1(:,4),100*(a1(:,ny)-zfit'),30,100*(a1(:,ny)-zfit'),'filled');xlabel('x');ylabel('y');zlabel('diff (in pixel)');
    for i=1:nTerms
        fprintf('%e ',c(i));
    end
    fprintf('\n');
    
elseif (strncmpi(mode,'fit',3)==1 && ~isempty(strfind(mode,'Extra')))

    if (~isempty(strfind(mode,'Poly8')))
        nTerms=(8+1)*(8+2)/2;
        c0=zeros(nTerms,1);
        c0(2)=-0.5;
    elseif (~isempty(strfind(mode,'Poly10')))
        nTerms=(10+1)*(10+2)/2;
        c0=zeros(nTerms,1);
        c0(2)=-0.5;
    elseif (~isempty(strfind(mode,'Z22')) ||  ~isempty(strfind(mode,'ZA22')) )
        nTerms=22;
        c0=zeros(nTerms,1);
    end
    
    xy = [a2(:,3)'; a2(:,4)'];
    z = a2(:,nx)'; %x prime
    if (~isempty(strfind(mode,'Poly8')))
        c=lsqcurvefit(@poly8_2D,c0,xy,z);%,'Display','off');
        zfit=poly8_2D(c,xy);
    elseif (~isempty(strfind(mode,'Poly10')))
        c=lsqcurvefit(@poly10_2D,c0,xy,z);%,'Display','off');
        zfit=poly10_2D(c,xy);  
    end    
    display(c);
    
    figure(3);clf;
    scatter3(a2(:,3),a2(:,4),a2(:,nx),'.r');xlabel('x');ylabel('y');
    hold on;
    scatter3(a2(:,3),a2(:,4),zfit,'+b');
    hold off;
    scatter3(a2(:,3),a2(:,4),100*(a2(:,nx)-zfit'),30,100*(a2(:,nx)-zfit'),'filled');xlabel('x');ylabel('y');zlabel('diff (in pixel)');
    for i=1:nTerms
        fprintf('%e ',c(i));
    end
    fprintf('\n');
    %

    %     c0=[0 0 0.5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    if (~isempty(strfind(mode,'Poly8')))
        nTerms=(8+1)*(8+2)/2;
        c0=zeros(nTerms,1);
        c0(2)=-0.5;
    elseif (~isempty(strfind(mode,'Poly10')))
        nTerms=(10+1)*(10+2)/2;
        c0=zeros(nTerms,1);
        c0(2)=-0.5;
    elseif (~isempty(strfind(mode,'Z22')) ||  ~isempty(strfind(mode,'ZA22')) )
        nTerms=22;
        c0=zeros(nTerms,1);
    end
    
    xy = [a2(:,3)'; a2(:,4)'];
    z = a2(:,ny)'; %y prime
    if (~isempty(strfind(mode,'Poly8')))
        c=lsqcurvefit(@poly8_2D,c0,xy,z);%,'Display','off');
        zfit=poly8_2D(c,xy);
    elseif (~isempty(strfind(mode,'Poly10')))
        c=lsqcurvefit(@poly10_2D,c0,xy,z);%,'Display','off');
        zfit=poly10_2D(c,xy);
    end
    display(c);
    
    figure(4);clf;
    scatter3(a2(:,3),a2(:,4),a2(:,ny),'.r');xlabel('x');ylabel('y');
    hold on;
    scatter3(a2(:,3),a2(:,4),zfit,'+b');
    hold off;
    scatter3(a2(:,3),a2(:,4),100*(a2(:,ny)-zfit'),30,100*(a2(:,ny)-zfit'),'filled');xlabel('x');ylabel('y');zlabel('diff (in pixel)');
    for i=1:nTerms
        fprintf('%e ',c(i));
    end
    fprintf('\n');
    
end
