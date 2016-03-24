%% Coded by
% Mohamed Mohamed El-Sayed Atyya
% mohamed.atyya94@eng-st.cu.edu.eg



% this program is used for showing a satalite about the earth
close all; clear all; clc;
%% constants
G=6.67428E-11; % gravitational constant
hz=15000; % simulation frequancy
step_sim=10; % simulation step
counter=0;  in=0;
%% Intial condition
M=[5.972*10^24,1000];       % [M1,M2]
R1=[0;0;0];                 % position of M(1)
R2=[8000e3;0;0];                 % position of M(2)
I=[0,0,0];              % location of initial axis
V1=[0;0;0];                 % velocity of M(1)
V2=[0;3;7.5]*1e3;                 % velocity of M(2)
% image (inter your URL of the image)
image_file = 'D:\4th year of Aerospace\1st\Orbital Mechanics\AER-427, Orbital Mechanics, Mohamed Mohamed Elsayed,SCE 2, BN 13  By MATLAB\week 11\satalite about earth with tracking/earth.jpg';
%% RK4 parameter
tf=24*3600*0+24*3600/8*1.13;   % final time of soution
dt=0.1*0+50;            % time step
X0=[R1;R2;V1;V2];
B=[0;0;0;0;0;0;0;0;0;0;0;0];
sol(1:12,1)=X0;
order=12;
%% RK4 solution
for n=1:length(0:dt:tf)
    b=G*M(2)/(norm(sol(1:3,n)-sol(4:6,n)))^3;
    c=-G*M(1)/(norm(sol(1:3,n)-sol(4:6,n)))^3;
    A=[0,0,0,0,0,0,1,0,0,0,0,0; ...
        0,0,0,0,0,0,0,1,0,0,0,0; ...
        0,0,0,0,0,0,0,0,1,0,0,0; ...
        0,0,0,0,0,0,0,0,0,1,0,0; ...
        0,0,0,0,0,0,0,0,0,0,1,0; ...
        0,0,0,0,0,0,0,0,0,0,0,1;...
        -b,0,0,b,0,0,0,0,0,0,0,0; ...
        0,-b,0,0,b,0,0,0,0,0,0,0; ...
        0,0,-b,0,0,b,0,0,0,0,0,0; ...
        -c,0,0,c,0,0,0,0,0,0,0,0; ...
        0,-c,0,0,c,0,0,0,0,0,0,0; ...
        0,0,-c,0,0,c,0,0,0,0,0,0 ];
    [ XX ] = RK4( A,B,sol(1:12,n),dt,n*dt,(n+1)*dt,order );
    sol(1:12,n+1)=XX(1:12,2);
end
R1_x=sol(1,:);
R1_y=sol(2,:);
R1_z=sol(3,:);
R2_x=sol(4,:);
R2_y=sol(5,:);
R2_z=sol(6,:);
V1_x=sol(7,:);
V1_y=sol(8,:);
V1_z=sol(9,:);
V2_x=sol(10,:);
V2_y=sol(11,:);
V2_z=sol(12,:);
%% projected path of the satalite on the earth
erot    = 7.2921158553e-5; % earth rotation rate (radians/sec)
for P=1:length(R1_x)
    points = intersectLineSphere([0,0,0,R2_x(P)-R1_x(P)+(P-1)*erot*dt,R2_y(P)-R1_y(P)+(P-1)*erot*dt,R2_z(P)-R1_z(P)+(P-1)*erot*dt], [0,0,0,6400e3]);
    XX1(P)=points(2,1);
    YY1(P)=points(2,2);
    ZZ1(P)=points(2,3);
end
%% calculating Right ascension and Declination
erot    = 7.2921158553e-5; % earth rotation rate (radians/sec)
for k=1:length(R1_x)
    [ RA(k), Dec(k) ] = lmn2RaDec( [R2_x(k), R2_y(k), R2_z(k)]-[R1_x(k), R1_y(k), R1_z(k)] );
    RA(k)=RA(k)+erot*(k-1)*dt*180/pi;
    while RA(k) > 360
        RA(k)=RA(k) - 360;
    end
end
%% Equatorial plane
    theta_vec=linspace(0,2*pi,30);
    r_vec=0:1e6:1.8e7;
    [theta_mat, r_mat]=meshgrid(theta_vec, r_vec);
    [x_mat, y_mat]=pol2cart(theta_mat, r_mat);
    z_mat=r_mat*0;
%% satellite tracking
figure(1);
% Texturemap the globe
% Load Earth image for texture map
cdata = imread(image_file);
J = imrotate(cdata,0,'bilinear');
v=1;
%--------------------------------------------------------------------------------------------------------------------------------------------------------
for p=1:step_sim:length(R1_x)
    set(0,'defaultfigureposition',[40 55 1300 600])
    subplot(1,2,1)
    cla;
    % Options
    space_color = 'k';
    npanels = 180;   % Number of globe panels around the equator deg/panel = 360/npanels
    alpha   = 1; % globe transparency level, 1 = opaque, through 0 = invisible
    % Earth texture image
    % Anything imread() will handle, but needs to be a 2:1 unprojected globe
    % Mean spherical earth
    erad    = 6371008.7714; % equatorial radius (meters)
    prad    = 6371008.7714; % polar radius (meters)
    %GMST0 = []; % Don't set up rotatable globe (ECEF)
    GMST0 = (4.89496121282306 - (p-1)*erot*dt); % Set up a rotatable globe at J2000.0
    set(gcf,'Color','w');
    % M2 trajectory
    plot3(I(1)*ones(1,length(R1_x))+R2_x-R1_x,I(2)*ones(1,length(R1_y))+R2_y-R1_y,I(3)*ones(1,length(R1_z))+R2_z-R1_z,'g','LineWidth',2);
    hold all;
    grid on;
    xlabel('X','Fontsize',18);
    ylabel('Y','Fontsize',18);
    zlabel('Z','Fontsize',18);
    title('Satellite about Earth','Fontsize',18);
    xlim auto;
    ylim auto;
    zlim auto;
    view(3);
    % M2 start
    plot3(R2_x(1,p)-R1_x(1,p),R2_y(1,p)-R1_y(1,p),R2_z(1,p)-R1_z(1,p),'o','color',[.9,0.2,0.5],'LineWidth',10);
    % Create wireframe globe
    % Create a 3D meshgrid of the sphere points using the ellipsoid function
    [x, y, z] = ellipsoid(0, 0, 0, erad, erad, prad, npanels);
    globe = surf(y+R1_x(1,p), -x+R1_y(1,p), -z+R1_z(1,p), 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]);
    if ~isempty(GMST0)
        hgx = hgtransform;
        set(hgx,'Matrix', makehgtform('zrotate',GMST0));
        set(globe,'Parent',hgx);
    end
    % Set image as color data (cdata) property, and set face color to indicate
    % a texturemap, which Matlab expects to be in cdata. Turn off the mesh edges.
    set(globe, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', alpha, 'EdgeColor', 'none');
    %Projected path of the satalite on the earth
    plot3(XX1,YY1,ZZ1,'r','LineWidth',2);
    % Equatorial plane
    if (R2_z(1,p)-R1_z(1,p)) >=0
        color='yellow';
    else
        color=[0.5,0.5,0.5];
    end
    surf(x_mat, y_mat, z_mat, 'FaceAlpha', 0.3, 'EdgeColor', color, 'FaceColor', color, 'EdgeAlpha', 0.3);
    legend('Satellite trajectory','Satellite','Projected path of the satalite on the earth','Equatorial plane','location','north');
    plot3([0,R2_x(p)],[0,R2_y(p)],[0,R2_z(p)],'b')
    
    
    
    % Satellite Tracking
    subplot(1,2,2)
    cla;
    image([0 360],[-90 90],cdata);
    title('Satellite Tracking','Fontsize',18);
    hold all;
    grid on;
    set(gca,'gridlinestyle','-')
    plot(RA(1:p), -Dec(1:p),'.','color','g','LineWidth',1)
    plot(RA(1), -Dec(1),'o','color','r','LineWidth',4)
    plot(RA(p), -Dec(p),'s','color',[.9,0.2,0.5],'LineWidth',2)
    %     scatter(RA(1:p), -Dec(1:p),'fill')
    %     gscatter (RA(1:p), -Dec(1:p))
    set(gca,'XTickLabel',{'-180','-129','-77','-26','26','77','129','180'})
    set(gca,'YTickLabel',{'80','60','40','20','0','-20','-40','-60','-80'})
    if RA(p) >= 204.2 && RA(p) <= 217.3 && Dec(p) <= 32.1 && Dec(p) >= 21.3 && in == 0
        counter = counter +1;
        in=1;
        visit=(p-1)*dt;
    elseif (RA(p) < 204.2 || RA(p) > 217.3) && (Dec(p) > 32.1 || Dec(p) < 21.3) && in == 1
        in=0;
        leave=(p-1)*dt;
    end
    xlim([0 360]);
    ylim([-90 90]);
    plot([204.2 217.3 217.3 204.2 204.2],[-21.3 -21.3 -32.1 -32.1 -21.3],'color','k','LineWidth',2)
    pause(1/hz);
end
%% effect of inclination angle (i)
figure(2);
[rz,Z]=max(R2_z-R1_z);
[ h, mag_h, i, omega, e, mag_e, w, theta ] = OrbitalElements ( [R2_x(Z)-R1_x(Z),R2_y(Z)-R1_y(Z),R2_z(Z)-R1_z(Z)],[V2_x(Z),V2_y(Z),V2_z(Z)],398600 );
set(gcf,'color','w');
% Texturemap the globe
% Load Earth image for texture map
cdata = imread(image_file);
J = imrotate(cdata,0,'bilinear');
v=1;
%--------------------------------------------------------------------------------------------------------------------------------------------------------
for p=0:5:180
    [ Q ] = RM ( omega*0, -(p-i), w*0, '313');
    EE2=Q*[R2_x;R2_y;R2_z];
    EE1=Q*[R1_x;R1_y;R1_z];
    EE0=Q*[XX1;YY1;ZZ1];
    set(0,'defaultfigureposition',[40 55 1300 600])
    for k=1:length(R1_x)
    [ RA(k), Dec(k) ] = lmn2RaDec( [EE2(1,k), EE2(2,k), EE2(3,k)]-[EE1(1,k), EE1(2,k), EE1(3,k)] );
    RA(k)=RA(k)+erot*(k-1)*dt*180/pi;
    while RA(k) > 360
        RA(k)=RA(k) - 360;
    end
    end
    % Satellite Tracking
    subplot(6,2,[2,4,6,8])
    cla;
    visit=[]; leave=[];
    v=1;
    image([0 360],[-90 90],cdata);
    title('Satellite Tracking','Fontsize',18);
    hold all;
    grid on;
    set(gca,'gridlinestyle','-')
    plot(RA, -Dec,'.','color','g','LineWidth',1)
    plot(RA(1), -Dec(1),'o','color','r','LineWidth',4)
    plot([204.2 217.3 217.3 204.2 204.2],[-21.3 -21.3 -32.1 -32.1 -21.3],'color','k','LineWidth',2)
    set(gca,'XTickLabel',{'-180','-129','-77','-26','26','77','129','180'})
    set(gca,'YTickLabel',{'100','80','60','40','20','0','-20','-40','-60','-80','-100'})
    for P=1:length(R2_x)
        if RA(P) >= 204.2 && RA(P) <= 217.3 && Dec(P) <= 32.1 && Dec(P) >= 21.3 && in == 0
            counter(v+1) = counter(v) +1;
            in=1;
            visit(v)=(P-1)*dt;
        elseif (RA(P) < 204.2 || RA(P) > 217.3) && (Dec(P) > 32.1 || Dec(P) < 21.3) && in == 1
            in=0;
            leave(v)=(P-1)*dt;
            v=v+1;
        end
    end
    xlim([0 360]);
    ylim([-100 100]);


    subplot(6,2,[1,3,5,7])
    cla;
    % Options
    space_color = 'k';
    npanels = 180;   % Number of globe panels around the equator deg/panel = 360/npanels
    alpha   = 1; % globe transparency level, 1 = opaque, through 0 = invisible
    % Earth texture image
    % Anything imread() will handle, but needs to be a 2:1 unprojected globe
    % Mean spherical earth
    erad    = 6371008.7714; % equatorial radius (meters)
    prad    = 6371008.7714; % polar radius (meters)
    %GMST0 = []; % Don't set up rotatable globe (ECEF)
    GMST0 = (4.89496121282306 - (p-1)*erot*dt); % Set up a rotatable globe at J2000.0
    set(gcf,'Color','w');
    % M2 trajectory
    plot3(I(1)*ones(1,length(R1_x))+EE2(1,:)-EE1(1,:),I(2)*ones(1,length(R1_y))+EE2(2,:)-EE1(2,:),I(3)*ones(1,length(R1_z))+EE2(3,:)-EE1(3,:),'g','LineWidth',2);
    hold all;
    grid on;
    xlabel('X','Fontsize',18);
    ylabel('Y','Fontsize',18);
    zlabel('Z','Fontsize',18);
    title(['Satellite about Earth with inclination angle = ' num2str(p)],'Fontsize',18);
    xlim auto;
    ylim auto;
    zlim auto;
    view(3);
%     view(-90,0);
    % M2 start
    plot3(EE2(1,1)-EE1(1,1),EE2(2,1)-EE1(2,1),EE2(3,1)-EE1(3,1),'o','color',[.9,0.2,0.5],'LineWidth',10);
    % Create wireframe globe
    % Create a 3D meshgrid of the sphere points using the ellipsoid function
    [x, y, z] = ellipsoid(0, 0, 0, erad, erad, prad, npanels);
    globe = surf(y, -x, -z, 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]);
    if ~isempty(GMST0)
        hgx = hgtransform;
        set(hgx,'Matrix', makehgtform('zrotate',GMST0));
        set(globe,'Parent',hgx);
    end
    % Set image as color data (cdata) property, and set face color to indicate
    % a texturemap, which Matlab expects to be in cdata. Turn off the mesh edges.
    set(globe, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', alpha, 'EdgeColor', 'none');
    %Projected path of the satalite on the earth
    plot3(EE0(1,:),EE0(2,:),EE0(3,:),'r','LineWidth',2);
    % Equatorial plane
    if (R2_z(1)-R1_z(1)) >=0
        color='yellow';
    else
        color=[0.5,0.5,0.5];
    end
    surf(x_mat, y_mat, z_mat, 'FaceAlpha', 0.3, 'EdgeColor', color, 'FaceColor', color, 'EdgeAlpha', 0.3);
    legend('Satellite trajectory','Satellite','Projected path of the satalite on the earth','Equatorial plane','location','north');
    plot3([0,EE2(1,1)],[0,EE2(2,1)],[0,EE2(3,1)],'b')

    subplot(6,2,[11,12])
    cla;
    axis off
    title('Satellite and Egypt related data','FontSize',16)
    text(0.05,0.9,'Satellite Visitation time (sec)','FontSize',12)
    text(0.05,0.5,'Satellite Leave time (sec)','FontSize',12)
    for T=1:length(leave)
        text(0.2+0.15*T,0.9,num2str(visit(T)),'FontSize',10)
        text(0.2+0.15*T,0.5,num2str(leave(T)),'FontSize',10)
    end
    pause(1/hz);
end
%% effect of right ascension of the ascending node (Omega)
figure(3);
[rz,Z]=max(R2_z-R1_z);
[ h, mag_h, i, omega, e, mag_e, w, theta ] = OrbitalElements ( [R2_x(Z)-R1_x(Z),R2_y(Z)-R1_y(Z),R2_z(Z)-R1_z(Z)],[V2_x(Z),V2_y(Z),V2_z(Z)],398600 );
set(gcf,'color','w');
% Texturemap the globe
% Load Earth image for texture map
cdata = imread(image_file);
J = imrotate(cdata,0,'bilinear');
v=1; V=1;
%--------------------------------------------------------------------------------------------------------------------------------------------------------
for p=0:10:360
    [ Q ] = RM ( -(p-omega), i*0, 0*w, '313');
    EE2=Q*[R2_x;R2_y;R2_z];
    EE1=Q*[R1_x;R1_y;R1_z];
    EE0=Q*[XX1;YY1;ZZ1];
    set(0,'defaultfigureposition',[40 55 1300 600])
    for k=1:length(R1_x)
        [ RA(k), Dec(k) ] = lmn2RaDec( [EE2(1,k), EE2(2,k), EE2(3,k)]-[EE1(1,k), EE1(2,k), EE1(3,k)] );
        RA(k)=RA(k)+erot*(k-1)*dt*180/pi;
        while RA(k) > 360
            RA(k)=RA(k) - 360;
        end
    end
    % Satellite Tracking
    subplot(6,2,[2,4,6,8])
    cla;
    visit=[]; leave=[];
    v=1;
    image([0 360],[-90 90],cdata);
    title('Satellite Tracking','Fontsize',18);
    hold all;
    grid on;
    set(gca,'gridlinestyle','-')
    plot(RA, -Dec,'.','color','g','LineWidth',1)
    plot(RA(1), -Dec(1),'o','color','r','LineWidth',4)
    plot([204.2 217.3 217.3 204.2 204.2],[-21.3 -21.3 -32.1 -32.1 -21.3],'color','k','LineWidth',2)
    set(gca,'XTickLabel',{'-180','-129','-77','-26','26','77','129','180'})
    set(gca,'YTickLabel',{'100','80','60','40','20','0','-20','-40','-60','-80','-100'})
    for P=1:length(R2_x)
        if RA(P) >= 204.2 && RA(P) <= 217.3 && Dec(P) <= 32.1 && Dec(P) >= 21.3 && in == 0
            counter(v+1) = counter(v) +1;
            in=1;
            visit(v)=(P-1)*dt;
        elseif (RA(P) < 204.2 || RA(P) > 217.3) && (Dec(P) > 32.1 || Dec(P) < 21.3) && in == 1
            in=0;
            leave(v)=(P-1)*dt;
            v=v+1;
        end
    end
    xlim([0 360]);
    ylim([-100 100]);


    subplot(6,2,[1,3,5,7])
    cla;
    % Options
    space_color = 'k';
    npanels = 180;   % Number of globe panels around the equator deg/panel = 360/npanels
    alpha   = 1; % globe transparency level, 1 = opaque, through 0 = invisible
    % Earth texture image
    % Anything imread() will handle, but needs to be a 2:1 unprojected globe
    % Mean spherical earth
    erad    = 6371008.7714; % equatorial radius (meters)
    prad    = 6371008.7714; % polar radius (meters)
    %GMST0 = []; % Don't set up rotatable globe (ECEF)
    GMST0 = (4.89496121282306 - (p-1)*erot*dt*0); % Set up a rotatable globe at J2000.0
    set(gcf,'Color','w');
    % M2 trajectory
    plot3(I(1)*ones(1,length(R1_x))+EE2(1,:)-EE1(1,:),I(2)*ones(1,length(R1_y))+EE2(2,:)-EE1(2,:),I(3)*ones(1,length(R1_z))+EE2(3,:)-EE1(3,:),'g','LineWidth',2);
    hold all;
    grid on;
    xlabel('X','Fontsize',18);
    ylabel('Y','Fontsize',18);
    zlabel('Z','Fontsize',18);
    title(['Satellite about Earth with right ascension of the ascending node = ' num2str(p) '^o'],'Fontsize',18);
    xlim auto;
    ylim auto;
    zlim auto;
    view(3);
    % M2 start
    plot3(EE2(1,1)-EE1(1,1),EE2(2,1)-EE1(2,1),EE2(3,1)-EE1(3,1),'o','color',[.9,0.2,0.5],'LineWidth',10);
    % Create wireframe globe
    % Create a 3D meshgrid of the sphere points using the ellipsoid function
    [x, y, z] = ellipsoid(0, 0, 0, erad, erad, prad, npanels);
    globe = surf(y, -x, -z, 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]);
    if ~isempty(GMST0)
        hgx = hgtransform;
        set(hgx,'Matrix', makehgtform('zrotate',GMST0));
        set(globe,'Parent',hgx);
    end
    % Set image as color data (cdata) property, and set face color to indicate
    % a texturemap, which Matlab expects to be in cdata. Turn off the mesh edges.
    set(globe, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', alpha, 'EdgeColor', 'none');
    %Projected path of the satalite on the earth
    plot3(EE0(1,:),EE0(2,:),EE0(3,:),'r','LineWidth',2);
    % Equatorial plane
    if (R2_z(1)-R1_z(1)) >=0
        color='yellow';
    else
        color=[0.5,0.5,0.5];
    end
    surf(x_mat, y_mat, z_mat, 'FaceAlpha', 0.3, 'EdgeColor', color, 'FaceColor', color, 'EdgeAlpha', 0.3);
    legend('Satellite trajectory','Satellite','Projected path of the satalite on the earth','Equatorial plane','location','north');
    plot3([0,EE2(1,1)],[0,EE2(2,1)],[0,EE2(3,1)],'b')

    subplot(6,2,[11,12])
    cla;
    axis off
    title('Satellite and Egypt related data','FontSize',16)
    text(0.05,0.9,'Satellite Visitation time (sec)','FontSize',12)
    text(0.05,0.5,'Satellite Leave time (sec)','FontSize',12)
    for T=1:length(leave)
        text(0.2+0.15*T,0.9,num2str(visit(T)),'FontSize',10)
        text(0.2+0.15*T,0.5,num2str(leave(T)),'FontSize',10)
    end
    pause(1/hz)
end
       %% getting photos for gif
%--------------------------------------------------------------------------------------------------------------------------------------------------------       
% satelite tracking
%--------------------------------------------------------------------------------------------------------------------------------------------------------       
       % % Texturemap the globe
       % % Load Earth image for texture map
       % cdata = imread(image_file);
       % J = imrotate(cdata,0,'bilinear');
       % v=1;
       % %--------------------------------------------------------------------------------------------------------------------------------------------------------
       % for p=1:step_sim:length(R1_x)
       %     figure(v);
       %     set(0,'defaultfigureposition',[40 55 1300 600])
       %     subplot(1,2,1)
       %     cla;
       %     % Options
       %     space_color = 'k';
       %     npanels = 180;   % Number of globe panels around the equator deg/panel = 360/npanels
       %     alpha   = 1; % globe transparency level, 1 = opaque, through 0 = invisible
       %     % Earth texture image
       %     % Anything imread() will handle, but needs to be a 2:1 unprojected globe
       %     % Mean spherical earth
       %     erad    = 6371008.7714; % equatorial radius (meters)
       %     prad    = 6371008.7714; % polar radius (meters)
       %     %GMST0 = []; % Don't set up rotatable globe (ECEF)
       %     GMST0 = (4.89496121282306 - (p-1)*erot*dt); % Set up a rotatable globe at J2000.0
       %     set(gcf,'Color','w');
       %     % M2 trajectory
       %     plot3(I(1)*ones(1,length(R1_x))+R2_x-R1_x,I(2)*ones(1,length(R1_y))+R2_y-R1_y,I(3)*ones(1,length(R1_z))+R2_z-R1_z,'g','LineWidth',2);
       %     hold all;
       %     grid on;
       %     xlabel('X','Fontsize',18);
       %     ylabel('Y','Fontsize',18);
       %     zlabel('Z','Fontsize',18);
       %     title('Satellite about Earth','Fontsize',18);
       %     xlim auto;
       %     ylim auto;
       %     zlim auto;
       %     view(3);
       %     % M2 start
       %     plot3(R2_x(1,p)-R1_x(1,p),R2_y(1,p)-R1_y(1,p),R2_z(1,p)-R1_z(1,p),'o','color',[.9,0.2,0.5],'LineWidth',10);
       %     % Create wireframe globe
       %     % Create a 3D meshgrid of the sphere points using the ellipsoid function
       %     [x, y, z] = ellipsoid(0, 0, 0, erad, erad, prad, npanels);
       %     globe = surf(y+R1_x(1,p), -x+R1_y(1,p), -z+R1_z(1,p), 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]);
       %     if ~isempty(GMST0)
       %         hgx = hgtransform;
       %         set(hgx,'Matrix', makehgtform('zrotate',GMST0));
       %         set(globe,'Parent',hgx);
       %     end
       %     % Set image as color data (cdata) property, and set face color to indicate
       %     % a texturemap, which Matlab expects to be in cdata. Turn off the mesh edges.
       %     set(globe, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', alpha, 'EdgeColor', 'none');
       %     %Projected path of the satalite on the earth
       %     plot3(XX1,YY1,ZZ1,'r','LineWidth',2);
       %     % Equatorial plane
       %     if (R2_z(1,p)-R1_z(1,p)) >=0
       %         color='yellow';
       %     else
       %         color=[0.5,0.5,0.5];
       %     end
       %     surf(x_mat, y_mat, z_mat, 'FaceAlpha', 0.3, 'EdgeColor', color, 'FaceColor', color, 'EdgeAlpha', 0.3);
       %     legend('Satellite trajectory','Satellite','Projected path of the satalite on the earth','Equatorial plane','location','north');
       %     plot3([0,R2_x(p)],[0,R2_y(p)],[0,R2_z(p)],'b')
       %
       %
       %
       %     % Satellite Tracking
       %     subplot(1,2,2)
       %     cla;
       %     image([0 360],[-90 90],cdata);
       %     title('Satellite Tracking','Fontsize',18);
       %     hold all;
       %     grid on;
       %     set(gca,'gridlinestyle','-')
       %     plot(RA(1:p), -Dec(1:p),'.','color','g','LineWidth',1)
       %     plot(RA(1), -Dec(1),'o','color','r','LineWidth',4)
       %     plot(RA(p), -Dec(p),'s','color',[.9,0.2,0.5],'LineWidth',2)
       % %     scatter(RA(1:p), -Dec(1:p),'fill')
       % %     gscatter (RA(1:p), -Dec(1:p))
       %     set(gca,'XTickLabel',{'-180','-129','-77','-26','26','77','129','180'})
       %     set(gca,'YTickLabel',{'80','60','40','20','0','-20','-40','-60','-80'})
       %     if RA(p) >= 204.2 && RA(p) <= 217.3 && Dec(p) <= 32.1 && Dec(p) >= 21.3 && in == 0
       %         counter = counter +1;
       %         in=1;
       %         visit=(p-1)*dt;
       %     elseif (RA(p) < 204.2 || RA(p) > 217.3) && (Dec(p) > 32.1 || Dec(p) < 21.3) && in == 1
       %         in=0;
       %         leave=(p-1)*dt;
       %     end
       %     xlim([0 360]);
       %     ylim([-90 90]);
       %     plot([204.2 217.3 217.3 204.2 204.2],[-21.3 -21.3 -32.1 -32.1 -21.3],'color','k','LineWidth',2)
       %     img= getframe(gcf);
       %     imwrite(img.cdata, [num2str(v), '.png']);
       %     v=v+1;
       %     close;
       % end
%--------------------------------------------------------------------------------------------------------------------------------------------------------       
% effect of inclination angle (i)
%--------------------------------------------------------------------------------------------------------------------------------------------------------       
       % [rz,Z]=max(R2_z-R1_z);
       % [ h, mag_h, i, omega, e, mag_e, w, theta ] = OrbitalElements ( [R2_x(Z)-R1_x(Z),R2_y(Z)-R1_y(Z),R2_z(Z)-R1_z(Z)],[V2_x(Z),V2_y(Z),V2_z(Z)],398600 );
       % set(gcf,'color','w');
       % % Texturemap the globe
       % % Load Earth image for texture map
       % cdata = imread(image_file);
       % J = imrotate(cdata,0,'bilinear');
       % v=1; V=1;
       % %--------------------------------------------------------------------------------------------------------------------------------------------------------
       % for p=0:5:180
       %     figure(V);
       % % for p=0:1
       %     [ Q ] = RM ( omega*0, -(p-i), w*0, '313');
       %     EE2=Q*[R2_x;R2_y;R2_z];
       %     EE1=Q*[R1_x;R1_y;R1_z];
       %     EE0=Q*[XX1;YY1;ZZ1];
       %     set(0,'defaultfigureposition',[40 55 1300 600])
       %     for k=1:length(R1_x)
       %     [ RA(k), Dec(k) ] = lmn2RaDec( [EE2(1,k), EE2(2,k), EE2(3,k)]-[EE1(1,k), EE1(2,k), EE1(3,k)] );
       %     RA(k)=RA(k)+erot*(k-1)*dt*180/pi;
       %     while RA(k) > 360
       %         RA(k)=RA(k) - 360;
       %     end
       %     end
       %     % Satellite Tracking
       %     subplot(6,2,[2,4,6,8])
       %     cla;
       %     visit=[]; leave=[];
       %     v=1;
       %     image([0 360],[-90 90],cdata);
       %     title('Satellite Tracking','Fontsize',18);
       %     hold all;
       %     grid on;
       %     set(gca,'gridlinestyle','-')
       %     plot(RA, -Dec,'.','color','g','LineWidth',1)
       %     plot(RA(1), -Dec(1),'o','color','r','LineWidth',4)
       %     plot([204.2 217.3 217.3 204.2 204.2],[-21.3 -21.3 -32.1 -32.1 -21.3],'color','k','LineWidth',2)
       %     set(gca,'XTickLabel',{'-180','-129','-77','-26','26','77','129','180'})
       %     set(gca,'YTickLabel',{'100','80','60','40','20','0','-20','-40','-60','-80','-100'})
       %     for P=1:length(R2_x)
       %         if RA(P) >= 204.2 && RA(P) <= 217.3 && Dec(P) <= 32.1 && Dec(P) >= 21.3 && in == 0
       %             counter(v+1) = counter(v) +1;
       %             in=1;
       %             visit(v)=(P-1)*dt;
       %         elseif (RA(P) < 204.2 || RA(P) > 217.3) && (Dec(P) > 32.1 || Dec(P) < 21.3) && in == 1
       %             in=0;
       %             leave(v)=(P-1)*dt;
       %             v=v+1;
       %         end
       %     end
       %     xlim([0 360]);
       %     ylim([-100 100]);
       %
       %
       %     subplot(6,2,[1,3,5,7])
       %     cla;
       %     % Options
       %     space_color = 'k';
       %     npanels = 180;   % Number of globe panels around the equator deg/panel = 360/npanels
       %     alpha   = 1; % globe transparency level, 1 = opaque, through 0 = invisible
       %     % Earth texture image
       %     % Anything imread() will handle, but needs to be a 2:1 unprojected globe
       %     % Mean spherical earth
       %     erad    = 6371008.7714; % equatorial radius (meters)
       %     prad    = 6371008.7714; % polar radius (meters)
       %     %GMST0 = []; % Don't set up rotatable globe (ECEF)
       %     GMST0 = (4.89496121282306 - (p-1)*erot*dt); % Set up a rotatable globe at J2000.0
       %     set(gcf,'Color','w');
       %     % M2 trajectory
       %     plot3(I(1)*ones(1,length(R1_x))+EE2(1,:)-EE1(1,:),I(2)*ones(1,length(R1_y))+EE2(2,:)-EE1(2,:),I(3)*ones(1,length(R1_z))+EE2(3,:)-EE1(3,:),'g','LineWidth',2);
       %     hold all;
       %     grid on;
       %     xlabel('X','Fontsize',18);
       %     ylabel('Y','Fontsize',18);
       %     zlabel('Z','Fontsize',18);
       %     title(['Satellite about Earth with inclination angle = ' num2str(p)],'Fontsize',18);
       %     xlim auto;
       %     ylim auto;
       %     zlim auto;
       %     view(3);
       % %     view(-90,0);
       %     % M2 start
       %     plot3(EE2(1,1)-EE1(1,1),EE2(2,1)-EE1(2,1),EE2(3,1)-EE1(3,1),'o','color',[.9,0.2,0.5],'LineWidth',10);
       %     % Create wireframe globe
       %     % Create a 3D meshgrid of the sphere points using the ellipsoid function
       %     [x, y, z] = ellipsoid(0, 0, 0, erad, erad, prad, npanels);
       %     globe = surf(y, -x, -z, 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]);
       %     if ~isempty(GMST0)
       %         hgx = hgtransform;
       %         set(hgx,'Matrix', makehgtform('zrotate',GMST0));
       %         set(globe,'Parent',hgx);
       %     end
       %     % Set image as color data (cdata) property, and set face color to indicate
       %     % a texturemap, which Matlab expects to be in cdata. Turn off the mesh edges.
       %     set(globe, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', alpha, 'EdgeColor', 'none');
       %     %Projected path of the satalite on the earth
       %     plot3(EE0(1,:),EE0(2,:),EE0(3,:),'r','LineWidth',2);
       %     % Equatorial plane
       %     if (R2_z(1)-R1_z(1)) >=0
       %         color='yellow';
       %     else
       %         color=[0.5,0.5,0.5];
       %     end
       %     surf(x_mat, y_mat, z_mat, 'FaceAlpha', 0.3, 'EdgeColor', color, 'FaceColor', color, 'EdgeAlpha', 0.3);
       %     legend('Satellite trajectory','Satellite','Projected path of the satalite on the earth','Equatorial plane','location','north');
       %     plot3([0,EE2(1,1)],[0,EE2(2,1)],[0,EE2(3,1)],'b')
       %
       %     subplot(6,2,[11,12])
       %     cla;
       %     axis off
       %     title('Satellite and Egypt related data','FontSize',16)
       %     text(0.05,0.9,'Satellite Visitation time (sec)','FontSize',12)
       %     text(0.05,0.5,'Satellite Leave time (sec)','FontSize',12)
       %     for T=1:length(leave)
       %         text(0.2+0.15*T,0.9,num2str(visit(T)),'FontSize',10)
       %         text(0.2+0.15*T,0.5,num2str(leave(T)),'FontSize',10)
       %     end
       %     img= getframe(gcf);
       %     imwrite(img.cdata, [num2str(V), '.png']);
       %     V=V+1;
       %     close;
       % end
%--------------------------------------------------------------------------------------------------------------------------------------------------------       
% effect of right ascension of the ascending node (Omega)
%--------------------------------------------------------------------------------------------------------------------------------------------------------       
        % [rz,Z]=max(R2_z-R1_z);
        % [ h, mag_h, i, omega, e, mag_e, w, theta ] = OrbitalElements ( [R2_x(Z)-R1_x(Z),R2_y(Z)-R1_y(Z),R2_z(Z)-R1_z(Z)],[V2_x(Z),V2_y(Z),V2_z(Z)],398600 );
        % set(gcf,'color','w');
        % % Texturemap the globe
        % % Load Earth image for texture map
        % cdata = imread(image_file);
        % J = imrotate(cdata,0,'bilinear');
        % v=1; V=1;
        % %--------------------------------------------------------------------------------------------------------------------------------------------------------
        % for p=0:10:360
        %     figure(V);
        %     [ Q ] = RM ( -(p-omega), i*0, 0*w, '313');
        %     EE2=Q*[R2_x;R2_y;R2_z];
        %     EE1=Q*[R1_x;R1_y;R1_z];
        %     EE0=Q*[XX1;YY1;ZZ1];
        %     set(0,'defaultfigureposition',[40 55 1300 600])
        %     for k=1:length(R1_x)
        %         [ RA(k), Dec(k) ] = lmn2RaDec( [EE2(1,k), EE2(2,k), EE2(3,k)]-[EE1(1,k), EE1(2,k), EE1(3,k)] );
        %         RA(k)=RA(k)+erot*(k-1)*dt*180/pi;
        %         while RA(k) > 360
        %             RA(k)=RA(k) - 360;
        %         end
        %     end
        %     % Satellite Tracking
        %     subplot(6,2,[2,4,6,8])
        %     cla;
        %     visit=[]; leave=[];
        %     v=1;
        %     image([0 360],[-90 90],cdata);
        %     title('Satellite Tracking','Fontsize',18);
        %     hold all;
        %     grid on;
        %     set(gca,'gridlinestyle','-')
        %     plot(RA, -Dec,'.','color','g','LineWidth',1)
        %     plot(RA(1), -Dec(1),'o','color','r','LineWidth',4)
        %     plot([204.2 217.3 217.3 204.2 204.2],[-21.3 -21.3 -32.1 -32.1 -21.3],'color','k','LineWidth',2)
        %     set(gca,'XTickLabel',{'-180','-129','-77','-26','26','77','129','180'})
        %     set(gca,'YTickLabel',{'100','80','60','40','20','0','-20','-40','-60','-80','-100'})
        %     for P=1:length(R2_x)
        %         if RA(P) >= 204.2 && RA(P) <= 217.3 && Dec(P) <= 32.1 && Dec(P) >= 21.3 && in == 0
        %             counter(v+1) = counter(v) +1;
        %             in=1;
        %             visit(v)=(P-1)*dt;
        %         elseif (RA(P) < 204.2 || RA(P) > 217.3) && (Dec(P) > 32.1 || Dec(P) < 21.3) && in == 1
        %             in=0;
        %             leave(v)=(P-1)*dt;
        %             v=v+1;
        %         end
        %     end
        %     xlim([0 360]);
        %     ylim([-100 100]);
        %
        %
        %     subplot(6,2,[1,3,5,7])
        %     cla;
        %     % Options
        %     space_color = 'k';
        %     npanels = 180;   % Number of globe panels around the equator deg/panel = 360/npanels
        %     alpha   = 1; % globe transparency level, 1 = opaque, through 0 = invisible
        %     % Earth texture image
        %     % Anything imread() will handle, but needs to be a 2:1 unprojected globe
        %     % Mean spherical earth
        %     erad    = 6371008.7714; % equatorial radius (meters)
        %     prad    = 6371008.7714; % polar radius (meters)
        %     %GMST0 = []; % Don't set up rotatable globe (ECEF)
        %     GMST0 = (4.89496121282306 - (p-1)*erot*dt*0); % Set up a rotatable globe at J2000.0
        %     set(gcf,'Color','w');
        %     % M2 trajectory
        %     plot3(I(1)*ones(1,length(R1_x))+EE2(1,:)-EE1(1,:),I(2)*ones(1,length(R1_y))+EE2(2,:)-EE1(2,:),I(3)*ones(1,length(R1_z))+EE2(3,:)-EE1(3,:),'g','LineWidth',2);
        %     hold all;
        %     grid on;
        %     xlabel('X','Fontsize',18);
        %     ylabel('Y','Fontsize',18);
        %     zlabel('Z','Fontsize',18);
        %     title(['Satellite about Earth with right ascension of the ascending node = ' num2str(p) '^o'],'Fontsize',18);
        %     xlim auto;
        %     ylim auto;
        %     zlim auto;
        %     view(3);
        %     % M2 start
        %     plot3(EE2(1,1)-EE1(1,1),EE2(2,1)-EE1(2,1),EE2(3,1)-EE1(3,1),'o','color',[.9,0.2,0.5],'LineWidth',10);
        %     % Create wireframe globe
        %     % Create a 3D meshgrid of the sphere points using the ellipsoid function
        %     [x, y, z] = ellipsoid(0, 0, 0, erad, erad, prad, npanels);
        %     globe = surf(y, -x, -z, 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]);
        %     if ~isempty(GMST0)
        %         hgx = hgtransform;
        %         set(hgx,'Matrix', makehgtform('zrotate',GMST0));
        %         set(globe,'Parent',hgx);
        %     end
        %     % Set image as color data (cdata) property, and set face color to indicate
        %     % a texturemap, which Matlab expects to be in cdata. Turn off the mesh edges.
        %     set(globe, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', alpha, 'EdgeColor', 'none');
        %     %Projected path of the satalite on the earth
        %     plot3(EE0(1,:),EE0(2,:),EE0(3,:),'r','LineWidth',2);
        %     % Equatorial plane
        %     if (R2_z(1)-R1_z(1)) >=0
        %         color='yellow';
        %     else
        %         color=[0.5,0.5,0.5];
        %     end
        %     surf(x_mat, y_mat, z_mat, 'FaceAlpha', 0.3, 'EdgeColor', color, 'FaceColor', color, 'EdgeAlpha', 0.3);
        %     legend('Satellite trajectory','Satellite','Projected path of the satalite on the earth','Equatorial plane','location','north');
        %     plot3([0,EE2(1,1)],[0,EE2(2,1)],[0,EE2(3,1)],'b')
        %
        %     subplot(6,2,[11,12])
        %     cla;
        %     axis off
        %     title('Satellite and Egypt related data','FontSize',16)
        %     text(0.05,0.9,'Satellite Visitation time (sec)','FontSize',12)
        %     text(0.05,0.5,'Satellite Leave time (sec)','FontSize',12)
        %     for T=1:length(leave)
        %         text(0.2+0.15*T,0.9,num2str(visit(T)),'FontSize',10)
        %         text(0.2+0.15*T,0.5,num2str(leave(T)),'FontSize',10)
        %     end
        %     img= getframe(gcf);
        %     imwrite(img.cdata, [num2str(V), '.png']);
        %     V=V+1;
        %     close;
        % end
        %
