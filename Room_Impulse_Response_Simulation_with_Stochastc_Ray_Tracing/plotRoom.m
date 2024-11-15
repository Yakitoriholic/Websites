function plotRoom(roomDimensions,receiverCoord,sourceCoord,figHandle)
    % PLOTROOM Helper function to plot 3D room with receiver/transmitter points
    figure(figHandle)
    X = [0;roomDimensions(1);roomDimensions(1);0;0];
    Y = [0;0;roomDimensions(2);roomDimensions(2);0];
    Z = [0;0;0;0;0];
    figure;
    hold on;
    plot3(X,Y,Z,"k",LineWidth=1.5);   
    plot3(X,Y,Z+roomDimensions(3),"k",LineWidth=1.5); 
    set(gca,"View",[-28,35]); 
    for k=1:length(X)-1
        plot3([X(k);X(k)],[Y(k);Y(k)],[0;roomDimensions(3)],"k",LineWidth=1.5);
    end
    grid on
    xlabel("X (m)")
    ylabel("Y (m)")
    zlabel("Z (m)")
    plot3(sourceCoord(1),sourceCoord(2),sourceCoord(3),"bx",LineWidth=2)
    plot3(receiverCoord(1),receiverCoord(2),receiverCoord(3),"ro",LineWidth=2)
    end
    
    function X=RandSampleSphere(N)
    % RANDSAMPLESPHERE Return random ray directions
    
    % Sample the unfolded right cylinder
    z = 2*rand(N,1)-1;
    lon = 2*pi*rand(N,1);
    
    % Convert z to latitude
    z(z<-1) = -1;
    z(z>1) = 1;
    lat = acos(z);
    
    % Convert spherical to rectangular co-ords
    s = sin(lat);
    x = cos(lon).*s;
    y = sin(lon).*s;
    
    X = [x y z];
    end
    
    function [surfaceofimpact,displacement] = getImpactWall(ray_xyz,ray_dxyz,roomDims)
    % GETIMPACTWALL Determine which wall the ray encounters
    surfaceofimpact = -1;
    displacement = 1000;
    %  Compute time to intersection with x-surfaces
    if (ray_dxyz(1) < 0)
        displacement = -ray_xyz(1) / ray_dxyz(1);
        if displacement==0
            displacement=1000;
        end
        surfaceofimpact = 0;
    elseif (ray_dxyz(1) > 0)
        displacement = (roomDims(1) - ray_xyz(1)) / ray_dxyz(1);
        if displacement==0
            displacement=1000;
        end
        surfaceofimpact = 1;
    end
    % Compute time to intersection with y-surfaces
    if ray_dxyz(2)<0
        t = -ray_xyz(2) / ray_dxyz(2);
        if (t<displacement) && t>0
            surfaceofimpact = 2;
            displacement = t;
        end
    elseif ray_dxyz(2)>0
        t = (roomDims(2) - ray_xyz(2)) / ray_dxyz(2);
        if (t<displacement) && t>0
            surfaceofimpact = 3;
            displacement = t;
        end
    end
    % Compute time to intersection with z-surfaces
    if ray_dxyz(3)<0
        t = -ray_xyz(3) / ray_dxyz(3);
        if (t<displacement) && t>0
            surfaceofimpact = 4;
            displacement = t;
        end
    elseif ray_dxyz(3)>0
        t = (roomDims(3) - ray_xyz(3)) / ray_dxyz(3);
        if (t<displacement) && t>0
            surfaceofimpact = 5;
            displacement = t;
        end
    end
    surfaceofimpact = surfaceofimpact + 1;
    
    displacement = displacement * ray_dxyz;
    
    end
    
    function N = getWallNormalVector(surfaceofimpact)
    % GETWALLNORMALVECTOR Get the normal vector of a surface
    switch surfaceofimpact
        case 1
            N = [1 0 0];
        case 2
            N = [-1 0 0];
        case 3
            N = [0 1 0];
        case 4
            N = [0 -1 0];
        case 5
            N = [0 0 1];
        case 6
            N = [0 0 -1];
    end
    
    end
    
    function bef = getBandedgeFrequencies(cf,fs)
    G = 2;
    BandsPerOctave = 1;
    fbpo = .5/BandsPerOctave;
    bef = [ cf(1)*G^-fbpo, ...
        sqrt(cf(1:end-1).*cf(2:end)), ...
        min(cf(end)*G^fbpo, fs/2) ];
    end