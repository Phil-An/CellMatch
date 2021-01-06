function Multireg = rotation_optimization(bloodvessels, Multireg, varargin)
% This function optimizes the overlap of blood vessels by rotating blood
% vessels detected in the 3D histology model relative to in vivo detected
% blood vessel structures.
%
% Inputs:
%          bloodvessels : a 3D matrix representing blood vessels detected
%                         in histology. Call detectvessels3D(Histostack, 
%                         maxvesselsize) for detecting blood vessels.
%
%          Multireg : Data structure for pre-processing and registration
%                      control pints for multimodal image alignment between
%                      histology and in vivo acquired Ca2+ imaging data.
%                      To obatin Multireg, call findworkingdistance(Histology,RECDSA)
% Outputs:
%          Multireg : Data structure for pre-processing and registration
%                      control pints for multimodal image alignment between
%                      histology and in vivo acquired Ca2+ imaging data.
%                      The Datastructure includes rotation parameters that
%                      optimize blood vessel overlap between in vivo and
%                      post hoc acquired images.
%
% Function is written by Philip Anner (2020)

% check for blood vessel image 
if isfield(Multireg,'cavesselsAffinereg')
    ca_vessels = Multireg.cavesselsAffinereg;
    if ~islogical(ca_vessels)
        ca_vessels =imbinarize(ca_vessels);
    end
else
    error('Initially aligned blood vessel structures detected in Ca2+ imaging not found. Please call findworkingdistance(Histology,RECDSA) first!')
end

% check if workingdistance was determined
if isfield(Multireg,'Histo')
    if isfield(Multireg.Histo,'vesselimgdepth')
        rotatedstack_viewpoint = im2double(bloodvessels(:,:,1:Multireg.Histo.vesselimgdepth));
    end
else
    rotatedstack_viewpoint = im2double(bloodvessels(:,:,:));
end

% Mask blood vessels in the center of the ROI
rotatedstack_viewpoint=imresize(rotatedstack_viewpoint,0.5);
ca_vessels=imresize(ca_vessels,0.5);
if size(varargin,1) == 1
    maskroi = varargin{1};
else
    figure, imshowpair(max(rotatedstack_viewpoint,[],3), ca_vessels);
    title('Please draw ROI of endoscope. Double click in ROI to finish drawing.')
    hold on
    h = imellipse;
    pos = wait(h);
    pos=h.getPosition;
    maskroi=h.createMask;
    rectangle('Position',pos,'Curvature',[1 1],'EdgeColor','r', 'LineWidth',3);
    pause(2)
    close
end


if ~isa(ca_vessels,'logical')
    ca_vessels_bw = imbinarize(ca_vessels);
else
    ca_vessels_bw =ca_vessels;
end

rotx = (-25 :0.5 :25);
roty = (-25 :0.5: 25);

ind = 1;
totcalc=length(rotx)*length(roty);
ca_vessels_orig = ca_vessels;
res = cell(length(rotx),length(roty));

h = waitbar(0,'Calculating rotations, please wait...');
for i=1:length(rotx)
    for j= 1:length(roty)
        ca_vessels = ca_vessels_orig;
        t = makehgtform('xrotate',deg2rad(rotx(i)),'yrotate',deg2rad(roty(j)));
        tform = affine3d(t);
        rotatedstack= imwarp(rotatedstack_viewpoint,tform, 'Interp', 'nearest', 'FillValues', 255);
        rotatedstack(rotatedstack>1)=0;
        
        t=max(rotatedstack(:,:,:),[],3);
        t=mat2gray(t);
        t_bw = imbinarize(t,'adaptive');
        t_bw =bwareaopen(t_bw,20);
        t = imadjust(t,stretchlim(t),[],2);
        
        if ~isa(ca_vessels,'logical')
            ca_vessels_bw = imbinarize(ca_vessels);
        else
            ca_vessels_bw = (ca_vessels);
        end
        
        if size(t_bw,1) < size(ca_vessels,1) | size(t_bw,2) < size(ca_vessels,2)
            [ ~, t_bw] = matchimagesize(ca_vessels, t_bw);
        end
        
        if size(t,1) < size(ca_vessels,1) | size(t,2) < size(ca_vessels,2)
            [ ~, t] = matchimagesize(ca_vessels, t);
        end
        
        
        [ca_vessels_bw_fit, ~ ] = matchimagesize(ca_vessels_bw, t_bw);
        [ca_vessels_fit, ~ ] = matchimagesize(ca_vessels, t);
        maskroit = maskroi;
        [maskroi_match, ~ ] = matchimagesize(maskroit, t_bw);
        
        if ~islogical(ca_vessels_bw_fit)
            ca_vessels_bw_fit= imbinarize(ca_vessels_bw_fit);
        end
        if ~islogical(t_bw )
            t_bw = imbinarize(t_bw);
        end
        t_bw = t_bw .*maskroi_match;
        t_bw = imbinarize(t_bw);
        
        ca_vessels_bw_fit = ca_vessels_bw_fit.*maskroi_match;
        ca_vessels_bw_fit =imbinarize(ca_vessels_bw_fit );
        
        
        res{ind}.corr= corr2(t,ca_vessels_fit);
        res{ind}.dice=dice(t_bw,ca_vessels_bw_fit);
        res{ind}.jaccard=jaccard(t_bw,ca_vessels_bw_fit);
        
        
        res{ind}.segimg = t;
        res{ind}.segimg_bw = t_bw;
        res{ind}.caimg = ca_vessels_bw_fit;
        res{ind}.rotx=rotx(i);
        res{ind}.roty=roty(j);
        
        waitbar(ind/totcalc);
        
        if ind == 1
            res{ind}.xrange = rotx;
            res{ind}.yrange = roty;
            res{ind}.maskroi = maskroi;
        end
        ind = ind+1;
    end
end
close(h)


simscore= cellfun(@(x) (x.jaccard), res);
simscore = simscore(:);
all_res_rotx= cellfun(@(x) (x.rotx), res);
all_res_roty= cellfun(@(x) (x.roty), res);

[x,y]= meshgrid(rotx,roty);
z=reshape(simscore,[size(rotx,2),size(roty,2)]);
figure, s=surf(x,y,z,'FaceColor' ,'Interp','FaceLighting','gouraud')
colormap parula

[~, I ] = max(simscore);

xlabel('Rotation (degrees) in x-axis','FontSize',12)
ylabel('Rotation (degrees) in y-axis','FontSize',12)
zlabel('Jaccard index','FontSize',12)
view(3)

if res{I}.rotx <= 10 && res{I}.roty <= 10
    Multireg.viewpoint.values.rotx = res{I}.rotx;
    Multireg.viewpoint.values.roty = res{I}.roty;
    title(['Accepting rotation optimization with jaccard index : ' num2str(res{I}.jaccard) ' at rotation x = ' num2str(res{I}.rotx)  ', y = ' num2str(res{I}.roty)  ])
else
    Multireg.viewpoint.values.rotx = 0;
    Multireg.viewpoint.values.roty = 0;
    title(['Refusing rotation optimization with jaccard index : ' num2str(res{I}.jaccard) ' at rotation x = ' num2str(res{I}.rotx)  ', y = ' num2str(res{I}.roty)  ])
end

end

