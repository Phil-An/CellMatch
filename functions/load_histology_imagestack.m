function imagestack = load_histology_imagestack
% load_histology_imagestack loads z-stack iamge sequences by user 
% interaction. Use the file picker tool and select all images belonging to 
% the first z-stack and click "Open". Another file picker tool will open
% for selecting images belonging to the next z-stack. When you have loaded
% the iamges of all histology sections press the escape button or hit the
% "Cancel" button to proceed with loading of the images. All images of a
% z-stack must have the same dimensions.
%
% Output :  
%           imagestack : cell array containing images of z-stacks
%
% Function is written by Philip Anner (2020)

% Select file names
loadsection = true;
% Counter for number of histology sections loaded
AllSections = 1;
while loadsection
    M = msgbox(['Please select all images of histology section ' num2str(AllSections) ' you would like to load. If you have loaded the images of all sections select Cancel.']);
    waitfor(M);
    [file{AllSections}, path{AllSections}] = uigetfile({'*.tif; *.tiff; *.jpg; *.jpeg; *.png', 'All Files (*.*)'}, ['Select images of section ' num2str(AllSections) '. If done - click [Cancel]'],'.\Data\Post hoc histology', 'MultiSelect', 'on');
    if ~ iscell(file{AllSections})
        loadsection = false;
        file(AllSections) = [];
        path(AllSections) = [];
        AllSections = AllSections-1;
    else
        AllSections = AllSections+1;
    end
end

% sum of all images
nImages = sum(cellfun(@numel,file));

%Load images
imagestack = cell(1,AllSections);
h = waitbar(0,'Loading images, please wait...');
count=1;
for section = 1: AllSections
    nimgStack = numel(file{section});
    img = cell(1,nimgStack);
    for i = 1:nimgStack
        img{i,1}= (imread([path{section} file{section}{i}]));
        % convert to grayscale
        if size(img{i,1},3) > 1
            img{i,1}=rgb2gray(img{i,1});
        end
        waitbar(count/nImages);
        count = count+1;
    end
    imagestack{section} =img;
    clear images
end
close(h)

try
    for i = 1:size(imagestack,2)
        imagestack{i}=cat(3,imagestack{i}{:});
    end
catch
    error('You must select at least one histology stack for processing. All images of a histology stack must have the same dimensions!')
end
