function main(inputdir, outputdir, fileEXT, noc, index_nucl_marker, index_antiX, scinfostruct, mpi_cmds, output_metrics, output_names)
if isfile('logfile')
    delete('logfile')
end

diary logfile
% Progress Bar init
f = uifigure;
d = uiprogressdlg(f,'Title','Cytoskaler Running',...
        'Message','Loading Data','Cancelable','on', 'Indeterminate','on');


%% Set Global Variables
% Note that these global variables are only set once throughout entire
% file, but referenced multiple times

% Struct of subceullar info (channel name, channel type, net to be used for
% segmenting)
global subcellinfoStruct
subcellinfoStruct = scinfostruct;

% Binary values indicating if user requested:
global BSC_req              % Binary Segmentation Coordinates
global MPI_req              % MPI values
global cc_req               % Correlation Coeficient values
global FBS_req              % Figure of Binary Segmentation
global FSS_req              % Figure of Subceullular Segmentation


BSC_req = output_metrics.BSC;
MPI_req = output_metrics.MPI;
cc_req = output_metrics.cc;
FBS_req = output_metrics.FBS;
FSS_req = output_metrics.FSS;

% assign tablenames for output excel file based on user output requests
global tablenames
if MPI_req == 1 && cc_req == 1
    tablenames = cellstr([getMPItitles(mpi_cmds) getCCtitles()]);  
elseif MPI_req == 1 && cc_req == 0
    tablenames = cellstr(getMPItitles(mpi_cmds));
elseif MPI_req == 0 && cc_req == 1
    tablenames = cellstr(getCCtitles(subcellinfoStruct));
else
    tablenames = 0;
end


% Java thresholding object for auto thresholds (ported from Image-J)
global javathresh_obj
javaaddpath("csThresholding");
javathresh_obj = AutoThresholder;


% Load neural nets from file
global net_skel
global net_plsm



%loadst = load('../networks/cytoskeletonNET.mat');
loadst = load('cytoskeletonNET.mat');
net_skel = loadst.cytoskeletonNET;
%loadst = load('../networks/cytoplasmNET.mat');
loadst = load('cytoplasmNET.mat');
net_plsm = loadst.cytoplasmNET;


% Add folder to "natsortfiles" files package for ordering filenames
%addpath("natsortfiles");




%% Load image filepaths
d.Message = 'Checking file paths...';
cd(inputdir)
fileList = dir(fullfile(inputdir, strcat("*", fileEXT))); 

% Check if all sets have same number of files
if length(fileList) / noc ~= floor(length(fileList) / noc)
    errordlg("All sets do not have an equal number of files. Aborting");
    return
end

% Store full filenames in cell array input_image_filepaths
num_files = length(fileList);
input_image_filepaths = cell(1, num_files);

for i = 1:num_files
    input_image_filepaths{i} = fullfile(fileList(i).folder, fileList(i).name);
end
% Sort input_image_filepaths
input_image_filepaths = natsortfiles(input_image_filepaths);


%% Begin Analysis
d.Message = 'Analyzing Set 1, Estimating total time remaining for all sets...';


% all MPI and correlation coefficient results stored in cell array named
% allAOI_results
allAOI_results = {};

% change to output directory
cd(outputdir)

%write output text names
if ~isempty(output_names)
    writecell(output_names, "tablenames.txt");
end



% set counter
counter = 0;

% for loop running through all sets
for i = 0:noc:num_files-1
    
    % check if user requested cancel;
    if d.CancelRequested
        errordlg("Operation Cancelled");
        break
    end
        
    % begin to track time for first set to get estimated time of entire
    % analysis
    if counter == 0
        tic;
    elseif counter == 1
        
        % get the amount of time first set took
        first_set_time = toc;
        
        % total time estimation
        total_time = first_set_time * num_files/noc;
        disp(total_time)
        
        % get correct units for total time estimation
        t_units = getEstTimeUNITS(total_time);
        
        % update loading bar
        d.Indeterminate = 'off';
        %d.Value = (i/noc)/(num_files/noc);
        d.Value = 0.5;
        
        time_remaining = getTimeRemaining(t_units, total_time - first_set_time * counter);
        %time_remaining = total_time - first_set_time * counter;
        %if time_remaining < 0
        %    time_remaining = 0;
        %end
        d.Message = strcat("Analyzing ","Set ",num2str(counter+1), "   ",  num2str(time_remaining, '%4.2f'), " ", t_units, " remaining for all sets...");

    else
        % update loading bar
        %d.Value = i/num_files-1;
        d.Value = 0.5;
        time_remaining = getTimeRemaining(t_units, total_time - first_set_time * counter);
        %time_remaining = total_time - first_set_time * counter;
        %if time_remaining < 0
        %    time_remaining = 0;
        %end
        d.Message = strcat("Analyzing Images,  ",  num2str(time_remaining, '%4.2f'), " ", t_units, " remaining...");
    end
    
    
    
    
    % read in images according to nuclear and anti x indexes user
    % specified
    disp("read in images according to nuclear and anti x indexes user");
    nucl_marker_img = imread(input_image_filepaths{i+index_nucl_marker});
    antiX_img = imread(input_image_filepaths{i+index_antiX});
    
    % create an Area of Interest struct that holds all channel names,
    % images, and segmented areas
    disp("create an Area of Interest struct that holds all channel names,");
    AOI_struct = createAOIstruct(nucl_marker_img, input_image_filepaths, i);
    assignin("base", "AOI_struct", AOI_struct);
    
    % get whole cell boundaries for each single cell.
    % whole cell is the union of every segmented compartment i.e 
    % (nucleus U cytoskeleton U cytoplasm)
    disp("get whole cell boundaries for each single cell.");
    singleCell_boundaries = getSCboundaries(AOI_struct);
    assignin("base", "singleCell_boundaries", singleCell_boundaries);
    
    %{
    %%colorbegin
    colored_img = zeros(1024,1024,3);

    redc = colored_img(:,:,1);
    greenc = colored_img(:,:,2);
    bluec = colored_img(:,:,3);

    for iii = 1:length(singleCell_boundaries)
    
        curr_whole_cell = singleCell_boundaries{iii};
        redc(curr_whole_cell) = rand;
        greenc(curr_whole_cell) = rand;
        bluec(curr_whole_cell) = rand;
    end

    final_img = cat(3, redc, greenc, bluec);
    imshow(final_img)
    %%color end
    %}
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cd(outputdir)
    
    % create a directory for current set
    disp("create a directory for current set");
    counter = counter + 1;
    mkdir(strcat("Set_", num2str(counter)));
    cd(strcat("Set_", num2str(counter)));
    
    
    % get single cell metrics that user specified
    disp("get single cell metrics that user specified");
    AOI_single_cell_results = singleCellAnalysis(singleCell_boundaries, mpi_cmds, antiX_img, AOI_struct);
    
    
    % convert cell array to table and write table to file
    disp("convert cell array to table and write table to file");
    AOI_single_cell_results_table = cell2table(AOI_single_cell_results, 'VariableNames', tablenames);
    writetable(AOI_single_cell_results_table, 'stats.csv');
    
    % combine previous set array results with new results
    disp("combine previous set array results with new results");
    allAOI_results = [allAOI_results; AOI_single_cell_results];
    
    
    cd ..
    
       
    
end


% convert cell array of ALL SETS results to table and write table to file
disp("convert cell array of ALL SETS results to table and write table to file");
allAOI_results_table = cell2table(allAOI_results, 'VariableNames', tablenames);
writetable(allAOI_results_table, 'stats.csv');

d.Message = "Done.";
close(f);

diary off
end

















%%
% Helper functions


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate Area of Interest Struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  AOI_struct = createAOIstruct(nucl_marker_img, input_image_filepaths, set_index)

global subcellinfoStruct
global net_skel
global net_plsm

AOI_struct(1+length(subcellinfoStruct)).origImage = [];
AOI_struct(1+length(subcellinfoStruct)).threshImg = [];
AOI_struct(1+length(subcellinfoStruct)).segmentedImgs = [];

AOI_struct(1).origImage = nucl_marker_img;
AOI_struct(1).threshImg = nucleus_thresh(AOI_struct(1).origImage);
AOI_struct(1).segmentedImgs = segemntNucleus(AOI_struct(1).threshImg);

for i = 2:length(subcellinfoStruct) + 1
    AOI_struct(i).origImage = imread(input_image_filepaths{set_index + subcellinfoStruct(i-1).index});
    
    if subcellinfoStruct(i-1).type == "skel"
        AOI_struct(i).threshImg = cytoskeleton_thresh(AOI_struct(i).origImage);
        AOI_struct(i).clroverlayed = genColorOvrly(AOI_struct(1).threshImg, AOI_struct(i).threshImg, AOI_struct(i).origImage);
        AOI_struct(i).segmentedImgs = segmentImages(AOI_struct(i).clroverlayed, net_skel, @post_processing_skel, "cytoskeleton");
    else
        AOI_struct(i).threshImg = cytoplasm_thresh(AOI_struct(i).origImage);
        AOI_struct(i).clroverlayed = genColorOvrly(AOI_struct(1).threshImg, AOI_struct(i).threshImg, AOI_struct(i).origImage);
        AOI_struct(i).segmentedImgs = segmentImages(AOI_struct(i).clroverlayed, net_plsm, @post_processing_plsm, "cytoplasm");   
    end    
end

end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Whole Cell Union Generator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function singleCell_boundaries = getSCboundaries(AOI_struct)
num_seperatedcells = length(AOI_struct(1).segmentedImgs);

singleCell_boundaries = cell(num_seperatedcells, length(AOI_struct)+1);

for i = 1:num_seperatedcells
    
    single_whole_cell = AOI_struct(1).segmentedImgs{i};
    
    for j = 2:length(AOI_struct)
        single_whole_cell = single_whole_cell | AOI_struct(j).segmentedImgs{i};
    end
    singleCell_boundaries{i,1} = imfill(single_whole_cell, 'holes');
    
end


for i = 1:num_seperatedcells
    for j = 1:length(AOI_struct)
        singleCell_boundaries{i,j+1} = singleCell_boundaries{i,1} & AOI_struct(j).threshImg;
    end
    
end


end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Single Cell Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function AOI_single_cell_results = singleCellAnalysis(singleCell_boundaries, mpi_cmds, antiX_orig,  AOI_struct)

global tablenames
global BSC_req              
global MPI_req              
global cc_req               
global FBS_req              


AOI_single_cell_results = cell(length(singleCell_boundaries), length(tablenames));
for i = 1:length(singleCell_boundaries)
    singleCell_boundaries_ROW = singleCell_boundaries(i, :);
    
    
    mkdir(strcat("Cell_",num2str(i)));
    cd(strcat("Cell_",num2str(i)));
    
    if BSC_req == 1
        [ycoords, xcoords] = find(singleCell_boundaries_ROW{1});
        writematrix(ycoords, 'ycoords.txt');
        writematrix(xcoords, 'xcoords.txt');
    end
    
    if MPI_req == 1 && cc_req == 1
        MPI_str_results = MPI_excecute(mpi_cmds, singleCell_boundaries_ROW, antiX_orig);
        correlationcoeff_str_results = corrcff_execute(singleCell_boundaries_ROW{1}, antiX_orig, AOI_struct);  
        single_cell_results = [MPI_str_results correlationcoeff_str_results]; 
        single_cell_results_table = cell2table(single_cell_results, 'VariableNames', tablenames);
        writetable(single_cell_results_table, 'stats.csv');
        AOI_single_cell_results(i,:) =  single_cell_results;
    elseif MPI_req == 1 && cc_req == 0
        MPI_str_results = MPI_excecute(mpi_cmds, singleCell_boundaries_ROW, antiX_orig);
        single_cell_results = MPI_str_results;
        single_cell_results_table = cell2table(single_cell_results, 'VariableNames', tablenames);
        writetable(single_cell_results_table, 'stats.csv');
        AOI_single_cell_results(i,:) =  single_cell_results;
    elseif MPI_req == 0 && cc_req == 1
       correlationcoeff_str_results = corrcff_execute(singleCell_boundaries_ROW{1}, antiX_orig, AOI_struct);  
       single_cell_results =  correlationcoeff_str_results;
       single_cell_results_table = cell2table(single_cell_results, 'VariableNames', tablenames);
        writetable(single_cell_results_table, 'stats.csv');
        AOI_single_cell_results(i,:) =  single_cell_results;
    end
    
   
    
    if FBS_req == 1
        sf = figure('visible','off');
        imshow(singleCell_boundaries_ROW{1})
        saveas(sf, 'WholeCell.png');
        close(sf);
    end
    
    cd ..
end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MPI Calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MPI_str_results = MPI_excecute(mpi_cmds, singleCell_boundaries_ROW, antiX_orig)

global subcellinfoStruct

MPI_cmds_size = size(mpi_cmds);
MPI_str_results = cell(1,MPI_cmds_size(1));

for i = 1:length(mpi_cmds)
    
    single_cmd_row = mpi_cmds{i};
    
    if length(single_cmd_row) == 1
        %index1 = find(strcmp([subcellinfoStruct.name], single_cmd_row{1})) + 2;
        index1 = find(strcmp({subcellinfoStruct(:).name}, single_cmd_row{1})) + 2;
        
        if isempty(index1)
            index1 = 2;
        end
        custom_boundary = singleCell_boundaries_ROW{index1};
    
    else
        for j = 1:length(single_cmd_row)
        
            if j == 1
                %index1 = find(strcmp([subcellinfoStruct.name], single_cmd_row{1})) + 2;
                index1 = find(strcmp({subcellinfoStruct(:).name}, single_cmd_row{1})) + 2;
                if isempty(index1)
                    index1 = 2;
                end
                %index2 = find(strcmp([subcellinfoStruct.name], single_cmd_row{3})) + 2;
                index2 = find(strcmp({subcellinfoStruct(:).name}, single_cmd_row{3})) + 2;
                if isempty(index2)
                    index2 = 2;
                end
                if single_cmd_row{2} == 'U' 
                    custom_boundary = singleCell_boundaries_ROW{index1} | singleCell_boundaries_ROW{index2};
                else
                    custom_boundary = singleCell_boundaries_ROW{index1};
                    custom_boundary(singleCell_boundaries_ROW{index2}) = 0;
                end
            elseif j ~= length(single_cmd_row) && single_cmd_row(j) ~= "U" && single_cmd_row(j) ~= "-"
                    %index1 = find(strcmp([subcellinfoStruct.name], single_cmd_row{j+2})) + 1;
                    index1 = find(strcmp({subcellinfoStruct(:).name}, single_cmd_row{j+2})) + 1;
                    if isempty(index1)
                        index1 = 2;
                    end
                    if single_cmd_row(j + 1) == 'U'
                        custom_boundary = custom_boundary | singleCell_boundaries_ROW{index1};
                    else
                        custom_boundary(singleCell_boundaries_ROW{index1}) = 0;
                    end
            end
        end
    end
    
   MPI_str_results{i} = num2str(mean(antiX_orig(custom_boundary)));
   
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Correlation Coefficient Calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function correlationcoeff_str_results = corrcff_execute(singleCell, antiX_orig, AOI_struct)


correlationcoeff_str_results = cell(1, length(AOI_struct));

wholeCellvector = singleCell(:);
indices_of_zero = find(~wholeCellvector);

anti_orig_vector = antiX_orig(:);
anti_orig_vector(indices_of_zero) = [];

for i = 1:length(AOI_struct)
   
    orig_img = AOI_struct(i).origImage;
    orig_img_vector = orig_img(:);
    orig_img_vector(indices_of_zero) = [];
    coefresult = corrcoef(double(orig_img_vector), double(anti_orig_vector));
    correlationcoeff_str_results{i} = num2str(coefresult(1,2));
    
end


end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% String table titles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mpi_titles = getMPItitles(mpi_cmds)
mpi_titles = cell(1,length(mpi_cmds));
for i = 1:length(mpi_cmds)
    mpi_titles{i} = strcat("MPICustomBoundary_", num2str(i));
end
end

function corrcoeftitles = getCCtitles()

global subcellinfoStruct

corrcoeftitles = cell(1, length(subcellinfoStruct));
for i = 1:length(subcellinfoStruct)+1
    corrcoeftitles{i} = strcat("CorrCoeff_", num2str(i));
end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thresholding Algorithms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Nucleus Thresholding:
function nucleus_threshed = nucleus_thresh(orig_nucleus_img)
    
global javathresh_obj

    min_size = 1000;
    nucl_hist = imhist(orig_nucleus_img);
    nucl_thresh_number = javathresh_obj.getThreshold("Otsu", nucl_hist);
    adjusted_nucl_thresh = nucl_thresh_number * 256;
    nucl_threshed = orig_nucleus_img > adjusted_nucl_thresh;
    filled_nucleus = imfill(nucl_threshed, 'holes');
    cleaned_nucleus = bwareaopen(filled_nucleus, min_size);
    border_cleared = imclearborder(cleaned_nucleus);
    nucleus_threshed = border_cleared;

end

% Cystoskeleton Thresholding:
function cytoskeleton_threshed = cytoskeleton_thresh(orig_cytoskeleton_img)

global javathresh_obj

    cytoskel_hist = imhist(orig_cytoskeleton_img);
    cytoskel_thresh_number = javathresh_obj.getThreshold("Moments", cytoskel_hist);
    adjusted_cytoskel_thresh = cytoskel_thresh_number * 256;
    cytoskeleton_threshed = orig_cytoskeleton_img > adjusted_cytoskel_thresh;
end

% Cytoplasm Thresholding:
function cytoplasm_threshed = cytoplasm_thresh(orig_cytoplasm_img)
 
global javathresh_obj
    cytoplas_hist = imhist(orig_cytoplasm_img);
    cytoplas_thresh_number = javathresh_obj.getThreshold("Percentile", cytoplas_hist);
    adjusted_cytoplas_thresh = cytoplas_thresh_number * 256;
    threshed_cytoplasm = orig_cytoplasm_img > adjusted_cytoplas_thresh;
    cytoplasm_threshed = medfilt2(threshed_cytoplasm);
 
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Binary Marking Algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function colored_images = genColorOvrly(marker_template, to_mark_threshed, to_mark_orig)


    stats = regionprops(marker_template, 'PixelList');
    stats_size_dim = size(stats);
    stats_size = stats_size_dim(1);
    

    masked_img = marker_template;
    masked_img(masked_img < 2) = 0;
    

    for pix_list_ind = 1:stats_size
 

        curr_list = stats(pix_list_ind);
        curr_list_array = struct2array(curr_list);
        size_dim_cla = size(curr_list_array);
        size_cla = size_dim_cla(1);
 

        for chng = 1:size_cla
            masked_img(curr_list_array(chng,2), curr_list_array(chng,1)) = 1;
        end 
 
    end
    

    single_nuclues_mask = masked_img;
    single_nuclues_mask(single_nuclues_mask < 2) = 0;
    

    stats = regionprops(marker_template, 'PixelList');
    stats_size_dim = size(stats);
    stats_size = stats_size_dim(1);
    


    if stats_size >= 1 
 
        images = cell(1, stats_size);
        
        for pix_list_ind = 1:stats_size
 

            single_nuclues_mask(single_nuclues_mask < 2) = 0;
 
            curr_list = stats(pix_list_ind);
            curr_list_array = struct2array(curr_list);
            size_dim_cla = size(curr_list_array);
            size_cla = size_dim_cla(1);
 

            for chng = 1:size_cla
                single_nuclues_mask(curr_list_array(chng,2), curr_list_array(chng,1)) = 1;
            end
            

            overlayed_cytoskeleton = bsxfun(@times, to_mark_orig, cast(to_mark_threshed, class(to_mark_orig)));
 

            cio_as_RGB = repmat(overlayed_cytoskeleton, [1 1 3]);
            redcio = cio_as_RGB(:,:,1);
            greencio = cio_as_RGB(:,:,2);
            bluecio = cio_as_RGB(:,:,3);
            

            redcio(marker_template) = 999999999999999999;
            greencio(marker_template) = 0;
            bluecio(marker_template) = 0;
            

            redcio(single_nuclues_mask) = 0;
            greencio(single_nuclues_mask) = 999999999999999999;
            bluecio(single_nuclues_mask) = 0;
            

            final_color_overlayed = cat(3, redcio, greencio, bluecio);
            final_color_overlayed_u8 = im2uint8(final_color_overlayed);
            images{pix_list_ind} = final_color_overlayed_u8;
            
        end

    end

    colored_images = images;
 

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nueral Network Segmentations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function segmented_BWs = segmentImages(orig_images, NeuralNet, postProcessing, classID)

segmented_BWs_init = cell(1, length(orig_images));

for i = 1:length(orig_images)
    
    currorig = orig_images{i};
    currseg = semanticseg(currorig, NeuralNet);
    post_processed_seg = postProcessing(currorig, currseg);
    finalsegmentation = post_processed_seg == classID;
    segmented_BWs_init{i} = finalsegmentation;
    
end
    segmented_BWs = segmented_BWs_init;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nucleus Segmentations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function nuclues_regions_segmented = segemntNucleus(nuclues_thresholded_img)
 
    stats = regionprops(nuclues_thresholded_img, 'PixelList');
    stats_size_dim = size(stats);
    stats_size = stats_size_dim(1);
    

    masked_img = nuclues_thresholded_img;
    masked_img(masked_img < 2) = 0;
    
    for pix_list_ind = 1:stats_size
 

        curr_list = stats(pix_list_ind);
        curr_list_array = struct2array(curr_list);
        size_dim_cla = size(curr_list_array);
        size_cla = size_dim_cla(1);

        for chng = 1:size_cla
            masked_img(curr_list_array(chng,2), curr_list_array(chng,1)) = 1;
        end 
 

    end
   

    single_nuclues_mask = masked_img;
    single_nuclues_mask(single_nuclues_mask < 2) = 0;
    

    stats = regionprops(nuclues_thresholded_img, 'PixelList');
    stats_size_dim = size(stats);
    stats_size = stats_size_dim(1);
    
    if stats_size >= 1 
        
        nuclues_regions_segmented_init = cell(1, stats_size);
        
         for pix_list_ind = 1:stats_size
             

             single_nuclues_mask(single_nuclues_mask < 2) = 0;
             

            curr_list = stats(pix_list_ind);
            curr_list_array = struct2array(curr_list);
            size_dim_cla = size(curr_list_array);
            size_cla = size_dim_cla(1);
 

            for chng = 1:size_cla
                single_nuclues_mask(curr_list_array(chng,2), curr_list_array(chng,1)) = 1;
            end
            
            nuclues_regions_segmented_init{pix_list_ind} = single_nuclues_mask;
         end
    end
    
    nuclues_regions_segmented = nuclues_regions_segmented_init;
    

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Post Processing Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Cytoskeleton Post Processing
function post_processed = post_processing_skel(orig_image, net_segged)
 
redc = orig_image(:,:,1);
greenc = orig_image(:,:,2);
bluec = orig_image(:,:,3);
 
nucleus = redc == 0 & greenc == 255 & bluec == 0;
other_nucleuses = redc == 255 & greenc == 0 & bluec == 0;
black_background = redc == 0 & greenc == 0 & bluec == 0;
 
biggest_cytoskeleton_predict = bwareafilt(net_segged == "cytoskeleton", 1);
 
extra_cytoskeleton_mspredict = net_segged == "cytoskeleton" & (~biggest_cytoskeleton_predict);
 
net_segged(extra_cytoskeleton_mspredict) = "othercytoskeleton";
 
net_segged(nucleus) = "nuclues";
net_segged(other_nucleuses) = "othernuclues";
 
net_segged(black_background) = "blackbackground";
 
 
post_processed = net_segged;
 
end

% Cytoplasm Post Processing
function post_processed = post_processing_plsm(orig_image, net_segged)
 
redc = orig_image(:,:,1);
greenc = orig_image(:,:,2);
bluec = orig_image(:,:,3);
 
nucleus = redc == 0 & greenc == 255 & bluec == 0;
other_nucleuses = redc == 255 & greenc == 0 & bluec == 0;
black_background = redc == 0 & greenc == 0 & bluec == 0;
 
biggest_cytoskeleton_predict = bwareafilt(net_segged == "cytoplasm", 1);
 
extra_cytoskeleton_mspredict = net_segged == "cytoplasm" & (~biggest_cytoskeleton_predict);
 
net_segged(extra_cytoskeleton_mspredict) = "othercytoplasm";
 
net_segged(nucleus) = "nucleus_cytoskeleton";
net_segged(other_nucleuses) = "othernucleus";
 
net_segged(black_background) = "blackbackground";
 
 
post_processed = net_segged;
 
end


% Get time units for time estimation
function t_units =getEstTimeUNITS(total_time)
t_units = "seconds";
if total_time > 60
    total_time = total_time / 60;
    t_units = "minutes";
    if total_time > 60
        t_units = "hours";
    end
end
end

function timeremaining = getTimeRemaining(t_units, total_time_rem)
if t_units == "seconds"
    timeremaining = total_time_rem;
elseif t_units == "minutes"
    timeremaining = total_time_rem/60;
else
    timeremaining = (total_time_rem/60)/60;    
end

if timeremaining < 0
    timeremaining = 0;
    
end

timeremaining = round(timeremaining, 1);


end







