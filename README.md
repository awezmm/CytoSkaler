# CytoSkaler 
[logo]: https://github.com/awezmm/CytoSkaler/blob/master/imagesForREADME/icon_resized.png "Logo Title Text 2"
CytoSkaler is an imaging tool that can quantify cellular distribution of antibody staining in sub-cellular areas. First, CytoSkaler uses deep-learning to segment densely populated sub-cellar regions in single image stains. Then, it analyzes pixel distribution patterns in specific regions to generate metrics that can quantify preferential binding at those specific regions.

CytoSkaler produces otherwise unknown quantitative data on relative preference of binding to different sub-cellular structures, not accounted for by assays like polyreactivity ELISAs.

It provides a method for cell segmentaion and an automated pipeline after segmentation that outputs specific metrics to measure antibody binding strength and specificity in any subcellular area that a user defines.

<p align="center">
  <img width="100" height="100" src="https://github.com/awezmm/CytoSkaler/blob/master/imagesForREADME/icon_resized.png">
</p>

# Requirements to Run
1. Mac OS 10.14.6 or later
2. MATLAB 2019A Runtime (https://www.mathworks.com/products/compiler/matlab-runtime.html) (install the 2019A MAC version)

A version for Windows and Linux will be released soon. The MATLAB Runtime is a free driver that allows you to run compiled MATLAB applications or components without installing MATLAB


# Type of Data to Use 
CytoSkaler uses one or multiple image sets for analysis. A single image set is a collection of stain images that focus on the same field of view (FOV). 
Here is an example of a single image set:

<img align="left"  src="https://github.com/awezmm/CytoSkaler/blob/master/imagesForREADME/ExampleSet.png">

These images all represent the same FOV. They specifically provide stains of the nucleus, cytoskeleton, and cytoplasm sub-cellular areas. The 4th image gives a stain of an antibody. In this case, CytoSkaler would quantify preferential binding of that antibody in the nucleus, cytoskeleton, and cytoplasm for each individual cell. 

However, CytoSkaler can quantify preferential binding in any type of sub-cellular region the user specifies.


# Methodology

## Part 1: Segment subcellular areas of densely populated cells using deep learning

Here is a schematic that roughly explains how CytoSkaler segments densely populated cells and their individual subcellular areas:

<p align="center">
  <img src="https://github.com/awezmm/CytoSkaler/blob/master/imagesForREADME/schematicfullpng.png">
</p>

In the first step, CytoSkaler takes in the input of raw stained images. Next, a custom thresholding algorithm is applied to each image to generate binary image that separates the subcellular pixels from the background. 

Unlike other semantic segmentaion models for single-cellular segmentation, CytoSkaler uses a system that reduces the complex problem of multi-instance class to single-instance classes by having a neural classify seperate subcellular pixels using the nuclear boundaries as markers. This system as resulted in a relatively high IOU accuracy for CytoSkaler's neural networks. 

Thus, in the third step, RGB images for each cell are generated with most nuclear markers colored red and a single nuclear marker colored green. In the fourth step, the neural net find subcellular pixels for a single cell - the cell with its nuclear marker colored green. This is repeated for each cell in the iamge. In the fifth step, a whole cell boundary is generated by creating a union of all subcellular boundaries.

Part f shows separation of individual sub-cellular areas.


## Part 2: Quantifying anitbody cellular distribution in specific subcellular regions
CytoSkaler primarily uses 2 metrics to quantify anitbody cellular distribution in specific subcellular regions - mean pixel intensity and Pearson correlation coefficeints.

CytoSkaler calculates the Pearson correlation coefficient according to the following formula:

<img  src="https://github.com/awezmm/CytoSkaler/blob/master/imagesForREADME/PCCFormula.png">

where A and B represent different vectors, containing pixel intensity values, of 2 different subcellular regions. Additionally, μ and σ represent the mean and standard deviation, respectively, of a data vector.

From Step 1, CytoSkaler generated individual boundaries of subcellular areas. Here is the continuing example for one cell:

<img align="left" src="https://github.com/awezmm/CytoSkaler/blob/master/imagesForREADME/boundaryexamples.png">

Now, to find antibody cellular distrubtion in a subcellular region like the Cytoplasm, we can simply overlay the mask. Recall, our orignal Antibody X stain image.

<img align="left"  src="https://github.com/awezmm/CytoSkaler/blob/master/imagesForREADME/boundaryoverlayexample.png">

Now, here are all of the boundaries we found from Step 1, being overlayed through the original Antibody X stain image.

<img align="left"  src="https://github.com/awezmm/CytoSkaler/blob/master/imagesForREADME/allboundaryflow.png">

CytoSkaler can be more precise by allowing custom boundary defintions. For instance, to be more precise, the "true" cytoplasm boundary would be found by subtracting the other subcellular areas(cytoskeleton and nucleus boundaries) from the cytoplasm boundary. 

CytoSkaler gives the option of subtracting all other boundaries from each boundary. More simply, this means that every pixel in a whole cell will only be assigned to one subcellular area. 

Here is an example of all the other subcellular boundaries (cytoskeleton and nucleus) are subtracted from the Cytoplasm boundary:

<img src="https://github.com/awezmm/CytoSkaler/blob/master/imagesForREADME/cytoplasmSubtraction.png">

A user can define any custom boundary using the '-' 'U' characters (to subtract or add boundaries to another boundary)

### Using Mean Pixel Intensity as a metric for preferential antibody binding strength
Let's take a look at the Mean Pixel Intesnity (MPI) values for the different "true" boundary subcellular areas we've found:

| Subcellular Area  | MPI   |
| ----------------- |:-----:|
| Nucleus           | 2114.1|
| Cytoskeleton      | 2369.7|
| Cytoplasm         | 1781.0|
| Whole Cell        | 1961.5|

A larger mean pixel intensity value indicates stronger preferential binding in that specific region. Remember the data here is just from 1 cell. With more cells we can determine the significance of these differences for this type of antibody.


### Using Pearson's Correlation as a metric for binding specificty

CytoSkaler calculates the Pearson's correlation coefficient after subcellular boundaries have been defined using a specific boundary overlayed onto it's corresponding orginal grayscale subcellular image and on the Antibody X grayscale image. Here are the correlation coefficient values for the same subcellular areas above:


| Colocalization Boundaries | Correlation Coefficient |
|---------------------------|:-----------------------:|
| Anti X vs Nucleus         |          0.0872         |
| Anti X vs Cytoskeleton    |          0.2263         |
| Anti X vs Cytoplasm       |          0.2651         |

### Note about these metrics
It is important to understand the nuanced difference between relative antibody binding strength and specificity. Moreover, these metrics are for a single cell. 
Most importantly, these images are for a single antibody. It is more valuable to compare these metrics against a different antibody to reveal important differences.

CytoSkaler conveniently ouputs these metrics in an excel file while grouping cells from accross different anitbody.

On the 2018 Apple iMac (8GB 2133MHz memory, 2.3GHz dual-core 7th-generation Intel Core i5 processor),
CytoSkaler can analyze 175 FOV images (each containing an average of 18 cells) in about 30 minutes. This means that it can roughly analyze 3,200 cells in 30 minutes or about 100 cells per minute.



# How to Install the Graphical Interface
1. Install the 2019A MAC version of MATLAB Runtime from https://www.mathworks.com/products/compiler/matlab-runtime.html 

2. Download this repository. If are you familiar with terminal you can run 
```git clone https://github.com/awezmm/CytoSkaler.git ``` 
to clone this repository. Otherwise, you can simply download a zip file from the 'clone or download' button on the main page.

3. Open Terminal and navigate to the directory of this downloaded repository.

4. Type ```chmod +x setup.sh``` and press enter.

5. Type ```./setup.sh``` and press enter.

6. There should now be a CytoSkaler application in the folder. Right click it and press open to begin.

# How to use the Graphical Interface

This video explains installation and running in greater detail:

Let's walk through an example using the demo images provided in the demo folder. 
First, note the ordering of the demo imges. Theres 2 images sets with this channel ordering:

*00.tif = nucleus stain

*01.tif = Antibody X stain

*02.tif = Cytoskeleton Stain

*03.tif = Cytoplasm Stain

CytoSkaler accepts any number of sets with 1-6 different channels. The files should all be named in the same order for each set and should all have the same extension type. All images from all sets should be in the same single folder, like the demo folder.

First, **select the number of channels per each set and the file extension for the images**. In the demo folder, there are 4 channels per set and the file extension is .tif.

Next, **select the directory that holds all image sets**. This would be the 'CytoSkaler/demo' directory for the demo images

Here's what your window should look like after these steps:
<img  src="https://github.com/awezmm/CytoSkaler/blob/master/imagesForREADME/Demo1.png">

On the right side, you will a list of the order in which CytoSkaler will process the images. Refer to this order when you assign image specification later. Press save to continue.

Now assign information for each type of image channel. First, **give each channel a name**. 
Possible names for the demo channels could simply be: "Nucleus", "AntiX", "Cytoskeleton", "Cytoplasm".

Next, **select the "type" for each channel.** 
A nuclues stain is required and should always be selected "Binary Marker (Nucleus Stain)" type. The antibody channel should always be selected as the "Probe X (to compare)" type. Subcellular areas can be selected as the "Other Subcellular Stain" type. Finally, the "skip" option type is avaiable to any extra files in a set you do not want included for analysis.

Next, **specify the segmentation strategy to use for each subcellular area** The neural networks were trained on a vast sample of images so look at the trainig set to see which one more closely relates with your data. 

Here is what you window should look like after filling out this information:
<img src="https://github.com/awezmm/CytoSkaler/blob/master/imagesForREADME/Demo2.png">

Next, **select the type of data you want to output**.

<img src="https://github.com/awezmm/CytoSkaler/blob/master/imagesForREADME/Demo3.png">
Binary segmentation Coordinates will provide the X and Y coordinates of each segmented subcellular area.
MPI provides mean pixel intensity values for a custom region overlayed onto the Antibody probe channel. 
Correlation Coefficients are simply the correlation values explained above between every whole subcellular area and the antibody probe channel.
The binary segmentation figure will provide images of each segmented subcellular area.

**Custom regions can be defined for calculating MPI by writing specifications, as separate lines, into the text box on the right, that appears when you click the MPI option**. 

Use the channel names you defined earlier to create custom boundaries and the 'U' and '-' characters with a space between every boundary and operator.

The 'U' character represents a union between two regions and the '-' represents a subtraction between two regions.

* If you wanted the MPI of the cytoskeleton region overlayed on the Antibody probe, you would simply write: Cytoskeleton 

* If you wanted the MPI of the cytoskeleton region (with the nucelus region subtracted), overlayed on the Antibody probe, you would simply write: Cytoskeleton - Nucleus

* If you wanted the MPI of the cytoskeleton region (with the cytoplasm region added), overlayed on the Antibody probe, you would simply write Cytoskeleton U Cytoplasm

Note that this can also apply to more complex defintions, with more than 2 regions, like: 

Nuclues U Cytoskeleton U Cytoplasm - Nucleus

This would indicate a whole cell region with a subtraction of the nucleus region (a whole cell region is defined by a union of all subcellular areas)

First this demo, let do two custom boundaries: Nucleus - Cytoskeleton and Cytoskeleton - Nuclues.
This is what the text box shoud look like for those custom boundaries:

<img src="https://github.com/awezmm/CytoSkaler/blob/master/imagesForREADME/Demo4.png">

Finally, **select an output directory** to write all output data files to. **Click run to begin**.


The output data will have metrics for each cell in csv files in separate folders. Each set will also have a csv file with metrics for all cells in the set and a final csv file will contain metrics for all cells in all sets. The tablenames for the stats csv files can found in the tablenames.txt file. Images and X/Y coordinate files will also be found in each subfolder.
