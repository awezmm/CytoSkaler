# CytoSkaler 
[logo]: https://github.com/awezmm/CytoSkaler/blob/master/imagesForREADME/icon_resized.png "Logo Title Text 2"
CytoSkaler is an imaging tool that can quantify cellular distribution of antibody staining in sub-cellular areas. First, CytoSkaler uses deep-learning to segment densely populated sub-cellar regions in single image stains. Then, it analyzes pixel distribution patterns in specific regions to generate metrics that can quantify preferential binding at those specific regions.

CytoSkaler provides otherwise unknown quantitative data on relative preference of binding to different sub-cellular structures, not accounted for by assays like polyreactivity ELISAs.


<p align="center">
  <img width="100" height="100" src="https://github.com/awezmm/CytoSkaler/blob/master/imagesForREADME/icon_resized.png">
</p>


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

<img align="left"  src="https://github.com/awezmm/CytoSkaler/blob/master/imagesForREADME/boundaryexamples.png">

Now, to find antibody cellular distrubtion in a subcellular region like the Cytoplasm, we can simply overlay the mask. Recall, our orignal Antibody X stain image.

<img align="left"  src="https://github.com/awezmm/CytoSkaler/blob/master/imagesForREADME/boundaryoverlayexample.png">

Now, here are all of the boundaries we found from Step 1, being overlayed through the original Antibody X stain image.

<img align="left"  src="https://github.com/awezmm/CytoSkaler/blob/master/imagesForREADME/allboundaryflow.png">

CytoSkaler can be more precise by allowing custom boundary defintions. For instance, to be more precise the "true" cytoplasm boundary would be found by subtracting the cytoskeleton and nucleus boundaries from the cytoplasm boundary:






