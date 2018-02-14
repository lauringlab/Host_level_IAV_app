# Host_level_IAV_app
A shiny app for exploring the iSNV in our data set especially those in transmission pairs.

### Overview

./data -  the relevant .csv file taken from the Host_level_IAV_evolution repository

./scripts  - script used to format the data from Host_level_IAV_evolution so that minimal processing is done in the app.

app.R - The script that controls the app


### Instructions
To run simply download the reposity, open app.R in Rstudio and click "run App".

The app contains two tabs. In the first it is possible to explore the frequency and 
identity of mutations relatvie to there position in the genome. Hovering over a mutaiton will display 
meta data regarding that mutations. It is possible to zoom in on a region by highlighting the desired
area and double clicking. A single click will reset the scale.

The purpose of the second tab is to explore the frequency trajectory of mutations within transmission pairs.
Selecting a house from the House Id panel will display possible transmission pairs.  It is possible
to display the identity of mutations and zoom in by hovering as before. Clicking on a mutation will highlight that
mutations trajectory in all samples. Mutations in donor samples are plotted in red. Those in recipient samples are in blue.
Solid lines connect longitudinal sample pairs while dashed lines connect samples used in the transmission bottleneck estimation.



### Dependencies
```
shiny
ggplot2
magrittr
dplyr
wesanderson
readr
```
