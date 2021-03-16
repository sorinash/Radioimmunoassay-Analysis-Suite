# Radioimmunoassay-Analysis-Suite

A video explaining the tool's purpose, inputs, outputs, and instructions for use can be found here: 

https://www.youtube.com/watch?v=XumJkcOr0Tk&feature=youtu.be

This program is intended to be an all-in-one preparation and analysis tool for radioimmunoassay (RIA) data, a technique that is used to quantify levels of serum corticosterone, a stress hormone. For details on RIA, please see the link provided at the end of this document. 

This program is able to calculate volumes of reagent needed, (ie, tritiated corticosterone and anti-corticosterone antibody), and generate a template file to copy-paste raw data into. 

When provided with an input of counts per minute (CPM) provided by a scintillation counter, as well as a reference file, the program will format triplicate data into a more readable format, generate a standard curve for values, interpolate mass of coticosterone in picograms based on said standard curve, and convert those masses to concentration, measured in micrograms per deciliter (ug/dL). It can also group, analyze, and plot data based on user-provided parameters. 

Users input commands via the command line. A number of said commands do not require an input file. For instance: 
* help provides some instructions on how to use the program. 
* --dir=<directory name> sets a default directory via Shelve
* -v calculates reagent volume needed, based on user input, provided by --num=<number of samples>
* -m <filename> creates a template file to enter data into. 

Full analysis of RIA data requires two separate files, namely an input and a reference. 

* The input file contains two columns, namely the labels for each sample, and the CPM measured in said sample. 
* The input file must contain samples for Total, Non-Specific Binding, and B0 (which has antibody, radiation, and charcoal, but no sample provided). It also must contain 9 triplicate measurements of pre-measured corticosterone standards, which are used to generate a standard curve. The user should also provide 15 quality control samples to measure experimenter error.
* The reference file contains four columns, namely the stress condition, sex, genotype, and treatment for each sample. 

An input file resembles the following:

![image](https://user-images.githubusercontent.com/50304901/111375626-799f3900-866c-11eb-95ca-f03599774dba.png)

A reference file resembles the following:

![image](https://user-images.githubusercontent.com/50304901/111375645-8328a100-866c-11eb-953d-348ffe77c270.png)

When performing analysis, the following commands are usable:
* -n <filename> designates an input file
* -r <filename> designates a reference file
* -s reshapes an input file into a more readable format, which displays triplicate values on the same row in Excel
* -i interpolates mass of corticosterone in picograms, based on a standard curve calculated from the user's standard samples. 
* -t calculates concentration of corticosterone
* -p creates plots and analysis from data
* --output <filename> names an output file
* -e saves the data. 

All information is saved into an Excel file, resembling the following: 
![image](https://user-images.githubusercontent.com/50304901/111376012-f16d6380-866c-11eb-80f4-45c123e63db0.png)

On a separate sheet, the standard error of quality controls will be displayed. Generally, a value less than .1 is desired. 
A two-way ANOVA is also calculated on this sheet, detailing p values for each individual condition, as well as interactions between said conditions. 

On a third sheet, basic plots of the data are shown, along with error bars, in bar graph form. 

Finally, two commands are used for user convenience
* -l <variable name> splits data based on parameters in the reference file. For instance, "-l Stress" would analyze data in stressed and non-stressed groups. 
* -a reshapes, interprets, and plots from a single file. However, it cannot be used with -l

Details on radioimmunoassay can be found here: 
https://www.antibodies-online.com/resources/17/1215/radioimmunoassay-ria
