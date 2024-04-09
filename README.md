# HDXBoxeR

HDXBoxeR is a tool designed to streamline various aspects of HDX (Hydrogen-Deuterium Exchange) data analysis:

1. **Data Reprocessing:** Reprocesses data to the format required for data publication.
2. **Parameter Calculation:** Calculates parameters for a general HDX summary table, including backexchange, peptide lengths, and statistical information.
3. **Output Conversion:** Converts output from HDXExaminer (formerly Sierra Analytics, now Trajan Scientific) to a more manageable and analyzable format.
4. **Statistical Analysis:** Identifies peptides significantly different between sets using Welch T-tests.
5. **Script Generation:** Generates scripts for Pymol visualization.
6. **Plot Generation:** Facilitates easy plot generation, including heat maps, robot plots (modified butterfly plots), significant peptides, volcano plots, and average deuteration.
7. **ExtReme Inputs:** Generates inputs for ExtReme, enabling comparison between different protein states and discrimination of significantly different peptides.

## Usage

HDXBoxeR can be used to enhance the efficiency and effectiveness of HDX data analysis, providing tools for data processing, analysis, and visualization.

## Installation & Loading

First, ensure that you have R (and RStudio) installed on your computer.
Then you can proceed with the installation of the HDXBoxeR package.
There are two methods to install the HDXBoxeR package:

1.  ## CRAN
To install, use the following command:

install.packages("HDXBoxeR") #execute only once

2.  ## GitHub
The HDXBoxeR package is available on GitHub and can be installed using the devtools package.
This method of installation is an alternative to the CRAN method.

install.packages("devtools") #if not installed
library(devtools) #run next two commends only once
devtools::install_github("mkajano/HDXBoxeR")

## Loading
library(HDXBoxeR)

## License

This project is licensed under the  GPL (>= 2) licence

