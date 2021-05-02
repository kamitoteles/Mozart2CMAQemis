# WRF-chemi to CMAQ emissions files

wrfchemi2camqemis.py is a python script designed to convert WRF-chemi hour emission files into 25 hours CMAQ emissions files.

## Getting the wrfchemi2cmaqemis conversor

1. Clone the repository into your system:

```bash
git clone https://github.com/kamitoteles/Mozart2CMAQemis.git
```

2. Install the required library versions specified in the requirements.txt file.

## Usage

### Input files

The intput WRF-chemi files to convert may be hourly files and have to be saved in the same input directory folder. The script would only convert the complete valid days that it founds.

**Comple valid day:** is defined as a day that have 24 WRF-chemi files for this day (from hour 00 to hour 23) and a WRF-chemi file that for the 00 hour of the netxt day. (eg. for the 2018_08_03 it shoul be 24 files whith this date and a file for the 2018_08_04_00 for the 2018_08_03 to be a complete valid day)

Alaso, all the WRF-chemi files must be named by the format "wrfchemi_d##_YYYY-MM-DD_hh_mm_ss" (eg. wrfchemi_d01_2018-09-01_03_00_00). Where:

- "\#\#" specifies the domain number of the grid (eg. "01" if domain is 1)
- "YYYY-MM-DD_hh_mm_ss" is the 4 digit year, two digit month, two digit day, two digit hour, two digit minute, and two digit second of the file.

**IMPORTANT:** if the day is not a complete valid day, the script would pass these date and convert the others that are complete in the input directory.

### Map conversor files

This excel are provided in order to give the script the information of how to convert between diferent WRF chemichal and aerosols mecahnisms into CMAQ cbo5_aero5_aq mechanism. Before running the script

- mozart2cbo5_conv_table.xlsx if your wrfchemi files are in Mozart
- racm2cbo5_conv_table.xlsx if your wrfchemi files are in RACM

**WARNNING:** these files are in beta testing and could be incomplete for certain species. Please check and edit the values before using the script to ensure better results for your cases.

## Proyection, coordinates and grid

You **MUST** change the cmaq_attrs dictionary values in the create_ncfile() function. There are de values for all the final CMAQ IO/API attributes that will be saved in the netCDF file. The most important are the ones related with the grid description. 

NCOLS, NROWS, NLAYS, XCELL, and YCELL are taken from de original WRF-chemi file values, but all grid attributes (P_ALP, P_BET, P_GAM, XCENT, YCENT, XORIG, and YORIG) are predefined for a test grid,  so you must change this for your specific case.

## Output

The output files would be saved whith the format "Emis_CMAQ_YYYYDDD.ncf", where the YYYYDDD references the date of the day whit the year "YYYY" and the day number for thar year "DDD".
