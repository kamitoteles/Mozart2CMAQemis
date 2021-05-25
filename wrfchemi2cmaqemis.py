
#+---------------------------------------------------------------------
#+ Python Script for converting wrfchemi files into CMAQ emission files
#+ Check the 'CHANGE' comments to define your directories for the run
#+ Author: Camilo Moreno
#+ Email = cama9709@gmail.com
#* WARNING: All input wrfchemi files must be hourly emission files
#* whith one layer on Z dimession.
#* The script will only convert the files if there are all 24 hours
#* of a day (hr 00 to 23) in wrfchemi files in the input directory 
#* and if each day has a 00 hr file for the netx day.
#+---------------------------------------------------------------------
#%%
import re
import glob 
import openpyxl
import numpy as np
import pandas as pd
from netCDF4 import Dataset
from collections import Counter
from os import listdir, scandir, getcwd
from datetime import datetime, date, timedelta


def get_complete_days(days):
    """Return the value of days that have 24 wrf hour files and that at
    least exist the 00 hr file for the netx day.

    Keyword arguments:
    days -- list containing the dates_hours in 'YYYYDDDHH' format
    """

    days_only = [day[:7] for day in days]
    days_only = list(dict.fromkeys(days_only))

    clean_days = []
    for day_o in days_only:
        sum_days = 0
        for day in days:
            if day_o in day:
                sum_days = sum_days + 1
        if sum_days == 24 and str(int(day_o) + 1) + '0' in days:
            clean_days.append(int(day_o))
    
    return clean_days

def get_wrffiles_days(wrf_dir):
    """Get all files that begins whith 'wrfchemi_' from the given 
    wrf_dir directory and extract the complete days found in there.

    Keyword arguments:
    wrf_dir -- string of the directory whre wrfchemi hour files are
    """

    all_files = [f for f in glob.glob(f'{wrf_dir}/wrfchemi_*')]
    all_files.sort()

    days = []
    for file in all_files:
        hr = int(file[-8:-6])
        day = int(file[-11:-9])
        month = int(file[-14:-12])
        year = int(file[-19:-15])

        date_indays = date(year, month, day).timetuple().tm_yday
        space = 3 - len(str(date_indays))
        string_date = str(year) + '0'*space + str(date_indays) + str(hr)
        days.append(string_date)

    days = get_complete_days(days)
    return all_files, days

def get_daywrffiles(all_files, day):
    """Return the 25 files of the day from the all_files parameter.

    Keyword arguments:
    all_files -- list of strings containing all the files found 
    day -- int representing the day of interest in fomrat YYYYDDD
    """

    day_str = str(day)
    year = int(day_str[0:4])
    day_year = int(day_str[4:])

    date = datetime(year, 1, 1) + timedelta(day_year - 1)
    daym = date.timetuple().tm_mday
    str_daym = '0'*(2 - len(str(daym))) + str(daym)
    month = date.timetuple().tm_mon
    str_month = '0'*(2 - len(str(month))) + str(month)

    next_date = datetime(year, 1, 1) + timedelta(day_year)
    next_daym = next_date.timetuple().tm_mday
    str_next_daym = '0'*(2 - len(str(next_daym))) + str(next_daym)
    next_month = next_date.timetuple().tm_mon
    str_next_month = '0'*(2 - len(str(next_month))) + str(next_month)
    next_year = next_date.timetuple().tm_year
    str_next_year = str(next_year)

    date_str = f'{year}-{str_month}-{str_daym}'
    next_date_str = f'{str_next_year}-{str_next_month}-{str_next_daym}_00'


    day_files = [file for file in all_files if date_str in file]
    next_day_00 = [file for file in all_files if next_date_str in file]

    day_files.extend(next_day_00)
    return day_files

def set_conv_map(map_file):
    """Return and set the maping species info in df_map.

    Keyword arguments:
    map_file -- string of the directory location of the mecanism conversion excel
    """

    df_map = pd.DataFrame(pd.read_excel(map_file))
    df_map = df_map[df_map['WRF_SPC'].notnull() & 
                    df_map['CMAQ_SPC'].notnull()]
    df_map.set_index(pd.Index(range(0,len(df_map))), inplace = True)
    return df_map

def create_cmaq_spc(df_map):
    """Return the empty dictionary of the CMAQ temporal especies 
    and the CMAQ species that will be created in the netCDF files.

    Keyword arguments:
    df_map -- data frame of the mecanism conversion map
    """

    cmaq_spc_names = list(df_map['CMAQ_SPC'])
    cmaq_spc_names = list(dict.fromkeys(cmaq_spc_names))
    dic_cmaq = {}

    for spc in cmaq_spc_names:
        dic_cmaq[spc] = np.array([])
    
    return dic_cmaq, cmaq_spc_names

def array_conv(day_files, df_map, dic_cmaq, day, hr):
    """Return filled cmaq arrays converted from wrfchemi species.

    Keyword arguments:
    day_files -- list of strings directions of the files for the day
    df_map -- data frame of the mecanism conversion map
    dic_cmaq -- empty dictionary of the CMAQ temporal especies
    day -- int of the day in format YYYYDDD
    hr -- initial hour if CMAQ files
    """

    count = 0
    for file in day_files:
        ds_wrf = Dataset(file, open = True, mode = 'r')
        wrf_dx = ds_wrf.DX
        wrf_dy = ds_wrf.DY

        for i in range(0,len(df_map)):
            cmaq_str = df_map.iloc[i]['CMAQ_SPC'] + '_temp'

            wrf_str = df_map.iloc[i]['WRF_SPC']
            conv_fact = df_map.iloc[i]['CONV_FACT']

            #CHANGE: In some specific cases, MCIP reduces the original wrfout
            # values by one in each side of the grid. You shoul notice if it is
            # necessary for your case to do the same for wrfarray value below.

            wrf_var = ds_wrf.variables[wrf_str]
            wrf_arr = np.ma.getdata(wrf_var[:][:][:][:])
            wrf_arr = wrf_arr[:,:,1:-1,1:-1]  # Uncoment if final grid is one cell smaller in all borders
         
            if wrf_var.units == 'mol km^-2 hr^-1':
                unit_conv = (wrf_dx/1000) * (wrf_dy/1000) / 3600
                if df_map.iloc[i]['UNITS_SDA'] == 'g/s':
                    mol_weigth = df_map.iloc[i]['MW']
                    unit_conv = unit_conv * mol_weigth
                if len(dic_cmaq[df_map.iloc[i]['CMAQ_SPC']]) == 0:
                    dic_cmaq[cmaq_str] = (wrf_arr * conv_fact * unit_conv)
                else:
                    dic_cmaq[cmaq_str] = dic_cmaq[cmaq_str] + (wrf_arr * conv_fact * unit_conv)

            elif wrf_var.units == 'ug m^-2 s^-1':
                unit_conv = (wrf_dx) * (wrf_dy) / 1000000
                if df_map.iloc[i]['UNITS_SDA'] == 'moles/s':
                    mol_weigth = df_map.iloc[i]['MW']
                    unit_conv = unit_conv / mol_weigth
                if len(dic_cmaq[df_map.iloc[i]['CMAQ_SPC']]) == 0:
                    dic_cmaq[cmaq_str] = (wrf_arr * conv_fact * unit_conv)
                else:
                    dic_cmaq[cmaq_str] = dic_cmaq[cmaq_str] + (wrf_arr * conv_fact * unit_conv )
        
        # TFLAG array
        dic_cmaq['TFLAG_temp'] = np.array([np.tile([day, hr], (len(cmaq_spc_names), 1))])

        if count == 0:
            count = 1
            dic_cmaq['TFLAG'] = dic_cmaq['TFLAG_temp']
            for spc in cmaq_spc_names:
                dic_cmaq[spc] = dic_cmaq[f'{spc}_temp']
                if np.isnan(dic_cmaq[f'{spc}_temp']).any():
                    print(f'ERROR: there are NaN values in the array result of interaction for these variables:\n')
                    print(f'Dictionary: \n {dic_cmaq[f"{spc}_temp"]} \n\nWrf array: \n {wrf_arr} \n\nConv Factor: {conv_fact}\n\nUnit Conversor: {pt_unit_conv}\n\nMol weight: {mol_weigth}')

        else: 
            dic_cmaq['TFLAG'] = np.concatenate((dic_cmaq['TFLAG'], dic_cmaq['TFLAG_temp']), axis = 0)
            for spc in cmaq_spc_names:
                dic_cmaq[spc] = np.concatenate((dic_cmaq[spc], dic_cmaq[f'{spc}_temp']), axis = 0)
        
        day = day + int((hr / 230000))
        hr = (hr + 10000) % 240000
    return dic_cmaq, ds_wrf

def monthly_date(day):
    """convert date from YYYYDDD to YYYYMMDD.

    Keyword arguments:
    day -- int of the day in format YYYYDDD
    """
    day_str = str(day)
    year = int(day_str[0:4])
    day_y = int(day_str[4:])

    date = datetime(year, 1, 1) + timedelta(day_y - 1)
    daym = date.timetuple().tm_mday
    str_daym = '0'*(2 - len(str(daym))) + str(daym)
    month = date.timetuple().tm_mon
    str_month = '0'*(2 - len(str(month))) + str(month)

    return str(year) + str_month + str_daym

def create_ncfile(save_dir, day, hr, ds_wrf, cmaq_spc_names, dic_cmaq, df_map):
    """Create Final NETCDF file.

    Keyword arguments:
    save_dir -- string of the location for saving the netCDF files
    day -- int of the day in format YYYYDDD
    hr -- initial hour if CMAQ files
    ds_wrf -- last wrf_dataset accessed
    cmaq_spc_names -- list of strings directions of the files for the day
    dic_cmaq -- empty dictionary of the CMAQ temporal especies
    df_map -- data frame of the mecanism conversion map
    """

    #* Define COLS, ROWS, and VAR values
    #CHANGE: In some specific cases, MCIP reduces the original wrfout
    # values by one in each side of the grid. You shoul notice if it is
    # necessary for your case to do the same for cols and rows values below.

    #cols = len(ds_wrf.dimensions['west_east'])
    #rows = len(ds_wrf.dimensions['south_north'])
    cols = len(ds_wrf.dimensions['west_east']) - 2 # Uncoment if MCIP reduced the wrfout grid laterals by one cell
    rows = len(ds_wrf.dimensions['south_north']) - 2 # Uncoment if MCIP reduced the wrfout grid laterals by one cell

    num_vars = len(cmaq_spc_names)

    #* Create new netCDF
    day_monthly = monthly_date(day)
    new_cmaq_file = f'{save_dir}/Emis_CMAQ_{day_monthly}.ncf'
    ds_new_cmaq = Dataset(new_cmaq_file, open = True, mode = 'w', format=  "NETCDF3_CLASSIC")

    #* Create dimenssions
    TSTEP = ds_new_cmaq.createDimension("TSTEP", None)
    DATE_TIME = ds_new_cmaq.createDimension("DATE-TIME", len(dic_cmaq['TFLAG'][0][0][:]))
    VAR = ds_new_cmaq.createDimension("VAR", num_vars)
    LAY = ds_new_cmaq.createDimension("LAY", len(ds_wrf.dimensions['emissions_zdim_stag']))
    ROW = ds_new_cmaq.createDimension("ROW", rows)
    COL = ds_new_cmaq.createDimension("COL", cols)

    #* Create variables
    tflag = ds_new_cmaq.createVariable("TFLAG","i4",("TSTEP","VAR", "DATE-TIME"))
    tflag.units = '<YYYYDDD,HHMMSS>'
    tflag.long_name = 'TFLAG'
    tflag.var_desc = 'Timestep-valid flags:  (1) YYYYDDD or (2) HHMMSS'

    for spc in cmaq_spc_names:
        unt = list(df_map[df_map['CMAQ_SPC'] == spc]['UNITS_SDA'])[0]
        var_temp = ds_new_cmaq.createVariable(spc,"f4",("TSTEP", "LAY", "ROW", "COL"))
        var_temp.units = unt
        var_temp.long_name = spc
        var_temp.var_desc = f'Model species {spc}'

    #* Fill variables
    ds_new_cmaq.variables['TFLAG'][:, :, :] = dic_cmaq['TFLAG']
    for spc in cmaq_spc_names:
        ds_new_cmaq.variables[spc][:, :, :] = dic_cmaq[spc]

    #* Creatae attributes
    varlist = ''
    for spc in cmaq_spc_names:
        space = 16 - len(spc)
        varlist = varlist + spc + ' '*space


    #CHANGE: This dictionary contains the atrribute values of the netCDF
    # file. It follows a IOAPI format and it may be edited to fit the 
    # wanted grid. See this link for guidance:
    # https://cmascenter.org/ioapi/documentation/all_versions/html/TUTORIAL.html
    cmaq_attrs = {'IOAPI_VERSION': '$Id: @(#) ioapi library version 3.1 $',
                'EXEC_ID': '????????????????',
                'FTYPE': np.int32(1),
                'CDATE': np.int32(2021100),
                'CTIME': np.int32(185404),
                'WDATE': np.int32(2021100),
                'WTIME': np.int32(185404),
                'SDATE': np.int32(day),
                'STIME': np.int32(hr),
                'TSTEP': np.int32(10000),
                'NTHIK': np.int32(1),
                'NCOLS': np.int32(cols),
                'NROWS': np.int32(rows),
                'NLAYS': np.int32(len(ds_wrf.dimensions['emissions_zdim_stag'])),
                'NVARS': np.int32(num_vars),
                'GDTYP': np.int32(2), # 2 is Lambert Conformal Conic
                'P_ALP': np.float64(-9.26763916015625),
                'P_BET': np.float64(19.6556396484375),
                'P_GAM': np.float64(-72.6330032348633),
                'XCENT': np.float64(-72.6330032348633),
                'YCENT': np.float64(5.19400024414062),
                'XORIG': np.float64(-1688312.),
                'YORIG': np.float64(-1592931.5),
                'XCELL': np.float64(ds_wrf.DX),
                'YCELL': np.float64(ds_wrf.DY),
                'VGTYP': np.int32(-1),
                'VGTOP': np.float32(0.0),
                'VGLVLS': np.array([np.float32(0.), np.float32(0.)]),
                'GDNAM': "2018_NSthAm_CROS",
                'UPNAM': "M3WNDW",
                'VAR-LIST': varlist,
                'FILEDESC': "Point source emissions data/ Camilo Moreno",
                'HISTORY': 'La historia es historia',}

    for attr in cmaq_attrs:
        ds_new_cmaq.setncattr(attr, cmaq_attrs[attr])

    #+ Close new netcdf file
    ds_new_cmaq.close()

    print(f"{day_monthly} emissions file DONE")

# %%
if __name__ == "__main__":
    #CHANGE: map_file is the path of the excel file where the map conversor is alocated
    map_file = '/Users/camilo/OneDrive - Universidad de los Andes/Estudio/Tesis_maestría/Code/Emissions_conversor/mozart2cbo5_conv_table.xlsx'

    #CHANGE: wrf_dir is the directory path where are all the wrfchemi files you want to comvert to CMAQ emission files
    wrf_dir = '/Volumes/Avispa/Emissions/wrfchemi_original/MOZART_sep'

    #CHANGE: save_dir is the directory path where you want to save all the new netCDF emission files
    save_dir = '/Volumes/Avispa/Emissions/CMAQ_emis/Sep_2018'
    
    df_map = set_conv_map(map_file)
    dic_cmaq, cmaq_spc_names = create_cmaq_spc(df_map)
    all_files, days = get_wrffiles_days(wrf_dir)
    hr = 0
    for day in days:
        day_files = get_daywrffiles(all_files, day)
        dic_cmaq_fill, ds_wrf = array_conv(day_files, df_map, dic_cmaq, day, hr)
        create_ncfile(save_dir, day, hr, ds_wrf, cmaq_spc_names, dic_cmaq_fill, df_map)
        ds_wrf.close()

        del day_files 
        del dic_cmaq_fill 
        del ds_wrf
    print('All files completed')
# %%