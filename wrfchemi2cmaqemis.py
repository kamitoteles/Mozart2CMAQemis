
#%%
import pandas as pd
import re
from netCDF4 import Dataset
import openpyxl
import numpy as np
import glob 
from os import listdir
from os import scandir, getcwd


#+ Read and set maping species info in df_map
def set_conv_map(map_file):
    df_map = pd.DataFrame(pd.read_excel(map_file))
    df_map = df_map[df_map['WRF_SPC'].notnull() & 
                    df_map['CMAQ_SPC'].notnull()]
    df_map.set_index(pd.Index(range(0,len(df_map))), inplace = True)
    return df_map


#+ Create empty arrays of the CMAQ temporal and final especies
def create_cmaq_spc(df_map):
    cmaq_spc_names = list(df_map['CMAQ_SPC'])
    cmaq_spc_names = list(dict.fromkeys(cmaq_spc_names))
    dic_cmaq = {}

    for spc in cmaq_spc_names:
        dic_cmaq[spc] = np.array([])
    
    return dic_cmaq, cmaq_spc_names


#+ Get all wrfout files from a directory
def get_wrffiles(wrf_dir):
    all_files = [f for f in glob.glob(f'{wrf_dir}/wrfchemi_*')]
    all_files.sort()

    day_0 = int(all_files[0][-11:-9])
    hr_0 = int(all_files[0][-8:-6])
    
    day_end = int(all_files[-1][-11:-9])
    hr_end = int(all_files[-1][-8:-6])

    return all_files, day_0, hr_0, day_end, hr_end


#+ Main loop for converting wrfchemi arrays into CMAQ arrays
#! It create a unique array for each variable that extends trougth the numbre of files in 'all_files'
#! asuming that each wrfchemi file contins a unique hour of data
def array_conv(all_files, df_map, dic_cmaq, day, hr):
    count = 0
    for file in all_files:
        ds_wrf = Dataset(file, open = True, mode = 'r')
        wrf_dx = ds_wrf.DX
        wrf_dy = ds_wrf.DY

        for i in range(0,len(df_map)):
            gas_unit_conv = ((wrf_dx * wrf_dy) / 1000) / 3600
            pt_unit_conv = (wrf_dx * wrf_dy) / 1000000

            cmaq_str = df_map.iloc[i]['CMAQ_SPC'] + '_temp'

            wrf_str = df_map.iloc[i]['WRF_SPC']
            conv_fact = df_map.iloc[i]['CONV_FACT']

            wrf_var = ds_wrf.variables[wrf_str]
            wrf_arr = np.ma.getdata(wrf_var[:][:][:][:])

            if len(dic_cmaq[df_map.iloc[i]['CMAQ_SPC']]) == 0:
                dic_cmaq[cmaq_str] = wrf_arr * conv_fact
    
            else:
                if wrf_var.units == 'mol km^-2 hr^-1':
                    if df_map.iloc[i]['UNITS_SDA'] == 'g/s':
                        mol_weigth = df_map.iloc[i]['MW']
                        gas_unit_conv = gas_unit_conv * mol_weigth
                    dic_cmaq[cmaq_str] = dic_cmaq[cmaq_str] + (wrf_arr * conv_fact * gas_unit_conv)

                elif wrf_var.units == 'ug m^-2 s^-1':
                    if df_map.iloc[i]['UNITS_SDA'] == 'moles/s':
                        mol_weigth = df_map.iloc[i]['MW']
                        pt_unit_conv = pt_unit_conv / mol_weigth
                    dic_cmaq[cmaq_str] = dic_cmaq[cmaq_str] + (wrf_arr * conv_fact * pt_unit_conv )
        
        # TFLAG array
        dic_cmaq['TFLAG_temp'] = np.array([np.tile([day, hr], (len(cmaq_spc_names), 1))])

        if count == 0:
            count = 1
            dic_cmaq['TFLAG'] = dic_cmaq['TFLAG_temp']
            for spc in cmaq_spc_names:
                dic_cmaq[spc] = dic_cmaq[f'{spc}_temp']
        else: 
            dic_cmaq['TFLAG'] = np.concatenate((dic_cmaq['TFLAG'], dic_cmaq['TFLAG_temp']), axis = 0)
            for spc in cmaq_spc_names:
                dic_cmaq[spc] = np.concatenate((dic_cmaq[spc], dic_cmaq[f'{spc}_temp']), axis = 0)
        
        day = day + int((hr / 230000))
        hr = (hr + 10000) % 240000
    return dic_cmaq, ds_wrf


#+ Create Final NETCDF file
def create_ncfile(ncdir, day, hr, ds_wrf, cmaq_spc_names, dic_cmaq, df_map):
    new_cmaq_file = f'{ncdir}/Emis_CMAQ_{day}.ncf'
    ds_new_cmaq = Dataset(new_cmaq_file, open = True, mode = 'w', format=  "NETCDF3_CLASSIC")

    #+ Create dimenssions
    TSTEP = ds_new_cmaq.createDimension("TSTEP", None)
    LAY = ds_new_cmaq.createDimension("LAY", len(ds_wrf.dimensions['emissions_zdim_stag']))
    ROW = ds_new_cmaq.createDimension("ROW", len(ds_wrf.dimensions['south_north']))
    COL = ds_new_cmaq.createDimension("COL", len(ds_wrf.dimensions['west_east']))
    VAR = ds_new_cmaq.createDimension("VAR", len(cmaq_spc_names))
    DATE_TIME = ds_new_cmaq.createDimension("DATE-TIME", len(dic_cmaq['TFLAG'][0][0][:]))

    #+ Create variables
    for spc in cmaq_spc_names:
        unt = list(df_map[df_map['CMAQ_SPC'] == spc]['UNITS_SDA'])[0]
        var_temp = ds_new_cmaq.createVariable(spc,"f4",("TSTEP", "LAY", "ROW", "COL"))
        var_temp.units = unt
        var_temp.long_name = spc
        var_temp.var_desc = f'Model species {spc}'

    tflag = ds_new_cmaq.createVariable("TFLAG","i4",("TSTEP","VAR", "DATE-TIME"))
    tflag.units = '<YYYYDDD,HHMMSS>'
    tflag.long_name = 'TFLAG'
    tflag.var_desc = 'Timestep-valid flags:  (1) YYYYDDD or (2) HHMMSS'

    #+ Fill variables
    for spc in cmaq_spc_names:
        ds_new_cmaq.variables[spc][:, :, :] = dic_cmaq[spc]
    ds_new_cmaq.variables['TFLAG'][:, :, :] = dic_cmaq['TFLAG']

    #+ Creatae attributes
    varlist = ''
    for spc in cmaq_spc_names:
        space = 16 - len(spc)
        varlist = varlist + spc + ' '*space

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
                'NCOLS': np.int32(len(ds_wrf.dimensions['west_east'])),
                'NROWS': np.int32(len(ds_wrf.dimensions['south_north'])),
                'NLAYS': np.int32(len(ds_wrf.dimensions['emissions_zdim_stag'])),
                'NVARS': np.int32(len(cmaq_spc_names)),
                'GDTYP': np.int32(2),
                'P_ALP': np.float64(-9.26763916015625),
                'P_BET': np.float64(19.6556396484375),
                'P_GAM': np.float64(-72.6330032348633),
                'XCENT': np.float64(-72.6330032348633),
                'YCENT': np.float64(5.19400024414062),
                'XORIG': np.float64(-1688312.),
                'YORIG': np.float64(-1592931.5),
                'XCELL': np.float64(27000.),
                'YCELL': np.float64(27000.),
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

    #+ Close netcdf file
    ds_new_cmaq.close()

    print("DONE")
# %%
if __name__ == "__main__":
    #! CHANGE, initial day and hour of the one day netCDF file
    day = 2018035
    hr = 0

    wrf_dir = '/Volumes/Avispa/wrfchemi_original/wrfchemi_d01_2018_09'
    all_files, day_0, hr_0, day_end, hr_end = get_wrffiles(wrf_dir)

    map_file = '/Users/camilo/OneDrive - Universidad de los Andes/Estudio/Tesis_maestría/Code/Emissions_conversor/map_conv_table.xlsx'
    ncdir = '/Users/camilo/Desktop'

    df_map = set_conv_map(map_file)
    dic_cmaq, cmaq_spc_names = create_cmaq_spc(df_map)
    
    dic_cmaq, ds_wrf = array_conv(all_files, df_map, dic_cmaq, day, hr)
    create_ncfile(ncdir, day, hr, ds_wrf, cmaq_spc_names, dic_cmaq, df_map)
    ds_wrf.close()

# %%