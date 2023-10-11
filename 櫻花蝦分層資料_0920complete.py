import pandas as pd
import netCDF4  # 讀取分析.nc檔
import numpy as np  # 轉換矩陣維度
import csv  # 檔案格式
import os
import math
from datetime import datetime
import warnings
print(min(20,9))
def addmiddle(alist):
    n = len(alist)
    new_list = []
    new_list.extend(alist)
    for i in range(n - 1):
        a = (alist[i] + alist[i + 1]) / 2
        new_list.insert(2 * i + 1, a)
    return new_list
##把一個二為數列的精度切成一半
def two_middle(b_array):
    lon_n=len(b_array[0,:])
    lat_n=len(b_array[:,0])
    new_arr=np.empty(shape=[2*lat_n-1,0])

    for i in range(lon_n):
        a = []
        a.extend(b_array[:, i])
        a=addmiddle(a)
        new_arr=np.column_stack((new_arr,a))

    new_arr2=np.empty(shape=[0,2*lon_n-1])
    for i in range(2*lat_n-1):
        a=[]
        a.extend(new_arr[i,:])
        a=addmiddle(a)
        new_arr2=np.vstack((new_arr2,a))
    return new_arr2

def fileNameList(filenameEnding):
    # list the file names in the folder
    filename = np.empty(shape=(1, 1), dtype=object)
    root = "./"
    for path, subdirs, files in os.walk(root):
        for x in files:
            # print(os.path.join(path, x))
            if x.endswith(filenameEnding):
                filename = np.column_stack((filename, os.path.join(path, x)))
    filename = filename.flatten()
    filename = sorted(filename[1:len(filename)])
    return filename

def num2dateTime(data):
    times = data.variables['MT']
    dates = netCDF4.num2date(times[:], times.units)
    dateAll = datetime.strptime(str(dates[0]), '%Y-%m-%d %H:%M:%S')
    time = str(dateAll.year) + str(dateAll.month).zfill(2) + str(dateAll.day).zfill(2)
    return time

def lonlatData(data,add_n_times):
    lats, lons = data.variables['Latitude'][:], data.variables['Longitude'][:]

    for i in range(add_n_times):
        lons=addmiddle(lons)
        lats=addmiddle(lats)
    csvData = np.empty(shape=(len(lons) * len(lats) * 9, 3), dtype=object)
    k = 0
    for i in range(len(lons)):
        for j in range(len(lats)):
            for l in range(9):
                csvData[k, :] = np.array([lons[i], lats[j], l])
                k = k + 1
    return csvData

##溫度分層
def tempData(data, csvData, data_from_2d,add_n_times):
    pot_temp = data.variables['pot_temp'][:]
    #temperature = np.array(pot_temp[0, 0:9, :, :])
    setdata = findMixlayer(data_2d=data_from_2d, data_3z=data)
    setdata = np.transpose(setdata)
    print(setdata.shape)
    a1 = len(setdata[:, 0])
    a2 = len(setdata[0, :])
    Tempdata = np.full(shape=(a1, a2, 9), fill_value=999, dtype=float)
    for i in range(a1):
        for j in range(a2):
            for k in range(min(setdata[i, j], 9)):
                Tempdata[i, j, k] =pot_temp[0,k,j,i]
    temperature = np.full(shape=(2 * a1 - 1, 2 * a2 - 1, 9), fill_value=999, dtype=float)
    for m in range(add_n_times):
        for k in range(9):
            temperature[:, :, k] = two_middle(Tempdata[:, :, k])
    temperature=np.where(temperature>100,999,temperature)
    temperature=temperature.flatten()
    csvData = np.column_stack((csvData, temperature))
    return csvData

def processTemp(dataList, outputFileName,datalist2,filterRange,add_times):
    warnings.filterwarnings('ignore')
    if filterRange == None:
        filterRange = [-999, 999, -999, 999]
        filterRange.append(filterRange)

    for i in range(len(dataList)):
        print('Read data from', dataList[i],'and',datalist2[i])
        data = netCDF4.Dataset(dataList[i])
        data_2d=netCDF4.Dataset(datalist2[i])
        time = num2dateTime(data)

        if i == 0:
            csvData = lonlatData(data,add_n_times=add_times)
            header = ['LON', 'LAT', 'layer']

        header = header[:] + time.split()
        csvData = tempData(data=data,csvData=csvData,data_from_2d=data_2d,add_n_times=add_times)

    outputCSV(outputFileName, csvData, header, filterRange)
    print("Save temperature data (in degrees) to CSV file")
    return csvData

def outputCSV(filename, csvData, header, filterRange):
    # Create the dataframe
    df = pd.DataFrame(csvData, columns=header)
    # Apply the function
    df = df.fillna(0)
    rslt_df = df.loc[(df['LAT'] >= filterRange[0]) &
                     (df['LAT'] <= filterRange[1]) &
                     (df['LON'] >= filterRange[2]) &
                     (df['LON'] <= filterRange[3])]
    # print the DataFrame
    rslt_df.to_csv(filename, index=False)
    return

##data1從_2d檔案中抓取
##data2=read from_3z
def findMixlayer(data_2d, data_3z):
    mlayer_thick = data_2d.variables['mixed_layer_thickness'][0, :, :]
    layer_depth = data_3z.variables['Depth'][:]
    a5 = len(mlayer_thick[:, 0])
    a6 = len(mlayer_thick[0, :])
    depth_data = np.zeros((a5, a6), dtype=int)
    for a in range(a5):
        for b in range(a6):
            for i in range(len(layer_depth)):
                if mlayer_thick[a, b] < layer_depth[i]:
                    depth_data[a, b] = i + 1
                    break
    return depth_data
def uvDirData(data, csvData,data_from_2d,add_n_times):
    u, v = data.variables['u'][:], data.variables['v'][:]
    setdata = findMixlayer(data_2d=data_from_2d, data_3z=data)
    setdata = np.transpose(setdata)
    print(setdata.shape)
    a1 = len(setdata[:, 0])
    a2 = len(setdata[0, :])
    UVdata = np.full(shape=(a1, a2, 9), fill_value=999, dtype=float)
    for i in range(a1):
        for j in range(a2):
            for k in range(min(setdata[i,j],9)):
                UVdata[i,j,k]=np.degrees(np.arctan2(v[0,k,j,i],u[0,k,j,i]))
    uvdir=np.full(shape=(2*a1-1,2*a2-1,9),fill_value=999,dtype=float)

    for m in range(add_n_times):
        for k in  range(9):
            uvdir[:,:,k]=two_middle(UVdata[:,:,k])
    uvdir=np.where(uvdir<0,uvdir+360,uvdir)
    uvdir=np.where(uvdir>360,999,uvdir)
    uvdir=uvdir.flatten()
    csvData = np.column_stack((csvData,uvdir))
    return csvData
def processUVDir(dataList, outputFileName, filterRange,datalist2,add_times):
    warnings.filterwarnings('ignore')
    if filterRange == None:
        filterRange = [-999, 999, -999, 999]
        filterRange.append(filterRange)

    for i in range(len(dataList)):
        print('Read data from', dataList[i],datalist2[i])
        data = netCDF4.Dataset(dataList[i])
        data_2d = netCDF4.Dataset(datalist2[i])
        time = num2dateTime(data)

        if i == 0:
            csvData = lonlatData(data,add_n_times=add_times)
            header = ['LON', 'LAT','depth']
        header = header[:] + time.split()
        csvData = uvDirData(data, csvData,data_2d,add_times)
    outputCSV(outputFileName, csvData, header, filterRange)
    print("Save current direction data (in degrees) to CSV file")
    return csvData
def uvvelData(data, csvData,data_from_2d,add_n_times):
    u, v = data.variables['u'][:], data.variables['v'][:]
    setdata = findMixlayer(data_2d=data_from_2d, data_3z=data)
    setdata = np.transpose(setdata)
    print(setdata.shape)
    a1 = len(setdata[:, 0])
    a2 = len(setdata[0, :])
    UVdata = np.full(shape=(a1, a2, 9), fill_value=999, dtype=float)
    for i in range(a1):
        for j in range(a2):
            for k in range(min(setdata[i,j],9)):
                UVdata[i,j,k] = np.sqrt(np.square(v[0,k,j,i]) + np.square(u[0,k,j,i]))
    uvvel=np.full(shape=(2*a1-1,2*a2-1,9),fill_value=999,dtype=float)

    for m in range(add_n_times):
        for k in  range(9):
            uvvel[:,:,k]=two_middle(UVdata[:,:,k])
    uvvel=np.where(uvvel>100,999,uvvel)
    uvvel=uvvel.flatten()
    csvData = np.column_stack((csvData,uvvel))
    return csvData
def processUVvel(dataList, outputFileName, filterRange,datalist2,add_times):
    warnings.filterwarnings('ignore')
    if filterRange == None:
        filterRange = [-999, 999, -999, 999]
        filterRange.append(filterRange)

    for i in range(len(dataList)):
        print('Read data from', dataList[i],datalist2[i])
        data = netCDF4.Dataset(dataList[i])
        data_2d = netCDF4.Dataset(datalist2[i])
        time = num2dateTime(data)

        if i == 0:
            csvData = lonlatData(data,add_n_times=add_times)
            header = ['LON', 'LAT','depth']
        header = header[:] + time.split()
        csvData = uvvelData(data, csvData,data_2d,add_times)
    outputCSV(outputFileName, csvData, header, filterRange)
    print("Save current velocity data (in m/s) to CSV file")

shrimp_range = [20, 28, 118, 122]
filename2D_day=fileNameList("12_2d.nc")
filename3D_day = fileNameList("12_3z.nc")

print("櫻花蝦_至混合層流向")
processUVDir(dataList=filename3D_day,outputFileName='櫻花_至混合層蝦流向.csv',filterRange=shrimp_range,datalist2=filename2D_day,add_times=1)

print("櫻花蝦_至混合層流速")
processUVvel(dataList=filename3D_day,outputFileName='櫻花_至混合層蝦流速.csv',filterRange=shrimp_range,datalist2=filename2D_day,add_times=1)

print("櫻花蝦_至混合層的水溫")
processTemp(dataList=filename3D_day,outputFileName='櫻花蝦至混合層的水溫.csv',filterRange=shrimp_range,datalist2=filename2D_day,add_times=1)
