import pandas as pd
import csv
import netCDF4
import numpy as np
import os
import math
import warnings
from datetime import datetime

### 在一個一維數列中內插中間值
def addmiddle(alist):
    n = len(alist)
    new_list = []
    new_list.extend(alist)
    for i in range(n - 1):
        a = (alist[i] + alist[i + 1]) / 2
        new_list.insert(2 * i + 1, a)
    return new_list
##對二維數列數列執行一次內插，size變成原來的(2n-1)^2倍
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
            if x.endswith(filenameEnding):
                filename = np.column_stack((filename, os.path.join(path, x)))

    filename = filename.flatten()
    filename = sorted(filename[1:len(filename)])
    return filename
#轉換時間參數
def num2dateTime(data):
    times = data.variables['MT']
    dates = netCDF4.num2date(times[:], times.units)
    dateAll = datetime.strptime(str(dates[0]), '%Y-%m-%d %H:%M:%S')
    time = str(dateAll.year) + str(dateAll.month).zfill(2) + str(dateAll.day).zfill(2)
    return time
def depthLayerIndex(data,layerLoc):
    L1=int(len(data[0,0,:,0]))
    L2=int(len(data[0, 0, 0, :]))

    totNumLayer=len(np.array(data[0, :, 0, 0]))
    layerNum = np.zeros((L1,L2),dtype=int)
    if layerLoc == 0:
        numLayer = 1
    elif layerLoc == -1:
        numLayer = 1
        for i in range(L1):
            for j in range(L2):
                depth_valid = np.array(data[0, :, i-1, j-1])
                tmp1=np.where(depth_valid[depth_valid>100])
                tmp2=np.array(tmp1)
                if tmp2.size == 0:
                    layerNum[i-1,j-1] = int(totNumLayer-1)
                else:
                    layerNum[i-1,j-1]=int(totNumLayer-2-max(tmp2[0]))
                    if layerNum[i-1,j-1] < 0:
                       layerNum[i-1,j-1] = 0
    else:
        print('Wrong water depth indexing')

    return layerNum, numLayer
def lonlatData(data,add_n_times):
    # get all variable names
    # print(data.variables.keys())
    lats, lons = data.variables['Latitude'][:], data.variables['Longitude'][:]
    #print(lons[268])

    #定義精度成0.4(使用定義的函式addmiddle)
    for i in range(add_n_times):
        lons=addmiddle(lons)
        lats=addmiddle(lats)

    csvData = np.empty(shape=(len(lons)*len(lats),2), dtype=object)
    k = 0
    for i in range(len(lons)):
        for j in range(len(lats)):
            csvData[k, :] = np.array([lons[i], lats[j]])
            k+=1
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

#海水溫度資料
def tempData(data, csvData, layerLoc,add_n_times):
    pot_temp = data.variables['pot_temp'][:]  # temperature variable
    [layerNum,numLayer] = depthLayerIndex(pot_temp, layerLoc)
    if layerLoc == 0:
        temperature = np.array(pot_temp[0, 0, :, :])
    elif layerLoc == -1:
        temperature = np.array(pot_temp[0, -1, :, :])
        for i in range(len(layerNum)):
           for j in range(len(layerNum[0])):
                if math.isnan(pot_temp[0, layerNum[i - 1, j - 1], i - 1, j - 1]):
                    temperature[i - 1, j - 1] = 999
                else:
                    temperature[i - 1, j - 1] = pot_temp[0, layerNum[i - 1, j - 1], i - 1, j - 1]
    for i in range(add_n_times):
        temperature = two_middle(temperature)
    temperature=np.transpose(temperature)
    temperature = temperature.flatten()
    temperature = [999.99 if ele > 100 else ele for ele in temperature]
    csvData = np.column_stack((csvData, temperature))
    return csvData
##add_times=內插執行次數
def processTemp(dataList, outputFileName, layerLoc, filterRange,add_times):
    warnings.filterwarnings('ignore')
    if filterRange == None:
        filterRange = [-999, 999, -999, 999]
        filterRange.append(filterRange)

    for i in range(len(dataList)):
        print('Read data from', dataList[i])
        data = netCDF4.Dataset(dataList[i])
        time = num2dateTime(data)

        if i == 0:
            csvData = lonlatData(data,add_n_times=add_times)
            header = ['LON', 'LAT']

        header = header[:] + time.split()
        csvData = tempData(data, csvData, layerLoc,add_n_times=add_times)

    outputCSV(outputFileName, csvData, header, filterRange)
    print("Save temperature data (in degrees) to CSV file")
    return csvData

#洋流流向(角度)
def uvDirData(data, csvData,layerLoc,add_n_times):
    u, v = data.variables['u'][:], data.variables['v'][:]
    if layerLoc == 0:
        u=np.where(u>100,1,u)
        v=np.where(v>100,0,v)
        currentDir = np.arctan2(v[0,0, :, :], u[0,0, :, :])
        currentDir = np.degrees(currentDir)
        currentDir = np.where(currentDir==0,999, currentDir)
    elif layerLoc == -1:
        currentDir = np.arctan2(v[0,-1, :, :], u[0,-1, :, :])
        for i in range(len(u[0,0,:,0])):
            for j in range(len(u[0,0,0,:])):
                for k in range(len(u[0,:,0,0,])):
                    if math.isnan(u[0,0,i,j]):
                        currentDir[i, j] = 999
                        break
                    elif math.isnan(u[0,k,i,j]):
                        currentDir[i, j] =np.degrees(np.arctan2(v[0, k-1, i, j],u[0,k-1,i,j]))
                        break
    currentDir = np.where(currentDir < 0, currentDir + 360, currentDir)

    for i in range(add_n_times):
        currentDir=two_middle(currentDir)
    currentDir = np.where((currentDir > 360), 999, currentDir)
    currentDir=np.transpose(currentDir)
    currentDir=currentDir.flatten()
    csvData = np.column_stack((csvData, currentDir))
    return csvData
def processUVDir(dataList, outputFileName, layerLoc, filterRange,add_times):
    warnings.filterwarnings('ignore')
    if filterRange == None:
        filterRange = [-999, 999, -999, 999]
        filterRange.append(filterRange)

    for i in range(len(dataList)):
        print('Read data from', dataList[i])
        data = netCDF4.Dataset(dataList[i])
        time = num2dateTime(data)
        if i == 0:
            csvData = lonlatData(data,add_n_times=add_times)
            header = ['LON', 'LAT']
        header = header[:] + time.split()
        csvData = uvDirData(data, csvData, layerLoc,add_n_times=add_times)
    outputCSV(outputFileName, csvData, header, filterRange)
    print("Save current direction data (in degrees) to CSV file")

#洋流流速
def uvVelData(data, csvData,layerLoc,add_n_times):
    u, v = data.variables['u'][:], data.variables['v'][:]
    vel_u = np.where(u > 100, 999, u)
    vel_v=np.where(v>100, 999,v)
    if layerLoc == 0:
        vel_v0=vel_v[0,0,:,:]
        vel_u0=vel_u[0,0,:,:]
        vel = np.sqrt(np.square(vel_v0) + np.square(vel_u0))
    elif layerLoc == -1:
        vel_v0 = vel_v[0, -1, :, :]
        vel_u0 = vel_u[0, -1, :, :]
        for i in range(len(u[0,0,:,0])):
            for j in range(len(u[0,0,0,:])):
                for k in range(len(u[0,:,0,0,])):
                    if math.isnan(u[0,0,i,j]):
                        vel_v0[i, j] = 999
                        vel_u0[i, j] = 999
                        break
                    elif math.isnan(u[0,k,i,j]):
                        vel_v0[i, j] =vel_v[0, k-1, i, j]
                        vel_u0[i, j] =vel_u[0, k-1,i,j]
                        break

        vel = np.sqrt(np.square(vel_v0) + np.square(vel_u0))
    for i in range(add_n_times):
        vel = two_middle(vel)
    vel = np.where(vel > 200, 999, vel)
    vel=np.transpose(vel)
    vel=vel.flatten()

    csvData = np.column_stack((csvData, vel))
    return csvData
def processUVvel(dataList, outputFileName, layerLoc, filterRange,addtimes):
    warnings.filterwarnings('ignore')
    if filterRange == None:
        filterRange = [-999, 999, -999, 999]
        filterRange.append(filterRange)

    for i in range(len(dataList)):
        print('Read data from', dataList[i])
        data = netCDF4.Dataset(dataList[i])
        time = num2dateTime(data)

        if i == 0:
            csvData = lonlatData(data,add_n_times=addtimes)
            header = ['LON', 'LAT']
        header = header[:] + time.split()
        csvData = uvVelData(data, csvData, layerLoc,add_n_times=addtimes)
    outputCSV(outputFileName, csvData, header, filterRange)
    print("Save current speed data (in m/s) to CSV file")

#海面高
def sshData(data, csvData,add_n_times):
    ssh = data.variables['ssh'][:]  # temperature variable
    ssh = np.array(ssh[0, :, :])
    for i in range(add_n_times):
        ssh=two_middle(ssh)
    ssh=np.transpose(ssh)
    ssh = ssh.flatten()
    ssh = [999.99 if ele > 100 else ele for ele in ssh]
    csvData = np.column_stack((csvData, ssh))
    return csvData
def processSSH(dataList, outputFileName, addtimes,filterRange=None,):
    if filterRange == None:
        filterRange = [-999, 999, -999, 999]
        filterRange.append(filterRange)
    for i in range(len(dataList)):
        print('Read data from', dataList[i])
        data = netCDF4.Dataset(dataList[i])
        time = num2dateTime(data)
        if i == 0:
            csvData = lonlatData(data,add_n_times=addtimes)
            header = ['LON', 'LAT']
        header = header[:] + time.split()
        csvData = sshData(data, csvData,add_n_times=addtimes)
    outputCSV(outputFileName, csvData, header, filterRange)
    print("Save surface elevation data (in meters) to CSV file")
    return csvData

#混合層深
def mixed_layer_thickData(data,csvData,add_n_times):
    layer_thick=data.variables['mixed_layer_thickness']
    layer_thick=layer_thick[0,:,:]
    for i in range(add_n_times):
        layer_thick=two_middle(layer_thick)
    layer_thick=np.transpose(layer_thick)
    layer_thick=layer_thick.flatten()
    layer_thick=[999 if ele==0 else ele for ele in layer_thick]
    csvData = np.column_stack((csvData, layer_thick))
    return csvData
def processmixthick(dataList, outputFileName, addtimes,filterRange=None,):
    if filterRange == None:
        filterRange = [-999, 999, -999, 999]
        filterRange.append(filterRange)
    for i in range(len(dataList)):
        print('Read data from', dataList[i])
        data = netCDF4.Dataset(dataList[i])
        time = num2dateTime(data)
        if i == 0:
            csvData = lonlatData(data,add_n_times=addtimes)
            header = ['LON', 'LAT']
        header = header[:] + time.split()
        csvData = mixed_layer_thickData(data, csvData,add_n_times=addtimes)
    outputCSV(outputFileName, csvData, header, filterRange)
    print("Save mixed_layer_thickness data (in meters) to CSV file")
    return csvData


#%%
#資料
    filename3D_day = fileNameList("12_3z.nc")
    filename2D_day = fileNameList("12_2d.nc")
    ##每月平均資料
    filename3D_month = fileNameList("_3z.timmean.nc")

######## Process 海洋漁業_蟳蟹 Data ###########
# 北緯25-26度，東經120.40-121.40度
filename3D_day = fileNameList("12_3z.nc")
filterRange_crab = [25, 26, 120.4, 121.4]
print("處理海洋漁業_蟳蟹_海底水溫")
processTemp(dataList=filename3D_day, outputFileName='csvOutput/蟳蟹海底水溫.csv', layerLoc=-1,filterRange=filterRange_crab,add_times=1)

##蟳蟹海底洋流流向
print("處理海洋漁業_蟳蟹_海底洋流流向")
filterRange_crab = [25, 26, 120.4, 121.4]
filename3D_month = fileNameList("_3z.timmean.nc")
processUVDir(dataList=filename3D_month, outputFileName='csvOutput/蟳蟹海底洋流流向.csv',layerLoc=-1,filterRange=filterRange_crab,add_times=1)


##蟳蟹海底洋流流速
print("處理海洋漁業_蟳蟹_海底洋流流速")
filterRange_crab = [25, 26, 120.4, 121.4]
filename3D_month = fileNameList("_3z.timmean.nc")
processUVvel(dataList=filename3D_month, outputFileName='csvOutput/蟳蟹海底洋流流速.csv',layerLoc=-1,filterRange=filterRange_crab,addtimes=1)


######## Process 海洋漁業_烏魚 Data ###########
print("處理海洋漁業_烏魚_海表水溫")
filename3D_day = fileNameList("12_3z.nc")
Black_fish=[22,28,119,123]
processTemp(dataList=filename3D_day, outputFileName='csvOutput/烏魚_海表水溫.csv', layerLoc=0,filterRange=Black_fish, add_times=1)

print("處理海洋漁業_烏魚_海面高度")
filename2D_day = fileNameList("12_2d.nc")
filterRange = [22, 28, 119, 123]
processSSH(dataList=filename2D_day, outputFileName='csvOutput/烏魚海面高度.csv', filterRange=filterRange,addtimes=1)

print("處理海洋漁業_烏魚_海面流向")
filename3D_day = fileNameList("12_3z.nc")
Black_fish=[22,28,119,123]
processUVDir(dataList=filename3D_day, outputFileName='csvOutput/烏魚海面流向.csv',layerLoc=0,filterRange=Black_fish,add_times=1)

print("處理海洋漁業_烏魚_海面流速")
filename3D_day = fileNameList("12_3z.nc")
Black_fish=[22,28,119,123]
processUVvel(dataList=filename3D_day, outputFileName='csvOutput/烏魚海面流速.csv',layerLoc=0,filterRange=Black_fish,addtimes=1)

####### Process 海洋漁業_鎖管 Data ##########
# 北緯24-30度，東經119-127度
filename3D_month = fileNameList("_3z.timmean.nc")
sq_filterRange=[24, 30, 119, 127]
print("處理海洋漁業_鎖管_海表水溫")
processTemp(dataList=filename3D_month,outputFileName='csvOutput/鎖管海表水溫.csv',layerLoc=0,filterRange=sq_filterRange,add_times=1)

# ######### Process 海洋漁業_櫻花蝦 Data ##########
# 北緯20-28度，東經118-122度
print("處理海洋漁業_櫻花蝦_海面高度")
filename2D_day = fileNameList("12_2d.nc")
shrimp_range=[20,28,118,122]
processSSH(dataList=filename2D_day, outputFileName='csvOutput/櫻花蝦海面高度.csv', filterRange=shrimp_range,addtimes=1)

print("處理海洋漁業_櫻花蝦_混合層深")
filename2D_day = fileNameList("12_2d.nc")
shrimp_range = [20, 28, 118, 122]
processmixthick(dataList=filename2D_day, outputFileName='csvOutput/櫻花蝦混合層深.csv', filterRange=shrimp_range,addtimes=1)
