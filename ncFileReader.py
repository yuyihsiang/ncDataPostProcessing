import pandas as pd
import netCDF4
import numpy as np
import csv
import os
from datetime import datetime

def fileNameList(filenameEnding):
# list the file names in the folder
    filename = np.empty(shape=(1, 1), dtype=object)
    for x in os.listdir():  
        if x.endswith(filenameEnding):
            filename = np.column_stack((filename, x))

    filename = filename.flatten()
    filename = sorted(filename[1:len(filename)])

    return filename



def num2dateTime(data):
    times = data.variables['MT']
    dates = netCDF4.num2date(times[:],times.units)
    
    dateAll=datetime.strptime(str(dates[0]),'%Y-%m-%d %H:%M:%S')
    time = str(dateAll.year)+ str(dateAll.month).zfill(2) +str(dateAll.day).zfill(2)

    return time    


def lonlatData(data):
    # get all variable names
    # print(data.variables.keys())

    lats, lons = data.variables['Latitude'][:], data.variables['Longitude'][:]

    csvData = np.empty(shape=(len(lons)*len(lats), 2), dtype=object)
    k = 0
    for i in range(len(lons)):
        for j in range(len(lats)):
            csvData[k, :] = np.array([lons[i], lats[j]])
            k = k+1

    return csvData


def tempData(data, csvData,layerNum):

    pot_temp = data.variables['pot_temp'][:]  # temperature variable
#    print(data.variables['pot_temp'])    

    temperature = np.array(pot_temp[0, layerNum, :, :])
    temperature = temperature.flatten()
    temperature = [999.99 if ele > 100 else ele for ele in temperature]

    csvData = np.column_stack((csvData, temperature))

    return csvData


def processTemp(dataList,outputFileName,layerLoc,filterRange=None):

    if filterRange == None:
        filterRange = [-999,999,-999,999]
        filterRange.append(filterRange)
                
    for i in range(len(dataList)):
        print('Read data from', dataList[i])
        data = netCDF4.Dataset(dataList[i])       
        time=num2dateTime(data)    
        layerNum=depthLayerIndex(data,layerLoc)
        
        if i == 0:
            csvData = lonlatData(data)
            header = ['LON', 'LAT']
        
        header = header[:] + time.split()
        
        csvData = tempData(data, csvData, layerNum)
        
    print("Save temperature data (in degrees) to CSV file")
    outputCSV(outputFileName, csvData, header,filterRange)
    
    return csvData


def uvDirData(data, csvData,layerNum):    
    u, v = data.variables['u'][:], data.variables['v'][:]

    # print(data.variables['u'])    
    # print(data.variables['v'])    
        
    currentDir = np.arctan2(v[0,layerNum,:,:],u[0,layerNum,:,:])

    currentDir = np.array(currentDir[:, :])
    currentDir = currentDir.flatten()    
    u = np.array(u[0, layerNum, :, :])
    u = u.flatten()    
    currentDir = [999.99 if ele > 100 else ele for ele in u]
    currentDir = 180/np.pi*np.array(currentDir[:])

    csvData = np.column_stack((csvData, currentDir))

    return csvData


def processUVDir(dataList,outputFileName,layerLoc,filterRange=None):

    if filterRange == None:
        filterRange = [-999,999,-999,999]
        filterRange.append(filterRange)
                
    for i in range(len(dataList)):
        print('Read data from', dataList[i])
        data = netCDF4.Dataset(dataList[i])
        time=num2dateTime(data)    
        layerNum=depthLayerIndex(data,layerLoc)

        if i == 0:
            csvData = lonlatData(data)
            header = ['LON', 'LAT']
        
        header = header[:] + time.split()
        
        csvData = uvDirData(data, csvData,layerNum)
        
    print("Save current direction data (in degrees) to CSV file")
    outputCSV(outputFileName, csvData, header,filterRange)
    
    return csvData


def sshData(data, csvData):

    ssh = data.variables['ssh'][:]  # temperature variable
#    print(data.variables['ssh'])

    ssh = np.array(ssh[0, :, :])
    ssh = ssh.flatten()
    ssh = [999.99 if ele > 100 else ele for ele in ssh]

    csvData = np.column_stack((csvData, ssh))

    return csvData


def processSSH(dataList,outputFileName,filterRange=None):

    if filterRange == None:
        filterRange = [-999,999,-999,999]
        filterRange.append(filterRange)
    
    for i in range(len(dataList)):
        print('Read data from', dataList[i])
        data = netCDF4.Dataset(dataList[i])
        time=num2dateTime(data)    

        if i == 0:
            csvData = lonlatData(data)
            header = ['LON', 'LAT']
        
        header = header[:] + time.split()
        
        csvData = sshData(data, csvData)
        
    print("Save surface elevation data (in meters) to CSV file")
    outputCSV(outputFileName, csvData, header,filterRange)
    
    return csvData


def depthLayerIndex(data,layer):
    depth = data.variables['Depth'][:]
    if layer == 0:
        layerNum = 0
    elif layer == 1:
        layerNum = len(depth[:])-1
    else:
        print('Wrong water depth indexing')

    return layerNum        


def outputCSV(filename,csvData,header,filterRange):    
    
    # Create the dataframe
    df = pd.DataFrame(csvData, columns=header)

    # Apply the function
    df = df.fillna(0)

    rslt_df = df.loc[(df['LAT'] >= filterRange[0]) & 
                     (df['LAT'] <= filterRange[1]) &
                     (df['LON'] >= filterRange[2]) & 
                     (df['LON'] <= filterRange[3]) ]

    # print the DataFrame
    rslt_df.to_csv(filename, index=False)
    
    return

# read nc file using netCDF4
filename3D = fileNameList("3z.timmean.nc")
filename2D = fileNameList("2d.timmean.nc")

# ######## Process 海洋漁業_蟳蟹 Data ###########
# 北緯25-26度，東經120.40-121.40度
filterRange=[25, 26, 120.4, 121.4]

print("處理海洋漁業_蟳蟹_海表水溫")
processTemp(dataList=filename3D,outputFileName='蟳蟹海表水溫.csv',layerLoc=0,filterRange=filterRange)

print("處理海洋漁業_蟳蟹_海底洋流流向")
processUVDir(dataList=filename3D,outputFileName='蟳蟹海表洋流流向.csv',layerLoc=0,filterRange=filterRange)
# #############################################


######## Process 海洋漁業_烏魚 Data ###########
# 北緯22-28度，東經119-123度
filterRange=[22, 28, 119, 123]

print("處理海洋漁業_烏魚_海表水溫")
processTemp(dataList=filename3D,outputFileName='烏魚海表水溫.csv',layerLoc=0,filterRange=filterRange)

print("處理海洋漁業_烏魚_海底洋流流向")
processUVDir(dataList=filename3D,outputFileName='烏魚海底洋流流向.csv',layerLoc=1,filterRange=filterRange)

print("處理海洋漁業_烏魚_海面高度")
processSSH(dataList=filename2D,outputFileName='烏魚海面高度.csv',filterRange=filterRange)

#############################################


# ######### Process 海洋漁業_鎖管 Data ##########
# 北緯24-30度，東經119-127度
filterRange=[24, 30, 119, 127]

print("處理海洋漁業_鎖管_海表水溫")
processTemp(dataList=filename3D,outputFileName='鎖管海表水溫.csv',layerLoc=0,filterRange=filterRange)
#############################################
