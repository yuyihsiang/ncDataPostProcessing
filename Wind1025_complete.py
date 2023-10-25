import pandas as pd
import netCDF4  # 讀取分析.nc檔
import numpy as np  # 轉換矩陣維度
import csv  # 檔案格式
import os
import math
from datetime import datetime
import warnings



def num2dateTime(data):
    times = data.variables['time']
    dates = netCDF4.num2date(times[:], times.units)
    dateAll = datetime.strptime(str(dates[0]), '%Y-%m-%d %H:%M:%S')
    time = str(dateAll.year) + str(dateAll.month).zfill(2) + str(dateAll.day).zfill(2)
    return time
def lonlatData(data):
    lats, lons = data.variables['latitude'][:], data.variables['longitude'][:]
    csvData = np.empty(shape=(len(lons)*len(lats),2), dtype=object)
    k = 0
    for i in range(len(lons)):
        for j in range(len(lats)):
            csvData[k, :] = np.array([lons[i], lats[j]])
            k+=1
    return csvData


def fileNameList(filenameBegin):
# list the file names in the folder

    filename = np.empty(shape=(1, 1), dtype=object)
    root = "./"
    for path, subdirs, files in os.walk(root):
        for x in files:
            if x.startswith(filenameBegin):
                filename = np.column_stack((filename, os.path.join(path, x)))

    filename = filename.flatten()
    filename = sorted(filename[1:len(filename)])
    return filename

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

def uvDirData(data, csvData):
    u, v = data.variables['UGRD_10maboveground'][:], data.variables['VGRD_10maboveground'][:]
    u = np.mean(u, axis=0)
    v = np.mean(v, axis=0)
    currentDir = np.degrees(np.arctan2(v[:,:], u[:,:]))
    currentDir = np.where(currentDir < 0, currentDir + 360, currentDir)
    currentDir=currentDir.flatten()
    csvData = np.column_stack((csvData, currentDir))
    return csvData
def processUVDir(dataList, outputFileName, filterRange):
    warnings.filterwarnings('ignore')
    if filterRange == None:
        filterRange = [-999, 999, -999, 999]
        filterRange.append(filterRange)

    for i in range(len(dataList)):
        print('Read data from', dataList[i])
        data = netCDF4.Dataset(dataList[i])
        time = num2dateTime(data)
        if i == 0:
            csvData = lonlatData(data)
            header = ['LON', 'LAT']
        header = header[:] + time.split()
        csvData = uvDirData(data, csvData)
    outputCSV(outputFileName, csvData, header, filterRange)
    print("Save current direction data (in degrees) to CSV file")
def windspeedData(data, csvData):
    u, v = data.variables['UGRD_10maboveground'][:], data.variables['VGRD_10maboveground'][:]
    u=np.mean(u,axis=0)
    v=np.mean(v,axis=0)
    vel = np.sqrt(np.square(u) + np.square(v))
    vel = np.where(vel < 0, vel + 360, vel)
    vel=vel.flatten()
    csvData=np.column_stack((csvData,vel))
    return csvData
def processwindspeed(dataList, outputFileName, filterRange):
    warnings.filterwarnings('ignore')
    if filterRange == None:
        filterRange = [-999, 999, -999, 999]
        filterRange.append(filterRange)

    for i in range(len(dataList)):
        print('Read data from', dataList[i])
        data = netCDF4.Dataset(dataList[i])
        time = num2dateTime(data)
        if i == 0:
            csvData = lonlatData(data)
            header = ['LON', 'LAT']
        header = header[:] + time.split()
        csvData = windspeedData(data, csvData)
    outputCSV(outputFileName, csvData, header, filterRange)
    print("Save current direction data (in degrees) to CSV file")


def main():
    readfrom=fileNameList('wnd10m')
    Black_fish = [22, 28, 119, 123]
    shrimp_range = [20, 28, 118, 122]

    #處理烏魚風向資料
    processUVDir(readfrom,outputFileName='csvOutput/烏魚風向資料.csv',filterRange=Black_fish)
    # 處理烏魚風速資料
    processwindspeed(readfrom, outputFileName='csvOutput/烏魚風速資料.csv', filterRange=Black_fish)

    #處理櫻花蝦風向資料
    processUVDir(readfrom, outputFileName='csvOutput/櫻花蝦風向資料.csv', filterRange=shrimp_range)
    #處理櫻花蝦風速資料
    processwindspeed(readfrom, outputFileName='csvOutput/櫻花蝦風速資料.csv', filterRange=shrimp_range)

if __name__ == '__main__':
    main()