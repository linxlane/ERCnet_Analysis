#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 10 19:53:19 2024

@author: linlane
"""
import os
import pandas
import csv
import matplotlib.pyplot as plt
import glob
import seaborn as sns
import numpy as np
import time
from upsetplot import plot, from_indicators
import matplotlib as mpl
import argparse

#Set up an argumanet parser
parser = argparse.ArgumentParser(description='Script for generating plots to analyze ERCnet run')

parser.add_argument('-p', '--path', type=str, metavar='', required=True, help='Full path to ERCnet output files to analyze.') 
parser.add_argument('-r', '--random', type=int, metavar='', required=True, help='Number of random replicates to generate and analyze.') 

#Define the parser
args = parser.parse_args()

#Store arguments
ercNetOutputFiles=args.path
numberOfRandom=args.random

def collectFiles(path, endStr = '.tsv'):
    files = []
    for file in os.listdir(path):
        if file.endswith(endStr):
            files.append(file)
    files.sort()
    return files

def navigateFinder(currentFolderPath, i):
    inds = [i for i, c in enumerate(currentFolderPath) if c=='/']
    sliceSpot = inds[i]
    return currentFolderPath[:sliceSpot]
    
def randomizeBCol(folderPath, file):
    ##Get contents of first file in list and create new tsv with two columns 
    colA = pandas.read_csv(folderPath + '/' + file, sep='\t', usecols=['GeneA_ID'])
    colB = pandas.read_csv(folderPath + '/' + file, sep='\t', usecols=['GeneB_ID'])

    randContents = colB.sample(frac=1)

    newIndex = randContents.reset_index()

    randColAB = pandas.concat([colA, newIndex], axis=1)

    return randColAB

def generateRandomizedFiles(dataFolder, randomSetsFolder, reps):
    for i in range(1, reps + 1):
        # Apply function to create a list of all files
        filesToRand = collectFiles(dataFolder, '.tsv')
        # Print the list
        print(filesToRand)
        
        randFolderPath = randomSetsFolder + '/randomizedFilesSet' + str(i)
        if not os.path.exists(randFolderPath):
            os.makedirs(randFolderPath)
        
        print("Randomize!")
        for file in filesToRand:
            randomized = randomizeBCol(dataFolder, file)
            newFilePath = randFolderPath + '/' + file.replace('.tsv', '') + '_rand.tsv'
            randomized.to_csv(newFilePath, sep='\t', index=False)

def randValueCounts(pres_ab_folderInput):
    avgRandFilePath = pres_ab_folderInput + '/randValueCounts.tsv'
    file1 = open(avgRandFilePath, "w")
    L = ["File_name\t", "1\t", "2\t", "3\t", "4\t", "5\t", "6\t", "7\t", "8\t", "9\t", "10\t", "11\t", "12\n",]
    file1.writelines(L)
    file1.close()
    
    files = glob.glob(pres_ab_folderInput + '/randomized*_presence_absence_counts.tsv')
    print('Randomized Files: ' + str(files))
    
    dfs = [pandas.read_csv(f, sep="\t") for f in files]
    
    i=0
    fileIndex = 0
    for table in dfs:
        i=0
        data = table['Total'].value_counts().to_frame()
     
        countsList = []
        length = len(list(data['count']))
        
        for i in range(0, length):
            countsList.append(list(data['count'])[i])
        
        delim = "\t"
 
        # using list comprehension to convert each element to string and adding delim
        res = delim.join([str(ele) for ele in countsList])
        
        # Append-adds at last
        file1 = open(avgRandFilePath, "a")  # append mode
        newLine = files[fileIndex] + "\t" + res + "\n"
        file1.writelines(newLine)
        file1.close()
        fileIndex+=1
    return avgRandFilePath

def writeAvgRandValCountsFile(dataPath, pres_ab_path):
  avgRandData = pandas.read_csv(dataPath, sep='\t')
  randDFCols = avgRandData.loc[:, ~avgRandData.columns.isin(['File_name', 'Proportion'])]
  meanRandData = randDFCols.mean()
  randStDev = list(randDFCols.std().values)
  meanData = meanRandData.to_frame(name='avg_counts').reset_index(drop=True)
  meanData = meanData.replace(np.nan, 0)
  totalCol = pandas.DataFrame({'Total': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]})
  newDF = pandas.concat([totalCol, meanData], axis=1)
  newDF.to_csv(pres_ab_path + '/avgRandValueCounts.tsv', sep='\t')


def standardizeAndConcat(folderPath, fileName):
    contents = pandas.read_csv(folderPath + '/' + fileName, sep='\t', usecols=['GeneA_ID','GeneB_ID'])
    removeNa = contents.dropna()
    castStringContents = removeNa.astype(str)
    standardizedIDs = castStringContents.apply(sorted, axis=1)
    
    pairedIDs = []
    for row in standardizedIDs:
        pairedIDs.append(str(row[0]) + '||' + str(row[1]))
    
    newFilePath = folderPath + '/' + fileName.replace('.tsv', '') + '_sorted.tsv'
    
    with open(newFilePath, 'w') as output:
        tsv_output = csv.writer(output, delimiter='\n')
        tsv_output.writerow(pairedIDs)
        
def genePairsByFile(folderPath, sortedFilesInput):
    pairsByFile = []
    # For loop to process each sorted file
    for file in sortedFilesInput:
        # Open and read the rows into a list without '\n'
        filePath = folderPath + '/' + file   
        f = open(filePath, "r")
        readNamesSet = set(f.read().splitlines())
 
        #Add list to outerlist 
        pairsByFile.append(readNamesSet)  
    return pairsByFile
            
def generateMasterList(pairsByFileInput):
    masterList = []
    #Iterate through list of lists and generate master list
    for lst in pairsByFileInput:
        masterList += lst
    return masterList

def presence_absence(inputGenePair, pairsByFile):
    presenceList = []
    for file in pairsByFile:
        if inputGenePair in file:
            presenceList.append(True)
        else:
            presenceList.append(False)
    total = sum(presenceList)
    presenceList.append(total)
    return presenceList

def writeTable(presenceTable, pres_abFolderPath, folderName):
    presenceTablePath = pres_abFolderPath + '/' + folderName +'_presence_absence_counts.tsv' 
    presenceTable.to_csv(presenceTablePath, sep='\t', index=False)
    return presenceTablePath

def generatePresAbTable(workingFolderPath):
    print(workingFolderPath)
    
    files = collectFiles(workingFolderPath)    
    print(files)
    
    ## Use the standardizeAndConcat method to create new files containing paired gene IDs     
    for file in files:
        standardizeAndConcat(workingFolderPath, file)
        
    ## Create a list of all sorted file containing paired gene IDs
    # Collect sorted files from folder path
    sortedFiles = collectFiles(workingFolderPath, '_sorted.tsv')
    print(sortedFiles)
    
    ## Generate a list of sets 
    ## where the inner sets contain the gene pairs from each sorted file 
    ## and the outerlists holds the lists of gene pairs in each sorted file
    print('Generate list of sets')
    pairsByFileOut = genePairsByFile(workingFolderPath, sortedFiles)   
    
    ## Generate master list of all possible gene pairs
    #Create empty master list
    print('Generate master list')
    masterList = generateMasterList(pairsByFileOut)
    
    ## Remove duplicates from the master list
    #Cast to set to remove duplciates and cast back to a list <-- re-evaluate for optimization?
    masterUnique = list(set(masterList))
    
    ## Produce a table the indicates presence or absence of gene pair in a file and sums the total occurrences
    # Create new dataframe with columns for gene pair, each file in sorted files, and total occurrences
    print('Generate presence table')
    presenceTable = pandas.DataFrame(columns=['Gene_Pair'] + sortedFiles + ['Total'])
    for genePair in masterUnique:
        # Determine what files contain the gene pair of interest
        newRow = [genePair] + presence_absence(genePair, pairsByFileOut)
        # Insert new row into dataframe
        presenceTable.loc[len(presenceTable)] = newRow
    return presenceTable

def writeDataValCountFile(dataPresAbPath, pres_abFolderPath):
    dataframe = pandas.read_csv(dataPresAbPath, sep='\t')
    valueCounts = dataframe['Total'].value_counts()
    print(valueCounts)
    newFolderPath = pres_abFolderPath + '/value_counts.tsv'
    valueCounts.to_csv(newFolderPath, sep='\t')
    return newFolderPath

def main(masterFolder, reps):
    start_time = time.time()
    
    pres_abFolderPath = masterFolder + '/presence_absence_files'
    
    if not os.path.exists(pres_abFolderPath):
        os.makedirs(pres_abFolderPath)

    value_countsFolderPath = masterFolder + '/value_count_files'
    
    if not os.path.exists(value_countsFolderPath):
        os.makedirs(value_countsFolderPath)
        
  ##Generate value counts file for random replicates
    randomFolderPath = masterFolder + '/randomSets'
    
    if not os.path.exists(randomFolderPath):
        os.makedirs(randomFolderPath)
    
    generateRandomizedFiles(masterFolder, randomFolderPath, reps)
    
    randSetFolderNames = collectFiles(randomFolderPath, '')
    
    for folderName in randSetFolderNames:
        randPresAbTable = generatePresAbTable(randomFolderPath + '/' + folderName) 
        writeTable(randPresAbTable, pres_abFolderPath, folderName)

    randValCountPath = randValueCounts(pres_abFolderPath)
    writeAvgRandValCountsFile(randValCountPath, value_countsFolderPath)

    dataPresenceTable = generatePresAbTable(masterFolder)
    dataPresenceTablePath = writeTable(dataPresenceTable, pres_abFolderPath, 'ERC_data')
    
    dataValCounts = writeDataValCountFile(dataPresenceTablePath, value_countsFolderPath)
    
    print("Program took", time.strftime("%Hh%Mm%Ss", time.gmtime(time.time() - start_time)), "to run")

main(ercNetOutputFiles, numberOfRandom)