#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 2024

@author: linlane
"""
import os
import pandas
import csv
import glob
import numpy as np
import time
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

def randomSetsOverlap(randomSets, writeLocation):
  files = glob.glob(randomSets + '/randomized*_presence_absence_counts.tsv')
  #print('Randomized Files: ' + str(files))

  dfs = [pandas.read_csv(f, sep="\t") for f in files]
  #print(dfs)

  overlapContents = pandas.DataFrame(columns=['File_Name', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12'])
  #print(overlapContents)

  for i in range(len(dfs)):
    newRow = [files[i]]
    data = dfs[i]['Total'].value_counts().to_list()
    newRow += data
    newRow += [0]*(13 - len(newRow))
    #print(newRow)
    overlapContents.loc[len(overlapContents.index)] = newRow

  overlapContents['Proportion_2+_Overlap'] = np.nan
  propList = []
  
  for row in overlapContents.index:
      numeratorList = list(overlapContents.iloc[row, 2:13])
      numerator = int(np.nansum(numeratorList))
      oneColumn = int(overlapContents.iloc[row, 1])
      denominator = numerator + oneColumn
      prop = numerator / denominator
      propList.append(prop)    
  overlapContents['Proportion_2+_Overlap'] = propList
  
  numericColumns = overlapContents.loc[:, ~overlapContents.columns.isin(['File_Name'])]
  avg = numericColumns.mean().to_list()
  overlapContents.loc[len(overlapContents.index)] = ['Average'] + avg
  stDev = numericColumns.std().to_list()
  overlapContents.loc[len(overlapContents.index)] = ['Standard Deviation'] + stDev
  overlapContents.to_csv(writeLocation + '/randSetSummary.tsv', sep='\t')

def standardizeAndConcat(folderPath, writeLocationPath, fileName):
    contents = pandas.read_csv(folderPath + '/' + fileName, sep='\t', usecols=['GeneA_ID','GeneB_ID'])
    removeNa = contents.dropna()
    castStringContents = removeNa.astype(str)
    standardizedIDs = castStringContents.apply(sorted, axis=1)
    
    pairedIDs = []
    for row in standardizedIDs:
        pairedIDs.append(str(row[0]) + '||' + str(row[1]))
    
    newFilePath = writeLocationPath + '/' + fileName.replace('.tsv', '') + '_sorted.tsv'
    
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

def generatePresAbTable(workingFolderPath, writeLocation):
    print(workingFolderPath)
    
    files = collectFiles(workingFolderPath)    
    print(files)
    
    ## Use the standardizeAndConcat method to create new files containing paired gene IDs     
    for file in files:
        standardizeAndConcat(workingFolderPath, writeLocation, file)
        
    ## Create a list of all sorted file containing paired gene IDs
    # Collect sorted files from folder path
    sortedFiles = collectFiles(writeLocation, '_sorted.tsv')
    print(sortedFiles)
    
    ## Generate a list of sets 
    ## where the inner sets contain the gene pairs from each sorted file 
    ## and the outerlists holds the lists of gene pairs in each sorted file
    print('Generate list of sets')
    pairsByFileOut = genePairsByFile(writeLocation, sortedFiles)   
    
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
    newFolderPath = pres_abFolderPath + '/ERCnetData_valueCounts.tsv'
    valueCounts.to_csv(newFolderPath, sep='\t')
    return newFolderPath

def writeOverlapFile(presenceTablePath, summaryFolderPath):
    fullPresenceTable = pandas.read_csv(presenceTablePath, sep='\t')
    genePairOverlap = fullPresenceTable[(fullPresenceTable['Total'] > 1)]
    genePairOverlap.drop(columns=['Total'], axis = 1)
    genePairOverlap.to_csv(summaryFolderPath + '/overlap.tsv', sep='\t', index=False)

def main(masterFolder, reps):
    start_time = time.time()

    outFolder = masterFolder + '/OUT_ERCnet_Analysis'
    
    pres_abFolderPath = outFolder + '/presence_absence_files'
    
    if not os.path.exists(pres_abFolderPath):
        os.makedirs(pres_abFolderPath)

    summaryFolderPath = outFolder + '/summary_files'
    
    if not os.path.exists(summaryFolderPath):
        os.makedirs(summaryFolderPath)
        
  ##Generate summar file for random replicates
    randomFolderPath = outFolder + '/randomSets'
    
    if not os.path.exists(randomFolderPath):
        os.makedirs(randomFolderPath)
    
    generateRandomizedFiles(masterFolder, randomFolderPath, reps)
    
    randSetFolderNames = collectFiles(randomFolderPath, '')
    
    for folderName in randSetFolderNames:
        randPresAbTable = generatePresAbTable(randomFolderPath + '/' + folderName, randomFolderPath + '/' + folderName) 
        writeTable(randPresAbTable, pres_abFolderPath, folderName)

    randomSetsOverlap(pres_abFolderPath, summaryFolderPath)

    sortedFolderPath = outFolder + '/sortedERCnetFiles'
    
    if not os.path.exists(sortedFolderPath):
        os.makedirs(sortedFolderPath)

    dataPresenceTable = generatePresAbTable(masterFolder, sortedFolderPath)
    dataPresenceTablePath = writeTable(dataPresenceTable, pres_abFolderPath, 'ERC_data')
    
    writeDataValCountFile(dataPresenceTablePath, summaryFolderPath)

    writeOverlapFile(dataPresenceTablePath, summaryFolderPath)
    
    print("Program took", time.strftime("%Hh%Mm%Ss", time.gmtime(time.time() - start_time)), "to run")

main(ercNetOutputFiles, numberOfRandom)