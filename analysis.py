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
import shutil
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
    print('Collecting ERCnet files to generate randomized version...')
    # Apply function to create a list of all files
    filesToRand = collectFiles(dataFolder, '.tsv')
    # Print the list
    print('Files collected:')
    print(filesToRand)

    for i in range(1, reps + 1):
        print('----------------------------------------------------------')
        print('Starting generation process for randomized set ' + str(i))
        randFolderPath = randomSetsFolder + '/randomizedFilesSet' + str(i)
        print('Making directory')
        if not os.path.exists(randFolderPath):
            os.makedirs(randFolderPath)
        
        print("Randomizing ERCnet data and writing in directory")
        for file in filesToRand:
            randomized = randomizeBCol(dataFolder, file)
            newFilePath = randFolderPath + '/' + file.replace('.tsv', '') + '_rand.tsv'
            randomized.to_csv(newFilePath, sep='\t', index=False)
        print('Finished randomization of randomized set ' + str(i))

def randomSetsOverlap(randomSets, writeLocation):
    print('\tCollect presence-absence counts from each randomized set... ', flush=True, end=' ')
    files = glob.glob(randomSets + '/randomized*_presence_absence_counts.tsv')
    #print('Randomized Files: ' + str(files))
    print('Successful')

    print('\tMake dataframes from each count file... ', flush=True, end=' ')
    dfs = [pandas.read_csv(f, sep="\t") for f in files]
    #print(dfs)
    print('Complete')

    print('\tCreate and populate summary dataframe... ', flush=True, end=' ')
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
    print('Complete')
    
    print('\tInclude summary stats avg and stDev... ', flush=True, end=' ')
    numericColumns = overlapContents.loc[:, ~overlapContents.columns.isin(['File_Name'])]
    avg = numericColumns.mean().to_list()
    overlapContents.loc[len(overlapContents.index)] = ['Average'] + avg
    stDev = numericColumns.std().to_list()
    overlapContents.loc[len(overlapContents.index)] = ['Standard Deviation'] + stDev
    print('Complete')

    print('Write completed summary dataframe... ', flush=True, end=' ')
    overlapContents.to_csv(writeLocation + '/randSetSummary.tsv', sep='\t')
    print('Complete')

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

def writePresAbTable(presenceTable, pres_abFolderPath, folderName):
    presenceTablePath = pres_abFolderPath + '/' + folderName +'_presence_absence_counts.tsv' 
    print('Beginning presence-absence table write...', flush=True, end=' ')
    presenceTable.to_csv(presenceTablePath, sep='\t', index=False)
    print('Complete')
    return presenceTablePath

def generatePresAbTable(workingFolderPath, writeLocation):
    #print(workingFolderPath)
    
    print('\tCollect source files...', flush=True, end=' ')
    files = collectFiles(workingFolderPath)    
    #print(files)
    print('Successful')
    
    print('\tSort and standardize gene pairs...', flush=True, end=' ')
    ## Use the standardizeAndConcat method to create new files containing paired gene IDs     
    for file in files:
        standardizeAndConcat(workingFolderPath, writeLocation, file)
    print('Complete')        

    ## Create a list of all sorted file containing paired gene IDs
    # Collect sorted files from folder path
    print('\tCollect sorted version of source files...', flush=True, end=' ')
    sortedFiles = collectFiles(writeLocation, '_sorted.tsv')
    #print(sortedFiles)
    print('Successful')
    
    ## Generate a list of sets 
    ## where the inner sets contain the gene pairs from each sorted file 
    ## and the outerlists holds the lists of gene pairs in each sorted file
    print('\tGenerate list of sets...', flush=True, end=' ')
    pairsByFileOut = genePairsByFile(writeLocation, sortedFiles)
    print('Complete')           
    
    ## Generate master list of all possible gene pairs
    #Create empty master list
    print('\tGenerate master list...', flush=True, end=' ')
    masterList = generateMasterList(pairsByFileOut)
    print('Complete')        
    
    ## Remove duplicates from the master list
    #Cast to set to remove duplciates and cast back to a list <-- re-evaluate for optimization?
    print('\tCast master list to set...', flush=True, end=' ')
    masterUnique = list(set(masterList))
    print('Complete')
    
    ## Produce a table the indicates presence or absence of gene pair in a file and sums the total occurrences
    # Create new dataframe with columns for gene pair, each file in sorted files, and total occurrences
    print('\tSearch for gene pairs in the master set, one by one...', flush=True, end=' ')
    presenceTable = pandas.DataFrame(columns=['Gene_Pair'] + sortedFiles + ['Total'])
    for genePair in masterUnique:
        # Determine what files contain the gene pair of interest
        newRow = [genePair] + presence_absence(genePair, pairsByFileOut)
        # Insert new row into dataframe
        presenceTable.loc[len(presenceTable)] = newRow
    print('Complete')
    print('Presence-absence table generation complete. Ready to write...')        
    return presenceTable

def writeDataValCountFile(dataPresAbPath, pres_abFolderPath):
    print('Fetching ERCnet data presence-absence table and summarizing count gene pair overlap... ', flush=True, end=' ')
    
    dataframe = pandas.read_csv(dataPresAbPath, sep='\t')
    valueCounts = dataframe['Total'].value_counts()
    print('Complete')
    
    print('Write collected data to file... ', flush=True, end=' ')

    newFolderPath = pres_abFolderPath + '/ERCnetData_valueCounts.tsv'
    valueCounts.to_csv(newFolderPath, sep='\t')
    
    print('Complete')

def writeOverlapFile(presenceTablePath, summaryFolderPath):
    print('Gather gene pairs in ERCnet data presence-absence table which are present in more than one ERCnet file... ', flush=True, end=' ')

    fullPresenceTable = pandas.read_csv(presenceTablePath, sep='\t')
    genePairOverlap = fullPresenceTable[(fullPresenceTable['Total'] > 1)]
    genePairOverlap.drop(columns=['Total'], axis = 1)
    
    print('Complete')
    
    print('Write collected data to file... ', flush=True, end=' ')

    genePairOverlap.to_csv(summaryFolderPath + '/overlap.tsv', sep='\t', index=False)

    print('Complete')

def main(masterFolder, reps):
    print('###################################')
    print('# Starting Analysis of ERCnet Run #')
    print('###################################')
    print()
    print('##########################################################')
    print('Begin by making directories for files...')
    print('##########################################################')

    #####OUTPUT DIRECTORY#####
    ##Establish directory for analysis output
    outFolder = masterFolder + '/OUT_ERCnet_Analysis'

    if os.path.exists(outFolder):
        print('Existing output folder found.')
        print('Deleting output folder to start fresh.')
        shutil.rmtree(outFolder)
        print("Directory deleted successfully.")
        print('----------------------------------------------------------')

    print('Output files will be written in: '  + outFolder)
    print('==========================================================')

    ##Check for presence-absence directory in output directory, supposed to generate if existing output directory was wiped 
    pres_abFolderPath = outFolder + '/presence_absence_files'
    
    if not os.path.exists(pres_abFolderPath):
        print('Making directory for writing presence-absence files...')
        os.makedirs(pres_abFolderPath)
        print('Successful')
    else:
        print('Directory for presence-absence files found. That\'s weird...')

    print('==========================================================')

    ##Check for summary directory in output directory, supposed to generate if existing output directory was wiped 
    summaryFolderPath = outFolder + '/summary_files'
    
    if not os.path.exists(summaryFolderPath):
        print('Making directory for writing summary files...')
        os.makedirs(summaryFolderPath)
        print('Successful')
    else:
        print('Directory for summary files found. That\'s weird...')

    print('==========================================================')

    ##Check for directory to write sorted version of ERCnet data, generate 
    sortedFolderPath = outFolder + '/sortedERCnetFiles'
    
    if not os.path.exists(sortedFolderPath):
        print('Making directory for writing ERCnet sorted files...')
        os.makedirs(sortedFolderPath)
        print('Successful')
    else:
        print('Directory for summary files found. That\'s weird...')

    print('==========================================================')

    #####RANDOM REPLICATES#####
    ##Check for randomSets directory in output directory, supposed to generate if existing output directory was wiped 
    randomFolderPath = outFolder + '/randomSets'
    
    if not os.path.exists(randomFolderPath):
        print('Making directory for writing random replicates...')
        os.makedirs(randomFolderPath)
        print('Successful')
    else:
        print('Directory for random replicates files found. That\'s weird...')
    
    print()
    print('##########################################################')
    print('Next, generate and analyze random replicates, if applicable')
    print('##########################################################')
    ##Generate random replicates of ERCnet data based on -r argument
    if reps != 0:
        print('-r argument was not 0. This analysis will include randomized replicates of ERCnet data.')
        print(str(reps) + ' replicates will be generated per -r argument.')
        print('----------------------------------------------------------')
        generateRandomizedFiles(masterFolder, randomFolderPath, reps)
        print('----------------------------------------------------------')
        print('Finished generating replicate sets.')
        print('==========================================================')
        print()

        print('Starting set analysis and generation of presence-absence table of each randomized replicates set')
        ##Collect directories for each randomized replicate
        randSetFolderNames = collectFiles(randomFolderPath, '')
        print('The following directories will be analyzed: ' + str(randSetFolderNames))
        print()
    
        ##Analyze random replicates of ERCnet data based and generate a presence-absence table for each set
        for folderName in randSetFolderNames:
            print('..........................................')
            print('Starting analysis of ' + folderName)
            print('..........................................')
            print('Generating presence-absence table')
            randPresAbTable = generatePresAbTable(randomFolderPath + '/' + folderName, randomFolderPath + '/' + folderName) 
            writePresAbTable(randPresAbTable, pres_abFolderPath, folderName)
            print()
        print('----------------------------------------------------------')

        print('Yay! Presence-absence tables for all randomized sets have been generated!')
        print('==========================================================')
        print()
        print('Producing summary file providing the gene overlap within randomized sets, prop 2+ overlap, and summary stats')
        randomSetsOverlap(pres_abFolderPath, summaryFolderPath)

    if reps == 0:
        print('-r argument was set to zero.') 
        print('Randomized replicates of ERCnet data will not be generated or analyzed.')

    print()
    print('##########################################################')
    print('Analyze real ERCnet data')
    print('##########################################################')

    dataPresenceTable = generatePresAbTable(masterFolder, sortedFolderPath)
    print()
    dataPresenceTablePath = writePresAbTable(dataPresenceTable, pres_abFolderPath, 'ERC_data')

    print()
    print('##########################################################')
    print('Finally, write summary files')
    print('##########################################################')
    
    writeDataValCountFile(dataPresenceTablePath, summaryFolderPath)
    print()
    writeOverlapFile(dataPresenceTablePath, summaryFolderPath)
    print()
    print('##########################################################')
    print('Full analysis complete!')
if __name__ == '__main__':
    main(ercNetOutputFiles, numberOfRandom)