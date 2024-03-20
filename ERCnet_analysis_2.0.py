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

def avgRandValueCounts(pres_ab_folderInput):
    avgRandFilePath = pres_ab_folderInput + '/avgRandValueCounts.tsv'
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

def totalsBarPlot(ercBarData, masterPath):
    xAxis = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
    ercDataDF = pandas.read_csv(ercBarData, sep='\t')
    ercDataValueCounts = ercDataDF['Total'].value_counts()
    mpl.rcParams['pdf.fonttype'] = 42
    plt.Figure(figsize=(3,2))
    plt.bar(ercDataValueCounts.index, ercDataValueCounts.values, bottom = 0.5)
    plt.title('Total Column Value Counts')
    plt.yscale('log')
    plt.xticks(xAxis)
    plt.xlabel('Gene Pair Overlap')
    plt.ylabel('Count (log)')
    plt.savefig(masterPath + '/dataBarplot.pdf', format = 'pdf', transparent = True) 
    plt.close()

def sidebysideBarplotAVG(ercBarData, avgRandData, masterPath):
    ercDataDF = pandas.read_csv(ercBarData, sep='\t')
    ercDataValueCounts = ercDataDF['Total'].value_counts()
    avgRandDF = pandas.read_csv(avgRandData, sep='\t')
    avgRandDF = avgRandDF[avgRandDF.columns.difference(['File_name'])]
    
    meanRandData = avgRandDF.mean()
    randStDev = list(avgRandDF.std().values)


    ercX = list(ercDataValueCounts.index)
    randX = list(meanRandData.index.astype('int64'))


    barWidth = 0.35
    xAxis = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
    mpl.rcParams['pdf.fonttype'] = 42
    plt.Figure(figsize=(16,12))
    plt.bar(ercX, ercDataValueCounts.values, bottom = 0.5, color='blue', width=-barWidth, align='edge', edgecolor ='black', label ='ERCnet Data')
    plt.bar(randX, meanRandData.values, bottom=0.5, yerr=randStDev, color='grey', width=barWidth, align='edge', edgecolor ='black', label ='Random')
    plt.title('ERC Data vs Rand Null Hyp')
    plt.yscale('log')
    plt.xticks(xAxis)
    plt.xlabel('Gene Overlap')
    plt.ylabel('Count (log)')
    plt.legend(loc='upper right')
    plt.savefig(masterPath + '/sidebysideBarplot.pdf', format = 'pdf', transparent = True) 
    plt.close()

def kde(avgRandValCountPath, dataPresAbPath, masterPath):
    repData = pandas.read_csv(avgRandValCountPath, sep='\t')
    repData['Proportion'] = np.nan
    propList = []
    
    for row in repData.index:
        numeratorList = list(repData.iloc[row, 2:13])
        numerator = int(np.nansum(numeratorList))
        oneColumn = int(repData.iloc[row, 1])
        denominator = numerator + oneColumn
        prop = numerator / denominator
        propList.append(prop)
        
    repData['Proportion'] = propList
    repData.to_csv(avgRandValCountPath, sep='\t')
     
    realDataValCounts = pandas.read_csv(dataPresAbPath, sep='\t')
    realDataValsList = realDataValCounts['count'].values.tolist()
    greaterThanOne = 0
    for ele in range(1, len(realDataValsList)):
        greaterThanOne = greaterThanOne + realDataValsList[ele]
    realDataProp = greaterThanOne/sum(realDataValsList)
    mpl.rcParams['pdf.fonttype'] = 42
    sns.kdeplot(repData['Proportion'])
    plt.axvline(x = realDataProp, color = 'red')
    plt.title('KDE')
    plt.savefig(masterPath + '/kde.pdf', format = 'pdf', transparent = True) 

    
def upsetPlot(presenceTablePath, masterPath):
    fullPresenceTable = pandas.read_csv(presenceTablePath, sep='\t')
    genePairOverlap = fullPresenceTable[(fullPresenceTable['Total'] > 1)]
    genePairOverlap.drop(columns=['Total'], axis = 1)
    genePairOverlap.to_csv(masterPath + '/overlap.tsv', sep='\t', index=False)
    upset = from_indicators(lambda genePairOverlap: genePairOverlap.select_dtypes(bool), data=genePairOverlap)
    mpl.rcParams['pdf.fonttype'] = 42
    fig = plt.figure(figsize=(40, 15))
    plot(upset, fig=fig, element_size=None, sort_categories_by='-input')
    plt.savefig(masterPath + '/upsetPlot.pdf', format = 'pdf', transparent = True) 
    plt.title('UpSet Plot')
    plt.close()
    
'''
def hypTest():    
    # Given information
    sample_mean = realDataProp
    population_mean = np.mean(propList)
    population_std = np.std(propList)
    sample_size = len(propList)
    alpha = 0.05# compute the z-score
    z_score = (sample_mean-population_mean)/(population_std/np.sqrt(len(propList)))
    print('Z-Score :',z_score)# Approach 1: Using Critical Z-Score# Critical Z-Score
    z_critical = stats.norm.ppf(1-alpha)
    print('Critical Z-Score :',z_critical)# Hypothesis
    if z_score >  z_critical:
        print("Reject Null Hypothesis")
    else:
      print("Fail to Reject Null Hypothesis")# Approach 2: Using P-value
        
    ### Report p_value as less than 1^-5!  
    # P-Value : Probability of getting less than a Z-score
    #p_value = 1-stats.norm.cdf(z_score)
    p_value = stats.norm.sf(abs(z_score))
    print('{0:.30f}'.format(p_value))
    #print('p-value :', p_value)# Hypothesis
    if p_value <  alpha:
        print("Reject Null Hypothesis")
    else:
      print("Fail to Reject Null Hypothesis")
'''

def main(masterFolder, reps):
    start_time = time.time()
    
    pres_abFolderPath = masterFolder + '/presence_absence_files'
    
    if not os.path.exists(pres_abFolderPath):
        os.makedirs(pres_abFolderPath)
        
    randomFolderPath = masterFolder + '/randomSets'
    
    if not os.path.exists(randomFolderPath):
        os.makedirs(randomFolderPath)
    
    generateRandomizedFiles(masterFolder, randomFolderPath, reps)
    
    randSetFolderNames = collectFiles(randomFolderPath, '')
    
    for folderName in randSetFolderNames:
        randPresAbTable = generatePresAbTable(randomFolderPath + '/' + folderName) 
        writeTable(randPresAbTable, pres_abFolderPath, folderName)

    avgRandValCountPath = avgRandValueCounts(pres_abFolderPath)
        
    dataPresenceTable = generatePresAbTable(masterFolder)
    dataPresenceTablePath = writeTable(dataPresenceTable, pres_abFolderPath, 'ERC_data')
    
    dataValCounts = writeDataValCountFile(dataPresenceTablePath, pres_abFolderPath)
    
    totalsBarPlot(dataPresenceTablePath, masterFolder)

    sidebysideBarplotAVG(dataPresenceTablePath, avgRandValCountPath, masterFolder)

    kde(avgRandValCountPath, dataValCounts, masterFolder)
    
    upsetPlot(dataPresenceTablePath, masterFolder)
    
    print("Program took", time.strftime("%Hh%Mm%Ss", time.gmtime(time.time() - start_time)), "to run")

main(ercNetOutputFiles, numberOfRandom)