#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 2024

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
parser = argparse.ArgumentParser(description='Script for generating plots only to analyze ERCnet run')

parser.add_argument('-p', '--path', type=str, required=True, help='Full path to folder containing value count files.') 
parser.add_argument('-w', '--write', type=str, required=True, help='Full path to folder where generated plot files will be written.') 
parser.add_argument('-t', '--totalsBarplot', action='store_true', required=False, help='Include this flag to generate a totals barplot with ERCnet output.')
parser.add_argument('-s', '--sbsBarplot', action='store_true', required=False, help='Include this flag to generate a side-by-side barplot with ERCnet output and random replicates.')
parser.add_argument('-k', '--kde', action='store_true', required=False, help='Include this flag to generate a Kernal Density Plot from random replicates.')
parser.add_argument('-u', '--upsetPlot', action='store_true', required=False, help='Include this flag to generate an upset plot.')

#Define the parser
args = parser.parse_args()

#Store arguments
valueCountFiles=args.path
writeLocation=args.write
barplotArg=args.totalsBarplot
sbsArg=args.sbsBarplot
kdeArg=args.kde
upsetArg=args.upsetPlot

def collectFiles(path, endStr = '.tsv'):
    files = []
    for file in os.listdir(path):
        if file.endswith(endStr):
            files.append(file)
    files.sort()
    return files

def totalsBarPlot(ERCnetData, writePath):
    xAxis = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
    mpl.rcParams['pdf.fonttype'] = 42
    plt.Figure(figsize=(3,2))
    plt.bar(ercDataDF['Total'], ercDataDF['count'], bottom = 0.5)
    plt.title('Total Column Value Counts')
    plt.yscale('log')
    plt.xticks(xAxis)
    plt.xlabel('Gene Pair Overlap')
    plt.ylabel('Count (log)')
    plt.savefig(writePath + '/dataBarplot.pdf', format = 'pdf', transparent = True) 
    plt.close()

def sidebysideBarplot(ercBarData, avgRandData, writePath):
    barWidth = 0.35
    xAxis = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
    mpl.rcParams['pdf.fonttype'] = 42
    plt.Figure(figsize=(16,12))
    plt.bar(ercBarData['Total'], ercBarData['count'], bottom = 0.5, color='blue', width=-barWidth, align='edge', edgecolor ='black', label ='ERCnet Data')
    plt.bar(avgRandData['Total'], avgRandData['avg_counts'], bottom=0.5, yerr=avgRandData['st_dev'], color='grey', width=barWidth, align='edge', edgecolor ='black', label ='Random')
    plt.title('ERC Data vs Rand Null Hyp')
    plt.yscale('log')
    plt.xticks(xAxis)
    plt.xlabel('Gene Overlap')
    plt.ylabel('Count (log)')
    plt.legend(loc='upper right')
    plt.savefig(writePath + '/sidebysideBarplot.pdf', format = 'pdf', transparent = True) 
    plt.close()

def kde(avgRandDataDF, dataPresAbPath, masterPath):
    repData = pandas.read_csv(avgRandDataDF, sep='\t')
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
    repData.to_csv(avgRandDataDF, sep='\t')
     
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
    

start_time = time.time()

argsList = [barplotArg, sbsArg, kdeArg, upsetArg]

if all(f == False for f in argsList):
    print("No flags passed. No plots will be generated.")

else:
    filesList = collectFiles(valueCountFiles)

    ercDataDF = pandas.read_csv(valueCountFiles + '/value_counts.tsv', sep='\t')
    avgRandDataDF = pandas.read_csv(valueCountFiles + '/avgRandValueCounts.tsv', sep='\t', index_col=0)
    print('ERCnet data totals')
    print(ercDataDF)
    print('Average totals from random replicates')
    print(avgRandDataDF)

    if barplotArg == True:
        totalsBarPlot(ercDataDF, writeLocation)
        print('Barplot showing total counts from ERCnet output sucessfully generated!')

    if sbsArg == True:
        sidebysideBarplot(ercDataDF, avgRandDataDF, writeLocation)
        print('Barplot showing total counts from ERCnet output and average random replicate counts sucessfully generated!')

    if kdeArg == True:
        kde(avgRandDataDF, dataValCounts, masterFolder)

    if upsetArg == True:
        upsetPlot(ercDataDF, masterFolder)

print("Program took", time.strftime("%Hh%Mm%Ss", time.gmtime(time.time() - start_time)), "to run")
