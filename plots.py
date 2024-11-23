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
parser.add_argument('-n', '--name', type=str, required=True, help='Name to add to plots') 
parser.add_argument('-w', '--write', type=str, required=True, help='Full path to folder where generated plot files will be written.') 
parser.add_argument('-b', '--barplot', action='store_true', required=False, help='Include this flag to generate a totals barplot with ERCnet output.')
parser.add_argument('-s', '--sbsBarplot', action='store_true', required=False, help='Include this flag to generate a side-by-side barplot with ERCnet output and random replicates.')
parser.add_argument('-k', '--kde', action='store_true', required=False, help='Include this flag to generate a Kernal Density Plot from random replicates.')
parser.add_argument('-u', '--upsetPlot', action='store_true', required=False, help='Include this flag to generate an upset plot.')
parser.add_argument('-a', '--all', action='store_true', required=False, help='Include this flag to generate all plots.')

#Define the parser
args = parser.parse_args()

#Store arguments
summaryFiles=args.path
dataName = args.name
writeLocation=args.write
barplotArg=args.barplot
sbsArg=args.sbsBarplot
kdeArg=args.kde
upsetArg=args.upsetPlot
allPlots = args.all

def collectFiles(path, endStr = '.tsv'):
    files = []
    for file in os.listdir(path):
        if file.endswith(endStr):
            files.append(file)
    files.sort()
    return files

def totalsBarPlot(ERCnetData, writePath):
    xAxis = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 , 13, 14, 15, 16, 17, 18, 19, 20]
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
    xAxis = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 , 13, 14, 15, 16, 17, 18, 19, 20]
    allCounts = avgRandData.iloc[1:len(avgRandData.index)-2, 1:13]
    print(type(allCounts))
    print(allCounts)
    newRandStDev = list(allCounts.std().values)
    avgRandDataColumns = avgRandData.iloc[len(avgRandData.index)-2, 1:13]
    #print(avgRandDataColumns)
    avgRandDataX= avgRandDataColumns.index.astype(int).to_list()
    #print(avgRandDataX)
    avgRandDataY = avgRandDataColumns.astype(int).to_list()
    #print(avgRandDataY)
    randStDevData = list(avgRandData.iloc[len(avgRandData.index)-1, 1:13])
    print(randStDevData)

    mpl.rcParams['pdf.fonttype'] = 42
    plt.Figure(figsize=(16,12))
    plt.bar(ercBarData['Total'], ercBarData['count'], bottom = 0.5, color='blue', width=-barWidth, align='edge', edgecolor ='black', label ='ERCnet Data')
    plt.bar(avgRandDataX, avgRandDataY, yerr=newRandStDev, bottom=0.5, color='grey', width=barWidth, align='edge', edgecolor ='black', label ='Random')
    plt.title('ERC Data vs Rand Null Hyp')
    plt.yscale('log')
    plt.xticks(xAxis)
    plt.xlabel('Gene Overlap')
    plt.ylabel('Count (log)')
    plt.legend(loc='upper right')
    plt.savefig(writePath + '/sidebysideBarplot.pdf', format = 'pdf', transparent = True) 
    plt.close()

def kde(avgRandDataDF, ercDataDF, writeLocation):
    propListFull = avgRandDataDF['Proportion_2+_Overlap'].to_list()
    propListRandOnly = propListFull[:len(propListFull)-2]
    #print(propListRandOnly)

    dataValueCounts = ercDataDF['count'].values.tolist()
    greaterThanOne = 0
    for ele in range(1, len(dataValueCounts)):
        greaterThanOne = greaterThanOne + dataValueCounts[ele]
    realDataProp = greaterThanOne/sum(dataValueCounts)

    mpl.rcParams['pdf.fonttype'] = 42
    sns.kdeplot(propListRandOnly)
    plt.axvline(x = realDataProp, color = 'red')
    plt.xlim(left=0)
    plt.title('KDE')
    plt.savefig(writeLocation + '/KDE.pdf', format = 'pdf', transparent = True) 
    
def upsetPlot(overlapFile, writeLocation):
    genePairOverlap = pandas.read_csv(overlapFile, sep='\t')
    
    upset = from_indicators(lambda genePairOverlap: genePairOverlap.select_dtypes(bool), data=genePairOverlap)
    mpl.rcParams['pdf.fonttype'] = 42
    fig = plt.figure(figsize=(40, 15))
    plot(upset, fig=fig, element_size=None, sort_categories_by='-input')
    plt.savefig(writeLocation + '/upsetPlot.pdf', format = 'pdf', transparent = True) 
    plt.title('UpSet Plot')
    plt.close()
    

argsList = [barplotArg, sbsArg, kdeArg, upsetArg, allPlots]

if all(f == False for f in argsList):
    print("No flags passed. No plots will be generated.")

else:
    if not os.path.exists(writeLocation):
        os.makedirs(writeLocation)
    
    filesList = collectFiles(summaryFiles)

    ercDataDF = pandas.read_csv(summaryFiles + '/ERCnetData_valueCounts.tsv', sep='\t')
    avgRandDataDF = pandas.read_csv(summaryFiles + '/randSetSummary.tsv', sep='\t', index_col=0)

    #print('ERCnet data totals')
    #print(ercDataDF)
    #print('Average totals from random replicates')
    #print(avgRandDataDF)

    if barplotArg == True or allPlots == True:
        totalsBarPlot(ercDataDF, writeLocation)
        print('Barplot showing total counts from ERCnet output sucessfully generated!')

    if sbsArg == True or allPlots == True:
        sidebysideBarplot(ercDataDF, avgRandDataDF, writeLocation)
        print('Barplot showing total counts from ERCnet output and average random replicate counts sucessfully generated!')

    if kdeArg == True or allPlots == True:
        kde(avgRandDataDF, ercDataDF, writeLocation)
        print('Kernel Density Estimate Plot sucessfully generated!')

    if upsetArg == True or allPlots == True:
        upsetPlot(summaryFiles + '/overlap.tsv', writeLocation)
        print('UpSet Plot sucessfully generated!')
    print("All Done!")