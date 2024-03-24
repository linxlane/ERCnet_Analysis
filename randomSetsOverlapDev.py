#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Sat Mar 23 2024

@author: linlane
"""

import os
import pandas
import csv
import glob
import numpy as np
import time
import argparse

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
  print(overlapContents)

  numericColumns = overlapContents.loc[:, ~overlapContents.columns.isin(['File_Name'])]
  avg = numericColumns.mean().to_list()
  overlapContents.loc[len(overlapContents.index)] = ['Average'] + avg
  overlapContents.to_csv(writeLocation, sep='\t')

randomSetsOverlap('/Users/linlane/Desktop/randSets_test', '/Users/linlane/Desktop/overlapcontents.tsv')
