#!/usr/bin/env python3

import pandas as pd
import glob
import numpy as np
from sys import argv
import os
import datetime
from pytz import timezone

print('\n**********************')
print('START TIME:', datetime.datetime.now(timezone('EST')))


# loading in files! --> We're not going to use the scores yet!! --> VERY DIFFERENT, WE'RE JUST COUNTING UP THE NUMBER SO FAR
res = pd.read_csv(str(argv[1]), sep='\t', header=None)
print(res.head())
res.columns = ['matches', 'misMatches', 'repMatches', 'nCount', 'qNumInsert', 'qBaseInsert', 'tNumInsert', 'tBaseInsert', 'strand', 'qName', 'qSize', 'qStart', 'qEnd', 'tName', 'tSize', 'tStart', 'tEnd', 'blockCount', 'blockSizes','qStarts', 'tStarts']
# sorting the blat responses
sorted_res = res.sort_values(by=['matches'], ascending=False)
sorted_res = sorted_res.drop_duplicates(subset='qName', keep="first")


# recall, no score csv --> using the 2nd argument for the svs
svs = pd.read_csv(str(argv[2]), comment='#', sep='\t')
print('PRE ANNOTATION')
print(svs.head())

# BY CHROMOSOME!!
chromo = str(argv[3])
svs = svs[svs.CHROM == chromo]
ogrow = len(svs)
print('og number of rows:', ogrow)

#adding a column for a merge..
for idx, row in svs.iterrows():
	svs.loc[idx, 'qName'] = ('POST_'+ row['CHROM']+ '_' +str(row['POS'])+ '_' + str(row['SVlen'])+ '_' + row['SV_Type'])

merged2 = pd.merge(svs, sorted_res, how = 'outer', on = 'qName', indicator = True)

for idx, row in merged2.iterrows():
	if row['_merge'] == 'left_only':
		merged2.loc[idx, 'Homology'] = False
	elif row['_merge'] == 'both':
		merged2.loc[idx, 'Homology'] = True

newrow = len(merged2)
print('new number of rows:', newrow)

# checking if the rows match (they should)
if ogrow == newrow:
	print('good merge')
	merged2.to_csv('/home/nboev/projects/def-sushant/nboev/preprocess/HGSVC2_v2_integrated_callset/sv/Blat/merges/'+chromo+'.Blatmerges.csv', sep='\t', index=False)
else:
	print('poor merged')


print('end of script')
print('END TIME:', datetime.datetime.now(timezone('EST')))
