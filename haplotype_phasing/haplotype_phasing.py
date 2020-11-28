
from os import path, makedirs
import sys
import math
import numpy as np
import pandas as pd
from fancyimpute import KNN

def clarks(ppllist,hlist):
	# put known haplotypes into knownlist
	knownlist = []
	for ppl in ppllist:
		matchlist = []
		#case where a person's genotype is 2 certain haplotypes
		if len(ppl)==1:
			for i in ppl:
				for j in i:
					hlist.append(j)
					if j not in knownlist:
						knownlist.append(j)
		else:
			# fine the maximum matches for a pair of haplotypes and define that as the final haplotype
			for i in ppl:
				matches=0
				for j in i:
					if j in knownlist:
						matches=matches+1
				matchlist.append(matches)
			m = matchlist.index(max(matchlist))
			for h in ppl[m]:
				hlist.append(h)
				if h not in knownlist:
					knownlist.append(h)
	return;



file_folder = path.dirname(path.abspath(__file__))
#file = open(sys.argv[1], "r+")
#glist = file.readlines()

file = pd.read_csv(sys.argv[1],sep=' ',names=range(0,50))


knn_imputer = KNN()
#svdimputer = IterativeSVD()
#imputer = IterativeImputer()
glist = file.copy(deep=True)
glist[glist=="*"]=np.nan
#gtest = glist.iloc[0:10,0:10]
#glist.iloc[:, :] = knn_imputer.fit_transform(glist)
glist.iloc[:,:] = np.round(knn_imputer.fit_transform(glist))
pd.options.display.float_format = '{:,.0f}'.format
np.savetxt(r"test_data_genotype_full_temp.txt", glist.values, fmt='%d')

file = open(r"test_data_genotype_full_temp.txt","r+")
glist = file.readlines()

# create haplotype phasing matrix
i = 0
# final haplotype list (used in clark's algorithm)
hlist = []
# list of people's potential haplotypes
# people = [] ## MOVED TO INSIDE THE LOOP so that clark's algorithm is only run on every 50 people
while len(glist)-i not in [4,8,12,16]:
	j = 0
	people=[]
	while j < 100:
		numbers = []
		if glist[i][j] == "0":
			numbers.append(0)
		if glist[i][j] == "1":
			numbers.append(1)
		if glist[i][j] == "2":
			numbers.append(2)
		if glist[i+1][j] == "0":
			numbers.append(0)
		if glist[i+1][j] == "1":
			numbers.append(1)
		if glist[i+1][j] == "2":
			numbers.append(2)
		if glist[i+2][j] == "0":
			numbers.append(0)
		if glist[i+2][j] == "1":
			numbers.append(1)
		if glist[i+2][j] == "2":
			numbers.append(2)
		if glist[i+3][j] == "0":
			numbers.append(0)
		if glist[i+3][j] == "1":
			numbers.append(1)
		if glist[i+3][j] == "2":
			numbers.append(2)
		if glist[i+4][j] == "0":
			numbers.append(0)
		if glist[i+4][j] == "1":
			numbers.append(1)
		if glist[i+4][j] == "2":
			numbers.append(2)
		onecount = 0
		for n in numbers:
			if n == 1:
				onecount += 1
		haplotypepairs = []
		if onecount == 0 or 1:
			solepair = [[], []]
			for n in numbers:
				if n == 0:
					solepair[0].append(0)
					solepair[1].append(0)
				if n == 2:
					solepair[0].append(1)
					solepair[1].append(1)
				if n == 1:
					solepair[0].append(0)
					solepair[1].append(1)
			haplotypepairs = [solepair]
		if onecount == 2:
			onesleft = 2
			twopairs = [ [[],[]], [[],[]] ]
			for n in numbers:
				if n == 0:
					for pair in twopairs:
						pair[0].append(0)
						pair[1].append(0)
				if n == 2:
					for pair in twopairs:
						pair[0].append(1)
						pair[1].append(1)
				if n == 1:
					if onesleft==2:
						twopairs[0][0].append(1)
						twopairs[0][1].append(0)
						twopairs[1][0].append(1)
						twopairs[1][1].append(0)
					if onesleft==1:
						twopairs[0][0].append(1)
						twopairs[0][1].append(0)
						twopairs[1][0].append(0)
						twopairs[1][1].append(1)
					oneseleft=onesleft-1
			haplotypepairs = twopairs
		if onecount == 3:
			onesleft = 3
			fourpairs = [ [[],[]], [[],[]], [[],[]], [[],[]]]
			for n in numbers:
				if n == 0:
					for pair in fourpairs:
						pair[0].append(0)
						pair[1].append(0)
				if n == 2:
					for pair in fourpairs:
						pair[0].append(1)
						pair[1].append(1)
				if n == 1:
					if onesleft==3:
						for pair in fourpairs:
							pair[0].append(1)
							pair[1].append(0)
					if onesleft==2:
						fourpairs[0][0].append(1)
						fourpairs[0][1].append(0)
						fourpairs[1][0].append(0)
						fourpairs[1][1].append(1)
						fourpairs[2][0].append(1)
						fourpairs[2][1].append(0)
						fourpairs[3][0].append(0)
						fourpairs[3][1].append(1)
					if onesleft==1:
						fourpairs[0][0].append(1)
						fourpairs[0][1].append(0)
						fourpairs[1][0].append(1)
						fourpairs[1][1].append(0)
						fourpairs[2][0].append(0)
						fourpairs[2][1].append(1)
						fourpairs[3][0].append(0)
						fourpairs[3][1].append(1)
					onesleft = onesleft-1
			haplotypepairs = fourpairs
		if onecount == 4:
			onesleft=4
			eightpairs = [ [[],[]], [[],[]], [[],[]], [[],[]], [[],[]], [[],[]], [[],[]], [[],[]]]
			for n in numbers:
				if n == 0:
					for pair in eightpairs:
						pair[0].append(0)
						pair[1].append(0)
				if n == 2:
					for pair in eightpairs:
						pair[0].append(1)
						pair[1].append(1)
				if n == 1:
					if onesleft==4:
						for pair in eightpairs:
							pair[0].append(1)
							pair[1].append(0)
					if onesleft==3:
						eightpairs[0][0].append(1)
						eightpairs[0][1].append(0)
						eightpairs[1][0].append(1)
						eightpairs[1][1].append(0)
						eightpairs[2][0].append(1)
						eightpairs[2][1].append(0)
						eightpairs[3][0].append(0)
						eightpairs[3][1].append(1)
						eightpairs[4][0].append(0)
						eightpairs[4][1].append(1)
						eightpairs[5][0].append(1)
						eightpairs[5][1].append(0)
						eightpairs[6][0].append(0)
						eightpairs[6][1].append(1)
						eightpairs[7][0].append(0)
						eightpairs[7][1].append(1)
					if onesleft==2:
						eightpairs[0][0].append(1)
						eightpairs[0][1].append(0)
						eightpairs[1][0].append(1)
						eightpairs[1][1].append(0)
						eightpairs[2][0].append(0)
						eightpairs[2][1].append(1)
						eightpairs[3][0].append(0)
						eightpairs[3][1].append(1)
						eightpairs[4][0].append(1)
						eightpairs[4][1].append(0)
						eightpairs[5][0].append(0)
						eightpairs[5][1].append(1)
						eightpairs[6][0].append(1)
						eightpairs[6][1].append(0)
						eightpairs[7][0].append(0)
						eightpairs[7][1].append(1)
					if onesleft==1:
						eightpairs[0][0].append(1)
						eightpairs[0][1].append(0)
						eightpairs[1][0].append(0)
						eightpairs[1][1].append(1)
						eightpairs[2][0].append(0)
						eightpairs[2][1].append(1)
						eightpairs[3][0].append(0)
						eightpairs[3][1].append(1)
						eightpairs[4][0].append(1)
						eightpairs[4][1].append(0)
						eightpairs[5][0].append(1)
						eightpairs[5][1].append(0)
						eightpairs[6][0].append(1)
						eightpairs[6][1].append(0)
						eightpairs[7][0].append(1)
						eightpairs[7][1].append(0)
					onesleft = onesleft-1
			haplotypepairs = eightpairs
		if onecount == 5:
			onesleft=5
			sixteenpairs = [ [[],[]], [[],[]], [[],[]], [[],[]], [[],[]], [[],[]], [[],[]], [[],[]], [[],[]], [[],[]], [[],[]], [[],[]], [[],[]], [[],[]], [[],[]], [[],[]]]
			for n in numbers:
				if n == 0:
					for pair in sixteenpairs:
						pair[0].append(0)
						pair[1].append(0)
				if n == 2:
					for pair in sixteenpairs:
						pair[0].append(1)
						pair[1].append(1)
				if n == 1:
					if onesleft==5:
						for pair in sixteenpairs:
							pair[0].append(1)
							pair[1].append(0)
					if onesleft==4:
						sixteenpairs[0][0].append(1)
						sixteenpairs[0][1].append(0)
						sixteenpairs[1][0].append(1)
						sixteenpairs[1][1].append(0)
						sixteenpairs[2][0].append(1)
						sixteenpairs[2][1].append(0)
						sixteenpairs[3][0].append(1)
						sixteenpairs[3][1].append(0)
						sixteenpairs[4][0].append(0)
						sixteenpairs[4][1].append(1)
						sixteenpairs[5][0].append(0)
						sixteenpairs[5][1].append(1)
						sixteenpairs[6][0].append(0)
						sixteenpairs[6][1].append(1)
						sixteenpairs[7][0].append(0)
						sixteenpairs[7][1].append(1)
						sixteenpairs[8][0].append(1)
						sixteenpairs[8][1].append(0)
						sixteenpairs[9][0].append(1)
						sixteenpairs[9][1].append(0)
						sixteenpairs[10][0].append(0)
						sixteenpairs[10][1].append(1)
						sixteenpairs[11][0].append(0)
						sixteenpairs[11][1].append(1)
						sixteenpairs[12][0].append(0)
						sixteenpairs[12][1].append(1)
						sixteenpairs[13][0].append(0)
						sixteenpairs[13][1].append(1)
						sixteenpairs[14][0].append(1)
						sixteenpairs[14][1].append(0)
						sixteenpairs[15][0].append(1)
						sixteenpairs[15][1].append(0)
					if onesleft==3:
						sixteenpairs[0][0].append(1)
						sixteenpairs[0][1].append(0)
						sixteenpairs[1][0].append(1)
						sixteenpairs[1][1].append(0)
						sixteenpairs[2][0].append(1)
						sixteenpairs[2][1].append(0)
						sixteenpairs[3][0].append(0)
						sixteenpairs[3][1].append(1)
						sixteenpairs[4][0].append(0)
						sixteenpairs[4][1].append(1)
						sixteenpairs[5][0].append(1)
						sixteenpairs[5][1].append(0)
						sixteenpairs[6][0].append(0)
						sixteenpairs[6][1].append(1)
						sixteenpairs[7][0].append(0)
						sixteenpairs[7][1].append(1)
						sixteenpairs[8][0].append(0)
						sixteenpairs[8][1].append(1)
						sixteenpairs[9][0].append(0)
						sixteenpairs[9][1].append(1)
						sixteenpairs[10][0].append(1)
						sixteenpairs[10][1].append(0)
						sixteenpairs[11][0].append(1)
						sixteenpairs[11][1].append(0)
						sixteenpairs[12][0].append(0)
						sixteenpairs[12][1].append(1)
						sixteenpairs[13][0].append(1)
						sixteenpairs[13][1].append(0)
						sixteenpairs[14][0].append(0)
						sixteenpairs[14][1].append(1)
						sixteenpairs[15][0].append(1)
						sixteenpairs[15][1].append(0)
					if onesleft==2:
						sixteenpairs[0][0].append(1)
						sixteenpairs[0][1].append(0)
						sixteenpairs[1][0].append(1)
						sixteenpairs[1][1].append(0)
						sixteenpairs[2][0].append(0)
						sixteenpairs[2][1].append(1)
						sixteenpairs[3][0].append(0)
						sixteenpairs[3][1].append(1)
						sixteenpairs[4][0].append(0)
						sixteenpairs[4][1].append(1)
						sixteenpairs[5][0].append(0)
						sixteenpairs[5][1].append(1)
						sixteenpairs[6][0].append(1)
						sixteenpairs[6][1].append(0)
						sixteenpairs[7][0].append(0)
						sixteenpairs[7][1].append(1)
						sixteenpairs[8][0].append(1)
						sixteenpairs[8][1].append(0)
						sixteenpairs[9][0].append(0)
						sixteenpairs[9][1].append(1)
						sixteenpairs[10][0].append(1)
						sixteenpairs[10][1].append(0)
						sixteenpairs[11][0].append(0)
						sixteenpairs[11][1].append(1)
						sixteenpairs[12][0].append(1)
						sixteenpairs[12][1].append(0)
						sixteenpairs[13][0].append(1)
						sixteenpairs[13][1].append(0)
						sixteenpairs[14][0].append(1)
						sixteenpairs[14][1].append(0)
						sixteenpairs[15][0].append(0)
						sixteenpairs[15][1].append(1)
					if onesleft==1:
						sixteenpairs[0][0].append(1)
						sixteenpairs[0][1].append(0)
						sixteenpairs[1][0].append(0)
						sixteenpairs[1][1].append(1)
						sixteenpairs[2][0].append(0)
						sixteenpairs[2][1].append(1)
						sixteenpairs[3][0].append(0)
						sixteenpairs[3][1].append(1)
						sixteenpairs[4][0].append(0)
						sixteenpairs[4][1].append(1)
						sixteenpairs[5][0].append(0)
						sixteenpairs[5][1].append(1)
						sixteenpairs[6][0].append(0)
						sixteenpairs[6][1].append(1)
						sixteenpairs[7][0].append(1)
						sixteenpairs[7][1].append(0)
						sixteenpairs[8][0].append(0)
						sixteenpairs[8][1].append(1)
						sixteenpairs[9][0].append(1)
						sixteenpairs[9][1].append(0)
						sixteenpairs[10][0].append(0)
						sixteenpairs[10][1].append(1)
						sixteenpairs[11][0].append(1)
						sixteenpairs[11][1].append(0)
						sixteenpairs[12][0].append(1)
						sixteenpairs[12][1].append(0)
						sixteenpairs[13][0].append(1)
						sixteenpairs[13][1].append(0)
						sixteenpairs[14][0].append(1)
						sixteenpairs[14][1].append(0)
						sixteenpairs[15][0].append(1)
						sixteenpairs[15][1].append(0)
					onesleft=onesleft-1
			haplotypepairs = sixteenpairs
		people.append(haplotypepairs)
		j += 2
	clarks(people,hlist)
	i += 5


while i < len(glist):
	people=[]
	j = 0
	while j < 100:
		numbers = []
		if glist[i][j] == "0":
			numbers.append(0)
		if glist[i][j] == "1":
			numbers.append(1)
		if glist[i][j] == "2":
			numbers.append(2)
		if glist[i+1][j] == "0":
			numbers.append(0)
		if glist[i+1][j] == "1":
			numbers.append(1)
		if glist[i+1][j] == "2":
			numbers.append(2)
		if glist[i+2][j] == "0":
			numbers.append(0)
		if glist[i+2][j] == "1":
			numbers.append(1)
		if glist[i+2][j] == "2":
			numbers.append(2)
		if glist[i+3][j] == "0":
			numbers.append(0)
		if glist[i+3][j] == "1":
			numbers.append(1)
		if glist[i+3][j] == "2":
			numbers.append(2)
		onecount = 0
		for n in numbers:
			if n == 1:
				onecount += 1
		haplotypepairs = []
		if onecount == 0 or 1:
			solepair = [[], []]
			for n in numbers:
				if n == 0:
					solepair[0].append(0)
					solepair[1].append(0)
				if n == 2:
					solepair[0].append(1)
					solepair[1].append(1)
				if n == 1:
					solepair[0].append(0)
					solepair[1].append(1)
			haplotypepairs = [solepair]
		if onecount == 2:
			onesleft=2
			twopairs = [ [[],[]], [[],[]] ]
			for n in numbers:
				if n == 0:
					for pair in twopairs:
						pair[0].append(0)
						pair[1].append(0)
				if n == 2:
					for pair in twopairs:
						pair[0].append(1)
						pair[1].append(1)
				if n == 1:
					if onesleft==2:
						twopairs[0][0].append(1)
						twopairs[0][1].append(0)
						twopairs[1][0].append(1)
						twopairs[1][1].append(0)
					if onesleft==1:
						twopairs[0][0].append(1)
						twopairs[0][1].append(0)
						twopairs[1][0].append(0)
						twopairs[1][1].append(1)
					onesleft=onesleft-1	
			haplotypepairs = twopairs
		if onecount == 3:
			onesleft=3
			fourpairs = [ [[],[]], [[],[]], [[],[]], [[],[]]]
			for n in numbers:
				if n == 0:
					for pair in fourpairs:
						pair[0].append(0)
						pair[1].append(0)
				if n == 2:
					for pair in fourpairs:
						pair[0].append(1)
						pair[1].append(1)
				if n == 1:
					if onesleft==3:
						for pair in fourpairs:
							pair[0].append(1)
							pair[1].append(0)
					if onesleft==2:
						fourpairs[0][0].append(1)
						fourpairs[0][1].append(0)
						fourpairs[1][0].append(0)
						fourpairs[1][1].append(1)
						fourpairs[2][0].append(1)
						fourpairs[2][1].append(0)
						fourpairs[3][0].append(0)
						fourpairs[3][1].append(1)
					if onesleft==1:
						fourpairs[0][0].append(1)
						fourpairs[0][1].append(0)
						fourpairs[1][0].append(1)
						fourpairs[1][1].append(0)
						fourpairs[2][0].append(0)
						fourpairs[2][1].append(1)
						fourpairs[3][0].append(0)
						fourpairs[3][1].append(1)
					onesleft = onesleft-1
			haplotypepairs = fourpairs
		if onecount == 4:
			onesleft=4
			eightpairs = [ [[],[]], [[],[]], [[],[]], [[],[]], [[],[]], [[],[]], [[],[]], [[],[]]]
			for n in numbers:
				if n == 0:
					for pair in eightpairs:
						pair[0].append(0)
						pair[1].append(0)
				if n == 2:
					for pair in eightpairs:
						pair[0].append(1)
						pair[1].append(1)
				if n == 1:
					if onesleft==4:
						for pair in eightpairs:
							pair[0].append(1)
							pair[1].append(0)
					if onesleft==3:
						eightpairs[0][0].append(1)
						eightpairs[0][1].append(0)
						eightpairs[1][0].append(1)
						eightpairs[1][1].append(0)
						eightpairs[2][0].append(1)
						eightpairs[2][1].append(0)
						eightpairs[3][0].append(0)
						eightpairs[3][1].append(1)
						eightpairs[4][0].append(0)
						eightpairs[4][1].append(1)
						eightpairs[5][0].append(1)
						eightpairs[5][1].append(0)
						eightpairs[6][0].append(0)
						eightpairs[6][1].append(1)
						eightpairs[7][0].append(0)
						eightpairs[7][1].append(1)
					if onesleft==2:
						eightpairs[0][0].append(1)
						eightpairs[0][1].append(0)
						eightpairs[1][0].append(1)
						eightpairs[1][1].append(0)
						eightpairs[2][0].append(0)
						eightpairs[2][1].append(1)
						eightpairs[3][0].append(0)
						eightpairs[3][1].append(1)
						eightpairs[4][0].append(1)
						eightpairs[4][1].append(0)
						eightpairs[5][0].append(0)
						eightpairs[5][1].append(1)
						eightpairs[6][0].append(1)
						eightpairs[6][1].append(0)
						eightpairs[7][0].append(0)
						eightpairs[7][1].append(1)
					if onesleft==1:
						eightpairs[0][0].append(1)
						eightpairs[0][1].append(0)
						eightpairs[1][0].append(0)
						eightpairs[1][1].append(1)
						eightpairs[2][0].append(0)
						eightpairs[2][1].append(1)
						eightpairs[3][0].append(0)
						eightpairs[3][1].append(1)
						eightpairs[4][0].append(1)
						eightpairs[4][1].append(0)
						eightpairs[5][0].append(1)
						eightpairs[5][1].append(0)
						eightpairs[6][0].append(1)
						eightpairs[6][1].append(0)
						eightpairs[7][0].append(1)
						eightpairs[7][1].append(0)
					onesleft = onesleft-1
			haplotypepairs = eightpairs
		people.append(haplotypepairs)
		j += 2
	clarks(people,hlist)
	i += 4



'''
i = 0
while i < 785400:
		print (hlist[i])
		i += 1

'''



firstfifty = []
secondfifty = []

i = 0
while i < len(glist):
	firstfifty.append("a b c d e f g h i j k l m n o p q r s t u v w x y z A B C D E F G H I J K L M N O P Q R S T U V W X ")
	secondfifty.append("a b c d e f g h i j k l m n o p q r s t u v w x y z A B C D E F G H I J K L M N O P Q R S T U V W X\n")
	i += 1


i = 0

while i < 785400:
	# 
	if i < 785100:
		# 
		k = 0
		# 
		while k < 5:
			# 
			j = i
			# 
			while j < (i + 100):
				# 
				if (j-i) < 50:
					if (j-i) == 0:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("a","0")
						else:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("a","1")
					if (j-i) == 1:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("b","0")
						else:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("b","1")							
					if (j-i) == 2:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("c","0")
						else:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("c","1")							
					if (j-i) == 3:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("d","0")
						else:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("d","1")							
					if (j-i) == 4:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("e","0")
						else:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("e","1")							
					if (j-i) == 5:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("f","0")
						else:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("f","1")							
					if (j-i) == 6:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("g","0")
						else:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("g","1")							
					if (j-i) == 7:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("h","0")
						else:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("h","1")							
					if (j-i) == 8:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("i","0")
						else:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("i","1")							
					if (j-i) == 9:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("j","0")
						else:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("j","1")							
					if (j-i) == 10:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("k","0")
						else:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("k","1")							
					if (j-i) == 11:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("l","0")
						else:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("l","1")							
					if (j-i) == 12:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("m","0")
						else:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("m","1")							
					if (j-i) == 13:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("n","0")
						else:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("n","1")							
					if (j-i) == 14:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("o","0")
						else:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("o","1")							
					if (j-i) == 15:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("p","0")
						else:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("p","1")							
					if (j-i) == 16:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("q","0")
						else:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("q","1")							
					if (j-i) == 17:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("r","0")
						else:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("r","1")							
					if (j-i) == 18:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("s","0")
						else:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("s","1")							
					if (j-i) == 19:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("t","0")
						else:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("t","1")							
					if (j-i) == 20:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("u","0")
						else:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("u","1")							
					if (j-i) == 21:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("v","0")
						else:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("v","1")							
					if (j-i) == 22:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("w","0")
						else:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("w","1")							
					if (j-i) == 23:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("x","0")
						else:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("x","1")							
					if (j-i) == 24:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("y","0")
						else:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("y","1")
					if (j-i) == 25:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("z","0")
						else:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("z","1")
					if (j-i) == 26:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("A","0")
						else:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("A","1")
					if (j-i) == 27:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("B","0")
						else:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("B","1")
					if (j-i) == 28:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("C","0")
						else:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("C","1")
					if (j-i) == 29:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("D","0")
						else:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("D","1")
					if (j-i) == 30:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("E","0")
						else:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("E","1")
					if (j-i) == 31:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("F","0")
						else:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("F","1")
					if (j-i) == 32:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("G","0")
						else:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("G","1")
					if (j-i) == 33:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("H","0")
						else:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("H","1")
					if (j-i) == 34:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("I","0")
						else:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("I","1")
					if (j-i) == 35:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("J","0")
						else:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("J","1")
					if (j-i) == 36:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("K","0")
						else:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("K","1")
					if (j-i) == 37:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("L","0")
						else:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("L","1")
					if (j-i) == 38:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("M","0")
						else:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("M","1")
					if (j-i) == 39:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("N","0")
						else:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("N","1")
					if (j-i) == 40:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("O","0")
						else:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("O","1")
					if (j-i) == 41:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("P","0")
						else:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("P","1")
					if (j-i) == 42:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("Q","0")
						else:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("Q","1")
					if (j-i) == 43:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("R","0")
						else:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("R","1")
					if (j-i) == 44:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("S","0")
						else:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("S","1")
					if (j-i) == 45:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("T","0")
						else:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("T","1")
					if (j-i) == 46:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("U","0")
						else:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("U","1")
					if (j-i) == 47:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("V","0")
						else:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("V","1")
					if (j-i) == 48:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("W","0")
						else:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("W","1")
					if (j-i) == 49:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("X","0")
						else:
							firstfifty[math.floor(5*(i/100) + k)] = firstfifty[math.floor(5*(i/100) + k)].replace("X","1")
							#
							#
				else:
					if (j-i) == 50:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("a","0")
						else:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("a","1")
					if (j-i) == 51:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("b","0")
						else:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("b","1")							
					if (j-i) == 52:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("c","0")
						else:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("c","1")							
					if (j-i) == 53:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("d","0")
						else:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("d","1")							
					if (j-i) == 54:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("e","0")
						else:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("e","1")							
					if (j-i) == 55:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("f","0")
						else:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("f","1")							
					if (j-i) == 56:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("g","0")
						else:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("g","1")							
					if (j-i) == 57:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("h","0")
						else:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("h","1")							
					if (j-i) == 58:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("i","0")
						else:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("i","1")							
					if (j-i) == 59:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("j","0")
						else:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("j","1")							
					if (j-i) == 60:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("k","0")
						else:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("k","1")							
					if (j-i) == 61:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("l","0")
						else:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("l","1")							
					if (j-i) == 62:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("m","0")
						else:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("m","1")							
					if (j-i) == 63:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("n","0")
						else:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("n","1")							
					if (j-i) == 64:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("o","0")
						else:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("o","1")							
					if (j-i) == 65:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("p","0")
						else:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("p","1")							
					if (j-i) == 66:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("q","0")
						else:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("q","1")							
					if (j-i) == 67:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("r","0")
						else:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("r","1")							
					if (j-i) == 68:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("s","0")
						else:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("s","1")							
					if (j-i) == 69:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("t","0")
						else:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("t","1")							
					if (j-i) == 70:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("u","0")
						else:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("u","1")							
					if (j-i) == 71:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("v","0")
						else:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("v","1")							
					if (j-i) == 72:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("w","0")
						else:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("w","1")							
					if (j-i) == 73:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("x","0")
						else:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("x","1")							
					if (j-i) == 74:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("y","0")
						else:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("y","1")
					if (j-i) == 75:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("z","0")
						else:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("z","1")
					if (j-i) == 76:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("A","0")
						else:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("A","1")
					if (j-i) == 77:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("B","0")
						else:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("B","1")
					if (j-i) == 78:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("C","0")
						else:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("C","1")
					if (j-i) == 79:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("D","0")
						else:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("D","1")
					if (j-i) == 80:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("E","0")
						else:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("E","1")
					if (j-i) == 81:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("F","0")
						else:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("F","1")
					if (j-i) == 82:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("G","0")
						else:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("G","1")
					if (j-i) == 83:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("H","0")
						else:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("H","1")
					if (j-i) == 84:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("I","0")
						else:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("I","1")
					if (j-i) == 85:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("J","0")
						else:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("J","1")
					if (j-i) == 86:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("K","0")
						else:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("K","1")
					if (j-i) == 87:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("L","0")
						else:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("L","1")
					if (j-i) == 88:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("M","0")
						else:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("M","1")
					if (j-i) == 89:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("N","0")
						else:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("N","1")
					if (j-i) == 90:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("O","0")
						else:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("O","1")
					if (j-i) == 91:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("P","0")
						else:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("P","1")
					if (j-i) == 92:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("Q","0")
						else:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("Q","1")
					if (j-i) == 93:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("R","0")
						else:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("R","1")
					if (j-i) == 94:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("S","0")
						else:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("S","1")
					if (j-i) == 95:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("T","0")
						else:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("T","1")
					if (j-i) == 96:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("U","0")
						else:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("U","1")
					if (j-i) == 97:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("V","0")
						else:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("V","1")
					if (j-i) == 98:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("W","0")
						else:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("W","1")
					if (j-i) == 99:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("X","0")
						else:
							secondfifty[math.floor(5*(i/100) + k)] = secondfifty[math.floor(5*(i/100) + k)].replace("X","1")
					#print(secondfifty[math.floor(5*(i/100) + k)])
				j += 1
				#
			k += 1
			#
	else:
		#
		k = 0
		#
		while k < 4:
			#
			j = i
			#
			while j < (i + 100):
				#
				if (j-i) < 50:
					if (j-i) == 0:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("a","0")
						else:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("a","1")
					if (j-i) == 1:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("b","0")
						else:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("b","1")							
					if (j-i) == 2:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("c","0")
						else:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("c","1")							
					if (j-i) == 3:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("d","0")
						else:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("d","1")							
					if (j-i) == 4:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("e","0")
						else:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("e","1")							
					if (j-i) == 5:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("f","0")
						else:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("f","1")							
					if (j-i) == 6:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("g","0")
						else:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("g","1")							
					if (j-i) == 7:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("h","0")
						else:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("h","1")							
					if (j-i) == 8:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("i","0")
						else:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("i","1")							
					if (j-i) == 9:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("j","0")
						else:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("j","1")							
					if (j-i) == 10:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("k","0")
						else:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("k","1")							
					if (j-i) == 11:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("l","0")
						else:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("l","1")							
					if (j-i) == 12:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("m","0")
						else:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("m","1")							
					if (j-i) == 13:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("n","0")
						else:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("n","1")							
					if (j-i) == 14:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("o","0")
						else:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("o","1")							
					if (j-i) == 15:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("p","0")
						else:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("p","1")							
					if (j-i) == 16:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("q","0")
						else:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("q","1")							
					if (j-i) == 17:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("r","0")
						else:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("r","1")							
					if (j-i) == 18:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("s","0")
						else:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("s","1")							
					if (j-i) == 19:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("t","0")
						else:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("t","1")							
					if (j-i) == 20:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("u","0")
						else:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("u","1")							
					if (j-i) == 21:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("v","0")
						else:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("v","1")							
					if (j-i) == 22:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("w","0")
						else:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("w","1")							
					if (j-i) == 23:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("x","0")
						else:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("x","1")							
					if (j-i) == 24:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("y","0")
						else:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("y","1")
					if (j-i) == 25:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("z","0")
						else:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("z","1")
					if (j-i) == 26:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("A","0")
						else:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("A","1")
					if (j-i) == 27:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("B","0")
						else:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("B","1")
					if (j-i) == 28:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("C","0")
						else:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("C","1")
					if (j-i) == 29:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("D","0")
						else:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("D","1")
					if (j-i) == 30:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("E","0")
						else:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("E","1")
					if (j-i) == 31:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("F","0")
						else:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("F","1")
					if (j-i) == 32:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("G","0")
						else:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("G","1")
					if (j-i) == 33:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("H","0")
						else:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("H","1")
					if (j-i) == 34:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("I","0")
						else:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("I","1")
					if (j-i) == 35:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("J","0")
						else:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("J","1")
					if (j-i) == 36:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("K","0")
						else:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("K","1")
					if (j-i) == 37:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("L","0")
						else:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("L","1")
					if (j-i) == 38:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("M","0")
						else:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("M","1")
					if (j-i) == 39:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("N","0")
						else:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("N","1")
					if (j-i) == 40:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("O","0")
						else:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("O","1")
					if (j-i) == 41:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("P","0")
						else:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("P","1")
					if (j-i) == 42:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("Q","0")
						else:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("Q","1")
					if (j-i) == 43:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("R","0")
						else:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("R","1")
					if (j-i) == 44:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("S","0")
						else:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("S","1")
					if (j-i) == 45:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("T","0")
						else:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("T","1")
					if (j-i) == 46:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("U","0")
						else:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("U","1")
					if (j-i) == 47:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("V","0")
						else:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("V","1")
					if (j-i) == 48:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("W","0")
						else:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("W","1")
					if (j-i) == 49:
						if hlist[j][k] == 0:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("X","0")
						else:
							firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = firstfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("X","1")
							#
							#
				else:
					if (j-i) == 50:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("a","0")
						else:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("a","1")
					if (j-i) == 51:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("b","0")
						else:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("b","1")							
					if (j-i) == 52:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("c","0")
						else:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("c","1")							
					if (j-i) == 53:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("d","0")
						else:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("d","1")							
					if (j-i) == 54:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("e","0")
						else:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("e","1")							
					if (j-i) == 55:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("f","0")
						else:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("f","1")							
					if (j-i) == 56:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("g","0")
						else:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("g","1")							
					if (j-i) == 57:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("h","0")
						else:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("h","1")							
					if (j-i) == 58:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("i","0")
						else:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("i","1")							
					if (j-i) == 59:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("j","0")
						else:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("j","1")							
					if (j-i) == 60:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("k","0")
						else:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("k","1")							
					if (j-i) == 61:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("l","0")
						else:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("l","1")							
					if (j-i) == 62:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("m","0")
						else:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("m","1")							
					if (j-i) == 63:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("n","0")
						else:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("n","1")							
					if (j-i) == 64:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("o","0")
						else:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("o","1")							
					if (j-i) == 65:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("p","0")
						else:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("p","1")							
					if (j-i) == 66:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("q","0")
						else:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("q","1")							
					if (j-i) == 67:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("r","0")
						else:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("r","1")							
					if (j-i) == 68:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("s","0")
						else:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("s","1")							
					if (j-i) == 69:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("t","0")
						else:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("t","1")							
					if (j-i) == 70:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("u","0")
						else:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("u","1")							
					if (j-i) == 71:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("v","0")
						else:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("v","1")							
					if (j-i) == 72:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("w","0")
						else:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("w","1")							
					if (j-i) == 73:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("x","0")
						else:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("x","1")							
					if (j-i) == 74:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("y","0")
						else:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("y","1")
					if (j-i) == 75:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("z","0")
						else:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("z","1")
					if (j-i) == 76:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("A","0")
						else:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("A","1")
					if (j-i) == 77:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("B","0")
						else:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("B","1")
					if (j-i) == 78:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("C","0")
						else:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("C","1")
					if (j-i) == 79:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("D","0")
						else:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("D","1")
					if (j-i) == 80:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("E","0")
						else:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("E","1")
					if (j-i) == 81:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("F","0")
						else:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("F","1")
					if (j-i) == 82:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("G","0")
						else:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("G","1")
					if (j-i) == 83:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("H","0")
						else:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("H","1")
					if (j-i) == 84:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("I","0")
						else:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("I","1")
					if (j-i) == 85:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("J","0")
						else:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("J","1")
					if (j-i) == 86:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("K","0")
						else:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("K","1")
					if (j-i) == 87:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("L","0")
						else:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("L","1")
					if (j-i) == 88:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("M","0")
						else:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("M","1")
					if (j-i) == 89:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("N","0")
						else:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("N","1")
					if (j-i) == 90:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("O","0")
						else:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("O","1")
					if (j-i) == 91:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("P","0")
						else:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("P","1")
					if (j-i) == 92:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("Q","0")
						else:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("Q","1")
					if (j-i) == 93:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("R","0")
						else:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("R","1")
					if (j-i) == 94:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("S","0")
						else:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("S","1")
					if (j-i) == 95:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("T","0")
						else:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("T","1")
					if (j-i) == 96:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("U","0")
						else:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("U","1")
					if (j-i) == 97:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("V","0")
						else:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("V","1")
					if (j-i) == 98:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("W","0")
						else:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("W","1")
					if (j-i) == 99:
						if hlist[j][k] == 0:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("X","0")
						else:
							secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))] = secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))].replace("X","1")
					#print(secondfifty[math.floor(5*(785100/100) + k + 4*( (i-785100) /100))])
					#
				j += 1
				#
			k += 1
			#
	i += 100



output = []

i = 0
while i < 39267:
	output.append( ( firstfifty[i] + secondfifty[i] ) )
	#print(output[i])
	i += 1


# print(len(output))


submission = open(r"test_data_sol.txt","w")
submission.writelines(output)
submission.close()



file.close()
print("phasing successful!")

