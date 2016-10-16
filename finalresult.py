from Bio import SeqIO
import pandas as pd
import re
import os 
import numpy as np
import pickle
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import minimum_spanning_tree
from scipy.sparse.csgraph import connected_components
from scipy.spatial.distance import cdist
import networkx as nx
import matplotlib.pyplot as plt
import pylab
from itertools import count

##Create membership matrix
###Need to change the directory which contains all of the binning result 
###Run metabat N times
N = 200
print(N+1)
files = os.listdir('/Users/bian0553/Downloads/finalresult/finalbin')
if len(files) > N:
	location = files.index('.DS_Store')
	del files[location]
num  = 1
headers = ["Contig",'Membership']
meta = pd.DataFrame()
for fn in files:
	handle = pd.read_csv("/Users/bian0553/Downloads/finalresult/finalbin/"+fn, names = headers,sep = '\t')
	meta['contigs']= handle['Contig']
	meta[num] = handle['Membership']
	num = num+1

meta = meta.replace(0,np.nan)

##check if NA shows up in the last row of matrix. 
lastrow = meta.iloc[meta.shape[0]-1][1:N+1]
##The last row cannot have NA
if(sum(pd.isnull(lastrow)) == N):
	##get to know the first elements without NA showing up
	for i in range(meta.shape[0]):
		row = meta.iloc[i][1:N+1]
		if(sum(pd.isnull(row)) == 0):
			num = i
			break
	##Move the row just found without NA showing to the last row 
	meta.loc[meta.shape[0]] = meta.iloc[num]
	meta = meta.drop(num)
	##get final membership 
	meta = meta.reset_index(drop = True)


##load the contig size document (note the directory changing)
contigsize = pd.read_csv("contigs.txt",sep = None,header =0,engine = 'python')

###################################do subsampling#####################################################
####If you get the membership data from the assembly.fa.subset, ignore this function. 
###If get the membership data from the whole dataset and you want to do the subsampling
###then run the subsampling function

###Input the membership data (meta) and the file contains contigs size table. Return subsampled membership matrix
def subsampling(result,congtigsize):

	#get how many genomes in total
	genome = []
	for name in result['contigs']:
		a = re.sub('_[0-9]*-[0-9]*',"",name)
		m = re.search('[A-Za-z\|0-9\.]*_\[([A-Za-z_0-9]*)\]',a)
		genome.append(m.group(1))

	genome_sum = np.unique(genome)
	result['genome'] = genome
	#get genome table
	columns = ['genome','number']
	index = np.arange(len(genome_sum)) # array of numbers for the number of samples
	gentable = pd.DataFrame(columns=columns, index = index)
	num = 0
	for unique in genome_sum:
		gentable['genome'][num]= unique
		gentable['number'][num]=genome.count(unique)
		num = num+1

	genome = []
	for num in range(len(gentable['genome'])):
		if gentable['number'][num] <= 1000 and gentable['number'][num] >= 300:
			genome.append(gentable['genome'][num])

	sample1 = {}
	for gen in genome:
		loc =np.random.choice(np.where(result['genome'] == gen)[0],100)
		contig = []
		for i in loc:
			contig.append(result['contigs'][i])
		sample1[gen] = contig

	##get size for each items
	sizedic = {}
	num = 0
	for key in sample1:
		size = []
		for item in sample1[key]:
			index = list(contigsize['Name']).index(item)
			size.append(contigsize['Size'][index])
		sizedic[num] = size
		num = num+1

	##get sampling
	sample = {}
	num = 0
	for key in sample1:
		sample[num] = sample1[key]
		num = num + 1

	total = []
	for key in sizedic:
		total.append(sum(sizedic[key]))

	new = {}
	for key in sizedic:
		count = sum(sizedic[key])
		while (count > 330000):
			deletloc = sizedic[key].index(max(sizedic[key]))
			del sample[key][deletloc]
			sizedic[key].remove(max(sizedic[key]))
			count = sum(sizedic[key])
		
	total = []
	for key in sizedic:
		total.append(sum(sizedic[key]))


	size = []
	newitem = []
	for key in sample:
		for item in sample[key]:
			newitem.append(item)

	loc = []
	for item in newitem:
		loc.append(list(result['contigs']).index(item))


	test = result.ix[loc]
	test = test.reset_index(drop = True)
	return test

#newdata = subsampling(meta,contigsize)
####Then use 'newdata' as the new generated membership matrix in the following steps
########################################################################################


### meta is the membership matrix
def coassociation_matrix(iter,meta):
	rows = []
	cols = []
	

	for label in meta[iter].unique()[~np.isnan(meta[iter].unique())]:
		indices = np.where(meta[iter] == label)[0]
		for index1 in indices:
			for index2 in indices:
				rows.append(index1)
				cols.append(index2)

	data = np.ones((len(rows),))
	coam = csr_matrix((data, (rows, cols)), dtype='float')
	return coam


def compresult(coassociation_matrix,t,membership_matrix):
	mst = minimum_spanning_tree(-coassociation_matrix)
	mst = mst.toarray()
	mst[mst>t] = 0
	number_of_clusters, labels = connected_components(mst)
	G = nx.from_scipy_sparse_matrix(csr_matrix(mst))
	comp = list(nx.connected_components(G))
	binning_test = {}
	for i in range(len(comp)):
		value = []
		for index in comp[i]:
			value.append(membership_matrix['contigs'][index])
		binning_test[i] = value
	return binning_test


def getsize(binning_test,contigsize):
	##get size for each items
	sizedic = {}
	num = 0
	for key in binning_test:
		if len(binning_test[key]) > 1:
			size = []
			for item in binning_test[key]:
				index = list(contigsize['Name']).index(item)
				size.append(contigsize['Size'][index])
			sizedic[num] = size
			num = num+1
	return sizedic


def cleanbin(binning_test):
	####remove the contig numbers
	binning = {}
	for key in binning_test:
		new = []
		for value in binning_test[key]:
			new.append(re.sub('_[0-9]*-[0-9]*',"",value))
		binning[key] = new

	###get the [..] genome name
	for key in binning:
		new = []
		for item in binning[key]:
			m = re.search('[A-Za-z\|0-9\.]*_\[([A-Za-z_0-9]*)\]',item)
			new.append(m.group(1))
		binning[key] = new

	newbin = {}
	num = 0
	for key in binning:
		if len(binning[key])>1:
			newbin[num] = binning[key]	
			num = num + 1
	return newbin


def genomesize(binning, sizedic):
	####calculate the total size which has the same genome name inside[..]
	total = {}
	for key in binning:
		value = np.unique(np.asarray(binning[key]))
		totalsize = []
		for item in value:
			loc = np.where(np.asarray(binning[key]) == item)[0]
			size = 0
			for i in loc:
				size = size + sizedic[key][i]	
			totalsize.append(size)
		total[key] = totalsize	
	return total


##input size dictionary and cleaned binning dictionary
def getrecall(total,binning):	
	### get the location with the largest sizes
	maxloc = []
	for key in total:
		maxloc.append(total[key].index(max(total[key])))
	###get the represented genome in each cluster
	genome = []
	for key in binning:
		allgen = np.unique(np.asarray(binning[key]))
		genome.append(allgen[maxloc[key]])


	######get the size of represented genome in each cluster
	gensize = []
	for i in range(len(genome)):
		size = 0
		for key in binning:
			allgen = np.unique(np.asarray(binning[key]))
			if(genome[i] in allgen):
				loc = np.where(allgen == genome[i])[0]
				size = size + total[key][loc]
		gensize.append(size)	


	recall = []
	for key in total:
		nume = max(total[key])
		denom = gensize[key]
		recall.append(float(nume)/float(denom))
	return recall

def getprecision(total):
	precision = []
	for key in total:
		nume = max(total[key])
		denom = sum(total[key])
		precision.append(float(nume)/float(denom))
	return precision


def gettable(genedic):
	df = pd.DataFrame()
	df['Reference'] = contigsize['New'].unique()
	num = 0
	for key in genedic:
		if len(genedic[key]) > 1:
			df[num] = np.zeros(len(contigsize['New'].unique()))
			for item in genedic[key]:
				loc = np.where(item == df['Reference'])[0]
				df[num][loc] = df[num][loc] + 1
			num = num+1
	return df

###get coassociation matrix
c = coassociation_matrix(1,meta)
for i in range(2,N+1):
	c = c + coassociation_matrix(i,meta)






########plot the relationship between the F1 mean score with each time metabat running.######## 
index = range(1,N+1)
fscore = []
finalf = []
for i in index:
	c = coassociation_matrix(i,meta)
	binning_test = compresult(c,0,meta)
	sizedic = getsize(binning_test,contigsize)
	binning = cleanbin(binning_test)
	total = genomesize(binning, sizedic)
	precision = np.asarray(getprecision(total))
	recall = np.asarray(getrecall(total,binning))
	f1score = 2*(precision*recall)/(precision+recall)
	fscore.append(f1score.mean())
	finalf.append(max(fscore))	
	print(i)

fig = plt.figure()
ax = fig.add_subplot(111)
plt.plot(index,finalf)
plt.xlabel('each time of metabat run')
plt.ylabel('F1 score')

plt.title('Plot for F1 score with each time metabat')
x = fscore.index(max(finalf))+1
y = max(finalf)
xy=(x,y)
ax.annotate('(%s, %s)' % xy, xy=xy, textcoords='data') 
plt.grid(True)
plt.savefig("onetimef1score.png")


#######################plot the curves
c = c/N
trange = np.arange(0.1,1.0,0.1)
fscore = []
roc = []
poc = []
for t in trange:
	binning_test = compresult(c,-t,meta)
	sizedic = getsize(binning_test,contigsize)
	binning = cleanbin(binning_test)	
	total = genomesize(binning, sizedic)
	precision = np.asarray(getprecision(total))
	recall = np.asarray(getrecall(total,binning))
	f1score = 2*(precision*recall)/(precision+recall)
	fscore.append(f1score.mean())
	roc.append(recall.mean())
	poc.append(precision.mean())
	print(t)

plt.plot(trange,roc,label = 'recall')
plt.plot(trange,poc, label = 'precision')
plt.plot(trange,fscore,label = 'f1 score')
plt.xlabel('threshold')
plt.ylabel('score')
plt.title('Score VS Threshold choosing')
plt.grid(True)
plt.legend(loc = 3)
plt.savefig("score.png")


####check the relationship between f1 score and the replication numbers
##This is computationally expensive. Will run couple days
c = coassociation_matrix(1,meta)
binning_test = compresult(c,0,meta)
sizedic = getsize(binning_test,contigsize)
binning = cleanbin(binning_test)
total = genomesize(binning, sizedic)
precision = np.asarray(getprecision(total))
recall = np.asarray(getrecall(total,binning))
f1score = 2*(precision*recall)/(precision+recall)   
score = []
score.append(f1score.mean())

for i in range(2,N+1):
	print(i)
	c = c + coassociation_matrix(i,meta)
	newc = c/i
	trange = np.arange(0.1,1.0,0.1)
	fscore = []
	for t in trange:
		print(t)
		binning_test = compresult(c,-t,meta)
		sizedic = getsize(binning_test,contigsize)
		binning = cleanbin(binning_test)
		total = genomesize(binning, sizedic)
		precision = np.asarray(getprecision(total))
		recall = np.asarray(getrecall(total,binning))
		f1score = 2*(precision*recall)/(precision+recall)
		fscore.append(f1score.mean())
	score.append(max(fscore))

plt.plot(range(1,N+1),score,label = 'f1 score')
plt.xlabel('replication')
plt.ylabel('F1 score')
plt.title('Score VS Replication number')
plt.grid(True)
plt.savefig("scoreandrep.png")










