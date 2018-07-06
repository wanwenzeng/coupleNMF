import numpy as np
from sklearn.decomposition import NMF
from sklearn.feature_selection import  SelectFdr,SelectPercentile,f_classif
from numpy import linalg as LA
import argparse
import itertools
import scipy.io as scio
import pandas as pd
import time
import scipy.stats as stats
from scipy.sparse import csr_matrix

def quantileNormalize(df_input):
    df = df_input.copy()
    #compute rank
    dic = {}
    for col in df:
        dic.update({col : sorted(df[col])})
    sorted_df = pd.DataFrame(dic)
    rank = sorted_df.mean(axis = 1).tolist()
    #sort
    for col in df:
        t = np.searchsorted(np.sort(df[col]), df[col])
        df[col] = [rank[i] for i in t]
    return df

def npmax(array):
    arrayindex = array.argmax(1)
    arrayvalue = array.max(1)
    i = arrayvalue.argmax()
    j = arrayindex[i]
    return i, j



parser = argparse.ArgumentParser(description='coupleNMF for joint clustering scRNA-seq and scATAC-seq.')
parser.add_argument('-k', dest='k', type=int, default=2, help='the number of clusters')
parser.add_argument('-E', type=argparse.FileType('r'), help='the location of singlecell expression E matrix')
parser.add_argument('-PeakO',type=argparse.FileType('r'), help='the location of singlecell ATAC-seq PeakO matrix')
parser.add_argument('-REO', type=argparse.FileType('r'), help='the location of REO matrix')
parser.add_argument('-E_symbol', type=argparse.FileType('r'), help='the location of E gene symbol matrix')
parser.add_argument('-s', type=str, help='the species (human or mouse)')
parser.add_argument('-ref', type=str, help='the reference genome (mm9, mm10, hg19 and hg38)')

args = parser.parse_args()

rep=50

print "Loading data..."
K=args.k
PeakO = np.loadtxt(args.PeakO)
REO   = np.loadtxt(args.REO)

E     = np.loadtxt(args.E)
E_symbol = []	
E_symbol = [line.strip() for line in args.E_symbol]

A        = np.load("RE_TG/"+args.s+"/"+args.ref +"/A.npy")
A_symbol = np.load("RE_TG/"+args.s+"/"+args.ref +"/A_symbol.npy")

print A
A = A.toarray()
E_symbol = np.asarray(E_symbol)
A_symbol = np.asarray(A_symbol)

E        = pd.DataFrame(E)
REO      = pd.DataFrame(REO)
PeakO    = pd.DataFrame(PeakO)
E     = quantileNormalize(E) 
REO   = quantileNormalize(REO)
PeakO = quantileNormalize(PeakO)
print "Initializing non-negative matrix factorization for E..."
E[E>10000] = 10000
X = np.log(1+E)

err1=np.zeros(rep)
for i in range(0,rep):
        model = NMF(n_components=K, init='random', random_state=i,solver='cd',max_iter=20)
        W20 = model.fit_transform(X)
        H20 = model.components_
        err1[i]=LA.norm(X-np.dot(W20,H20),ord = 'fro')

model = NMF(n_components=K, init='random', random_state=np.argmin(err1),solver='cd',max_iter=20)
W20 = model.fit_transform(X)
H20 = model.components_

model = NMF(n_components=K, init='custom',solver='cd',max_iter=1000)
W20 = model.fit_transform(X,W=W20,H=H20)
H20 = model.components_
S20=np.argmax(H20,0)

print "Selecting dynamic genes..."
dynamic = np.var(X,axis = 1,ddof = 1)
temp = int(len(dynamic)/5)+2
dymindexs = dynamic.argsort()[-temp:]

dym = set(E_symbol[dymindexs]).intersection(set(A_symbol))
test = E_symbol.tolist()
indexs1  = []
for item in dym:
	indexs1.append(test.index(item))

print "Selecting differentially expressed genes..."
statistic, pvalue = f_classif(np.transpose(X.ix[indexs1,:]),S20)
pvalue[np.isnan(pvalue) ] = 1
scores = -np.log10(pvalue)
temp = int(len(E_symbol)/100)
indexs2 = scores.argsort()[-temp:][::-1]
temp = E_symbol[indexs1]
E_subsymbol = temp[indexs2]

print "Selecting coupled A matrix..."
A_subindex = []
E_subindex = []
for item in E_subsymbol:
	E_subindex.append(E_symbol.tolist().index(item))
        A_subindex.append(A_symbol.tolist().index(item))

A = A[A_subindex]
REO = REO.ix[np.sum(np.abs(A),axis = 0)>0,:]
A = A[:,np.sum(np.abs(A),axis = 0)>0]

print "Initializing non-negative matrix factorization for PeakO..."
PeakO = np.log(PeakO+1)
err=np.zeros(rep)
for i in range(0,rep):
	model = NMF(n_components=K, init='random', random_state=i,solver='cd',max_iter=20)
	W10 = model.fit_transform(PeakO)
	H10 = model.components_
	err[i]=LA.norm(PeakO-np.dot(W10,H10),ord = 'fro')

model = NMF(n_components=K, init='random', random_state=np.argmin(err),solver='cd',max_iter=20)
W10 = model.fit_transform(PeakO)
H10 = model.components_

model = NMF(n_components=K, init='custom',random_state=np.argmin(err),solver='cd',max_iter=1000)
W10 = model.fit_transform(PeakO,W=W10,H=H10)
H10 = model.components_
S10=np.argmax(H10,0)

print "Initializing non-negative matrix factorization for REO..."
REO = np.log(1+REO)
SW10 = np.dot(REO,LA.pinv(H10))
SW10[SW10<0] = 0

print S10
print S20
time.sleep(5)

print "Initializing hyperparameters lambda1, lambda2 and mu..."
set1=[1,10,100,1000,10000]
set2=[0.0001,0.001,0.01,0.1]
#set1 = [1000]
#set2 = [1]
mu      = 1
eps = 0.001
detr = np.zeros((len(set1),len(set2)))
S1_all = np.zeros((len(set1)*len(set2),REO.shape[1]))
S2_all = np.zeros((len(set1)*len(set2),E.shape[1]))
print "Starting coupleNMF..."
count = 0
for x in range(len(set1)):
	for y in range(len(set2)):
		lambda1 = set1[x]
		lambda2 = set2[y]
		W1 = W10
		W2 = W20
		H1 = H10
		H2 = H20
		SW1 =SW10
		A = A
		print lambda1,lambda2
	
		perm = list(itertools.permutations(range(K)))
		score = np.zeros(len(perm))
		for i in range(len(perm)):
			temp = W2[:,perm[i]]
			score[i] = np.trace(np.dot(np.dot(np.transpose(temp[E_subindex,:]),A),SW1))
	
		match = np.argmax(score)
		W2 = W2[:,perm[match]]
		H2 = H2[perm[match],:]
		
		print "Iterating coupleNMF..."
		maxiter   = 1000
		err       = 1
		terms     = np.zeros(maxiter)
		it        = 0
		rate      = 1
		terms[it] = max(score)
		while it < maxiter-1 and err >1e-6:
			it  = it +1
			T1 = 0.5*lambda2*np.dot(np.transpose(A),W2[E_subindex,:])
			T1[T1<0] = 0
			W1  = W1*np.dot(PeakO,np.transpose(H1))/(eps+np.dot(W1,np.dot(H1,np.transpose(H1)))+0.5*mu*W1)
			H1  = H1*(np.dot(np.transpose(W1),PeakO)+rate*np.dot(np.transpose(SW1),REO))/(eps+np.dot(np.dot(np.transpose(W1),W1),H1)+rate*np.dot(np.dot(np.transpose(SW1),SW1),H1))
			SW1 = SW1*np.asarray(np.sqrt((np.dot(REO,np.transpose(H1))+T1)/(eps+np.dot(SW1,np.dot(np.dot(np.transpose(SW1),REO),np.transpose(H1))))))
			T2  = np.zeros((W2.shape))
                        T2[E_subindex,:] = 0.5*(lambda2/lambda1+eps)*np.dot(A,SW1)
                        T2[T2<0] = 0
			W2  = W2*(np.dot(X,np.transpose(H2))+T2)/(eps+np.dot(W2,np.dot(H2,np.transpose(H2)))+0.5*mu*W2)
                        H2  = H2*(np.dot(np.transpose(W2),X)/(eps+np.dot(np.dot(np.transpose(W2),W2)+mu*np.ones((K,K)),H2)))
			m1  = np.zeros((K,K))
			m2  = np.zeros((K,K))
			for z in range(K):
				m1[z,z] = LA.norm(H1[z,:])
				m2[z,z] = LA.norm(H2[z,:])
			
			SW1 = np.dot(SW1,m1)
			W2  = np.dot(W2,m2)
			W1  = np.dot(W1,m1)
			H1  = np.dot(LA.inv(m1),H1)
			H2  = np.dot(LA.inv(m2),H2)
			
			m1  = np.zeros((K,K))
                        m2  = np.zeros((K,K))
			for z in range(K):
			        m1[z,z] = LA.norm(H1[z,:])
			        m2[z,z] = LA.norm(H2[z,:])
			
			SW1_n = np.dot(SW1,m1)
			W2_n  = np.dot(W2,m2)
			
			score = np.zeros(len(perm))
			for i in range(len(perm)):
			        temp = W2_n[:,perm[i]]
				temp1 = np.transpose(temp[E_subindex,:])
				temp2 = np.transpose(np.dot(A,SW1_n))
				for z in range(K):	
					temp3 = np.corrcoef(temp1[z,:],temp2[z,:])
					score[i] = score[i] + temp3[0,1]
			
			match = np.argmax(score)
			W2 = W2[:,perm[match]]
			H2 = H2[perm[match],:]
			terms[it]  = np.max(score)
			err = abs(terms[it]-terms[it-1])/abs(terms[it-1])
				
				
		for i in range(5):
			H1 = H1*np.dot(np.transpose(W1),PeakO)/(np.dot(np.dot(np.transpose(W1),W1),H1)+eps)
		
		m1  = np.zeros((K,K))
		m2  = np.zeros((K,K))
		for z in range(K):
			m1[z,z] = LA.norm(H1[z,:])
			m2[z,z] = LA.norm(H2[z,:])
		
		W1  = np.dot(W1,m1)
		SW1 = np.dot(SW1,m1)
		W2  = np.dot(W2,m2)
		H1  = np.dot(LA.inv(m1),H1)
		H2  = np.dot(LA.inv(m2),H2)
		S2=np.argmax(H2,0)
		S1=np.argmax(H1,0)
		
		S1_all[count] = S1
		S2_all[count] = S2
		count = count + 1
		
		SW2 = np.zeros((len(E_subindex),K))
		SW22 = np.zeros((len(E_subindex),K))
		for i in range(K):
			SW2[:,i] = np.transpose(np.mean(X.ix[E_subindex,S2 == i],axis =1))
			SW22[:,i] = np.transpose(np.dot(A,np.mean(REO.ix[:,S1 == i],axis = 1)))
		
		cmat = np.zeros((K,K))
		for p in range(K):
			for q in range(K):
				temp = np.corrcoef(SW2[:,p],SW22[:,q])
				cmat[p,q] = temp[0,1]

		detr[x,y] = LA.det(cmat)
		
[i,j] = npmax(detr)

index = detr.argmax()
S1_final = S1_all[index,:]
S2_final = S2_all[index,:]


fout1 = open("scATAC-result.txt","w")
fout2 = open("scRNA-result.txt","w")

print S1_final
print S2_final

for item in S1_final:
	fout1.write(str(item)+"\t")
fout1.write("\n")


for item in S2_final:
        fout2.write(str(item)+"\t")
fout2.write("\n")

