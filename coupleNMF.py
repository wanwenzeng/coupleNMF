import numpy as np
from sklearn.decomposition import NMF
from sklearn.feature_selection import  SelectFdr,SelectPercentile,f_classif
from numpy import linalg as LA
import math
import argparse
import itertools
import scipy.io as scio
import pandas as pd
import time
import scipy.stats as stats
from statsmodels.stats.weightstats import ttest_ind
from scipy  import sparse

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
parser.add_argument('-E_symbol', type=argparse.FileType('r'), help='the location of E gene symbol matrix')
parser.add_argument('-P_symbol', type=argparse.FileType('r'), help='the location of Peak symbol matrix')
parser.add_argument('-pe', type=argparse.FileType('r'), help='the location of pre-calculated peak-gene interactions')
parser.add_argument('-lambda1', dest='lambda1', type=float, help='lambda1, hyperparameters to control the term NMF for E')
parser.add_argument('-lambda2', dest='lambda2', type=float, help='lambda2, hyperparameters to control the coupled term')

args = parser.parse_args()

rep=50

print "Loading data..."

K=args.k
PeakO = np.loadtxt(args.PeakO)

E     = np.loadtxt(args.E)
E_symbol = []	
E_symbol = [line.strip() for line in args.E_symbol]

P_symbol = []
P_symbol = [line.strip() for line in args.P_symbol]

A = np.zeros((E.shape[0],PeakO.shape[0]))
for line in args.pe:
	data = line.strip().split()
	pindex = P_symbol.index(data[0])
	eindex = E_symbol.index(data[1])
	temp1 = float(data[3])
	if temp1 <0:
		temp1 = 0
	temp2 = float(data[2])
	A[eindex,pindex] = math.exp(-temp2/30000)*temp1

E_symbol = np.asarray(E_symbol)
P_symbol = np.asarray(P_symbol)
E        = pd.DataFrame(E)
PeakO    = pd.DataFrame(PeakO)
E     = quantileNormalize(E) 
PeakO = quantileNormalize(PeakO)

print "Initializing non-negative matrix factorization for E..."
E[E>10000] = 10000
X = np.log(1+E)

err1=np.zeros(rep)
for i in range(0,rep):
        model = NMF(n_components=K, init='random', random_state=i,solver='cd',max_iter=50)
        W20 = model.fit_transform(X)
        H20 = model.components_
        err1[i]=LA.norm(X-np.dot(W20,H20),ord = 'fro')

model = NMF(n_components=K, init='random', random_state=np.argmin(err1),solver='cd',max_iter=1000)
W20 = model.fit_transform(X)
H20 = model.components_
S20=np.argmax(H20,0)

print "Initializing non-negative matrix factorization for PeakO..."
PeakO = np.log(PeakO+1)
err=np.zeros(rep)
for i in range(0,rep):
        model = NMF(n_components=K, init='random', random_state=i,solver='cd',max_iter=50)
        W10 = model.fit_transform(PeakO)
        H10 = model.components_
        err[i]=LA.norm(PeakO-np.dot(W10,H10),ord = 'fro')

model = NMF(n_components=K, init='random', random_state=np.argmin(err),solver='cd',max_iter=1000)
W10 = model.fit_transform(PeakO)
H10 = model.components_
S10=np.argmax(H10,0)

print "Selecting differentially expressed genes..."
p2 = np.zeros((X.shape[0],K))
for i in range(K):
	for j in range(X.shape[0]):
		statistic, p2[j,i],df  = ttest_ind(X.ix[j,S20==i], X.ix[j,S20!=i] ,alternative='smaller')

WP2 = np.zeros((W20.shape))
p2[np.isnan(p2) ] = 1
scores = -np.log10(p2)
temp = int(len(E_symbol)/20)
for i in range(K):
	indexs = scores[:,i].argsort()[-temp:][::-1]
	WP2[indexs,i] = 1

print "Selecting differentially open peaks..."
p1 = np.zeros((PeakO.shape[0],K))
for i in range(K):
        for j in range(PeakO.shape[0]):
                statistic, p1[j,i],df  = ttest_ind(PeakO.ix[j,S10==i], PeakO.ix[j,S10!=i] ,alternative='smaller')

WP1 = np.zeros((W10.shape))
p1[np.isnan(p1) ] = 1
scores = -np.log10(p1)
temp = int(len(P_symbol)/20)
for i in range(K):
        indexs = scores[:,i].argsort()[-temp:][::-1]
        WP1[indexs,i] = 1

perm = list(itertools.permutations(range(K)))
score = np.zeros(len(perm))
for i in range(len(perm)):
        score[i] = np.trace(np.dot(np.dot(np.transpose(WP2),A),WP1))

match = np.argmax(score)
W20 = W20[:,perm[match]]
H20 = H20[perm[match],:]
S20=np.argmax(H20,0)

print "Initializing hyperparameters lambda1, lambda2 and mu..."
lambda10 = pow(LA.norm(X-np.dot(W20,H20),ord = 'fro'),2)/pow(LA.norm(PeakO-np.dot(W10,H10),ord = 'fro'),2)
lambda20 = pow(np.trace(np.dot(np.dot(np.transpose(W20),A),W10)),2)/pow(LA.norm(PeakO-np.dot(W10,H10),ord = 'fro'),2)
if type(args.lambda1) == type(None) and type(args.lambda2) == type(None):
	set1=[lambda10*pow(5,0),lambda10*pow(5,1),lambda10*pow(5,2),lambda10*pow(5,3),lambda10*pow(5,4)]
	set2=[lambda20*pow(5,-4),lambda20*pow(5,-3),lambda20*pow(5,-2),lambda20*pow(5,-1),lambda20*pow(5,0)]
elif type(args.lambda1) == type(None):
	set1=[lambda10*pow(5,0),lambda10*pow(5,1),lambda10*pow(5,2),lambda10*pow(5,3),lambda10*pow(5,4)]
	set2=[args.lambda2]
elif type(args.lambda2) == type(None):
        set1=[args.lambda1]
        set2=[lambda20*pow(5,-4),lambda20*pow(5,-3),lambda20*pow(5,-2),lambda20*pow(5,-1),lambda20*pow(5,0)]

else:
	set1=[args.lambda1*lambda10]
	set2=[args.lambda2*lambda20]


mu      = 1
eps = 0.001
detr = np.zeros((len(set1),len(set2)))
detr1 = np.zeros((len(set1),len(set2)))
S1_all = np.zeros((len(set1)*len(set2),PeakO.shape[1]))
S2_all = np.zeros((len(set1)*len(set2),E.shape[1]))
P_all = np.zeros((len(set1)*len(set2),K,PeakO.shape[0]))
E_all = np.zeros((len(set1)*len(set2),K,E.shape[0]))
P_p_all = np.zeros((len(set1)*len(set2),K,PeakO.shape[0]))
E_p_all = np.zeros((len(set1)*len(set2),K,E.shape[0]))
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
		print lambda1,lambda2
	
		
		print "Iterating coupleNMF..."
		maxiter   = 500
		err       = 1
		terms     = np.zeros(maxiter)
		it        = 0
		terms[it] = lambda1*pow(LA.norm(X-np.dot(W2,H2),ord = 'fro'),2)+pow(LA.norm(PeakO-np.dot(W1,H1),ord = 'fro'),2)+lambda2*pow(np.trace(np.dot(np.dot(np.transpose(W2),A),W1)),2)+mu*(pow(LA.norm(W1,ord = 'fro'),2)+pow(LA.norm(W2,ord = 'fro'),2))
		while it < maxiter-1 and err >1e-6:
			it  = it +1
			T1 = 0.5*lambda2*np.dot(np.transpose(A),W2)
			T1[T1<0] = 0
			W1  = W1*np.dot(PeakO,np.transpose(H1))/(eps+np.dot(W1,np.dot(H1,np.transpose(H1)))+0.5*mu*W1)
			H1  = H1*(np.dot(np.transpose(W1),PeakO))/(eps+np.dot(np.dot(np.transpose(W1),W1),H1))
                        T2 = 0.5*(lambda2/lambda1+eps)*np.dot(A,W1)
                        T2[T2<0] = 0
			W2  = W2*(np.dot(X,np.transpose(H2))+T2)/(eps+np.dot(W2,np.dot(H2,np.transpose(H2)))+0.5*mu*W2)
			H2  = H2*(np.dot(np.transpose(W2),X)/(eps+np.dot(np.dot(np.transpose(W2),W2),H2)))
			m1  = np.zeros((K,K))
			m2  = np.zeros((K,K))
			for z in range(K):
				m1[z,z] = LA.norm(H1[z,:])
				m2[z,z] = LA.norm(H2[z,:])
			
			W2  = np.dot(W2,m2)
			W1  = np.dot(W1,m1)
			H1  = np.dot(LA.inv(m1),H1)
			H2  = np.dot(LA.inv(m2),H2)
			
			terms[it] = lambda1*pow(LA.norm(X-np.dot(W2,H2),ord = 'fro'),2)+pow(LA.norm(PeakO-np.dot(W1,H1),ord = 'fro'),2)+lambda2*pow(np.trace(np.dot(np.dot(np.transpose(W2),A),W1)),2)+mu*(pow(LA.norm(W1,ord = 'fro'),2)+pow(LA.norm(W2,ord = 'fro'),2))
			err = abs(terms[it]-terms[it-1])/abs(terms[it-1])
 
		S2=np.argmax(H2,0)
		S1=np.argmax(H1,0)

		p2 = np.zeros((X.shape[0],K))
		for i in range(K):
        		for j in range(X.shape[0]):
                		statistic, p2[j,i],df  = ttest_ind(X.ix[j,S2==i], X.ix[j,S2!=i] ,alternative='smaller')

		WP2 = np.zeros((W2.shape))
		p2[np.isnan(p2) ] = 1
		scores = -np.log10(p2)
		temp = int(len(E_symbol)/20)
		for i in range(K):
        		indexs = scores[:,i].argsort()[-temp:][::-1]
        		WP2[indexs,i] = 1
			E_all[count,i,indexs] = 1
			E_p_all[count,i,indexs] = p2[indexs,i]

		p1 = np.zeros((PeakO.shape[0],K))
		for i in range(K):
        		for j in range(PeakO.shape[0]):
                		statistic, p1[j,i],df  = ttest_ind(PeakO.ix[j,S1==i], PeakO.ix[j,S1!=i] ,alternative='smaller')

		WP1 = np.zeros((W1.shape))
		p1[np.isnan(p1) ] = 1
		scores = -np.log10(p1)
		temp = int(len(P_symbol)/20)
		for i in range(K):
        		indexs = scores[:,i].argsort()[-temp:][::-1]
        		WP1[indexs,i] = 1
			P_all[count,i,indexs] = 1
			P_p_all[count,i,indexs] = p1[indexs,i]

        	T = np.dot(np.dot(np.transpose(WP2),A),WP1)
		temp = np.sum(np.sum(T))*np.diag(1/np.sum(T,axis=0))*T*np.diag(1/np.sum(T,axis=1))
		detr1[x,y] = np.trace(temp)
		detr[x,y] = np.trace(T)
		S1_all[count] = S1
		S2_all[count] = S2
		count = count + 1
			
[i,j] = npmax(detr)
print "Score is :"+ str(detr1[i,j]/K)
print "If the score >=1, the clustering matching for scRNA-seq and scATAC-seq is well. Otherwise, we sugguest to tune the parameters."

index = detr.argmax()
S1_final = S1_all[index,:]+1
S2_final = S2_all[index,:]+1
E_final  = E_all[index,:,:]
P_final  = P_all[index,:,:]
E_p_final  = E_p_all[index,:,:]
P_p_final  = P_p_all[index,:,:]

fout1 = open("scATAC-result.txt","w")
fout2 = open("scRNA-result.txt","w")
fout3 = open("cluster-specific-peaks-genes-pairs.txt","w")

print S1_final
print S2_final

for item in S1_final:
	fout1.write(str(item)+"\t")
fout1.write("\n")


for item in S2_final:
        fout2.write(str(item)+"\t")
fout2.write("\n")

for i in range(K):
	temp = np.dot(np.reshape(E_final[i,:],(E.shape[0],1)),np.reshape(P_final[i,:],(1,PeakO.shape[0])))*A
	p, q = np.nonzero(temp)
	for j in range(len(p)):
		fout3.write("cluster "+str(i+1)+": "+E_symbol[p[j]]+"\t"+P_symbol[q[j]]+"\t"+str(E_p_final[i,p[j]])+"\t"+str(P_p_final[i,q[j]])+"\n")

