import matplotlib
matplotlib.use('Agg')
import os
import sys
import getopt
import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt
import warnings
from scipy.stats import *
import math
warnings.filterwarnings('ignore')
from functools import reduce
from sklearn.ensemble import IsolationForest
from sklearn.neighbors import KernelDensity
from sklearn.cluster import DBSCAN
from multiprocessing import Pool, cpu_count



plt.close()
inp_= sys.argv[1] 
outp_= sys.argv[2]
print('input:', inp_)
print('output:', outp_)

input_ = inp_#'/hpcnfs/scratch/PGP/SCMseq/AML4/merged_enriched/BC_demux'
output_ =  outp_#'/hpcnfs/scratch/PGP/SCMseq/AML4/merged_enriched/genotype_imputation'
if not os.path.exists(output_):
    os.mkdir(output_)

output=output_

df_all = pd.read_csv(os.path.join(input_, 'df_bc_all.csv'))

x1=df_all.columns[2:]

gene_list=[]
for x in x1:
    if x.split("_")[0] not in gene_list:
        gene_list.append(x.split("_")[0])
gene_list

def mutated_cell_confidence(n_mutated_reads, accuracy=0.95):
        return 1 - (1-accuracy)**n_mutated_reads

def process_bkpoints(bkpoints, stat_list): 
    processed = []
    for p, s in zip(bkpoints, stat_list): 
        if s =='neither': 
            processed.append(p)
        elif s=='max': 
            processed = processed + [p, p]
    return processed

def find_stationary(f, unique_vals, th): 
    stat_list = []
    for val in unique_vals: 
        if (f(val-th)-f(val))*(f(val+th)-f(val)) <=0: 
            stat_list.append('neither')
        else:
            if f(val-th) < f(val): 
                stat_list.append('max')
            else: 
                stat_list.append('min')
    return stat_list
            

def compute_likelihood(value, kde, grid, kde_grid):
    pdf_val = kde.evaluate(value)[0]
    a = list(filter(lambda x: abs(x-pdf_val)<0.001, kde_grid))
    b = [grid[kde_grid.index(y)] for y in a]
    consec_diff = [y-x for x, y in zip(b[:-1], b[1:])]
    big_vals = list(filter(lambda x: x>=0.001, consec_diff))
    unique_vals = [b[0]]+[b[consec_diff.index(x)+1] for x in big_vals]
    stat_list = find_stationary(kde.evaluate, unique_vals, th=0.005)
    proc_bkpts = process_bkpoints(unique_vals, stat_list)
    bkpoints = [-10]+ proc_bkpts+[10]
    integration_intervals = [(bkpoints[2*i], bkpoints[2*i+1]) for i in range(int(len(bkpoints)/2))]
    likelihood = sum([kde.integrate_box_1d(x[0], x[1]) for x in integration_intervals])
    return likelihood

def compute_confidence1(n_wt, kde, grid, kde_grid, approx=True, alpha=0.95, verbose=False): 
    if n_wt == 1: 
        return 0 
    else:
        li = list(map(lambda x: compute_likelihood(round(x/n_wt, 3), kde, grid, kde_grid), 
                      range(1, n_wt)))
        prob = list(map(lambda x: binom(n_wt, 1-alpha).pmf(x), 
                         range(1, n_wt)))
        if prob !=[]: 

            if verbose: 
                for x, y, z in zip(li, prob, range(1, n_wt)): 
                    print(z, x, y)
                print([x*y for x, y in zip(li, prob)])
                print(sum(prob))


            return   1-sum([x*y for x, y in zip(li, prob)])/sum(prob)
        else: 
            return 0

def chunks(l, n):
    """Yield n number of sequential chunks from l."""
    d, r = divmod(len(l), n)
    for i in range(n):
        si = (d+1)*(i if i < r else r) + d*(0 if i < r else i - r)
        yield l[si:si+(d+1 if i < r else d)]

def compute_confidence_map(li): 
    return [compute_confidence(x, kde,grid,kde_grid,alpha=accuracy) for x in li]

def confidence_observable(array, nmax, accuracy=0.95):
    kde = gaussian_kde(array)
    grid = np.linspace(0, 1, 10000)
    kde_grid = list(kde.evaluate(grid))
    conf_sample  = [compute_confidence1(x, kde,grid,kde_grid,alpha=accuracy) for x in range(1, nmax)]
    
    return conf_sample

def how_many_errors_needed(n_reads, clusters_mean_af, clusters_mean_std): 
    d = {}
    for l, v in clusters_mean_af.items(): 
        for i in range(1, n_reads+1): 
            if abs((n_reads-i)/n_reads-v) <= clusters_mean_std[l]: 
                break
        if abs((n_reads-i)/n_reads-v) <= clusters_mean_std[l]:
            d[l] = i
        
    return d

def compute_confidence(n_wt, clusters_mean_af, clusters_mean_std, scores, alpha=0.95, verbose=False):
    if n_wt==1: 
        return 0
    else: 
        how_many_errors = how_many_errors_needed(n_wt, clusters_mean_af, clusters_mean_std)
        probs = {l:binom(n_wt, 1-alpha).pmf(v) for l, v in how_many_errors.items()}
        confidence = 1 - sum(scores[i]*probs[i] for i in probs.keys())
        return confidence


for MUT in gene_list:
    plt.ioff()
    df = df_all.drop(list(filter(lambda x: MUT not in x and 'barcode' not in x, df_all.columns)), axis=1)
    df.columns = ['barcode', 'alt', 'ref', 'mis']
    print(df)
    df['Tot'] = df.alt + df.ref + df.mis
    df = df.sort_values(by='Tot')
    df = df[df.Tot!=0]
    df['accuracy'] = (df.alt + df.ref)/df.Tot
    accuracy = np.mean(df.accuracy[df.accuracy>0])
    df = df[df.accuracy>=0.9]
    #####wt
    wt = df[df.alt == 0]
    ####MUT
    mut = df[df.alt > 0]
    #####compute mutated cell imputation confidence
    mut = mut.sort_values(by='alt')
    accuracy = 0.9
    mut['Confidence'] = mut.alt.apply(lambda x: mutated_cell_confidence(x, accuracy=accuracy))
    ####plot alt_confidence
    plt.scatter(mut.alt, mut.Confidence)
    plt.title('Mutation: '+ str(MUT))
    plt.xlabel('Number of mutated reads')
    plt.ylabel('Confidence')
    plt.savefig(os.path.join(output, str(MUT)+'__mut_Tot.pdf'))
    plt.close()
    #plt.show()
    #####save to csv
    mut.to_csv(os.path.join(output, str(MUT)+'_mutated_NoFilter.csv'))
    #####compute wildtype cell imputation conficence
    mut_and_wt = df[np.logical_and(df.alt > 0, df.ref > 0)]
    if len(mut_and_wt)< 10 :
        print("len is less than 10 BC")
        continue
    #sns.distplot(mut_and_wt.Tot, bins=10)
    #plt.savefig(os.path.join(output, f'{MUT}__mut_and_wt_Tot.pdf'))
    #plt.show()
    ### compute AF distributions only based on the cells which have at least minimum_n_reads reads
    minimum_n_reads = np.quantile(mut_and_wt.Tot,0.5)
    print(minimum_n_reads)
    mut_and_wt = mut_and_wt[mut_and_wt.Tot>minimum_n_reads]
    mut_and_wt['AF'] = mut_and_wt.ref/mut_and_wt.Tot
    if mut_and_wt['AF'].size < 2:
        continue
    kde = gaussian_kde(mut_and_wt['AF'])
    grid = np.linspace(0, 1, 10000)
    kde_grid = list(kde.evaluate(grid))

    ######plot AF
    sns.distplot(mut_and_wt['AF'], bins=10,kde_kws = {'bw' : 0.5})
    plt.plot(grid, kde_grid)
    plt.suptitle(str(MUT))
    plt.savefig(os.path.join(output, str(MUT)+'__mut_and_wt_AF.pdf'))
    plt.close()
    fitted_log_dens = KernelDensity(kernel='cosine').fit(np.array(mut_and_wt['AF']).reshape(-1, 1))
    log_dens = fitted_log_dens.score_samples(np.linspace(0, 1, 1000).reshape(-1, 1))
    ###plot AF_scatter
    plt.scatter(np.linspace(0, 1, 1000), np.exp(log_dens), s=0.1)
    plt.scatter(mut_and_wt['AF'], np.zeros(len(mut_and_wt['AF'])), marker='*')
    plt.savefig(os.path.join(output, str(MUT)+'__mut_and_wt_AF_scatter.pdf'))
    #plt.show()
    plt.close()
    #plt.show()
   # nmax = max(wt.Tot)+1
   # conf_sample_ref = confidence_observable(mut_and_wt.ref/mut_and_wt.Tot, nmax, accuracy)
   # conf_sample_alt = confidence_observable(mut_and_wt.alt/mut_and_wt.Tot, nmax, accuracy)
   # conf_sample = np.maximum(np.array(conf_sample_ref), np.array(conf_sample_alt))
  #  isof = IsolationForest(random_state=0, n_estimators=100)
  #  outliers = isof.fit_predict(np.reshape(conf_sample, (-1, 1)))
  #  plt.scatter(range(1, nmax), conf_sample, c=outliers)
  #  plt.savefig(os.path.join(output, str(MUT)+'__wt.conf_sample.pdf'))
  #  plt.close()
    ##
  #  confidence = dict(zip(range(1, nmax+1), conf_sample))
  #  wt['Confidence'] = wt.Tot.apply(lambda x: confidence[x])
    ######new approach
    labels = DBSCAN(eps=0.03).fit_predict(np.array(mut_and_wt['AF']).reshape(-1, 1))
    ###plot AF by DBSCAN
    #plt.scatter(mut_and_wt['AF'], labels)
    #plt.show()
    #plt.close()
    ####
    unique_labels = set(labels)
    af_array = np.array(mut_and_wt['AF'])
    scores = {l: sum(list(map(lambda x: x==l, labels)))/len(labels) for l in unique_labels}
    clusters = {l: [af_array[i] for i in range(len(mut_and_wt['AF'])) if labels[i]==l]for l in unique_labels}
    clusters_mean_af = {k:np.mean(v) for k, v in clusters.items()}
    clusters_mean_std = {k:np.std(v) for k, v in clusters.items()}
    how_many_errors_needed(20, clusters_mean_af, clusters_mean_std)
    nmax = max(wt.Tot)+1
    confidences = {i: compute_confidence(i, clusters_mean_af,clusters_mean_std, scores, alpha=0.5) for i in range(1, nmax)}
    wt['Confidence'] = wt.Tot.apply(lambda x: confidences[x])
    wt.to_csv(os.path.join(output, str(MUT)+'_wt_90acc_median.csv'))
    sns.distplot(wt.Confidence,bins=20,kde_kws = {'bw' : 0.5})
    plt.savefig(os.path.join(output, str(MUT)+'__wt.Confidence.pdf'))
    plt.close()
    gc.collect()

