#!/usr/bin/env python

# plot original algorithms for each query vs k (k=1,10,100,1000)
# toma el archivo results.csv para hacer las estadisticas

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import tikzplotlib


# In[3]:


alg_label = ["G. LOUDS Backtracking",
             "G. LOUDS Optimal Order",
             "G. DFUDS Backtracking",
              "Gradual DFUDS Optimal Order"]#,
             # "Ranked LOUDS Backtracking",
             # "Ranked LOUDS Optimal Order",
             # "Ranked DFUDS Backtracking",
             # "Ranked DFUDS Optimal Order"]
queries_label = ["j3","j4","p2","p3","p4","s1","s2","s3","s4","t2","t3","t4","ti2","ti3","ti4","tr1","tr2"]
queries_title = ["J3","J4","P2","P3","P4","S1","S2","S3","S4","T2","T3","T4","Ti2","Ti3","Ti4","Tr1","Tr2"]

# box plots for each query
j3,j4,p2,p3,p4,s1,s2,s3,s4,t2,t3,t4,ti2,ti3,ti4,tr1,tr2 = [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]
data = [j3,j4,p2,p3,p4,s1,s2,s3,s4,t2,t3,t4,ti2,ti3,ti4,tr1,tr2]
datasets =[[],[],[],[],[],[],[],[]]
for i,type_fun in enumerate([0]):#,1]):
    file = f"results-f{type_fun}-v2.csv"
    print(file)
    partialLoudsBack = pd.read_csv(f'../partial/louds/backtracking/{file}', delimiter=';')
    partialLoudsNon = pd.read_csv(f'../partial/louds/nonFixedQueue/{file}', delimiter=';')
    partialDfudsBack = pd.read_csv(f'../partial/dfuds/backtracking/{file}', delimiter=';')
    partialDfudsNon = pd.read_csv(f'../partial/dfuds/nonFixedQueue/{file}', delimiter=';')

    # rankedLoudsBack = pd.read_csv(f'../ranked/louds/backtracking/{file}', delimiter=';')
    # rankedLoudsNon = pd.read_csv(f'../ranked/louds/nonFixedQueue/{file}', delimiter=';')
    # rankedDfudsBack = pd.read_csv(f'../ranked/dfuds/backtracking/{file}', delimiter=';')
    # rankedDfudsNon = pd.read_csv(f'../ranked/dfuds/nonFixedQueue/{file}', delimiter=';')
    
    datasets[i] = [partialLoudsBack,
                   partialLoudsNon,
                   partialDfudsBack,
                   partialDfudsNon]#,
                   # rankedLoudsBack,
                   # rankedLoudsNon,
                   # rankedDfudsBack,
                   # rankedDfudsNon]
    #traditional = pd.read_csv(f'../original/results.csv',delimiter=';')
    # print(i)
    for j,query in enumerate(queries_label):
        data[j] = [partialLoudsBack[query], partialLoudsNon[query], partialDfudsBack[query], partialDfudsNon[query] ]#, rankedLoudsBack[query], rankedLoudsNon[query], rankedDfudsBack[query], rankedDfudsNon[query]]


variables = datasets[0][0].columns[1:]  # Excluir la columna 'k'


#colors = ['maroon', 'red', 'lightsalmon']#, 'gold', 'dodgerblue', 'darkturquoise', 'mediumspringgreen', 'lime']
colors = ['red', 'peru','lightsalmon', 'gold']#, 'dodgerblue', 'darkturquoise', 'mediumspringgreen', 'lime']

#fig, (j3,j4,p2,p3,p4) = plt.subplots(1,5, sharey=True)
#fig, ((j3,j4,p2,p3,p4),(s1,s2,s3,s4,t2,t3),(t4,ti2,ti3,ti4,tr1,tr2)) = plt.subplots(3, 6, layout='constrained', sharey=True)

# Crear la figura
fig = plt.figure(figsize=(10, 13))

# Crear un GridSpec con la distribución deseada
gs = gridspec.GridSpec(3, 6, figure=fig, hspace=0.3, wspace=0.8)

# Añadir los subplots a la figura
j3 = fig.add_subplot(gs[0, 0])
j4 = fig.add_subplot(gs[0, 1], sharex=j3, sharey=j3)
p2 = fig.add_subplot(gs[0, 2], sharex=j3, sharey=j3)
p3 = fig.add_subplot(gs[0, 3], sharex=j3, sharey=j3)
p4 = fig.add_subplot(gs[0, 4], sharex=j3, sharey=j3)
color_legend = fig.add_subplot(gs[0, 5]) # Space for the legend
color_legend.axis('off')  # Hide the axis

s1 = fig.add_subplot(gs[1, 0], sharex=j3, sharey=j3)
s2 = fig.add_subplot(gs[1, 1], sharex=j3, sharey=j3)
s3 = fig.add_subplot(gs[1, 2], sharex=j3, sharey=j3)
s4 = fig.add_subplot(gs[1, 3], sharex=j3, sharey=j3)
t2 = fig.add_subplot(gs[1, 4], sharex=j3, sharey=j3)
t3 = fig.add_subplot(gs[1, 5], sharex=j3, sharey=j3)

t4 = fig.add_subplot(gs[2, 0], sharex=j3, sharey=j3)
ti2 = fig.add_subplot(gs[2, 1], sharex=j3, sharey=j3)
ti3 = fig.add_subplot(gs[2, 2], sharex=j3, sharey=j3)
ti4 = fig.add_subplot(gs[2, 3], sharex=j3, sharey=j3)
tr1 = fig.add_subplot(gs[2, 4], sharex=j3, sharey=j3)
tr2 = fig.add_subplot(gs[2, 5], sharex=j3, sharey=j3)

j3.set_ylabel('Time (ms)')
s1.set_ylabel('Time (ms)')
t4.set_ylabel('Time (ms)')

t4.set_xlabel('k results')
ti2.set_xlabel('k results')
ti3.set_xlabel('k results')
ti4.set_xlabel('k results')
tr1.set_xlabel('k results')
tr2.set_xlabel('k results')

#type_fun=0



for i,query in enumerate([j3,j4,p2,p3,p4,s1,s2,s3,s4,t2,t3,t4,ti2,ti3,ti4,tr1,tr2]):
    for j, alg in enumerate(alg_label):
        query.plot(datasets[0][0]['k'], data[i][j], marker='', color=colors[j], label=alg_label[j])
    query.set_xscale('log')
    query.set_yscale('log')
    query.set_title(queries_title[i])
    query.axvline(x=10, color='gray', linestyle='--', linewidth=0.7)
    query.axvline(x=100, color='gray', linestyle='--', linewidth=0.7)
    # query.axvline(x=1000, color='gray', linestyle='--', linewidth=0.7)
    if(query != j3 and query != s1 and query != t4):
        query.axes.get_yaxis().set_visible(False)
    if(query != t4 and query != ti2 and query != ti3 and query != ti4 and query != tr1 and query != tr2):
        query.axes.get_xaxis().set_visible(False)

handles, labels = j3.get_legend_handles_labels()
color_legend.legend(handles, labels, loc='center', fontsize='11')


plt.savefig('/Users/asugomez/Desktop/Magister/Tesis/tesisQdags/imagenes/all_queries_k_1_10_100_1000.pdf')

plt.show()

