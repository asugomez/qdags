#!/usr/bin/env python

# toma cada archivo que termina en k1000 para hacer las estadisticas

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import seaborn as sns

alg_label = ["G. LOUDS Backtrack.",
             "G. LOUDS Op. Order",
             "G. DFUDS Backtrack."]#,
            #"G. DFUDS Op. Order",
            #  "R. LOUDS Backtrack.", "R. LOUDS Op. Order",
            # "R. DFUDS Backtrack.", "R. DFUDS Op. Order"]
queries_label = ["j3","j4","p2","p3","p4","s1","s2","s3","s4","t2","t3","t4","ti2","ti3","ti4","tr1","tr2"]
queries_title = ["J3","J4","P2","P3","P4","S1","S2","S3","S4","T2","T3","T4","Ti2","Ti3","Ti4","Tr1","Tr2"]

# box plots for each query
data = []
rows, cols = 17, 9
row_template = [0] * cols
data = [row_template[:] for _ in range(rows)]
# box plots for each query
# TDO hacer un for para cada query y función y tamaño k
data = []
# Replication with *
rows, cols = 17, 3 #9 TODO: change this
row_template = [0] * cols
data = [row_template[:] for _ in range(rows)]

for type_fun in [0]:#,1]:
    for k in [1000]: # TODO: maybe only k=1000 ?
        for i,query in enumerate(queries_label): #[j3,j4,p2,p3,p4,s1,s2,s3,s4,t2,t3,t4,ti2,ti3,ti4,tr1,tr2]:
            file = f"{query}-f{type_fun}-k{k}.txt"
            print(file)
            partialLoudsBack = np.loadtxt(f'../partial/louds/backtracking/{file}', dtype=float)
            partialLoudsNon = np.loadtxt(f'../partial/louds/backtracking/{file}', dtype=float)# np.loadtxt(f'../partial/louds/nonFixedQueue/{file}', dtype=float)
            partialDfudsBack = np.loadtxt(f'../partial/louds/backtracking/{file}', dtype=float)#np.loadtxt(f'../partial/dfuds/backtracking/{file}', dtype=float)
            #partialDfudsNon = np.loadtxt(f'../partial/louds/backtracking/{file}', dtype=float)#np.loadtxt(f'../partial/dfuds/nonFixedQueue/{file}', dtype=float)

            #rankedLoudsBack = np.loadtxt(f'../partial/louds/backtracking/{file}', dtype=float)#np.loadtxt(f'../ranked/louds/backtracking/{file}', dtype=float)
            #rankedLoudsNon = np.loadtxt(f'../partial/louds/backtracking/{file}', dtype=float)#np.loadtxt(f'../ranked/louds/nonFixedQueue/{file}', dtype=float)
            #rankedDfudsBack = np.loadtxt(f'../partial/louds/backtracking/{file}', dtype=float)#np.loadtxt(f'../ranked/dfuds/backtracking/{file}', dtype=float)
            #rankedDfudsNon = np.loadtxt(f'../partial/louds/backtracking/{file}', dtype=float)#np.loadtxt(f'../ranked/dfuds/nonFixedQueue/{file}', dtype=float)

            #traditional = np.loadtxt(f'../partial/louds/backtracking/{file}', dtype=float)#np.loadtxt(f'../all/{query}.txt',dtype=float)

            # Convertir a floats explícitamente (aunque np.loadtxt ya debería cargar como float)
            partialLoudsBack = partialLoudsBack.astype(float)
            partialLoudsNon = partialLoudsNon.astype(float)
            partialDfudsBack = partialDfudsBack.astype(float)

            data[i] = [partialLoudsBack, partialLoudsNon, partialDfudsBack]#,
            #[traditional,
            #           partialLoudsBack, partialLoudsNon, partialDfudsBack]#,
            #partialDfudsNon, rankedLoudsBack, rankedLoudsNon, rankedDfudsBack, rankedDfudsNon]

            #for j,alg in enumerate(data[i]):
            #   if len(alg) < 50: # Rellenar con NaN si hay menos elementos
            #      data[i][j] = np.pad(alg, (0, 50 - len(alg)), constant_values=np.nan)


# colors:
# https://matplotlib.org/stable/gallery/color/named_colors.html
#colors = ['hotpink','red','peru']#, 'lightsalmon', 'gold', 'dodgerblue', 'darkturquoise', 'mediumspringgreen', 'lime']
colors = ['red', 'peru','lightsalmon']#, 'gold', 'dodgerblue', 'darkturquoise', 'mediumspringgreen', 'lime']
# Crear la figura
fig = plt.figure(figsize=(14,10))

# Crear un GridSpec con la distribución deseada
gs = gridspec.GridSpec(3, 6, figure=fig, hspace=0.1, wspace=0.2)

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


for i,query in enumerate([j3,j4,p2,p3,p4,s1,s2,s3,s4,t2,t3,t4,ti2,ti3,ti4,tr1,tr2]):
    for j, alg in enumerate(alg_label):
        bp = query.boxplot(data[i][j], positions=[j],whis=1, widths=0.4, patch_artist=True,
                           boxprops=dict(facecolor=colors[j],color=colors[j], label=alg_label[j]),
                           medianprops=dict(color='blueviolet'))
    query.set_yscale('log')
    query.set_title(queries_title[i])
    query.axes.get_xaxis().set_visible(False)
    if(query != j3 and query != s1 and query != t4):
        query.axes.get_yaxis().set_visible(False)


handles, labels = j3.get_legend_handles_labels()
color_legend.legend(handles, labels, loc='center', fontsize='x-small')


plt.savefig('/Users/asugomez/Desktop/Magister/Tesis/tesisQdags/imagenes/all_queries_k_1000.pdf')

plt.show()