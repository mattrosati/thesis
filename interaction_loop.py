import pandas as pd
import numpy as np
import networkx as nx
import pdb
import matplotlib as mpl
import matplotlib.pyplot as plt
import nxviz as nz
from pyvis.network import Network
import markov_clustering as mc

#import data
suppressed = pd.read_csv("/Users/matteorosati/Desktop/MacMicking Lab/IFN_microarray/Data Analysis/s_entrez.csv", dtype={'entrezgene': object})
activated = pd.read_csv("/Users/matteorosati/Desktop/MacMicking Lab/IFN_microarray/Data Analysis/a_entrez.csv", dtype={'entrezgene': object})
bioGrid = pd.read_csv("/Users/matteorosati/Desktop/MacMicking Lab/IFN_microarray/BioPlex/BioPlex_interactionList_v4a.tsv", sep='\t')

#remove row names
suppressed = suppressed.iloc[:,2:]
activated = activated.iloc[:, 2:]
suppressed.set_index('entrezgene', inplace=True)
activated.set_index('entrezgene', inplace=True)

print(bioGrid.head(5))


interactions = bioGrid[["GeneA", "GeneB"]] #make edgelist
interactions = interactions.applymap(str)
initial = nx.from_pandas_edgelist(interactions, 'GeneA', 'GeneB', create_using = nx.DiGraph()) #make network from edgelist

print(len(list(initial.nodes())))
colnames = activated.columns.values


#add metadata to nodes
for i in list(initial.nodes()):
	# if i == '57172':
	# 	pdb.set_trace()
	if i in list(activated.index):
		behavior = dict(zip(list(colnames), list(activated.loc[i].values+(-1)*suppressed.loc[i].values)))
		initial.node[i]['behavior'] = behavior
	else:
		initial.node[i]['behavior'] = dict(zip(list(colnames), [0]*8))

initial_save = initial.copy()
for i in list(initial_save.nodes()):
	for j in range(0,8):
		if (i in list(activated.index)) or (i in list(suppressed.index)):
			activated_gene = list(activated.loc[i].values)
			suppressed_gene = list(suppressed.loc[i].values)
			if activated_gene[j] == 1:
				initial_save.node[i][colnames[j]] = 1
			elif suppressed_gene[j] == 1:
				initial_save.node[i][colnames[j]] = -1
			else:
			    initial_save.node[i][colnames[j]] = 0
		else:
			initial_save.node[i][colnames[j]] = 0
for i in list(initial_save.nodes()):
	del initial_save.node[i]['behavior']

print(initial.nodes['57172'])

nx.write_gml(initial_save, '/Users/matteorosati/Desktop/MacMicking Lab/IFN_microarray/Processed Data/network_raw.gml')

#think about pruning the network, 
# i.e. taking only the network genes in the matrices.
# this means make the last else statement a node removal

#pruned network
networkp = initial.copy()
for i in list(initial.nodes()):
	if (i not in list(activated.index)) and (i not in list(suppressed.index)):
		networkp.remove_node(i)

print(len(list(networkp.nodes())) - len(list(initial.nodes())) )

#207 genes removed

#degree histogram
deg_count = nx.degree(initial)
deg = [v for n, v in deg_count]
# plt.hist(deg)
# plt.show()


#pagerank and plot results
page = nx.pagerank_scipy(initial)
# plt.hist(list(page.values()))
# plt.show()


#degree centrality and plot
deg_cent = nx.degree_centrality(initial)
# plt.hist(list(deg_cent.values()))
# plt.show()

#all of these histograms are uninteresting


'''
find a way to quantify ISG-ness (but consider also that a gene can be both 
activated and suppressed among cell types). Need an ordering of the possible
vectors?
'''

network = initial.copy()

for n, d in network.nodes(data=True):
	network.node[n]['mISG'] = sum(list(network.nodes[n]['behavior'].values()))
	if (1 in list(network.nodes[n]['behavior'].values())) and (-1 in list(network.nodes[n]['behavior'].values())):
			network.node[n]['note'] = 'weird'

#add pr values to node metadata
for i in network.nodes():
	network.node[i]['PageRank'] = page[i]
	network.node[i]['deg_cen'] = deg_cent[i]

#build a label function to get gene symbol, add to network

symbol_entrez_a = bioGrid.loc[:,('GeneA', 'SymbolA')]
symbol_entrez_b = bioGrid.loc[:,('GeneB', 'SymbolB')]

new_columns = ["Gene", "Symbol"]
symbol_entrez_a.columns = new_columns
symbol_entrez_b.columns = new_columns

symbol_entrez = pd.concat([symbol_entrez_a, symbol_entrez_b], ignore_index=True)
symbol_entrez.applymap(str)
symbol_entrez.drop_duplicates(keep='first', inplace=True)
symbol_entrez.set_index('Gene', inplace=True)
symbol_entrez.index = symbol_entrez.index.astype(str)

for i in list(network.nodes()):
	network.node[i]['symbol'] = list(symbol_entrez.loc[i])[0]


#nx.write_gml(network, '/Users/matteorosati/Desktop/MacMicking Lab/IFN_microarray/Processed Data/network_pr_isg.gml')

'''
find closeness centrality of activated and suppressed genes for every cell type. 
Make box plots!

also, need better computer for cytoscape
'''

#cliques comparison between pruned and not pruned. 
#first need to generate undirected equivalents

undir = network.to_undirected()
undirp = networkp.to_undirected()

clength_p = [undirp.subgraph(c) for c in nx.connected_components(undirp)]
clength_np = [undir.subgraph(c) for c in nx.connected_components(undir)]
print(list(len(list(c.nodes())) for c in clength_p), list(len(list(c.nodes())) for c in clength_np))
#unpruned has a lot less unconnected components

#closeness centrality boxplots
close = nx.closeness_centrality(undir)

for i in undir.nodes():
	undir.node[i]['closeness'] = close[i]

#first define stimulation function
def abs_sum(g):
	sum = 0
	for n in g:
		sum += abs(n)
	return sum

df = pd.DataFrame(columns=['gene', 'c_cen', 'p_rank', 'cell_type', 'd_cen', 'label'])

for n in list(undir.nodes()):
	for j in colnames:
		row = [n, undir.nodes[n]['closeness'], undir.nodes[n]['PageRank'], j, undir.nodes[n]['deg_cen'], '']
		if undir.nodes[n]['behavior'][j] == 1:
			row[-1] = 'activated'
			df = df.append(dict(zip(df.columns.values, row)), ignore_index=True)
		elif undir.nodes[n]['behavior'][j] == -1:
			row[-1] = 'suppressed'
			df = df.append(dict(zip(df.columns.values, row)), ignore_index=True)

print(df.head(5))

df.to_csv('/Users/matteorosati/Desktop/MacMicking Lab/IFN_microarray/Processed Data/boxplots_data.csv', index=False)

df_analysis = pd.DataFrame(columns=['gene', 'c_cen', 'p_rank', 'd_cen', 'symbol', 'stimul', 'stimul_a', 'stimul_s', 'mISG'])

for n in list(undir.nodes()):
	s = abs_sum(list(undir.node[n]['behavior'].values()))
	s_a = 0
	s_s = 0
	for i in list(undir.node[n]['behavior'].values()):
		if i == 1:
			s_a += 1
		elif i == -1:
			s_s += 1
	row = [n, undir.nodes[n]['closeness'], undir.nodes[n]['PageRank'], undir.nodes[n]['deg_cen'], undir.nodes[n]['symbol'], s, s_a, s_s, undir.nodes[n]['mISG']]
	df_analysis = df_analysis.append(dict(zip(df_analysis.columns.values, row)), ignore_index=True)

df_analysis.to_csv('/Users/matteorosati/Desktop/MacMicking Lab/IFN_microarray/Data Analysis/bigdata.csv', index=False)


df_sub = pd.DataFrame(columns=['gene', 'cc_sub', 'cell_type', 'dc_sub', 'label'])


# for i in colnames:
# 	subgraph_a = undir.subgraph(df['cell_type' == i and 'label' == 'activated']['gene'])
# 	subgraph_s = undir.subgraph(df['cell_type' == i and 'label' == 'suppressed']['gene'])
# 	cc_a = nx.closeness_centrality(subgraph_a)
# 	cc_s = nx.closeness_centrality(subgraph_s)
# 	print('cc done', end='')
# 	dc_a = nx.degree_centrality(subgraph_a)
# 	dc_s = nx.degree_centrality(subgraph_s)
# 	print('dc done')


#find all subgraphs of neighbors and rank them by 
#how stimulated the components are
'''
discussion on ave_stimul: 2 possible ways:
1. Take average of stimulation over all nodes in cluster.
2. Take average of stimulation over just stimulated nodes in cluster.
While 1 will decrease the more unstimulated nodes there are, 2 will only
consider the stimulated ones, therefore making the unstimulated ones irrelevant.
2 will thus when ranked generate possible candidates of ISGs because the 
subgraphs with a high ave_stimul might have unstimulated genes in the graph.
Why do these genes stay unstimulated?
'''

#get stimulation values for all nodes and add to network as attribute
#also store all subgraphs in list
sub_list = list([])
stimul = undir.copy()
for i in list(stimul.nodes()):
	stimul.node[i]['stimulation'] = abs_sum(list(stimul.node[i]['behavior'].values()))
	sub_list += [stimul.subgraph(list(stimul.neighbors(i))+[i])]



#in all subgraphs, find average stimulation value
#and create list of lists of subgraph and ave_stimul
ave_s = []
for u in range(len(sub_list)):
	sub = sub_list[u]
	for j in sub.nodes():
		stimul_list = [sub.node[n]['stimulation'] for n in sub.nodes()]
		ave = sum(stimul_list)/len(stimul_list)
		ave_s += [ave]

sub_list_s = list(sorted(zip(sub_list, ave_s), key=lambda x:x[1], reverse=True))

'''
actually, instead of the subgraphs, possibly the markov clusters are more useful.
find the markov clusters and again do the sorting based on ave_simul
that is done above.
However, to avoid weirdness, take only the major connected component graph.
'''
clength_np_dir = nx.weakly_connected_component_subgraphs(network)
small = max(clength_np_dir,key=len).copy() #taking only connected graph
small_undir = small.to_undirected()

matrix = nx.to_scipy_sparse_matrix(small_undir) #doing mcl, but optimizing Q with loop.


result = mc.run_mcl(matrix, inflation=1.8)
clusters = mc.get_clusters(result)
print("number of clusters: ", len(clusters))

len([clusters[i] for i in range(len(clusters)) if len(clusters[i]) >= 3])

# for inflation in [i / 200 for i in range(280, 320)]:
#     result = mc.run_mcl(matrix, inflation=inflation)
#     clusters = mc.get_clusters(result)
#     print('clusters done, now calculating Q...')
#     Q = mc.modularity(matrix=result, clusters=clusters)
#     print("inflation:", inflation, "modularity:", Q)
#     print("number of clusters: ", len(clusters))
'''
optimization results:
1.5 Q: 0.9
2 Q: 0.7
But increases as one approaches inflation=1.
Pick 1.8 because it gives number of clusters close to literature result.
'''

#clusters is list of lists, therefore just fetch subgraph. Directed? Undirected?
#Let's look at directed.

#first only pick clusters with >=3 proteins.
int_clusters = [clusters[i] for i in range(len(clusters)) if len(clusters[i]) >= 3]

#save clusters!
array = np.array(int_clusters)
np.save('/Users/matteorosati/Desktop/MacMicking Lab/IFN_microarray/Processed Data/clusters.csv', array)

#fetch subgraphs
clusters_subgraphs = list([])
for i in range(len(int_clusters)):
	sub_nodes = []
	for j in int_clusters[i]:
		node = list(small.nodes())[j]
		small.node[node]['stimulation'] = abs_sum(list(small.node[node]['behavior'].values()))
		sub_nodes += [node]
	clusters_subgraphs += [small.subgraph(sub_nodes)]


#finding values of ave_s based on 1) [ave_s] and 2) [cool_s]
ave_s = []
cool_s = []
for u in clusters_subgraphs:
	stimul_list_val = [u.node[n]['stimulation'] for n in u.nodes()]
	stimul_list = [m for m in list(u.nodes()) if u.node[m]['stimulation'] != 0]
	ave = sum(stimul_list_val)/len(stimul_list_val)
	if len(stimul_list) == 0:
		cool = 0
	else:
		cool = sum(stimul_list_val)/len(stimul_list)
	ave_s += [ave]
	cool_s += [cool]

subgraphs_ave = list(sorted(zip(clusters_subgraphs, ave_s), key=lambda x:x[1], reverse=True))
subgraphs_cool = list(sorted(zip(clusters_subgraphs, cool_s), key=lambda x:x[1], reverse=True))

quantities_df = pd.DataFrame({'Method 1': ave_s, 'Method 2': cool_s})

quantities_df.to_csv('/Users/matteorosati/Desktop/MacMicking Lab/IFN_microarray/Processed Data/subgraph_quant_data_wund.csv', index=False)

#graph top subgraphs for both statistics (actually, graph top two)
pos_ave = [nx.layout.spring_layout(subgraphs_ave[0][0]), nx.layout.spring_layout(subgraphs_ave[1][0])]
pos_cool = [nx.layout.spring_layout(subgraphs_cool[0][0]), nx.layout.spring_layout(subgraphs_cool[1][0])]

def Graph_sub(layout, subgraph, title):
	sizes = [subgraph.node[n]['stimulation']*15 + 5 for n in list(subgraph.nodes())]
	nx.draw_networkx_nodes(subgraph, layout, 
		node_size=sizes, node_color='blue')
	nx.draw_networkx_edges(subgraph, layout, node_size=sizes, 
		arrowstyle='->', arrowsize=10, edge_color='teal', width=2)
	labels = {n:subgraph.node[n]['symbol'] for n in list(subgraph.nodes())}
	nx.draw_networkx_labels(subgraph, layout, labels=labels, font_size=12, alpha = 0.8)
	plt.title(title)
	plt.axis('off')
	plt.show()

Graph_sub(pos_ave[0], subgraphs_ave[0][0], 'First Subgraph, Method 1')
Graph_sub(pos_ave[1], subgraphs_ave[1][0], 'Second Subgraph, Method 1')

Graph_sub(pos_cool[0], subgraphs_cool[0][0], 'First Subgraph, Method 2')
Graph_sub(pos_cool[1], subgraphs_cool[1][0], 'Second Subgraph, Method 2')

#graphing all clusters
cluster_network = nx.DiGraph()

#define a function which determines given a node which cluster it is in
def findC(node):
	for i in clusters_subgraphs:
		if node in list(i.nodes()):
			return i
		else:
			pass


#generate edge list
edgelist = []
for i in clusters_subgraphs:
	for j in list(i.nodes()):
		neighbors = small.neighbors(j)
		for k in neighbors:
			if k in list(i.nodes()):
				pass
			else:
				edgelist += [(i, findC(k))]

edgelist = list(set(edgelist))

#add edges
cluster_network.add_edges_from(edgelist)
cluster_network.remove_node(None)

#visualize
pos = nx.layout.shell_layout(cluster_network)
# ave_dict = {str(n[0]):n[1] for n in subgraphs_ave}
# cool_dict = {str(n[0]):n[1] for n in subgraphs_cool}
for n in list(cluster_network.nodes()):
	#pdb.set_trace()
	stimul_list_val = [n.node[m]['stimulation'] for m in list(n.nodes())]
	stimul_list = [p for p in list(n.nodes()) if n.node[p]['stimulation'] != 0]
	cluster_network.node[n]['ave'] = sum(stimul_list_val)/len(stimul_list_val)
	if len(stimul_list) == 0:
		cluster_network.node[n]['cool'] = 0
	else:
		cluster_network.node[n]['cool'] = sum(stimul_list_val)/len(stimul_list)


sizes_ave = [cluster_network.node[n]['ave']*10 + 5 for n in list(cluster_network.nodes())]
sizes_cool = [cluster_network.node[n]['cool']*10 + 5 for n in list(cluster_network.nodes())]
nx.draw_networkx_nodes(cluster_network, pos, node_size=sizes_cool, node_color='blue')
nx.draw_networkx_edges(cluster_network, pos, node_size=sizes_cool, 
	arrowstyle='->', arrowsize=5, edge_color='black', width=0.5, alpha = 0.5)
plt.title('Clusters in Network, Size by Method 2')
plt.axis('off')
plt.show()






