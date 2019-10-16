import networkx as nx
import matplotlib.pyplot as plt
import collections as col
from numpy import zeros, argmax


N_RL = 1000  # How many lines the program reads from the file

N = 5  # Number of time frames
PGD = 5
PCN = 5
PJC = 5
PA = 5
PPA = 5

parameter_list = [PGD, PCN, PJC, PA, PPA]


def min_max():  # Calculates the minimum and the maximum timestamp
    try:
        f = open('min_max' + str(N_RL) + '.txt')
        f = f.read().split(',')
        max = int(f[0])
        min = int(f[1])
    except FileNotFoundError:
        with open('sx-stackoverflow.txt') as f:
            line = f.readline()
            max = int(line.split(' ')[2].replace('\n', ''))
            min = int(line.split(' ')[2].replace('\n', ''))
            for j in range(N_RL):
                line = f.readline().split(' ')
                line[2] = int(line[2].replace('\n', ''))
                if line[2] > max:
                    max = line[2]
                elif line[2] < min:
                    min = line[2]
        print('Max is: ', str(max),
              '\nMin is: ', str(min))
        f = open('min_max' + str(N_RL) + '.txt', 'w')
        f.write(str(max) + ',' + str(min))
    return max, min


def graph(t, centrality='degree'):  # Creates a graph
    g = nx.DiGraph()
    with open('sx-stackoverflow.txt') as f:
        for j in range(N_RL):
            line = f.readline().split(' ')
            line[2] = int(line[2].replace('\n', ''))
            if qN[t][0] <= line[2] < qN[t][1] or (t == N - 1 and qN[t][0] <= line[2] <= qN[t][1]) \
                    and line[0] != line[1]:
                g.add_edge(line[0], line[1])
    nx.draw_networkx(g, node_size=150, font_size=10)
    plt.show()
    centrality_histogram(g, centrality)


def centrality_histogram(g, c):  # Creates the centrality histogram specified
    if c == 'degree':
        degree_sequence = sorted([val for key, val in nx.degree_centrality(g).items()])
    elif c == 'in_degree':
        degree_sequence = sorted([val for key, val in nx.in_degree_centrality(g).items()])
    elif c == 'out_degree':
        degree_sequence = sorted([val for key, val in nx.out_degree_centrality(g).items()])
    elif c == 'closeness':
        degree_sequence = sorted([val for key, val in nx.closeness_centrality(g).items()])
    elif c == 'betweenness':
        degree_sequence = sorted([val for key, val in nx.betweenness_centrality(g).items()])
    elif c == 'eigenvector':
        degree_sequence = sorted([val for key, val in nx.eigenvector_centrality(g).items()])
    elif c == 'katz':
        degree_sequence = sorted([val for key, val in nx.katz_centrality(g).items()])
    degree_count = col.Counter(degree_sequence)
    deg, cnt = zip(*degree_count.items())
    plt.bar(deg, cnt, width=0.01, color='b')
    plt.title('Degree Histogram')
    plt.ylabel('Count')
    plt.xlabel('Degree')
    plt.show()


def graph_star(t):  # Finds V* and E*
    if t < N - 1:
        g1 = nx.DiGraph()
        g2 = nx.DiGraph()
        with open('sx-stackoverflow.txt') as f:
            for j in range(N_RL):
                line = f.readline().split(' ')
                line[2] = int(line[2].replace('\n', ''))
                if qN[t][0] <= line[2] < qN[t][1]:
                    if line[0] != line[1]:
                        g1.add_edge(line[0], line[1])
                elif qN[t + 1][0] <= line[2] < qN[t + 1][1]:
                    if line[0] != line[1]:
                        g2.add_edge(line[0], line[1])
        v_star = []
        e1_star = []
        e2_star = []
        for i in g1.nodes:
            for j in g2.nodes:
                if i == j:
                    v_star.append(i)
                    break
        for i in g1.edges:
            found1 = False
            found2 = False
            for j in v_star:
                if i[0] == j:
                    found1 = True
                elif i[1] == j:
                    found2 = True
            if found1 and found2:
                e1_star.append(i)

        for i in g2.edges:
            found1 = False
            found2 = False
            for j in v_star:
                if i[0] == j:
                    found1 = True
                elif i[1] == j:
                    found2 = True
            if found1 and found2:
                e2_star.append(i)
        print('t =', t)
        print('V*[t_', t, ',', 't_', t + 2, '} = ', v_star)
        print('E*[t_', t, ',', 't_', t + 1, '] =', e1_star)
        print('E*[t_', t + 1, ',', 't_', t + 2, '] =', e2_star)
        if similarity_matrices(e1_star, v_star) == -1:
            return -1


def similarity_matrices(edges, nodes):  # Calculates the similarity matrices
    if len(nodes) == 0:
        print('V* is empty, skipping to next t.')
        return -1
    g = nx.DiGraph()
    g.add_edges_from(edges)
    ung = nx.Graph(g)

    gd = zeros((len(nodes), len(nodes)))
    cn = zeros((len(nodes), len(nodes)))
    jc = zeros((len(nodes), len(nodes)))
    a = zeros((len(nodes), len(nodes)))
    pa = zeros((len(nodes), len(nodes)))
    for i in range(len(nodes)):
        for j in range(len(nodes)):
            try:
                gd[i][j] = nx.shortest_path_length(g, nodes[i], nodes[j])
            except nx.NetworkXNoPath:
                gd[i][j] = -1
                pass
            except nx.NodeNotFound:
                gd[i][j] = -1
                pass
            try:
                cn[i][j] = len(sorted(nx.common_neighbors(ung, nodes[i], nodes[j])))
            except nx.NetworkXError:
                cn[i][j] = -1
                pass
            try:
                for u, v, p in nx.jaccard_coefficient(ung, [(nodes[i], nodes[j])]):
                    jc[i][j] = p
            except:
                jc[i][j] = -1
                pass
            try:
                for u, v, p in nx.adamic_adar_index(ung, [(nodes[i], nodes[j])]):
                    a[i][j] = p
            except ZeroDivisionError:
                a[i][j] = -1
                pass
            except nx.NetworkXError:
                a[i][j] = -1
                pass
            try:
                for u, v, p in nx.preferential_attachment(ung, [(nodes[i], nodes[j])]):
                    pa[i][j] = p
            except:
                pa[i][j] = -1
                pass
    ####
    k = 0
    for par in parameter_list:
        ind_list = []
        if k == 0:
            ref = gd
            t = 'Pgd'
        elif k == 1:
            ref = cn
            t = 'Pcn'
        elif k == 2:
            ref = jc
            t = 'Pjc'
        elif k == 3:
            ref = a
            t = 'Pa'
        else:
            ref = pa
            t = 'Ppa'

        for i in range(par):
            flat_ind = argmax(ref)
            dim_ind = tuple((flat_ind // len(nodes), flat_ind % len(nodes)))
            ref[dim_ind[0]][dim_ind[1]] = -1
            ind_list.append(dim_ind)
        cnt = 0
        for j in ind_list:
            if tuple((nodes[j[0]], nodes[j[1]])) in edges or tuple((nodes[j[1]], nodes[j[0]])) in edges:
                cnt += 1
        k += 1
        print(t, cnt / par)


max_time, min_time = min_max()
dt = max_time - min_time
tN = []
qN = []
dT = dt / N
for j in range(N + 1):
    tN.append(int(min_time + j * dT))
for j in range(len(tN) - 1):
    qN.append([tN[j], tN[j + 1]])
del tN

for i in range(N):
    graph(i, centrality='degree')
    if graph_star(i) == -1:
        continue
