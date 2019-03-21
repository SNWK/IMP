import sys
import time
import heapq
import numpy as np
import math

NETWORK = []
MIIA = []
MIOA = []
nodetable = []
dijk = []
Vnumber = 0
Enumber = 0
MAX = 9999


def MIA(k):
    S = []
    incinf = [0] * Vnumber
    for v in range(Vnumber):
        for i in range(len(MIIA[v])):
            #print(MIIA[v][i][0])
            incinf[MIIA[v][i][0]] += getAlpha_initial(v, MIIA[v][i][0])
    ###main loop
    for i in range(0, k):
        u = maxp(incinf, S)
        for j in range(len(MIOA[u])):
            if MIOA[u][j][0] in S:
                continue
            else:
                incinf[j] -= getAlpha(MIOA[u][j][0],u,S)*(1-ap(u,S,MIOA[u][j][0]))
        S.append(u)
        for j in range(len(MIOA[u])):
            if MIOA[u][j][0] in S:
                continue
            else:
                incinf[j] += getAlpha(MIOA[u][j][0],u,S)*(1-ap(u,S,MIOA[u][j][0]))
    return S


def maxp(incinf, S):
    maxp = 0
    maxv = 0
    for i in range(Vnumber):
        if i in S:
            continue
        elif incinf[i] > maxv:
            maxv = incinf[i]
            maxp = i
    return maxp


def getAlpha(v, u, S):
    if v == u:
        return 1
    else:
        w = MIOA[u][0][0]  ### 0 ???
        if w in S:
            return 0
        else:
            listV = []
            for i in range(len(MIIA[v])):
                listV.append(MIIA[v][i][0])
            Nin = []
            for i in range(len(MIIA[w])):
                if MIIA[w][i][0] in listV:
                    Nin.append(MIIA[w][i])

            pun = 1
            for i in range(len(Nin)):
                if Nin[i][0] == u:
                    continue
                else:
                    pun = pun * (1 - ap(Nin[i][0], S, v) * Nin[i][1])
            return getAlpha_initial(v, w) * MIOA[u][0][1] * (1 - pun)




def getAlpha_initial(v, u):
    if v == u:
        return 1
    else:
        #print(u, MIOA[u])
        if MIOA[u] == []:
            return 0
        w = MIOA[u][0][0]  ### 0 ???
        return getAlpha_initial(v, w) * MIOA[u][0][1]


def ap(u, S, v):
    listV = []
    for i in range(len(MIIA[v])):
        listV.append(MIIA[v][i][0])
    Nin = []
    for i in range(len(MIIA[u])):
        if MIIA[u][i][0] in listV:
            Nin.append(MIIA[u][i])

    if u in S:
        return 1
    elif Nin == []:
        return 0
    else:
        pun = 1
        for i in range(len(Nin)):
            pun = pun * (1 - ap(Nin[i][0], S, v) * Nin[i][1])  ##CUN YI !!!  ap(Nin[i][0], S, v)
        return 1 - pun


def Load_data_IC(network_file):
    global NETWORK, nodetable, Vnumber, Enumber, MIIA, MIOA
    file = open(network_file, "r", encoding="utf-8")
    Vnumber, Enumber = file.readline().split()
    Vnumber = int(Vnumber)
    Enumber = int(Enumber)
    nodetable = [[MAX for j in range(Vnumber)] for j in range(Vnumber)]
    # (weight, dijk)
    for i in range(0, Vnumber):
        NETWORK.append([])
        MIIA.append([])
        MIOA.append([])
    for i in range(0, Enumber):
        start, end, weight = file.readline().split()
        start = int(start)
        end = int(end)
        weight = float(weight)
        NETWORK[start].append([end, weight])
        nodetable[start][end] = math.log(weight)
    return nodetable
    # 行数 = v
    # 每一行[:]，[end, weight]


def getdijk_IC():
    global dijk
    for i in range(0, Vnumber):
        ##print(i)
        dijk.append(Dijkstra(nodetable, i))
    return dijk


def Dijkstra(table, node):
    final = [0] * (Vnumber)
    distance = [0] * (Vnumber)
    path = [0] * (Vnumber)
    for i in range(0, Vnumber):
        distance[i] = table[node][i]
        if distance[i] != MAX:
            path[i] = node
        else:
            path[i] = MAX
    final[node] = 1
    path[node] = node
    k = 0
    for i in range(0, Vnumber):
        min = MAX
        for j in range(0, Vnumber):
            if distance[j] < min and final[j] == 0:
                min = distance[j]
                k = j
        final[k] = True
        for j in range(0, Vnumber):
            if (distance[j] > min + table[k][j]) and final[j] == 0:
                distance[j] = min + table[k][j]
                path[j] = k
    distance[node] = 0
    return distance


def getMMM(theta):
    global MIIA, MIOA
    for i in range(Vnumber):
        for j in range(Vnumber):
            if dijk[i][j] > 1 or i == j:
                continue
            elif dijk[i][j] >= theta:
                MIIA[j].append((i,dijk[i][j]))
                MIOA[i].append ((j, dijk[i][j]))


if __name__ == '__main__':
    # python IMP.py –i network.txt -k 5 -m IC -t 60
    time_start = time.time()
    #network_file = sys.argv[2]
    #seed_size = int(sys.argv[4])
    #diffusionmodel = sys.argv[6]
    #TIME = int(sys.argv[8])
    network_file = "network.txt"
    nodetable = Load_data_IC(network_file)
    print(nodetable)
    print(MIOA[10])
    print("begindijk", time.time() - time_start)
    dijk = getdijk_IC()
    getMMM(0.05)
    print ("beginmia", time.time() - time_start)
    S = MIA(5)
    time_end = time.time()
    print(S)
    print(time_end - time_start)
#[33, 51, 60, 46, 53]