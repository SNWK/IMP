import heapq
import multiprocessing
import operator
import sys
import time
import math
from functools import reduce

import numpy as np

NETWORK = {}
Vnumber = 0
Enumber = 0
process = 8

def Load_data(network_file):
    global NETWORK, Vnumber, Enumber
    file = open(network_file, "r", encoding="utf-8")
    Vnumber, Enumber = file.readline().split()
    Vnumber = int(Vnumber) + 1
    Enumber = int(Enumber)
    # (weight, dijk)
    for i in range(0, Vnumber):
        NETWORK[i] = []
    for i in range(0, Enumber):
        start, end, weight = file.readline().split()
        start = int(start)
        end = int(end)
        weight = float(weight)
        NETWORK[end].append((start, weight))
    # dict行数 = v
    # 每一行[:]，[end, weight]


def C(n, k):
    return reduce(operator.mul, range(n - k + 1, n + 1)) / reduce(operator.mul, range(1, k + 1))


def getLamb(EpsilonNew, n, k, l):
    a = 2 + (2 / 3) * EpsilonNew
    b = math.log(C(n, k)) + l * math.log(n) + math.log(math.log2(n))
    c = EpsilonNew **2
    return (a * b * n )/ c


def generate(v, NETWORK):

    Rnew = set()
    len0 = len(NETWORK[v])
    for i in range(len0):
        Rnew.add(NETWORK[v][i][0])
    R = Rnew.copy()
    while Rnew:
        a_temp = set()
        for i in Rnew:
            len1 = len(NETWORK[i])
            for j in range(len1):
                if np.random.uniform(0,1) <= NETWORK[i][j][1]:
                    if NETWORK[i][j][0] not in R:
                        a_temp.add(NETWORK[i][j][0])
                        R.add(NETWORK[i][j][0])
        Rnew = a_temp.copy()
    return R

def generate_LT(v,NETWORK):

    Rnew = set()
    len0 = len(NETWORK[v])
    R = Rnew.copy()
    for i in range(len0):
        Rnew.add(NETWORK[v][i][0])
    while Rnew:
        a_temp = set()
        for i in Rnew:
            if NETWORK[i] == []:
                continue
            len1 = len(NETWORK[i])
            pp = np.random.randint(0,len1)
            if NETWORK[i][pp][0] not in R:
                a_temp.add(NETWORK[i][pp][0])
                R.add(NETWORK[i][pp][0])
        Rnew = a_temp.copy()
    return R


def cover(S, Rset):
    len0 = len(Rset)
    S = set(S)
    cont = 0
    for i in range(len0):
        for iterm in Rset[i]:
            if iterm in S:
                cont += 1
                break
    return cont


def nodeSelection(Rset, k):
    # start = time.time()
    Sk = []
    Rdict = dict()
    num = []
    total = len(Rset)
    active = 0
    for i in range(Vnumber):
        num.append(0)
    for i in range(total):
        for j in Rset[i]:
            num[j] += 1
            if j not in Rdict.keys():
                Rdict[j] = [i]
            else:
                Rdict[j].append(i)
    while len(Sk) < k:
        maxx = max(num)
        s = num.index(maxx)
        active += maxx
        Sk.append(s)
        rr = Rdict[s].copy()
        for i in rr:
            for j in Rset[i]:
                num[j] -= 1
                Rdict[j].remove(i)

    # print("noed",time.time()-start)
    return Sk, active / total


def getLambStar(k, l, n, e):
    alpha = (l * math.log(n) + math.log(2)) ** 0.5
    beta = ((1 - 1 / math.e) * (math.log(C(n, k)) + alpha * alpha)) ** 0.5
    return 2 * n * (((1 - 1 / math.e) * alpha + beta) ** 2) * (e ** (-2))


def addRset(que, times, Vnumber, Network):
    Rset = []
    for i in range(times):
        vv = np.random.randint(0, Vnumber)
        Rset.append(generate(vv,Network))
    que.put(Rset)

def addRset_LT(que, times, Vnumber, Network):
    Rset = []
    for i in range(times):
        vv = np.random.randint(0, Vnumber)
        Rset.append(generate_LT(vv,Network))
    que.put(Rset)

def sampling(NETWORK, k, Epsilon, l):
    # start = time.time()
    Rset = []
    LB = 1
    EpsilonNew = (2 ** 0.5) * Epsilon
    n = Vnumber
    log2n = int(math.log2(n))
    lambda_a = getLamb(EpsilonNew, n, k, l)
    lambda_b = getLambStar(k, l, Vnumber, Epsilon)
    for i in range(1, log2n):
        Rset = []
        x = n / (2**i)
        Thetai = lambda_a / x

        p = multiprocessing.Pool(process)
        times = int((Thetai - len(Rset)) // (process)) + 1
        que = multiprocessing.Manager().Queue()
        for t in range(process-1):
            p.apply_async(addRset, args=(que, times, Vnumber, NETWORK))
        addRset(que, times, Vnumber, NETWORK)
        while not que.empty():
            Rset += que.get()
        p.close()
        p.join()
        while not que.empty():
            Rset += que.get()

        Si, fraction = nodeSelection(Rset, k)
        if n * fraction >= (1 + EpsilonNew) * x:
            LB = n * fraction / (1 + EpsilonNew)
            break
    # print("sed1", time.time() - start)
    theta = lambda_b / LB
    least = theta - len(Rset)
    times = int(least//(process))+1
    p = multiprocessing.Pool(process)
    que = multiprocessing.Manager().Queue()
    for t in range(process - 1):
        p.apply_async(addRset, args=(que, times, Vnumber, NETWORK))
    addRset(que, times, Vnumber, NETWORK)
    while not que.empty():
        Rset += que.get()
    p.close()
    p.join()
    while not que.empty():
        Rset += que.get()
    return Rset

def sampling_LT(NETWORK, k, Epsilon, l):
    # start = time.time()
    Rset = []
    LB = 1
    EpsilonNew = (2 ** 0.5) * Epsilon
    n = Vnumber
    log2n = int(math.log2(n))
    lambda_a = getLamb(EpsilonNew, n, k, l)
    lambda_b = getLambStar(k, l, Vnumber, Epsilon)
    for i in range(1, log2n):
        Rset = []
        x = n / (2**i)
        Thetai = lambda_a / x

        p = multiprocessing.Pool(process)
        times = int((Thetai - len(Rset)) // (process)) + 1
        que = multiprocessing.Manager().Queue()
        for t in range(process-1):
            p.apply_async(addRset_LT, args=(que, times, Vnumber, NETWORK))
        addRset(que, times, Vnumber, NETWORK)
        while not que.empty():
            Rset += que.get()
        p.close()
        p.join()
        while not que.empty():
            Rset += que.get()

        Si, fraction = nodeSelection(Rset, k)
        if n * fraction >= (1 + EpsilonNew) * x:
            LB = n * fraction / (1 + EpsilonNew)
            break
    # print("sed1", time.time() - start)
    theta = lambda_b / LB
    least = theta - len(Rset)
    times = int(least//(process))+1
    p = multiprocessing.Pool(process)
    que = multiprocessing.Manager().Queue()
    for t in range(process - 1):
        p.apply_async(addRset_LT, args=(que, times, Vnumber, NETWORK))
    addRset(que, times, Vnumber, NETWORK)
    while not que.empty():
        Rset += que.get()
    p.close()
    p.join()
    while not que.empty():
        Rset += que.get()
    return Rset


if __name__ == '__main__':
    #python IMP.py –i network.txt -k 5 -m IC -t 60
    time_start = time.time()
    network_file = sys.argv[2]
    seed_size = int(sys.argv[4])
    diffusionmodel = sys.argv[6]
    TIME = int(sys.argv[8])
    # network_file = "network5.txt"
    # seed_size = 5
    Load_data(network_file)
    # print(NETWORK)
    # print("s",time.time() - time_start)
    l = 1 + math.log(2)/math.log(Vnumber)
    if diffusionmodel == "IC":
        Rset = sampling(NETWORK, seed_size, 0.07, l)
    else:
        Rset = sampling_LT(NETWORK, seed_size, 0.07, l)
    # print("n",time.time() - time_start)
    S, _ = nodeSelection(Rset, seed_size)
    for i in S:
        print(i)
    # print(S)
    print(time.time() - time_start)
