from random import random, randint, sample, uniform
import functools

class solution():
    def __init__(self, numvar, l, chrom=None, bounds=None, decoder=None, evaluator=None):
        '''
        numvar - number of variables
        l - length (in bits) of each variable
        bounds - [ [lo, hi], [lo, hi], [lo, hi]...], length is n
        '''
        self.n = 0 # domination counter, the # of solutions that dominate this one
        self.distance = 0 
        self.fitness = 0
        self.numvar = numvar 
        self.l = l
        if not chrom:
            self.chromosome = ''.join(['1' if flip(0.5) else '0' for i in range(numvar*l)])
        else:
            self.chromosome = chrom
        if len(bounds) != numvar:
            raise Exception('Bounds must have the same length as n')
        if len(self.chromosome) != numvar*l:
            raise Exception('Lengh of chromosome error')
        self.decoder = decoder
        self.evaluator = evaluator
        self.obj = []
        self.bounds = bounds
        self.deltas = [] # used in decode to convert binary 
        for b in bounds:
            d = (b[1]-b[0]) / (2**l-1)
            self.deltas.append(d) # store delta
        # self.obj = evaluator(decoder(self.chromosome))
        
    def decode(self):
        if not self.evaluator:
            raise Exception('Evaluator not specified')
        return self.decoder(self.chromosome, self.numvar, self.l, self.bounds, self.deltas)
        
    def evaluate(self):
        if not self.evaluator:
            raise Exception('Evaluator not specified')
        self.obj = self.evaluator(self.decode())
    
    def __lt__(self, other): # return true if there is at least one < and others <= or <
        num_proper_smaller_than = 0
        for i in range(len(self.obj)):
            if self.obj[i] > other.obj[i]:
                return False
            if self.obj[i] < other.obj[i]:
                num_proper_smaller_than += 1
        return num_proper_smaller_than > 0

    def __gt__(self, other):
        num_proper_greater_than = 0
        for i in range(len(self.obj)):
            if self.obj[i] < other.obj[i]:
                return False
            if self.obj[i] > other.obj[i]:
                num_proper_greater_than += 1
        return num_proper_greater_than > 0
    
    def dominate(self, other, direction='minimize'): # return true if this object dominates the other
        if direction=='minimize':
            return self < other
        else:
            return self > other


def fast_non_domnated_sort(P): # P is a list of individual objects
    F = {1:[]}                        #. stores index in R (P and Q combined)
    S = {}
    for i in range(len(P)):
        p = P[i]
        S[i] = []
        p.n = 0
        for j in range(len(P)):
            q = P[j]
            if p < q:                 # If p dominates q (#. assuming minimization)
                S[i].append(j)        # Add q to the set of solutions dominated by p
            elif q < p:
                p.n += 1              # Increament the domination counter of p
        if p.n == 0:                  # p belongs to the first front
            p.rank = 1
            F[1].append(i)
    i = 1                             # Initialize the front counter
    while F[i] != []:
        Q = []                        # Used to store the members of the next front
        for pi in F[i]:
            for qi in S[pi]:
                P[qi].n -= 1
                if P[qi].n == 0:      # q belongs to the next front
                    P[qi].rank = i+1
                    Q.append(qi)
        i += 1
        F[i] = Q

    F1 = {}
    for f in F:
        F1[f] = [P[i] for i in F[f]]  #. change indexes to objects for return
    return F1

def crowding_distance_assignment(I):
    # I is a list of nondominated solutions
    l = len(I)
    if l == 0:
        return
    for m in range(len(I[0].obj)):
        # sort by obj m
        I.sort(key=lambda x: x.obj[m])
        I[0].distance = I[-1].distance = float('inf')
        fm_min = I[0].obj[m]
        fm_max = I[-1].obj[m]
        if fm_min == fm_max:
            for i in range(1, l-1):
                I[i].distance = 0
        else:
            for i in range(1, l-1):
                I[i].distance += (I[i+1].obj[m] - I[i-1].obj[m]) / (fm_max - fm_min)

def crowded_comparison_operator(i, j):
    # i and j are two individual objects
    # i is BETTER than j if:
    #     i_rank < j_rank 
    #     or 
    #     i_rank = j_rank and i_distance > j_distance
    #
    if i.rank < j.rank:
        return 1
    if i.rank == j.rank and i.distance > j.distance:
        return 1
    if i.rank == j.rank and i.distance == j.distance:
        return 0
    return -1

#
# Functions for a general purpose GA using binary strings
#

def decoder(chrom, numvar, l, bounds, deltas):
    '''
    This is a general purpose decoder for binary strings.
    
    chrom  - a binary string with a total length of numvar*l
    numvar - the number of variables encoded in the string
    l      - length (in bits) of each variable
    bounds - [ [lo, hi], [lo, hi], [lo, hi]...], length is n
    deltas - a list of constants used to decode the variables 
             To decode, we use var = lo + decimal(var_binary) * (hi - lo) / (2^l - 1) for each variable
             where lo and hi are the lower and upper bounds of the variable, and l length in bits.
             Each of the delta value is (hi-lo)/(2^l-1).

    Return - list of variable values
    '''
    vars = []
    for i in range(numvar):
        var = bounds[i][0] + int(chrom[i*l:(i+1)*l], 2)*deltas[i]
        vars.append(var)
    return vars


def selection(population):
    # tournament
    # population is a list of individual objects
    i, j = sample(range(len(population)), 2)
    i, j = population[i], population[j]
    comp = crowded_comparison_operator(i, j)
    if comp > 0:
        return i
    elif comp < 0:
        return j
    else:
        return i # i and j doesn't matter

def recombination(i, j, p):
    c1 = i.chromosome
    c2 = j.chromosome
    if not flip(p):
        return
    point = randint(0, i.numvar*i.l-1)
    c1 = i.chromosome[:point] + j.chromosome[point:]
    c2 = j.chromosome[:point] + i.chromosome[point:]
    i1 = solution(i.numvar, i.l, c1, i.bounds, i.decoder, i.evaluator)
    j1 = solution(i.numvar, i.l, c2, i.bounds, i.decoder, i.evaluator)
    return i1, j1

def flip(prob):
    if random() < prob:
        return True
    return False

def mutation(i, p):
    # reverse_chars = ['1', '0']
    # chrom = ''.join([c if not flip(p) else reverse_chars[int(c)] for c in i.chromosome])
    chrom = ''.join([c if not flip(p) else str(1-int(c)) for c in i.chromosome])
    i.chromosome = chrom

def make_new_pop(P, crossover_prob, mutation_prob):
    Q = []
    while len(Q) < len(P):
        i = selection(P)
        j = selection(P)
        recomb = recombination(i, j, crossover_prob)
        if recomb: # not None
            i1, j1 = recomb
            mutation(i1, mutation_prob)
            mutation(j1, mutation_prob)
            i1.evaluate()
            j1.evaluate()
            Q.append(i1)
            Q.append(j1)
    return Q


def nsga2_x(N, T, crossover_prob, mutation_prob, numvar, varlen, bounds, decoder, evaluator):
    '''
    N - population size
    T - number of generations
    bounds
    decoder
    evaluator
    '''
    R = {}
    P = {}
    Q = {}

    # created P0 randomly
    t = 0
    P[t] = [solution(numvar, varlen, bounds=bounds, decoder=decoder, evaluator=evaluator) for _ in range(N)] # evaluate by default?
    for i in P[t]:
        i.evaluate()

    fast_non_domnated_sort(P[t]) # need the fitness/rank for tournament selection

    Q[t] = make_new_pop(P[t], crossover_prob, mutation_prob)

    # main loop
    while t <= T:
        R[t] = P[t] + Q[t]
        F = fast_non_domnated_sort(R[t])
        P[t+1] = []
        i = 1
        while len(P[t+1])+len(F[i]) <= N:
            crowding_distance_assignment(F[i])
            P[t+1] = P[t+1] + F[i]
            i += 1
        # i is not the front that hasn't been included
        # will need some solutions from it if the previous loop hasn't got N solutions
        # F[i] is not assigned fitness/rank yet (overlooked in original paper)
        
        crowding_distance_assignment(F[i]) # NOTE: this step is missed in the original paper!
        F[i].sort(key=functools.cmp_to_key(crowded_comparison_operator), reverse=True) 
        P[t+1] = P[t+1] + F[i][:N-len(P[t+1])]
        Q[t+1] = make_new_pop(P[t+1], crossover_prob, mutation_prob)
        t += 1
    return P

