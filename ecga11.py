#!/usr/bin/env python
# coding: utf-8

# **Implementação do ECGA**

# In[1]:


import numpy as np
import operator
import math
import itertools


# In[2]:


'''
# Funçao de avaliacao: Problema armadilha
# k é o tamanho do BB e o u é o vetor de bits deste BB
def trapk(k,u):
    if u.count(1) == k:
        return(k)
    else:
        return((k-1)-u.count(1))
        # saida: int
'''

# Funçao de avaliacao: Problema armadilha f3deceptive k = 3
# k é o tamanho do BB e o u é o vetor de bits deste BB
def trapk(k,u):
    if u.count(1) == 0:
        return(0.9)
    elif u.count(1) == 1:
        return(0.8)
    elif u.count(1) == 2:
        return(0)
    else:
        return(1)
        # saida: int
    
# Soma dos resultados da função de avaliação de cada BB do cromomosso x
def ftrapk(x,k):
    stop = 0
    aux = 0
    #print(int(len(x)/k))
    for i in range(int(len(x)/k)):
        aux = aux + trapk(k,x[stop:k+stop])
        stop = stop + k
    return(aux)
    # saida: int


# In[3]:


# Gera o cromossomo da partir da lista de probabilidade para cada partição
def gera_cromossomo(l_prob):
    # l_pro: lista de probalidade para cada partição do cromossomo
    cromossomo = []
    for t in range(len(l_prob)):
        probabilidade = [item[1] for item in l_prob[t]]
        possibilidades = [item[0] for item in l_prob[t]]
        lista_possibilidades = [i for i in range(len(possibilidades))]
        escolheu = np.random.choice(lista_possibilidades, 1, p = probabilidade)
        cromossomo = cromossomo + possibilidades[escolheu[0]]
    return(cromossomo)
    # [int, int,...,int]


# In[4]:


# Gerar uma populacao com uma lista de probabilidade para cada partição
def gera_populacao(tam_pop, l_prob, l_part, k,n_a):
    # tam_cro: tamanho do cromossomo
    # tam_pop: tamanho da populacao
    # l_pro: lista de probalidade para cada partição
    populacao = []
    for i in range(tam_pop):
        cromossomo = gera_cromossomo(l_prob)
        populacao.append([cromossomo,ftrapk(cromossomo,k)])
        n_a = n_a + 1
    return(populacao,n_a)
    # saida: [[int, int, ..., int],[int, int, ..., int],...,[int, int, ..., int]]


# In[6]:


# Torneio: Seleciona 50% cromossomo com maior fitness
def torneio(populacao):
    # ordenar pelo fiteness maior para o menor
    populacao = sorted(populacao, key=operator.itemgetter(1), reverse=True)
    # escolhemos metade da populacao
    # caso seja impar arredondamos para cima
    del populacao[math.ceil(len(populacao)/2):len(populacao)]        
    return(populacao)
    # saida: [[[int, int, ..., int],int],[[int, int, ..., int],int],...,[[int, int, ..., int],int]


# In[7]:


# Calcula a porcentagem de cada particao dada
def cal_porcent(populacao, l_ind_part):
    tamanho = len(populacao)
    l_prob = []
    
    for particao in l_ind_part:
        l_part_perm = list(itertools.product([0, 1], repeat=len(particao)))
        pop_part = [[cromossomo[0][i] for i in particao] for cromossomo in populacao]
        l_aux = []
        for item in l_part_perm:
            l_aux.append([list(item),pop_part.count(list(item))/tamanho])
        l_prob.append(l_aux)
        
    return(l_prob)    
    # saida: [[lista_part[0],%],...,[lista_part[tam_particao],%]]      


# In[8]:


# Calcula Cp: Complexidade da Populacao Comprimida
# Cp = tam_pop*sum(1 a m)Cp_i
# m: numero de particoes no modelo
def Compl_p(tam_pop,l_prob):
    # Calcula o Cpi = sum(i a 2^k_i)-p_ijlog2(pij)
    
    # Calcula Cp
    Cp = 0
    aux = 0
    m = len(l_prob)
    for i in range(m):
        # Calcula Cpi
        Cpi = 0
        for j in range(len(l_prob[i])):
            pi = l_prob[i][j][1]
            if pi != 0:
                Cpi = Cpi - l_prob[i][j][1]*math.log2(l_prob[i][j][1])
        aux = aux + Cpi
    Cp = tam_pop*aux
    return(Cp)


# In[9]:


def Compl_m(tam_pop,l_prob):
    Cm = 0
    m = len(l_prob)
    for j in range(m):
        Cm = Cm + (2**len(l_prob[j][0][0])-1)
    Cm = math.log2(tam_pop)*Cm
    return(Cm)


# In[10]:


def calcula_Cc(pop,n,l_part):
    # calcula l_prob para as partições
    l_prob = []
    l_prob = cal_porcent(pop, l_part)

    # calcula Cc para as partições
    return(Compl_p(n,l_prob) + Compl_m(n,l_prob))


# In[11]:


def preenchendo_lista(l_part, comb):
    # l_part: lista de partição atual
    # comb: combinação escolhida, no caso é uma tupla
    
    # lista todos elementos que nao estao na lista
    lista = [x for x in l_part if x not in comb]
    
    # transforma a comb em lista
    aux = []
    for j in range(len(comb)):
        aux = aux + comb[j]
    aux.sort()
    # adiciona na lista e ordena
    lista.append(aux)
    lista.sort()
    return(lista) 


# In[12]:


def combinatoria(tam_comb,l_part):
    # tam_comb: tamanho da combinacao
    # l_part: lista de particao atual

    l_l_comb = []
    combinacoes = itertools.combinations(l_part, tam_comb)
    for comb in list(combinacoes):
        aux = preenchendo_lista(l_part, comb)
        l_l_comb.append(aux)
    
    return(l_l_comb)


# In[13]:


def gulosa(pop,n,tam_comb,Cc_atual,l_part_atual):
    # gera a lista de combinações possiveis com a l_part_atual
    l_l_comb = combinatoria(tam_comb,l_part_atual)
    # calcula o Cc para cada um e verifica se tem um menor que o Cc atual
    # caso encontre, atualiza
    #print(len(l_l_comb))
    encontrado = False
    for l_comb in l_l_comb:
        Cc = calcula_Cc(pop,n,l_comb)
        #print(Cc)
        if Cc < Cc_atual:
            encontrado = True
            l_part_atual = l_comb
            Cc_atual = Cc
            #print("escolhido",Cc_atual)
    if encontrado == True:
        return(gulosa(pop,n,tam_comb,Cc_atual,l_part_atual))
    else:
        return(l_part_atual,Cc_atual)


# In[14]:


def ecga(n,l,tam_comb,ciclo_Max,k):

    # Numero de avaliações
    n_a = 0
    # primeira particao realizada: [[0],[1],...,[l-1]]
    l_part_atual = [[i] for i in range(l)]
    # lista de probabilidade inicial
    l_prob_ini = []
    for i in range(l):
        l_prob_ini.append([[[0], 0.5], [[1], 0.5]])

    # gera a populacao
    pop,n_a = gera_populacao(n, l_prob_ini, l_part_atual,k,n_a)
    # realiza o torneio
    pop = torneio(pop)
    
    # calcula Cc para as partições de tamanho 1
    Cc_atual = calcula_Cc(pop,n,l_part_atual)

    # lista de particoes atual
    #print(Cc_atual, l_part_atual)

    ciclo = 0
    estagnou = False
    
    # verifica se o modelo encontrou o otimo
    otimo = True
    for cromossomo in pop:
        if otimo == True:
            if cromossomo[1] < l-1:
                otimo = False
                break
                
    while ciclo < ciclo_Max and estagnou == False and otimo == False:
        # Modelando
        l_part_atual,Cc_atual = gulosa(pop,n,tam_comb,Cc_atual,l_part_atual)
        #print("escolhido",Cc_atual,l_part_atual)
        # calcula as probabilidades
        l_prob = cal_porcent(pop, l_part_atual)
        
        # verifica se o modelo estagnou
        estagnou = True
        for particoes in l_prob:
            if estagnou == True:
                for part in particoes:
                    if part[1] == 0 or part[1] == 1:
                        estagnou = True
                    else:
                        estagnou = False
                break
            else:
                break
                        
        # verifica se o modelo encontrou o otimo
        n_otimo = 0
        for part in l_prob:
            for j in part:
                procurando = False
                if (j[0].count(1) == len(j[0])) and (j[1] == 1):
                    procurando = True
                    n_otimo = n_otimo + 1
                    break
            
        otimo = False
        if n_otimo == 1 and len(l_prob) == 1:
            otimo = True    
        if n_otimo >= len(l_prob) - 1 and len(l_prob) > 1:
            otimo = True
                                        
        if estagnou != True and otimo != True:
            #l_part_atual = [[i] for i in range(l)]
            ciclo = ciclo + 1
            # gerando a populacao a partir do modelo e fazendo um torneio
            pop,n_a = gera_populacao(n, l_prob, l_part_atual,k,n_a)

            pop = torneio(pop)
            
        #print(l_prob)
        #print('\n')
    #fitness = [row[1] for row in pop]
    #print("escolhido",Cc_atual,l_part_atual)
    return(n_a,ciclo, estagnou,l_part_atual,l_prob,otimo,n_otimo, len(l_prob))

