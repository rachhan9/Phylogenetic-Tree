import os
from Bio import SeqIO
import numpy as np
import pandas as pd
from anytree import AnyNode, Node, RenderTree, exporter

def preprocessing():
    # Preprocess datasets, turn character to upper, remove character not A T C  G
    directory = os.getcwd()+"/nCoV2019"
    sequences =  []
    labels = []
    for filename in os.listdir(directory):

        record = list(SeqIO.parse("nCoV2019/" + filename, "fasta"))[0]
        seq = record.seq.upper()
        seq = "".join(list(filter(lambda x: x in  ('A','T','G','C'),seq)))
        sequences.append(seq)
        labels.append(record.id)
    return sequences, labels


def cgr(seq, order, k):
    # cgr function copy from Assignment3.pdf
    ln = len(seq)
    pw = 2**k
    out = [[0 for i in range(pw)] for j in range(pw)]
    x = 2**(k-1)
    y = 2**(k-1)
    for i in range(0,ln):
        x=x//2
        y=y//2
        if(seq[i] == order[2] or seq[i] == order[3]):
            x = x + (2**(k-1))
        if(seq[i] == order[0] or seq[i] == order[3]):
            y = y + (2**(k-1))
        if(i>=k-1):
            out[y][x] = out[y][x]+1
    return out


def convert_seq_to_cgr(sequences,k):
    # covert all Sequences into Cgr arrays
    cgrs = []
    for seq in sequences:
       c =  cgr(seq,"ACGT",k)
       cgrs.append(c)
    return cgrs


def distance(A,B):
    # Calculate euclidian distance between A and B and normalize it
    # Require: A and B are cgr
    A = np.array(A)
    B = np.array(B)
    eucli_dist = np.linalg.norm(A-B)
    normal_eucli_dist = eucli_dist / ((np.sum(A**2) + np.sum(B ** 2)) ** 0.5)
    return normal_eucli_dist



def pairwise_distance_matrix(cgrs):
    # Construct a pairwise distance matrix base on all cgrs
    num_seq = len(cgrs)
    d = np.zeros([num_seq,num_seq])
    for i in range(num_seq):
        for j in range(num_seq):
            d[i][j] = distance(cgrs[i],cgrs[j])

    for i in range(num_seq):
        for j in range(num_seq):
            if i == j:
                assert(d[i][j] == 0)
            assert(d[i][j] == d[j][i])
            assert(d[i][j] >= 0 )
            assert(d[i][j] <= 1)
    return d

def UPGMA_alg(distance_matrix,labels):
    # UPGMA algorithm  in the lecture slides.
    mat = pd.DataFrame(distance_matrix, columns=labels, index=labels)

    node_dic = {}
    # PlaceHolder = [100]

    def find_smallest_dist(mat):
        # find smallest distance for the distance matrix mat
        smallest  = [100,-1,-1]
        for i in mat.index:
            for j in mat.columns:
                if i == j: continue
                if  mat[i][j] < smallest[0]:
                    smallest = [mat[i][j], i, j]
        return smallest
    
    def construct_new_mat(i,j,mat):
        # construct a new distance matrix through combining node i and j

        if i not in node_dic:
            node_dic[i] = Node(i,leaf=True,label=i)
        if j not in node_dic:
            node_dic[j] = Node(j,leaf=True,label=j)

        old_mat = mat
        mat = mat.drop(columns=[i, j],axis=1)
        mat = mat.drop(index =[i, j])
        new_dist = []
        for r in mat.index:
            new_dist.append(0.5 * (old_mat[r][i]+ old_mat[r][j]))
        new_name = "(" + i + "," + j + ")"
        mat.loc[:,new_name] = new_dist
        new_dist.append(0) 
        mat.loc[new_name,:] = new_dist

        # combining node
        node_dic[new_name] = Node(new_name,leaf=False,label="")
        node_dic[i].parent = node_dic[new_name]
        node_dic[j].parent = node_dic[new_name]

        del node_dic[i]
        del node_dic[j]

        return mat

    while len(mat) != 1:
        smallest = find_smallest_dist(mat)
        i = smallest[1]
        j = smallest[2]
        mat = construct_new_mat(i,j,mat)

    key = list(node_dic.keys())[0]

    return mat.index[0], node_dic[key]



def generate_tree(k):
    sequences, labels = preprocessing()
    cgrs = convert_seq_to_cgr(sequences,k)
    distance_matrix = pairwise_distance_matrix(cgrs)
    print("Distance matrix for k = "  + str(k))
    print(distance_matrix)
    print()
    # labels = list(map(lambda x: str(x),range(0,len(sequences))))
    newick, tree = UPGMA_alg(distance_matrix,labels)
    exporter.DotExporter(tree,nodeattrfunc=lambda n: 'label="{}"'.format(n.label)).to_picture("tree{}.png".format(k))
    print("Newick form of tree for k = "+ str(k))
    print(newick)
    print()


def  main():
    generate_tree(3)
    generate_tree(9)

main()
