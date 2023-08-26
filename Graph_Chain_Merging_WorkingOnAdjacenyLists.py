import networkx, matplotlib.pyplot as plot 

def kmers(seq,k):
    return [seq[i:i+k] for i in range(len(seq)-k+1)]

def DeBruijnGraph(reads,k):
    dicti={}
    E=[]
    for read in reads:
        KMers=kmers(read,k)
        edges=[kmers(mer,k-1) for mer in KMers]
        for edge in edges:
            if edge[0] not in dicti.keys(): dicti[edge[0]]=[]
            if edge[1] not in dicti.keys(): dicti[edge[1]]=[]
            dicti[edge[0]].append(edge[1])
            E.append(edge)
    V=list(dicti.keys())
    graph={'nodes':V,'edges':E}
    return (graph,dicti)

def visualizeDBGraph(graph):

    #Draw with labels of edges on graph 
    dbGraph=networkx.DiGraph()

    #Incase compacted graph has only a single node left with no edges as in test case 3 
    if(graph['edges'] == []):
         dbGraph.add_node(graph['nodes'][0])
         pos = networkx.shell_layout(dbGraph)
         networkx.draw_networkx(dbGraph, pos,with_labels=True, node_size=10000)
         #plot.savefig("Draw_DeBruijnGraph.png")
         plot.show()

    else:
        edges = []
        nodes = []
        for i in graph['edges']:
           if(len(i[0]) <= len(i[1])):
              edges.append( i[0]+ i[1][ len(i[0]) - 1 : ] )
           else:
              edges.append( i[0]+ i[1][-1] )
           nodes.append(tuple(i))
           dbGraph.add_edge(i[0],i[1])

        pos = networkx.shell_layout(dbGraph)

        # 2 Differnet Styles
        #networkx.draw_networkx(dbGraph,pos,with_labels=True,node_size=3000)
        networkx.draw_networkx(dbGraph,pos,with_labels=True,connectionstyle='arc3, rad = 0.1',node_size=3000)

        graphNodes = dict(zip(tuple(nodes), edges))
        networkx.draw_networkx_edge_labels(G = dbGraph,pos= pos, label_pos = 0.6,font_color = 'k',clip_on = True,edge_labels=graphNodes)
        #plot.savefig("Draw_DeBruijnGraph.png")
        plot.show()

#Compacting/Merging Nodes to decrease graph size

def DeBruijnGraph_Compacted(reads,k):

    graph,dicti=DeBruijnGraph(reads,k)    

    Flag = True

    while (Flag) :

        iteration = 0

        for i,j in dicti.items():

             iteration += 1

             #Sort the dictionary by value in a descending order, so that end node of adjacency list [] allows the iteration variable
             #to be equal to the size of the final dictionary to exit the while loop
             dicti = {k: v for k, v in sorted(dicti.items(), key=lambda item: item[1],reverse=True)}

             #Check Output Degree of this node is equal to 1
             if(len(j) == 1):  
                 curr = dicti[i][0]
                 count = 0

                 #If yes, then check that its only output node is of input degree of one
                 for val in dicti.values():
                     if curr in val:
                        count = count + 1
                     if(count > 1):
                         break

                 #Then, start compacting by overrididng existing dictionary
                 if(count == 1):

                    newnode = i + curr[ k - 2 : ] #New node string

                    dicti[newnode] = dicti[j[0]] #Add the children of the second node to the new node and add the new node to the dictionary

                    for key,val in dicti.items(): #Replace first node from all adjacency lists by the new node 
                        newChildren = []
                        for element in val:
                            if(element != i):
                               newChildren.append(element)
                            else:
                               newChildren.append(newnode)
                        dicti[key] = newChildren

                    #Delete the 2 nodes from the dictionary
                    del dicti[i]
                    del dicti[curr]

                    #Start iterarting again on the newly edited dictionary, to use previous edits
                    break

             else:

                 if(iteration == len(dicti)): #If reached the end of the dictionary and no new compacting happened, exit 
                    Flag = False
                    break

    #Form new graph from new adjaceny lists dictionary for visulaization
    newGraph = {'nodes':[],'edges':[]}

    for i,j in dicti.items():
        newGraph['nodes'].append(i)
        for val in j:
            newGraph['edges'].append([i,val])

    return newGraph,dicti

#Testing Functions

graph,dicti=DeBruijnGraph(['TTACGTT','CCGTTA','GTTAC','GTTCGA','CGTTC'],5)
newGraph,newDicti = DeBruijnGraph_Compacted(['TTACGTT','CCGTTA','GTTAC','GTTCGA','CGTTC'],5)

#graph,dicti=DeBruijnGraph(['ATGG', 'TGCC', 'TAAT', 'CCAT', 'GGG', 'GGATG', 'ATGTT'],3)
#newGraph,newDicti = DeBruijnGraph_Compacted(['ATGG', 'TGCC', 'TAAT', 'CCAT', 'GGG', 'GGATG', 'ATGTT'],3)

#graph,dicti=DeBruijnGraph(["CGATTCTAAGT"],4)
#newGraph,newDicti = DeBruijnGraph_Compacted(["CGATTCTAAGT"],4)

print("Adjacency List Before Compacting Nodes: \n",dicti)

print("\n")

print("Old Graph is: \n",graph)

print("\n")

visualizeDBGraph(graph)

print("Adjacency List After Compacting Nodes: \n",newDicti)

print("\n")

print("New Graph is: \n",newGraph)

print("\n")

visualizeDBGraph(newGraph)


       