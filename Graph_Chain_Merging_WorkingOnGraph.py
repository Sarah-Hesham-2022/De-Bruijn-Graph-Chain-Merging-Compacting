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
    #print(graph)
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

        # 2 Different Styles
        #networkx.draw_networkx(dbGraph,pos,with_labels=True,node_size=3000)
        networkx.draw_networkx(dbGraph,pos,with_labels=True,connectionstyle='arc3, rad = 0.1',node_size=3000)

        graphNodes = dict(zip(tuple(nodes), edges))
        networkx.draw_networkx_edge_labels(G = dbGraph,pos= pos, label_pos = 0.6,font_color = 'k',clip_on = True,edge_labels=graphNodes)
        #plot.savefig("Draw_DeBruijnGraph.png")
        plot.show()

#Compacting/Merging Nodes(Chains) to decrease graph size by working on graph dictionary

def DeBruijnGraph_Compacted(graph):

    kmer = len(graph['nodes'][0]) + 1
    Flag = True

    while(Flag and graph['edges']):

       iteration = 0

       #Separate nodes of edges out/in
       outs = []
       ins  = []
       for i in graph['edges']:
         outs.append(i[0])
         ins.append(i[1])

       for i in range(len(graph['edges'])):
    
         iteration += 1

         #Check Output Degree of this node is equal to 1
         if(outs.count(graph['edges'][i][0]) == 1):

             #If yes, then check that its only output node is of input degree of one
             if(ins.count(graph['edges'][i][1]) == 1):

                 #New node string
                 newNode = graph['edges'][i][0]+graph['edges'][i][1][ kmer - 2 : ] 

                 #Collect all existing children of the second node 
                 children = [] 
                 for j in range(len(outs)):
                     if(outs[j] == graph['edges'][i][1]):
                         children.append(ins[j])
                
                 #Then, start compacting by overriding existing edges nodes
                 for j in range(len(graph['edges'])):

                     if(i != j):

                         if(graph['edges'][j][1] == graph['edges'][i][0]): #Replace first node by new node
                             graph['edges'][j][1] = newNode

                         if(graph['edges'][j][0] == graph['edges'][i][1]): #Replace second node by new node and making all its children assigned to new node
                            graph['edges'][j][0] = newNode
                            graph['edges'][j][1] = children[-1]
                            children.remove(children[-1])

                 #Remove the edge from the graph dictionary 
                 graph['edges'].remove(graph['edges'][i])

                 #Start iterating again on the newly edited graph dictionary, to use previous edits of edges list
                 break
         else:
            if(iteration == len(graph['edges'])): #If reached the end of the edges list and no new compacting happened, exit 
               Flag = False
               break
 
    #Form new nodes for the graph from newly compacted nodes
    nodes = []

    if(graph['edges'] == []): #Incase we have one node only after compacting graph, as in test case 3

       node = graph['nodes'][0]

       for i in range(1,len(graph['nodes'])):
         node += graph['nodes'][i][-1]

       graph['nodes'] = [node]

    else:

      for i in graph['edges']:
         nodes.append(i[0])
         nodes.append(i[1])

      nodes = set(nodes)
      graph['nodes'] = list(nodes) 

#Testing Functions

graph,dicti=DeBruijnGraph(['TTACGTT','CCGTTA','GTTAC','GTTCGA','CGTTC'],5)
#graph,dicti=DeBruijnGraph(['ATGG', 'TGCC', 'TAAT', 'CCAT', 'GGG', 'GGATG', 'ATGTT'],3)
#graph,dicti=DeBruijnGraph(["CGATTCTAAGT"],4)

print("Old Graph is: \n",graph)

print("\n")

visualizeDBGraph(graph)

DeBruijnGraph_Compacted(graph)

print("New Graph is: \n",graph)

print("\n")

visualizeDBGraph(graph)

       
