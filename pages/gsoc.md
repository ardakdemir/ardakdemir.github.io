---
layout: page
title: Google Summer of Code 2019 BioJulia
---



## De-Bruijn Graph Constructor Package for De-novo Genome Assembly

Mentor : [Ben J. Ward](http://www.earlham.ac.uk/ben-ward)

Relevant Google Summer of Code 2019 [webpage](https://summerofcode.withgoogle.com/projects/#5248824832950272)


### Project Abstract

De-novo sequence assembly is the process of constructing a contiguous long sequence out of shorter sub-sequences produced by sequencing platforms, without referring to a reference genome. It is an essential task in many biological studies today, including population and medical studies. The initial stages of de-novo assembly require the construction of a de-bruijn graph (DBG) from sequencing reads, the compression of a de-bruijn graph into a unitig graph, and the compression of multiple unitigs and nodes into contigs, supported by evidence from mapping paired-end reads. A coherent ecosystem of computational tools and packages allow researchers to quickly implement and test their ideas. For bioinformatics, Julia already offers such an ecosystem in the form of the BioJulia & EcoJulia projects, and additional independent packages. This project will add sequence assembly tools to the BioJulia ecosystem, specifically: 1) DBG construction from reads, 2) UG construction from a DBG and 3) Constructing contigs using unitigs. These tools will allow researchers to quickly construct and analyze the contigs obtained from a set of reads.


### Community Bonding Period (6 - 27 May)

I have started connecting with other GSoC students and exploring the open source community.
It is a fascinating experience to be in a developing team for the Julia language.
Going over the source code of the relevant package I already gained a deeper understanding of the underlying concepts
behind the Julia language. Even though I had previous experience coding in Julia, I realized that I never truly understood some key
concepts especially related to Types. I believe that through this GSoC project I will gain valuable insights about the fundamental concepts of  Julia programming language and programming languages in general.

My primary aim during this period is to get familiar with the relevant packages and get a head start to the project if it is posssible.
I am spending most of my time reading the source code, looking at issue related to the packages (trying to fix them if possible) and applying what I have learnt to write some simple code.

Apart from my contributions to the BioSequenceGraph package, I am also storing my experiments with using the BioSequences package in a separate jupyter notebook. This notebook can be thought of as a coding  diary. As I go through the source code and add new functionalities I test them and create example codes. I believe at the end  of the  project this notebook will also be useful for users who would like to see some working examples. It can be considered as an extended cheat sheet for BioSequences package which  to the best of my knowledge does not exist online.

I started contributing to the package two weeks before  the end of the community bonding period. I started by adding some simple functionalities and going over some bugs (or concepts which I thought are bugs but turned out not to be :)).
Initial commits to the BioSequenceGraph repository:

- Fixed several minor bugs in BioSequences repository.

- Added the subsequencing functionality to BioSequence types. This functionality allows querying a dna sequences with a vector:

```
dna_1 = dna"ATGC"
dna_1[[1,3]]
```


- Random Kmer generator is implemented. This functionality is mainly useful for speeding up the testing of the code.
Similar functions are already available in BioSequences package in the randseq.jl file.

- Added several functionalities which will form the basis of constructing De Bruijn Graphs from SequenceGraphNodes.

- Initialized the DeBruijnGraph type and started writing the core functionalities.


#### Weekly meetings with Mentor

We are constantly in communication with my mentor Dr. Ben J. Ward through emails and issues on GitHub.
Apart from that, we are making weekly meetings where we discuss about how to proceed and I go over the work I have done.
These meetings are important for me especially at the initial phases of the project as they help me a lot to understand basic concepts of the BioSequence packages  and all the underlying design concepts. Also it is a valuable time for me to ask my questions.


#### DeBruijnGraph



Before waiting for the end of community bonding period we have started implementing the DeBruijnGraph package along with some functionalities.

In our formalization of the DeBruijnGraph we represent each dna sequence (of arbitrary length) as a Node on the Graph. There exists a Link between two Nodes which represent the overlaps between nodes/sequences. Each Node is of type SequenceGraphNode and each Link is of type SequenceGraphLink. These types have their special constructors and functionalities. Below is an example of a DeBruijnGraph where sequences are represented as nodes.

<a href="../pages/publpics/debru3.png">
    <img src="../pages/publpics/debru3.png"
          title="DeBruijnGraph" alt="dbg"  height="420" width="420"/></a>

DeBruijnGraph is a special type of SequenceGraph with its own constraints ( e.g. links between arbitrary nodes are not allowed).
It is made up of two fields:

- A vector of Nodes Vector

- A vector of vector of Links

```

struct DeBruijnGraph
    nodes::Vector{SequenceGraphNode}
    links::Vector{Vector{SequenceGraphLink}}
end

```

Implemented Functionalities:

**DeBruijnGraph Constructor:**

At the beginning, we have designed a constructor to represent a static DeBruijnGraph where we assume no operation such as node merging will later be performed. This constructor receives as input a list of kmers of type Kmer{T,K} where T denotes the NucleicAcidType (DNA or RNA) and K denotes the length k. Using these kmers the constructor checks for overlaps of length $k-1$ and creates directed Links from source to destination. Source is the node which has the overlap as a suffix and destination is the node which has the overlap as a prefix.


**Query Functions:**

Before moving into the next milestone which is to build a Unitig Graph (UG) from the kmers in the DeBruijn Graph, we will implement some more core functionalities to ease implementation of the more complicated stages of the package. Some of these core functionalities are query functions which are necessary for finding paths on dbg suitable for merging. Below is the list of these queries:

- count_indegree : counts the number of incoming edges to a vertex
- count_outdegree : counts the number of outgoing edges from a vertex
- is_a_path  : given a sequence of nucleotides checks whether a path exists that yields the given sequences
- is_in_node : this function will be useful for detecting substrings of a node label which will be useful during path checking as explained below.

These query functions will be useful during UG construction.

For now the counters are calculated assuming that we have a typical dbg as a directed graph. So for each node we look at the source and sink ends to count indegree and outdegree respectively.

**is_a_path** : If we formulate the dbg as static in that we initialize it with some kmers and do not update it, this query is trivial to implement. We can exhaustively check all the nodes that have the first k letters of the path as its label and move down on its children. Since we do not have multiple nodes for the same kmer this can take at most O(\|V\|^2) or O(\|V\|  x \|E\| ) even if we implement the query naively. However, if we allow edge collapsing/vertex merging, a more involved solution is necessary to handle all possible graphs. Imagine the following intermediate representation that occurs after merging two 3mers "ATG" and "TGT":


<img src="publpics/is_path_1.png?" alt="is_path" width="300" height="300">


If we check all the edges and look for an exact match than the algorithm will return a false negative in the cases where we search for a sequence that starts with "TGT". For this reason we implement the node_search as a query over the substrings of the label of a node rather than an exact match.

### Official Coding Period (28 May -)

We have already implemented the core functionalities and defined a basic de Bruijn Graph type for BioJulia during the community bonding period.



After our  discussions with Dr. Ben Ward we decided to revise the constructor function for the dbg. The new design represents each kmer and its reverse complement using the same node and uses (+) and (-) end of a node for labeling edges between nodes.

Below is the  pseudocode for the  new  constructor :  

```
Make an empty graph

Make two empty Vector{Tuple{DNAKmer{K-1}, Int64}}. Call one `kmer_ovl_bw_nodes`, call the other one `kmer_ovl_fw_nodes`.

For each kmer in kmerset...
    Make the canonical form of the kmer. (there should be a canonical method in BioSequences).
    Add the canonical kmer to the nodes of the graph, and note it's ID.
    Take the prefix (k-1) of the canonical kmer.
    If the prefix is canonical, push the tuple (canonical(prefix), +nodeid) to `kmer_ovl_fw_nodes`.
    Else push the tuple(canonical(prefix), +nodeid) to the `kmer_ovl_bw_nodes`.
    Take the suffix (k-1) of the canonical kmer.
    If the suffix is canonical, push the tuple (canonical(suffix), -nodeid) to `kmer_ovl_bw_nodes`.
    Else push the tuple (canonical(suffix), - nodeid) to `kmer_ovl_fw_nodes`.
end

Sort `kmer_ovl_fw_nodes`.
Sort `kmer_ovl_bw_nodes`.

for kbn in kmerovl_bw_nodes
        for kfn in kmerovl_fw_nodes
            if first(kbn) == first(kfn)
                add link to graph with source=last(kbn), destination=last(kfn) and distance=-k+1
            end
        end
end

Return the graph.
```

This constructor returns  a dbg which also consists of two vectors.  Yet the  main difference is that we only represent the kmers  in their canonical form and represent each kmer and its reverse complement with the  same node, both decreases the memory footprint of  the graph. Canonical form of  a kmer  is the lexicographically lesser  one of a  kmer  and  its reverse complement.

Example:

```
julia> kmer  = Kmer{DNA,4}("ACTT")
DNA 4-mer:
ACTT

julia> canonical(kmer)
DNA 4-mer:
AAGT

julia> kmer  == canonical(kmer)
false
```

kmer ACTT is not in the canonical form as the reverse complement AAGT is lexicographically  smaller. Thus both kmers ACTT, and AAGT is represented   using only a single node as (+) AAGT (-). The plus  end  denotes the prefix in canonical form and also denotes suffix in the non-canonical form.


We have finalized the implementation of the  new constructor and did some  initial tests. For random 5 3mer listed  below:


```
5-element Array{Kmer{DNA,3},1}:
 AAC
 TCA
 GAG
 AAT
 GCG
```

The constructor generates the following dbg :

```

DeBruijnGraph(SequenceGraphNode[SequenceGraphNode{Kmer{DNA,3}}(AAC, true), SequenceGraphNode{Kmer{DNA,3}}(TCA, true), SequenceGraphNode{Kmer{DNA,3}}(CTC, true), SequenceGraphNode{Kmer{DNA,3}}(AAT, true), SequenceGraphNode{Kmer{DNA,3}}(CGC, true)], Array{SequenceGraphLink,1}[[], [SequenceGraphLink(2, -3, -2)], [], [], [], []], 3)
```

which consists of 5 nodes with labels AAC, TCA, CTC, AAT and CGC. As one  can see 3mers GAG and GCG are replaced with their canonical  forms in the  new graph.


Next step in our schedule is to revise the  is_a_path and other related dbg functionalities  to be  compatible with the new constructor.
Then, we are planning to implement some functions for detecting the simple paths in a dbg. is_simple_path function is already initialized for the DeBruijnGraph type.

Below is an example of a case where we need to consider during addition of the edges. If the suffix/prefix are equivalent to their reverse_complement both directions  should be taken  into account otherwise some paths will be missing in the dbg.s

<a href="../pages/publpics/GSOC.png">
    <img src="../pages/publpics/GSOC.png"
          title="DeBruijnGraph" alt="dbg"  height="420" width="420"/></a>



Our plan is to combine the nodes on a simple path to create unitigs from kmers. This will increase the nucleotide length of some nodes.
A simple path is defined as a series of nodes where the internal nodes have incoming and outgoing degrees of one and the start node and the end node have incoming and outgoing degrees >1 respectively. Thus merging the nodes in this path does not create ambiguity.



#### Combining Simple paths

Next step in our milestone is to collapse the nodes and merge all the intermediate  edges.
However  this step requires special care as our current design traverses the links  and edges by their node_id.
The vector structure  should be  converted to either a  hash table or we must make updates to the all elements of the  nodes  and links vectors  which is  very costly.

That is why I decided to change the design of the dbg. Previously we represented nodes and links as vectors  but now they will be  represented as dictionaries. This way we  do not have to update the index information. Dictionary data structure allows us to cleanly access each node and link after unitigging. 
