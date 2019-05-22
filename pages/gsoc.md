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

'''
dna_1 = dna"ATGC"
dna_1[[1,3]]
'''

- Random Kmer generator is implemented. This functionality is mainly useful for speeding up the testing of the code.
Similar functions are already available in BioSequences package in the randseq.jl file.

- Added several functionalities which will form the basis of constructing De Bruijn Graphs from SequenceGraphNodes.

- Initialized the DeBruijnGraph type and started writing the core functionalities.


#### Weekly meetings with Mentor

We are constantly in communication with my mentor Dr. Ben J. Ward through emails and issues on GitHub.
Apart from that, we are making weekly meetings where we discuss about how to proceed and I go over the work I have done.
These meetings are important for me especially at the initial phases of the project as they help me a lot to understand basic concepts of the BioSequence packages  and all the underlying design concepts. Also it is a valuable time for me to ask my questions.


### DeBruijnGraph

Before waiting for the end of community bonding period we have started implementing the DeBruijnGraph package along with some functionalities.

[[ardakdemir.github.io/assets/publpics/debru2.png|alt=debru]]
