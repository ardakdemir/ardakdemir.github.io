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

Initial commits to the BioSequenceGraph repository:

- Fixed several minor bugs in BioSequences repository.

- Added the subsequencing functionality to BioSequencess

- Added several function which will form the basis of constructing De Bruijn Graphs from SequenceGraphNodes.
