---
title: "Flume: A Lightweight Framework for Multidisciplinary Design Optimization based on Directed Acyclic Graphs"
tags:
  - Multidisciplinary design optimization
  - Optimization framework
  - Python
authors:
  - name: Cameron S. Smith
    orcid: 0000-0002-3397-2223
    affiliation: 1
    corresponding: true
  - name: Graeme J. Kennedy
    orcid:
    affiliation: 1
affiliations:
  - name: Georgia Institute of Technology, United States
    index: 1
    ror: 01zkghx44
date: 9 October 2025
bibliography: paper.bib
---

# Summary

In an effort to solve a variety of multidisciplinary optimization (MDO) problems with gradient-based techniques, we have developed a tool to facilitate the system-level construction and execution for both the forward and reverse analyses.
This framework, entitled _Flume_, is designed around systems that can be described by directed acyclic graphs (DAG).
While this specific architecture prevents systems that have coupled relationships, it is applicable for a wide variety of problems, including topology and trajectory optimization.
Following the DAG structure, each node represents an individual analysis that needs to be performed for the MDO problem, and the edges represent connections between analyses and denote the flow of information.
To ensure that the framework is extensible, lightweight, and minimalist, three base classes have been implemented in Python to capture the necessary functionality: _State_, _Analysis_, and _System_.
The first class, _State_, simply provides an object that stores numerical data along with some additional metadata, including type, shape, and object source information.
_Analysis_ is the foundation of Flume, and its primary task is to execute the forward and adjoint procedures to obtain the outputs and derivatives needed for optimization.
Finally, _System_ provides a set of methods to declare the objective function, constraints, and design variables that will be used within the optimization problem, as well as a means to visualize the DAG.
To utilize Flume, a user's primary responsibility is to construct the individual analyses, inherited from the _Analysis_ base class, that are needed for their _System_.
By adhering to a few architectural requirements when scripting these analyses, the framework's backend will automatically connect outputs to variables that share the same name.
This provides the user with a streamlined workflow, enabling them to focus on implementing new features and procedures instead of managing the integration.

**potentially insert image here that would help visualize the DAG structure, specifically as it relates to Flume**

# Statement of Need

Numerical optimization is a field with a diverse range of applications, and, in engineering, is often used as a means of formulating problems in a formal, mathematical manner to identify optimal designs.
The development of frameworks for organizing and solving these problems, which are often multidisciplinary in nature, is a topic that has been addressed by others in the past.
Investigations of these frameworks has exposed a set of requirements for MDO frameworks, including modularity, intuitive user interfaces, object-oriented principles, and minimal overhead [@salas1998framework, @padula2006multidisciplinary].
These attributes, among others, are critical components to ensuring that an MDO framework is extensible to a multitude of disciplines and accessible by a variety of individuals.
While there are many MDO frameworks that have been developed since the late 20th century, a few are listed below to highlight the current state of the art.

- _ASTROS_: one of the first examples of MDO frameworks, _ASTROS_ performs preliminary structural design using numerical optimization based on the finite-element method [@astros]
- _DAKOTA_: developed by Sandia National Laboratories, _DAKOTA_ is a software suite written in C++ that provides methods for a variety of analyses, including gradient-based and gradient-free optimization [@dakota]
- _Isight_: a commercial tool by Dassault Systèmes that utilizes and object-oriented approach to connect a variety of simulation-based models [@isight]
- _pyMDO_: an object-oriented solution written in Python, providing a means of defining MDO problems independently of the optimization algorithm used with the use of inheritance and operator overloading [@martins2009pymdo]
- _OpenMDAO_: another open-source framework constructed within Python that utilizes gradient-based techniques for optimization of systems constructed with distinct components [@gray2019openmdao]

As evident from this list, users are presented with many viable options to perform numerical optimization for their system of interest.
_Flume_ provides another resource for those solving MDO problems, but it was particularly architected to facilitate gradient-based design optimization with the adjoint method.
For each _Analysis_ that a user wants to integrate into their model, they are responsible for providing the forward and adjoint analysis techniques that are specific to the computations they want to perform.
The backend of _Flume_ will then orchestrate the assembly of total derivatives through the adjoint method by accumulating contributions from the objective function and any constraints for the system.
Thus, each individual _Analysis_ must only consider its respective variable-output combinations, and the adjoint variables will propagate information through the model to compute the necessary total derivative information.
This assists the user by allowing them to prioritize the development of the individual _Analysis_ objects, and the framework, paired with the proper _System_ construction, will address the data distribution for both computational directions.

# Applications of _Flume_

To date, the framework has primarily been tested for two applications beyond the simpler examples provided within the repository to detail the construction of _Analysis_ disciplines and _System_ architectures.
In a nascent stage, it was utilized to perform optimization under uncertainty for a next-generation Mars helicopter.
While the core functionality of the framework remains the same, the changes to the _Analysis_ and addition of the _System_ base classes cause this example to be out of date.
Now, _Flume_ has primarily been tested in the field of topology optimization, with specific applications for initial post-buckling behavior and inverse design problems.
The network complexity of these systems, specifically regarding the flow of information between distinct analyses, emphasizes the importance of a tool like _Flume_.
Both of these demonstrations were instrumental in designing the framework and ultimately has resulted in its present state, and publications on both of these topics are in preparation for submission.

# _Flume_ by Example: Constrained Rosenbrock

TODO: include the code here for the constrained Rosenbrock example, highlighting the key aspects that contribute to the Analysis/System construction and the optimization with SciPy

<!-- NOTE: this was used for the grant proposal for Flume (unsure if this got submitted); not deleting this yet, but will likely not add to the paper. To further its development and ensure its applicability for an increased number of disciplines, we propose to extend Flume in two main areas: parallelism and remote procedure calls.
While there are no inherent limitations that prevent a user from enabling parallelism within specific analyses, this has not been a development focus to date.
However, to support large computing effort, we aim to implement and test parallel execution features.
Additionally, since one of the primary benefits of Flume is the automatic management and communication between _Analysis_ instances in a _System_, it also lends itself to remote procedure calls.
In scenarios where certain information or _Analysis_ instances can only be performed on secure machines but cannot be conducted entirely in this location, a remote procedure call is necessary and could be facilitated with a framework that automates data distribution.
With these additional applications, the framework would be further equipped to solve a breadth of MDO problems that are discipline independent. -->

<!-- # To-do List

- DONE: Draft two paragraphs and potentially an image explaining Flume. Intent is to use this for the grant proposal, as well as the _Summary_ or _Statement of Need_ section for JOSS submission
  - Should mention stuff about parallelism and remote procedures for the grant proposal version
- TODO: Start to populate this file with the metadata for the JOSS submission
- DONE: Document the code that has missing docstrings; ensure to follow the proper format with arguments, kwargs, etc.
- DONE: Construct an interface between Flume and Scipy, which can be used upon installation of Flume (unlike ParOpt, which requires additional dependencies)
- DONE: Implement a few more unittests that will perform optimization checks (for Rosenbrock and something else) and integration with Scipy minimize
- TODO: update examples directory to reflect the current structure of the code; i.e. remove any examples that do not use updated Analysis/System classes; make sure that there are a few different examples covering a variety of disciplines
    - examples to include: Hohmann transfer from OpenMDAO (opt), Rosenbrock (opt), multi-objective (demonstrating more complex system construction), methods demo Jupyter notebook, COPS sphere problem (opt)
- TODO: add a script/jupyter notebook/something that details the useful helper methods and outlines how to start constructing a Flume example (want to incorporate adjoint test methods, network visualization, etc.)
- TODO: update README.md with installation instructions, and point the user to various helpful portions of the documentation (e.g. new users check out the Jupyter notebook, for descriptive details point to the paper, etc.)
-->
