# Summary

In an effort to solve a variety of MDO problems with gradient-based techniques, we have developed a tool to facilitate the system-level construction and execution for both the forward and reverse analyses.
This framework, Flume, is designed around systems that can be described by directed acyclic graphs (DAG).
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

To date, the framework has been tested for two applications.
In a nascent stage, it was utilized to perform optimization under uncertainty for a next-generation Mars helicopter.
Now, Flume has primarily been tested in the field of topology optimization, with specific applications for initial post-buckling behavior and inverse design problems.
Both of these demonstrations were instrumental in designing the framework and ultimately has resulted in its present state.
To further its development and ensure its applicability for an increased number of disciplines, we propose to extend Flume in two main areas: parallelism and remote procedure calls.
While there are no inherent limitations that prevent a user from enabling parallelism within specific analyses, this has not been a development focus to date.
However, to support large computing effort, we aim to implement and test parallel execution features.
Additionally, since one of the primary benefits of Flume is the automatic management and communication between _Analysis_ instances in a _System_, it also lends itself to remote procedure calls.
In scenarios where certain information or _Analysis_ instances can only be performed on secure machines but cannot be conducted entirely in this location, a remote procedure call is necessary and could be facilitated with a framework that automates data distribution.
With these additional applications, the framework would be further equipped to solve a breadth of MDO problems that are discipline independent.

# To-do List

- Draft two paragraphs and potentially an image explaining Flume. Intent is to use this for the grant proposal, as well as the _Summary_ or _Statement of Need_ section for JOSS submission
  - Should mention stuff about parallelism and remote procedures for the grant proposal version
- Start to populate this file with the metadata for the JOSS submission
- Document the code that has missing docstrings; ensure to follow the proper format with arguments, kwargs, etc.
- Construct an interface between Flume and Scipy, which can be used upon installation of Flume (unlike ParOpt, which requires additional dependencies)
- Implement a few more unittests that will perform optimization checks (for Rosenbrock and something else) and integration with Scipy minimize
