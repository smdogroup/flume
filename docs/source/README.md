# Flume

Flume is an open-source framework for multidisciplinary analysis and adjoint evaluations based on directed acyclic graphs. It is intended to be a lightweight framework that affords users with a degree of flexibility in implementing their own analysis and optimization procedures while providing a common underlying structure. _Flume_ also has a set of interfaces between some common optimizers, which allow the user to easily setup an optimization problem that corresponds to their specific needs.
 
## Getting Started

For new users, it is recommended that you read the _Flume_ [overview notebook](https://github.com/smdogroup/flume/blob/main/examples/flume_overview.ipynb), which outlines in detail the critical aspects of the framework in context of an example. This discusses the nomenclature that is used within the framework and illustrates how to construct _State_, _Analysis_, and _System_ objects for use with a _Flume_ interface. Also, for a review of other useful methods that are included for _Analysis_ classes, check out the [methods demonstration notebook](https://github.com/smdogroup/flume/blob/main/examples/analysis_methods_demo.ipynb). The [examples gallery](https://github.com/smdogroup/flume/tree/main/examples) also includes several different types of problems, which showcase the construction of many types of _Analysis_ classes, the assembly of a _System_, and how to interface with an optimizer interface. 

## Installation

To install _Flume_, the only necessary dependency is a working version of Python that is at least version 3.10. It is generally encouraged that you use a virtual environment. If you do not have a virtual environment already, navigate to the directory where you want to use _Flume_ and execute the following.
```
python -m venv venv
source venv/bin/activate
```
Then, with the virtual environment activated, simply install the package from the GitHub repository with
```
pip install git+https://github.com/smdogroup/flume.git
```
This will install the base classes, optimizer interfaces, and utilities into a Python package named "flume". For those who want to make changes to the code or access the unit tests and examples, you can also clone the repository and then install in editable mode.
```
cd {your_chosen_git_directory}
git clone https://github.com/smdogroup/flume.git
cd flume
python -m venv venv
source venv/bin/activate
pip install -e .
```
After following either of these processes, you should be able to access the base classes and get started with developing your own scripts within the framework! It should be noted that ParOpt is a separate package that is not included with _Flume_ during the build process. For those who want to utilize this interface, we refer you to the [ParOpt installation instructions](https://github.com/smdogroup/paropt/blob/master/INSTALL.txt).
