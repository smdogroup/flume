# FLUME OVERVIEW

This repository contains the base classes for the analysis and optimization framework used within the SMDO group. It is intended to be minimalistic and provide the user with lots of flexibility for implementing their own analysis and optimization procedures while still providing a common structure. The key terminology and a description of each class is given below. Some examples are provided to establish a foundation for how to utilize this framework, but it is ultimatley structured to be flexible and fit the user's individual needs. 

## Class Descriptions

*This section is a work in progress as the classes and their intentions are finalized for the framework.*

## Key Terminology

Within this framework, there exist a few naming conventions to distinguish between different components of the framework. For an analysis object, the following convention is utilized.

- **Parameter:** Numerical, Boolean, string, or other values that are nominally set during object instantiation and act as inputs to an analysis procedure. These are intended to be set once and not changed after instantiation, so they are generally not considered as design variables during optimiztaion.
- **Variable:** Numerical quantities (either floats or Numpy arrays) that are candidates for design variables during optimization and are inputs to an analysis procedure. Variables are numerical quantities that *may* change during optimization, but do not have to if the variable is not declared as a design variable. All variables are instances of a State object, whose attributes contain information about the value, data type, shape, description, derivative value, and source information. If sub-analyses exist for another analysis class, then variables may be sourced from another object if its variable name shares the same name as an output of a connected sub-analysis. 
- **Output:** The construct of an output is similar to that of a variable, except is is an output of an analysis procedure. Outputs are also instance of State objects, so information about an output is accessed in the same way as a variable. During optimization, seed values are set on outputs for the adjoint analysis. Outputs of sub-analysis classes can be linked to variables for top-level analyses if they share the same name. 

