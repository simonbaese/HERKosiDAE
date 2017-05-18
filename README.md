## Introduction

The numerical solution of differential-algebraic equations can be expensive, if the underlying model has high dimensions. When using Runge-Kutta methods, large nonlinear systems of equations have to be solved during the iteration. The main focus for this algorithm is to utilize half-explicit Runge-Kutta methods for a special class of differential-algebraic equations to make these computations more efficient. We exploit the structure and properties of the given differential-algebraic equations to construct an algorithm that allows to solve reduced nonlinear systems of equations. The algorithm was tested with multiple examples from the fields of mechanical, electrical and chemical sciences to validate the approach.

## Setup and Execution

Please read Manual.pdf contained in the Manual folder for a short introdution to the setup of the program. For a thorough introduction to the topic and the mathematical background we refer to HERKosiDAE.pdf in the Manual folder. Please inform yourself carefully before using the software - your model may not satisfy all prerequisites.

## Examples

The Examples folder contains different systems and setups which are also documented in HERKosiDAE.pdf, i.e.

**Mathematical Pendulum**

![Mathematical Pendulum](Examples/Mathematical%20Pendulum/MathematicalPendulum.PNG)

**Linear Circuit with one CV Loop**

![Linear Circuit with one CV Loop](Examples/Linear%20Circuit%20with%20one%20CV%20Loop/LinearCircuitOneCVLoop.PNG)

**Spring Mass Chain**

![Spring Mass Chain](Examples/Spring%20Mass%20Chain/SpringMassChain.PNG)

## Feedback and Support

Please feel free to open an issue if you find a bug or seek support. Also, I encourage you to fork the project.

## Release Notes

Version 1.0 - Initial Release
