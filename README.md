## Introduction

The numerical solution of differential-algebraic equations can be expensive, if the underlying model has high dimensions. When using Runge-Kutta methods, large nonlinear systems of equations have to be solved during the iteration. The main focus for this algorithm is to utilize half-explicit Runge-Kutta methods for a special class of differential-algebraic equations to make these computations more efficient. We exploit the structure and properties of the given differential-algebraic equations to construct an algorithm that allows to solve reduced nonlinear systems of equations. The algorithm was tested with multiple examples from the fields of mechanical, electrical and chemical sciences to validate the approach.

## Setup and Execution

Please read Manual.pdf contained in the Manual folder for a short introdution to the setup of the program. For a thorough introduction to the topic and the mathematical background we refer to HERKosiDAE.pdf in the Manual folder. Please inform yourself carefully before using the software - your model may not satisfy all prerequisites.

## Feedback and Support

Please feel free to open an issue if you find a bug or seek support. Also, I encourage you to fork the project.

## Release Notes

Version 1.0 - Initial Release

## License

Copyright 2017 Simon Michael Baese

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.