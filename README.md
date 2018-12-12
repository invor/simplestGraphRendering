The Simplest Graph Renderer there is

Dependencies:

    * GLEW
    * OpenGL
    * GLFW

To build super lazy on Linux:

    ./build.sh 


To build:
    
    mkdir build
    cd build
    cmake ../

To run (in build):

    ./simple -gf resources/test.gl

gl format:

```
number of nodes
number of edges
X Y of node1
...
X Y of nodeN
nodeFrom nodeTo width color
...
nodeFrom nodeTo width color
```

Example:

Two nodes: node0 (1, 1) and node1 (2, 2). One arc (node0, node1) of 1px width with color 1

```
2
1
1 1
2 2
0 1 1 2
```
