import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import sys

def showMeshPlot(nodes, elements, prefix):

    y = nodes[:,0]
    z = nodes[:,1]

    triang = tri.Triangulation(y, z, triangles)
    
    # plot nodes
    #plt.plot(y,z, marker="o", ls="", color="crimson")

    plt.triplot(triang, 'go-', lw=1.0, label=prefix)


def getNodesAndTriangles(prefix):
    nodes = None

    with open(prefix + "_vertex.dat") as vert_file:
        nodes = np.array([ [ float(coord) for coord in line.split(' ')] for line in vert_file] )

    triangles = None
    with open(prefix + "_triangles.dat") as triang_file:
        triangles = np.array([  [int(index) for index in line.split(' ')] for line in triang_file] )
    
    return nodes, triangles

if len(sys.argv) < 2:
    print("Fourning prefix")
    sys.exit()

prefixes = sys.argv[1:]

plt.figure()
plt.gca().set_aspect('equal')

plt.title('This is the plot for meshes')
plt.xlabel('X Axis')
plt.ylabel('Y Axis')
    
for prefix in prefixes:
    nodes, triangles = getNodesAndTriangles(prefix)
    showMeshPlot(nodes, triangles, prefix)

plt.legend()

plt.show()
