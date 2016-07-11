from visual import *

class solid:
    """
    This class allows to draw a rectangle 3D solid, giving the center position,
    the size of each side and the axis orientation.
    """

    def __init__(self, size, center, axis, num):
        """
        Method which initialize a solid instance.
        Input:
            size: array = [sx,sy,sz], dimentions of the solid.
            center: array = [cx,cy,cz], center of mass.
            axis: 3x3 array = [[ix,iy,iz],[jx,jy,jz],[kx,ky,kz]], orientation.
        """

        self.size = vector(size)
        self.center = vector(center)
        self.i = vector(axis[0]).norm()
        self.j = vector(axis[1]).norm()
        self.k = vector(axis[2]).norm()
        self.color = (0.62,0.16,0.11)
        self.list_of_vertices = []
        self.list_of_edges = [[],[],[],[]]
        self.list_of_curves = []
        self.coef = [[1,-1,-1],[1,1,-1],[1,1,1],[1,-1,1],[-1,-1,-1],[-1,1,-1],[-1,1,1],[-1,-1,1]]
        self.tray = False
        self.dots = points(size = 2, color = color.red, shape = 'round')
        self.dotv = points(size = 2, color = color.blue, shape = 'round')
        #self.num = label(pos = self.center, color = color.blue, text = str(num))


    def draw(self):
        """
        Method which draws the solid.
        """
        self.obj = box(pos = self.center,axis = self.i,up = self.j,size = self.size,color = self.color, make_trail = True)
        #self.draw_edges()        
        if self.tray:
            self.dots.append(self.center)
            self.dotv.append(self.list_of_vertices[2])

    def get_vertices(self):
        """
        Method which gets the initial vertices of the solid.
        Each vertex is computed as follows:
            V_n = OM + c_nx*i + c_ny*j + c_nz*k, n=0:7
            where
            OM : segment between the origin and the center point of the brick.
            c_nx : coefficient of the vertex n in the direction i.
            c_ny : coefficient of the vertex n in the direction j.
            c_nz : coefficient of the vertex n in the direction k.
        """
        for i in range(8):            
            vertex_i = self.coef[i][0]*(self.size[0]/2)*self.i
            vertex_j = self.coef[i][1]*(self.size[1]/2)*self.j
            vertex_k = self.coef[i][2]*(self.size[2]/2)*self.k
            vertex =  self.center + vertex_i + vertex_j + vertex_k
            self.list_of_vertices.append(vertex)

    def draw_edges(self):
        """
        Method which draws each edge of the solid.
        """
        
        self.get_vertices()

        vertices_list = [[0,1,2,3,0],[4,5,6,7,4],[1,5,6,2,1],[0,4,7,3,0]]

        for i in range(len(vertices_list)):
            for j in range(len(vertices_list[i])):
                self.list_of_edges[i].append(self.list_of_vertices[vertices_list[i][j]])
                
        for edge in self.list_of_edges:
            self.list_of_curves.append(curve(pos = edge, color = color.white))


    def set_position(self, center, axis):
        """
        Method which set the position and the axis of the solid.
        """

        self.center = vector(center)
        self.obj.pos = vector(center)        
        if self.tray:
            self.dots.append(self.center)
            self.dotv.append(self.list_of_vertices[2])       

        self.i = vector(axis[0]).norm()
        self.j = vector(axis[1]).norm()
        self.k = vector(axis[2]).norm()

        self.obj.axis = self.size[0]*self.i
        self.obj.up = self.j
        
    def set_vertices(self):
        """
        Method which set the vertices of the solid.
        """
        for i in range(8):            
            vertex_i = self.coef[i][0]*(self.size[0]/2)*self.i
            vertex_j = self.coef[i][1]*(self.size[1]/2)*self.j
            vertex_k = self.coef[i][2]*(self.size[2]/2)*self.k
            vertex =  self.center + vertex_i + vertex_j + vertex_k
            self.list_of_vertices[i] = (vertex)

    def set_edges(self):
        """
        Method which set the edges of the solid.
        """
        
        self.set_vertices()

        vertices_list = [[0,1,2,3,0],[4,5,6,7,4],[1,5,6,2,1],[0,4,7,3,0]]

        for i in range(len(vertices_list)):
            for j in range(len(vertices_list[i])):
                self.list_of_edges[i][j] = self.list_of_vertices[vertices_list[i][j]]
        
        for i in range(len(self.list_of_curves)):
            self.list_of_curves[i].pos = self.list_of_edges[i]


