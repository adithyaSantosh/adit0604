#Write a Python script that models a Mobius strip using parametric equations and computes key geometric properties.

#1. Importing of libraries

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial.distance import euclidean
from scipy.integrate import simps

#2. Class Definition

#R: radius of the center circle of the strip.
#w: width of the strip.
#n: resolution — number of discrete steps in the grid.
class MobiusStrip:
    def __init__(self, R=1.0, w=0.2, n=200):
        self.R = R
        self.w = w
        self.n = n
#3.Create a mesh grid in parameters 𝑢∈[0,2𝜋],u∈[0,2π] and 𝑣∈[−𝑤/2,𝑤/2],v∈[−w/2,w/2].
        self.u = np.linspace(0, 2 * np.pi, n)
        self.v = np.linspace(-w / 2, w / 2, n)
        self.U, self.V = np.meshgrid(self.u, self.v)
        self.X, self.Y, self.Z = self._compute_mesh()
        #Generates the 3D coordinates using the parametric equations.
        
#4. Compute Mesh
    def _compute_mesh(self):
        U, V = self.U, self.V
        X = (self.R + V * np.cos(U / 2)) * np.cos(U)
        Y = (self.R + V * np.cos(U / 2)) * np.sin(U)
        Z = V * np.sin(U / 2)
        return X, Y, Z
        #Applies the Möbius strip parametric equations to convert (u, v) → (x, y, z).
        
#5. Surface Area Calculation
    def compute_surface_area(self):
        #Computes the numerical surface area by integrating over the mesh.
        du = 2 * np.pi / (self.n - 1)
        dv = self.w / (self.n - 1)

        Xu = np.gradient(self.X, du, axis=1)
        Yu = np.gradient(self.Y, du, axis=1)
        Zu = np.gradient(self.Z, du, axis=1)

        Xv = np.gradient(self.X, dv, axis=0)
        Yv = np.gradient(self.Y, dv, axis=0)
        Zv = np.gradient(self.Z, dv, axis=0)

        Nx = Yu * Zv - Zu * Yv
        Ny = Zu * Xv - Xu * Zv
        Nz = Xu * Yv - Yu * Xv

        dA = np.sqrt(Nx**2 + Ny**2 + Nz**2)

        area = simps(simps(dA, self.v), self.u)
        return area
#6. Edge Length Calculation
    def compute_edge_length(self):
        edge_pts = [
            (self.X[i, :], self.Y[i, :], self.Z[i, :]) 
            for i in [0, -1]
        ]

        total_length = 0.0
        for X_edge, Y_edge, Z_edge in edge_pts:
            for i in range(1, len(X_edge)):
                p1 = (X_edge[i - 1], Y_edge[i - 1], Z_edge[i - 1])
                p2 = (X_edge[i], Y_edge[i], Z_edge[i])
                total_length += euclidean(p1, p2)
        return total_length
#7. Plotting the Strip
    def plot(self):
        fig = plt.figure(figsize=(10, 6))
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_surface(self.X, self.Y, self.Z, cmap='plasma', rstride=1, cstride=1, alpha=0.8)
        ax.set_title('Möbius Strip')
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        plt.tight_layout()
        plt.show()

# Example usage
if __name__ == "__main__":
    mobius = MobiusStrip(R=1, w=0.3, n=300)
    mobius.plot()
    print(f"Surface Area ≈ {mobius.compute_surface_area():.4f}")
    print(f"Edge Length ≈ {mobius.compute_edge_length():.4f}")
    
#Initializes an instance of MobiusStrip.
#Plots the 3D strip.
#Prints calculated surface area and edge length.

