import numpy as np
import matplotlib.pyplot as plt

def discrete_circle( R, Nrad, join=False ) -> tuple:
    theta1 = np.linspace( 0, 2*np.pi, Nrad)
    X,Y = R * np.sin(theta1), R*np.cos(theta1)
    return X[:None if join else -1],Y[:None if join else -1]

def discrete_track(  L, R, N=None, d=None, join=False  ) -> tuple:
    if N is None and d is None:
       raise ValueError("Please specify either d or N")
    if N is not None and d is not None:
       raise ValueError("Please specify either d or N")
    if d is None:
       Ltotal = 2*np.pi*R+2*L
       d = Ltotal / float(N)
    Nrad = round( np.pi*R / d)
    Nline = round( L / d )
    theta1 = np.linspace(0, np.pi, Nrad)
    theta2 = np.linspace(np.pi, 2*np.pi, Nrad)
    X1,Y1 = R*np.sin(theta1) + L/2, R*np.cos( theta1 )
    X2,Y2 = np.linspace( L/2,-L/2,Nline), -R*np.ones(Nline)
    X2,Y2 = X2[1:-1], Y2[1:-1]#Lazy way to prevent double points
    X3,Y3 = R*np.sin(theta2) - L/2, R*np.cos( theta2 )
    X4,Y4 = np.linspace(-L/2,L/2,Nline), R*np.ones(Nline)
    X4,Y4 = X4[1:],Y4[1:]
    X, Y = np.r_[ X1, X2, X3, X4 ], np.r_[Y1, Y2, Y3, Y4 ]
    return X[:None if join else -1],Y[:None if join else -1]

if __name__ == "__main__":
    R = 0.2
    L = 0.6
    X,Y = discrete_track(L,  R, d=0.025 , join=False)
    plt.plot(X,Y,"-ok")
    plt.axis("equal")
    plt.show()
