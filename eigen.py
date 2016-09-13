from __future__ import print_function
import numpy as np
import pylab as pl

def load_K(filename):
    return np.loadtxt(filename, delimiter=',')

def eig_from_file(filename):
    K = load_K(filename)
    print(np.shape(K))
    eigenvalues, eigenvectors = np.linalg.eig(K)
    print(eigenvalues)
    return eigenvalues, eigenvectors

def plot_timescales_eigenvalues(eigenvalues, color='k', eigenvectors=None):
    """plot complex eigenvalues as scatter."""
    keepers = []
    for i in range(len(eigenvalues)):
        if abs(eigenvalues[i].imag) > 0.:
            keepers.append(i)
    keepers = np.array(keepers)
    reals = np.array([-1/eig.real for eig in eigenvalues])
    imags = np.array([abs(-1/(eig.imag)) for eig in eigenvalues])
    #print(reals)
    #print(imags)
    if len(keepers)>0:
        pl.scatter(reals[keepers],imags[keepers], color=color)
        pl.xlabel('Real component timescale (years)',fontsize=20)
        pl.ylabel('Imaginary component timescale (years)', fontsize=20)
        pl.loglog()

        return reals[keepers],eigenvectors

    else:
        return None, None

#pl.figure()
#for filename,color in zip(['B_default.txt','B_median.txt','B_low.txt'],
#                          ['k','r','b']):
#    eigenvalues,x = eig_from_file(filename)
#    plot_timescales_eigenvalues(eigenvalues,color=color)


K1 = load_K('B_default.txt')
K2 = load_K('B_low.txt')
K3 = load_K('B_median.txt')

diff_matrix = K3-K2
dK = np.linspace(0,1,501)
eigs = []
ds = []
vecs = []

for change in dK:
    eigenvals, eigenvecs = np.linalg.eig(K3-change*diff_matrix)
    eigs += list(eigenvals)
    ds += [change]*len(eigenvals)
    for x in eigenvecs:
        vecs.append(x)

pl.figure()

modes = {1:[],2:[],3:[],4:[]}
colormodes = {1:[],2:[],3:[],4:[]}
def get_mode(realeig):
    if realeig < 10:
        return 1
    elif realeig < 50:
        return 2
    elif realeig < 100:
        return 3
    elif realeig < 1000:
        return 4

for eig,change,vec in zip(eigs,ds,vecs):
    realeig, eigvec = plot_timescales_eigenvalues([eig],color=(change,0,0),
                                                  eigenvectors=vec)
    if not eigvec==None:
        #print(realeig[0])
        modes[get_mode(realeig[0])].append(eigvec)
        colormodes[get_mode(realeig[0])].append(change)

print([np.shape(modes[i]) for i in range(1,5)])
print(modes[1][0])



#H, xedges, yedges = np.histogram2d(x, y, bins=(128,128))
#H.shape, xedges.shape, yedges.shape

#extent = [yedges[0], yedges[-1], xedges[-1], xedges[0]]

#plt.imshow(H, extent=extent, interpolation='nearest')

#plt.colorbar()



pl.figure()
for mode in modes:
    H1 = np.zeros((18,31))
    pl.subplot(2,2,mode)
    m = modes[mode]
    c = colormodes[mode]
    for vec,col in zip(m,c):
        vec = [v.real for v in vec]
        #print(vec)
        #pl.scatter(range(1,19),vec, alpha=0.02, color=(col,0,0))
        H, xs, ys = np.histogram2d(range(1,19),vec,bins=(18,31),range=((0.5,18.5),(-1,1)))
        H1 += H
    pl.pcolor(xs,ys,H1.T)
    pl.colorbar()

pl.figure()
for mode in modes:
    H1 = np.zeros((18,31))
    pl.subplot(2,2,mode)
    m = modes[mode]
    c = colormodes[mode]
    for vec,col in zip(m,c):
        vec = [v.imag for v in vec]
        #print(vec)
        #pl.scatter(range(1,19),vec, alpha=0.02, color=(col,0,0))
        H, xs, ys = np.histogram2d(range(1,19),vec,bins=(18,31),range=((0.5,18.5),(-0.5,0.5)))
        H1 += H
    pl.pcolor(xs,ys,H1.T)
    pl.colorbar()


#average eigenvector (abs) for each mode in eigenvalues

pl.show()
