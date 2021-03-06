import fast3tree
import numpy as np
from mpi4py import MPI
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
N=20
def calNN(data,boxsize):
    if rank ==0:
        print "calculate NN"
    min = boxsize/306/2.
    max = boxsize/2.-1.
    npoint = len(data)
    rho = npoint/boxsize**3
    count = np.zeros((N+1,),dtype=np.float64)
    xi = np.zeros((N+1,),dtype=np.float64)
    r = np.zeros(N+1)
    dx = (np.log10(max)-np.log10(min))/(N)

    start_n = rank*npoint/size
    stop_n = np.array([(rank+1)*npoint/size-1,npoint-1]).min()
    comm.Barrier()
    if(stop_n > start_n):
        with fast3tree.fast3tree(data) as tree:
            tree.set_boundaries(0.0,boxsize)
            tree.rebuild_boundaries()
            for j in range(start_n,stop_n+1):
                if rank == 0:
                    if j%((stop_n-start_n)/10) == 0:
                        print "process",j*100./(stop_n-start_n),"%"
                for i in range(N+1):
                    upper_r = 10.**(np.log10(min)+i*dx)
                    idx = tree.query_radius(data[j],upper_r, periodic=True,output='count') - 1
                    r[i] = upper_r
                    count[i] += idx
    
    comm.Barrier()
    # the 'totals' array will hold the sum of each 'data' array
    if comm.rank==0:
        print "Reducing data"
        # only processor 0 will actually get the data
        totals = np.zeros_like(count)
    else:
        totals = None

    # use MPI to get the totals 
    comm.Reduce(
        [count, MPI.DOUBLE],
        [totals, MPI.DOUBLE],
        op = MPI.SUM,
        root = 0
    )

    if rank == 0:
        xi[0] = totals[0]/(4./3*np.pi*r[0]**3)/rho/npoint 
        for i in range(1,N+1):
            dV = (4./3*np.pi*r[i]**3-4./3*np.pi*r[i-1]**3)
            xi[i] = (totals[i]-totals[i-1])/rho/dV/npoint 
        print xi
    
        #pylab.plot(r,xi)
        #pylab.xscale('log')
        #pylab.show()
    #xi = comm.bcast(xi, root=0)
    return (r,xi)
def cal_error(data,boxsize,nsub):
    sublength = boxsize/nsub
    xi0 = np.zeros(N+1,dtype=np.float64)
    xi2 =  np.zeros(N+1,dtype=np.float64)
    for i in range(nsub):
        for j in range(nsub):
            for k in range(nsub):
                if rank == 0:
                    print "calculate error",i,j,k 
                cond = np.where(~((data[:,0] > i*sublength) & (data[:,0] < (i+1)*sublength) \
                       & (data[:,1] > j*sublength) & (data[:,1] < (j+1)*sublength) \
                       & (data[:,2] > k*sublength) & (data[:,2] < (k+1)*sublength)))[0]
                ddata = data[cond]
                if rank == 0:
                    print cond
                    print ddata
                    (r,xi) = calNN(ddata,boxsize) 
                xi0 += xi
                xi2 += (xi)**2.
    delta = np.sqrt((nsub**3-1)*(xi2/nsub**3-(xi0/nsub**3)**2))
    return delta
# def main():
#     boxsize = 47.0
#     npoint = 10000
#     if rank == 0:
#         data = boxsize*np.random.rand(npoint, 3)
#     else:
#         data = None
#     data = comm.bcast(data, root=0)
#     (r,xi) = calNN(data,boxsize)
#     print xi
# if __name__ == "__main__":
#     main()    
