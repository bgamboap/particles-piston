import numpy as np

# Uniformize the time discretization of the data
def uniformize(arr, tarray):
    # arr: data input
    # tarray: desired times
    size = len(arr)
    
    # Determine the time discretization and number of points
    dt = tarray[1] - tarray[0]
    NT = len(tarray)
    
    
    # for each time t, find what is next closest time in arr 
    first_in_bin = np.zeros(NT,dtype='int') - 1
    for i in range(size-1,-1,-1):
        t = arr[i,0]
        b = min(int(t/dt),NT-1)
        first_in_bin[b] = i
        
    for i in range(NT-1,-1,-1):
        if first_in_bin[i] == -1:
            first_in_bin[i] = first_in_bin[i+1]
            
    
    # Perform the uniformization
    uniformed = np.zeros([NT, len(arr[0])])
    uniformed[0] = arr[0]
    t = 0
    for i in range(1,NT):
        n = first_in_bin[i]
        prev_n = first_in_bin[i-1]
        # print(i,n, prev_n)
        
        xprev = arr[n-1,1]
        xnext = arr[n,1]
        
        tprev = arr[n-1,0]
        tnext = arr[n,0]
        
        t = i*dt
        x = xprev + (xnext - xprev)*(t-tprev)/(tnext - tprev) # Assume linear motion

        # Assume the other variables do not change during the time
        vL2 = arr[n-1,2]
        vR2 = arr[n-1,3]
        dP = np.sum([arr[m-1,4] for m in range(prev_n, n)])
        
        uniformed[i] = [t,x,vL2, vR2, dP]
        
        
    return uniformed

if __name__ == "__main__":


    # Test the uniformize function
    import matplotlib.pyplot as plt

    # The input data is of this form: [t, xm, vL2, vR2, dP], where
    # t is the time
    # xm is the piston position
    # vL2 is the average of the square of the left particles' velocities
    # vR2 is the average of the square of the right particles' velocities
    # dP is the momentum transfer

    data = np.array([
        [0.000000,0.500000,0.167626,0.196591,0.001000],
        [0.012884,0.499170,0.170723,0.196591,-0.001162],
        [0.022077,0.497679,0.176708,0.196591,-0.000770],
        [0.029314,0.496019,0.176708,0.188299,-0.000876],
        [0.033055,0.494899,0.183650,0.188299,-0.000263],
        [0.042134,0.491699,0.198621,0.188299,-0.000715],
        [0.046211,0.490072,0.212542,0.188299,-0.000238],
        [0.059652,0.483760,0.238946,0.188299,-0.001092],
        [0.067975,0.479155,0.238946,0.188299,-0.000832],
        [0.072235,0.476530,0.238946,0.188299,-0.000426],
        [0.083044,0.469056,0.238946,0.188299,-0.001081],
        [0.087670,0.465500,0.263997,0.188299,-0.000303],
        [0.089980,0.463681,0.292059,0.188299,-0.000053],
        [0.091495,0.462486,0.292059,0.188299,-0.000152],
        [0.094558,0.460000,0.307121,0.188299,-0.000215],
        [0.100064,0.455346,0.337140,0.188299,-0.000377],
        [0.102594,0.453149,0.389598,0.188299,0.000050],
        [0.107685,0.448690,0.437821,0.188299,-0.000238],
        [0.109659,0.446944,0.470943,0.188299,-0.000010],
        [0.111380,0.445424,0.470943,0.188299,-0.000172]])



    Tmax = data[-1,0]

    fig,axs = plt.subplots(2,2)
    fig.set_figheight(7)
    fig.set_figwidth(13)


    NTd    = 32
    lastT  = 0.08
    tarray = np.linspace(0,lastT,NTd)
    unif   = uniformize(data, tarray)

    # Plot
    axs[0,0].plot(data[:,0], data[:,1],'o')
    for t in tarray:
        axs[0,0].axvline(t,color='k')
    axs[0,0].plot(unif[:,0], unif[:,1],'o')


    NTd    = 4
    lastT  = 0.08
    tarray = np.linspace(0,lastT,NTd)
    unif   = uniformize(data, tarray)

    # Plot
    axs[0,1].plot(data[:,0], data[:,1],'o')
    for t in tarray:
        axs[0,1].axvline(t,color='k')
    axs[0,1].plot(unif[:,0], unif[:,1],'o')



    NTd    = 6
    lastT  = 0.08
    tarray = np.linspace(0,lastT,NTd)
    unif   = uniformize(data, tarray)

    # Plot
    axs[1,0].plot(data[:,0], data[:,1],'o')
    for t in tarray:
        axs[1,0].axvline(t,color='k')
    axs[1,0].plot(unif[:,0], unif[:,1],'o')




    NTd    = 9
    lastT  = 0.1
    tarray = np.linspace(0,lastT,NTd)
    unif   = uniformize(data, tarray)

    # Plot
    axs[1,1].plot(data[:,0], data[:,1],'o')
    for t in tarray:
        axs[1,1].axvline(t,color='k')
    axs[1,1].plot(unif[:,0], unif[:,1],'o')
    plt.show()


    # Pressure
    NTd    = 4
    lastT  = 0.1
    tarray = np.linspace(0,lastT,NTd)
    print("TEST")
    unif   = uniformize(data, tarray)
    print(unif[:,4])


