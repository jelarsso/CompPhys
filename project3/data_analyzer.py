import numpy as np
import matplotlib.pyplot as plt
import subprocess as sb
import time
from IPython import embed

plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble'] = [r'\usepackage{lmodern}']

def read_output(filename):
    """
    Input:
    filename (str): the filename to read the data from. In the format specified by the headers from SolarSystem::write_to_file in the c++ code.
    
        Reads the data dumped to (filename). In the format specified by the headers from SolarSystem::write_to_file in the c++ code.
    
    Returns:
    positions (array), an array of shape timesteps * number of bodies * dims, the positions of every body at every timestep
    velocities (array), an array of shape timesteps * number of bodies * dims, the velocities of every body at every timestep
    timesteps (int), the number of timesteps read
    dt (float), the timestep-length used.
    nbodies (int), the number of bodies.
    dims (int), the number of dimensions.
    
    """
    with open(filename, "r") as infile:
        timesteps, dt, nbodies, dims = infile.readline().split(" ")[2::2]
        timesteps = int(timesteps)
        dt = float(dt)
        nbodies = int(nbodies)
        dims = int(dims)
    
    data = np.loadtxt(filename)
    data = data.reshape(timesteps,nbodies,2*dims)
    positions = np.zeros((timesteps,nbodies,dims))
    velocities = np.zeros_like(positions)

    positions = data[:,:,::2]
    velocities = data[:,:,1::2]
    
    return positions,velocities,timesteps,dt,nbodies,dims


def t3c_different_dt():
    """
    This functions runs the earth_sun_system simulations according to 3c. 
    Simulates the earth-sun system for many different timestep-lengths.
    Generates the plots for postions, energies and simulation wall-time.
    """
    dts = [1,1e-1,1e-2,1e-3,1e-6,1e-5,1e-4]
    sim_length = 10

    verlet_timings = []
    euler_timings = []
    last_pos_verlet = []
    last_pos_euler = []

    for dt in dts:
        start_t = time.time_ns()
        sb.run(["./euler_earth_sun",str(sim_length),str(dt)])
        euler_timings.append((time.time_ns()-start_t)*1e-9)
        euler_pos,euler_vel,*_ = read_output("output.data")
        start_t = time.time_ns()
        sb.run(["./verlet_earth_sun",str(sim_length),str(dt)])
        verlet_timings.append((time.time_ns()-start_t)*1e-9)
        verlet_pos,verlet_vel,*_ = read_output("output.data")

        last_pos_euler.append(np.linalg.norm(euler_pos[0,:,:]) - np.linalg.norm(euler_pos[-1,:,:]))
        last_pos_verlet.append(np.linalg.norm(verlet_pos[0,:,:]) - np.linalg.norm(verlet_pos[-1,:,:]))
    
    
        plt.plot(euler_pos[:,0,0],euler_pos[:,0,1],label="Euler")
        plt.plot(verlet_pos[:,0,0],verlet_pos[:,0,1],label="Verlet")
        plt.title(f"Simulation of circular orbit, $T={sim_length}$ year, $dt={dt}$ year")
        plt.xlabel("position x [AU]")
        plt.ylabel("position y [AU]")
        plt.legend()
        ax = plt.gca()
        ax.set_aspect("equal")
        plt.savefig(f"figs/pos_3c_dt{int(np.log(dt))}.pdf")
        plt.clf()
    
    
    verlet_kin_eng = 1/2*np.linalg.norm(verlet_vel,axis=2)**2
    euler_kin_eng = 1/2*np.linalg.norm(euler_vel,axis=2)**2
    verlet_pot_eng = 4*np.pi**2/np.linalg.norm(verlet_pos,axis=2)
    euler_pot_eng = 4*np.pi**2/np.linalg.norm(euler_pos,axis=2)

    mean_kin = np.mean(verlet_kin_eng+euler_kin_eng)/2
    mean_pot = np.mean(verlet_pot_eng+euler_pot_eng)/2

    t = np.arange(0,sim_length,dt)[::1000]
    
    plt.plot(t,(mean_kin - verlet_kin_eng)[1::1000]/mean_kin, label="Verlet Kinetic Energy")
    plt.plot(t,(mean_pot - verlet_pot_eng)[1::1000]/mean_pot, label="Verlet Potential Energy")
    plt.legend()
    plt.title(f"Verlet: Kinetic and potential energy relative difference from mean dt={dt} year")
    plt.xlabel("time [year]")
    plt.ylabel(r"relative energy")#[$M_\otimes\frac{\mathrm{AU}^2}{\mathrm{year}^2}$]")
    plt.savefig(f"figs/verlet_energy_dt{np.log10(dt)}.pdf")
    plt.clf()

    plt.plot(t,(mean_kin - euler_kin_eng)[1::1000]/mean_kin, label="Euler Kinetic Energy")
    plt.plot(t,(mean_pot - euler_pot_eng)[1::1000]/mean_pot, label="Euler Potential Energy")
    plt.legend()
    plt.title(f"Euler: Kinetic and potential energy relative difference from mean dt={dt} year")
    plt.xlabel("time [year]")
    plt.ylabel(r"relative energy")# [$M_\otimes\frac{\mathrm{AU}^2}{\mathrm{year}^2}$]")
    plt.savefig(f"figs/euler_energy_dt{np.log10(dt)}.pdf")
    plt.clf()

    plt.plot(np.log10(dts),np.log10(verlet_timings),"rx",label="Verlet")
    plt.plot(np.log10(dts),np.log10(euler_timings),"bx",label="Euler")
    plt.legend()
    plt.xlabel("dt [log(year)]")
    plt.ylabel("time [log(s)]")
    plt.title("Wall-time duration of simulations")
    plt.savefig("figs/elapsed_time.pdf")
    plt.clf()


def t3d_kepler_second():
    """
    Simulates the Earth-Sun system with three different initial velocities. 
    Breaks the simulations into time-intervals, and sums the covered area between
    Earth and the Sun in those intervals. 

    Output: One plot of the covered area for each time interval
    """
    dt = 1e-5
    sim_length = 2
    area_timesteps = 1000


    start_vels = [2*np.pi,np.pi,2.5*np.pi]
    for start_vel in start_vels:
        sb.run(["./verlet_earth_sun_start_vel",str(sim_length),str(dt),str(start_vel)])
        pos,vel,timesteps,*d = read_output("output.data")

        triangles = 1/2*(pos[:-1,0,0]*pos[1:,0,1] - pos[:-1,0,1]*pos[1:,0,0])
        
        sums = np.zeros((timesteps//area_timesteps))

        for i in range(timesteps//area_timesteps):
            sums[i] = np.sum(triangles[i*area_timesteps:(i+1)*area_timesteps-1])

        plt.plot(sums)
    plt.title('Covered area for each time-interval')
    plt.xlabel('Interval number / n'); plt.ylabel('Covered area over last time-interval / AU$^2$')
    plt.legend(['Initial velocity = 2$\pi$AU/yr','Initial velocity = $\pi$AU/yr','Initial velocity = 2.5$\pi$AU/yr'])
    plt.ylim(0, 5e-2)
    plt.show()


def t3e_beta():
    dt = 1e-5
    sim_length = 10
    start_vel = 5
    beta = [2,2.1,2.5,2.8,3]

    for b in beta:
        sb.run(["./verlet_es_v0_beta",str(sim_length),str(dt),str(start_vel),str(b)])
        pos,vel,timesteps,dt, *d = read_output("output.data")
        """
        plt.plot(pos[::100,0,0],pos[::100,0,1],label="position")
        plt.gca().set_aspect("equal")
        plt.suptitle(f"Simulation of orbit")
        plt.title(f"$v_0 = {start_vel:.2f}$AU/year, $\\beta={b}$, $T={sim_length}$ year, $dt={dt}$ year")
        plt.xlabel("position x [AU]")
        plt.ylabel("position y [AU]")
        plt.legend()
        ax = plt.gca()
        ax.set_aspect("equal")
        plt.savefig(f"figs/pos_3e_beta_{b}_dt{int(np.log(dt))}_vel{start_vel:.2f}.pdf")
        #plt.show()
        plt.clf()
        """

        verlet_kin_eng = 1/2*np.linalg.norm(vel,axis=2)**2
        verlet_pot_eng = -4*np.pi**2/(np.linalg.norm(pos,axis=2)**(b-1) )*1/(b-1)

        mean_kin = np.mean(verlet_kin_eng)
        mean_pot = np.mean(verlet_pot_eng)

        t = np.arange(0,sim_length,dt)[::100]
        
        #plt.plot(t,(mean_kin - verlet_kin_eng)[1::100]/mean_kin, label="Verlet Kinetic Energy")
        #plt.plot(t,(mean_pot - verlet_pot_eng)[1::100]/mean_pot, label="Verlet Potential Energy")
        #plt.plot(t,((verlet_kin_eng+verlet_pot_eng)/np.mean(verlet_pot_eng+verlet_kin_eng))[1::100],label="Mechanical Energy")
        """
        plt.plot(t,verlet_kin_eng[1::100],label="Kinetic energy")
        plt.plot(t,verlet_pot_eng[1::100],label="Potential energy")
        plt.plot(t,verlet_kin_eng[1::100]+verlet_pot_eng[1::100],label="Mechanical energy")
        plt.legend()
        plt.suptitle(f"Mechanical, Kinetic and potential energy.")
        plt.title(f"$v_0 = {start_vel:.2f}$AU/year, $\\beta={b}$, $T={sim_length}$ year, $dt={dt}$ year")
        plt.xlabel("time [year]")
        plt.ylabel(r"energy [$M_\otimes\frac{\mathrm{AU}^2}{\mathrm{year}^2}$]")
        plt.savefig(f"figs/energy_3e_beta{b}dt{np.log10(dt)}_vel{start_vel:.2f}.pdf")
        #plt.show()
        plt.clf()"""
        
        area_timesteps=5000

        triangles = 1/2*(pos[:-1,0,0]*pos[1:,0,1] - pos[:-1,0,1]*pos[1:,0,0])
        
        sums = np.zeros((timesteps//area_timesteps))

        for i in range(timesteps//area_timesteps):
            sums[i] = np.sum(triangles[i*area_timesteps:(i+1)*area_timesteps-1])

        plt.plot(sums)
        plt.xlabel("Time interval number")
        plt.ylabel("Summed area $[AU^2]$")
        plt.suptitle("The sweeped area of the orbit.")
        plt.title(f"$v_0 = {start_vel:.2f}$AU/year, $\\beta={b}$, $T={sim_length}$ year, $dt={dt}$ year",y=1.03)
        plt.savefig(f"figs/area_3e_beta{b}dt{np.log10(dt)}_vel{start_vel:.2f}.pdf")
        plt.clf()
        print(b)
        #kepler 2nd law ..




def t3f_esacpe():
    """
    Sets initial conditions so that Earth escapes the gravitational potential. Continously 
    simulates Earth-Sun systems until the Earth no longer escapes. 

    Output: Prints calculated velocity, exact escape velocity and relative error.  
    """
    dt = 1e-2
    sim_length = 100
    v = [2*np.pi,2.5*np.pi,2.8*np.pi,3*np.pi]
    beta = 2

    v0 = 2.9*np.pi
    check = False
    h = 1e-2
    while check==False:
        sb.run(["./verlet_es_v0_beta",str(sim_length),str(dt),str(v0),str(beta)])
        pos,vel,timesteps,*d = read_output("output.data")
        v0 = v0-h
        vr = np.dot(vel[:, 0, :], pos[:, 0, :].T) / np.linalg.norm(pos[:, 0, :], axis=1)  
        check = np.any(vr[0]<0)
        
    print('Calculated escape-velocity: ', v0, 'AU/yr')
    print('Excact escape-velocity: ', np.sqrt(8*np.pi**2), 'AU/yr')
    print('Relative error', (np.sqrt(8*np.pi**2) - v0)/np.sqrt(8*np.pi**2))

def t3g_tbp():
    dt = 1e-8
    sim_length = 2
    jupiter_mass_factors = [1000]

    for jm in jupiter_mass_factors: 
        sb.run(["./tbp_jupiter",str(sim_length),str(dt),str(jm)])
        pos,vel,timesteps,*d = read_output("output.data")
        
        for i in range(2):
            plt.plot(pos[:,i,0],pos[:,i,1])
        plt.legend(["Earth orbit","Jupiter orbit"])
        plt.suptitle(f"Simulation of orbit of Jupiter and Earth")
        plt.title(f"Jupiter mass factor = {jm}, $T={sim_length}$ year, $dt={dt}$ year")
        plt.xlabel("position x [AU]")
        plt.ylabel("position y [AU]")
        ax = plt.gca()
        ax.set_aspect("equal")
        plt.savefig(f"figs/pos_3g_jm{jm}_dt{int(np.log(dt))}.pdf")
        plt.show()
        plt.clf()

        diff = pos[:,1,:] - pos[:,0,:]
        plt.plot(np.linspace(0,sim_length,timesteps),np.linalg.norm(diff,axis=-1))
        plt.suptitle("Distance between Earth and Jupiter")
        plt.title(f"Jupiter mass factor = {jm}, $T={sim_length}$ year, $dt={dt}$ year")
        plt.xlabel("Time [year]")
        plt.ylabel("Distance [AU]")
        plt.savefig(f"figs/distance_je_3g_jm{jm}_dt{int(np.log(dt))}.pdf")
        plt.clf()

        plt.plot(np.linspace(0,sim_length,timesteps),np.linalg.norm(pos[:,0,:],axis=1))
        plt.suptitle("Distance between Earth and Sun")
        plt.title(f"Jupiter mass factor = {jm}, $T={sim_length}$ year, $dt={dt}$ year")
        plt.xlabel("Time [year]")
        plt.ylabel("Distance [AU]")
        plt.savefig(f"figs/distance_js_3g_jm{jm}_dt{int(np.log(dt))}.pdf")
        plt.clf()

def t3i_mercury():
    dt = 1e-6
    sim_length = 400 #number of
    each_sim = 0.25
    x0 = 0.3075
    y0 = 0
    vx0 = 0
    vy0 = 12.44
    apoapsis = []
    time = []
    xy_coord = []

    for i in range(sim_length):
        print(i)
        out = sb.run(["./gr_merc",str(each_sim),str(dt),str(x0),str(y0),str(vx0),str(vy0)],capture_output=True)
        if out.returncode!=0:
            raise ValueError("Code crashed")
        pos,vel,timesteps,*d = read_output("output.data")
        t = np.linspace(0,each_sim,timesteps)
        r = np.linalg.norm(pos[:,0,:],axis=1)
        ap = np.argmin(r)
        apoapsis.append(r[ap])
        time.append(t[ap]+i*each_sim)
        print(t[ap] + i*each_sim)
        xy_coord.append(pos[ap,0,:])
        x0 = pos[-1,0,0]
        y0 = pos[-1,0,1]
        vx0 = vel[-1,0,0]
        vy0 = vel[-1,0,1]

    xy_coord = np.asarray(xy_coord)
    precession = 648000/np.pi*np.arctan2(xy_coord[:,1],xy_coord[:,0])
    plt.plot(time,precession)
    plt.title("Precession of Mercury's orbit with GR-correction")
    plt.xlabel("time [year]")
    plt.ylabel("precessesion of perihelion [arcseconds]") 

def t3h_finmod():
    """
    Simulates a solar system with the Sun, all planets and Pluto. 

    Output: 
    Plot of all orbits.
    Plot of distance from Earth to Sun over the time of the simulation. 
    Plot of distance from Earth to mass center over the time of the simulation.
    """
    dt = 5e-5
    sim_length = 60

    sb.run(["./finalmod",str(sim_length),str(dt)])
    pos,vel,*d = read_output("output.data")

    #r = 40
    for i in range(10):
        plt.plot(pos[:,i,0],pos[:,i,1])
    #plt.xlim(-r, r); plt.ylim(-r, r)

    dist = pos[:,0,:]-pos[:,3,:]
    dist = np.linalg.norm(dist, axis=1)
    plt.figure()
    plt.plot(dist)

    plt.figure()
    plt.gca().set_aspect("equal")
    plt.grid()
    plt.title('Final model of the solar system')
    plt.xlabel('Distance / AU'); plt.ylabel('Distance / AU')
    plt.legend(['Sun', 'Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune', 'Pluto'])
    plt.show()






if __name__=="__main__":
    #uncomment one or more:
    
    #t3c_different_dt()
    t3d_kepler_second()
    #t3e_beta()
    #t3f_esacpe()
    #t3g_tbp()
    #t3i_mercury()
    #t3h_finmod()
