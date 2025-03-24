import scipy.constants
import numpy as np
import matplotlib.pyplot as plt
import json
from matplotlib.animation import FuncAnimation



G = 1.18638e-4

class Planet:
    '''Class for defining a planet in the system'''
    
    def __init__(self,name,colour,mass,iniVel,iniPos):
        self.name = name
        self.colour = colour
        self.mass = mass
        self.iniVel = np.array(iniVel, dtype=np.float64)
        self.iniPos = np.array(iniPos, dtype=np.float64)
    
    def plColour(self):
        return self.colour
    
    def plName(self):
         return self.name

def initialize(planets,timestamps):
        num_planets = len(planets)

        pos = np.zeros((timestamps, num_planets,2))

        acc = np.zeros((timestamps, num_planets,2))

        vel = np.zeros((timestamps, num_planets,2))

        

        for i, planet in enumerate(planets):

                pos[0,i] = planet.iniPos

                vel[0,i] = planet.iniVel
    

        for j in range(num_planets):
            sum_acc = np.zeros(2)

            for k in range(num_planets):
                    if k != j:
                        r = pos[0, j] - pos[0, k]
                        sum_acc += (planets[k].mass / np.linalg.norm(r)**3) * r
            acc[0, j] = -G * sum_acc

        return num_planets, pos, acc, vel


def beeman_positions(planets,timestamps,dt):

    num_planets, pos, acc, vel = initialize(planets,timestamps)
                    
                    
    for i in range(timestamps-1):
        
        pos[i+1] = pos[i] + vel[i] * dt + 0.5 * acc[i] * dt**2
            
            
        for j in range(num_planets):
            sum_acc = np.zeros(2)

            for k in range(num_planets):

                if k != j:
                    r = pos[i+1, j] - pos[i+1, k]
                    sum_acc += (planets[k].mass / np.linalg.norm(r)**3) * r

            acc[i+1, j] = -G * sum_acc
            
            
        vel[i+1] = vel[i] + 0.5 * (acc[i] + acc[i+1]) * dt

    return pos, vel

def eulercromer_positions(planets,timestamps,dt):

    num_planets, pos, acc, vel = initialize(planets,timestamps)
                    
                    
    for i in range(timestamps-1):
        
        vel[i+1] = vel[i] + acc[i] * dt
            
        pos[i+1] = pos[i] + vel[i+1] * dt

        for j in range(num_planets):
            sum_acc = np.zeros(2)

            for k in range(num_planets):

                if k != j:
                    r = pos[i+1, j] - pos[i+1, k]
                    sum_acc += (planets[k].mass / np.linalg.norm(r)**3) * r

            acc[i+1, j] = -G * sum_acc

    return pos, vel           
        
def euler_positions(planets,timestamps,dt):
    num_planets, pos, acc, vel = initialize(planets,timestamps)
                    
                    
    for i in range(timestamps-1):
            
        pos[i+1] = pos[i] + vel[i] * dt

        vel[i+1] = vel[i] + acc[i] * dt
            
            

        for j in range(num_planets):
            sum_acc = np.zeros(2)

            for k in range(num_planets):

                if k != j:
                    r = pos[i+1, j] - pos[i+1, k]
                    sum_acc += (planets[k].mass / np.linalg.norm(r)**3) * r

            acc[i+1, j] = -G * sum_acc
                
    return pos, vel      

class Animation:

    def __init__(self,method,planets,timestamps,dt):

        self.planets = planets
        self.timestamps = timestamps
        self.dt = dt
        self.fig, self.ax = plt.subplots()
        self.method = str(method)

        self.pos, _ = method(planets,timestamps,dt)
        _, self.vel = method(planets,timestamps,dt)
        
        #x_min, x_max = np.min(self.pos[:, :, 0]),  np.max(self.pos[:, :, 0])
        y_min, y_max = np.min(self.pos[:, :, 1]),  np.max(self.pos[:, :, 1])
        self.ax.set_xlim(-30.3,30.3)
        self.ax.set_ylim(y_min+3,y_max+3)
        self.ax.set_xlabel("x (Mm)")
        self.ax.set_ylabel("y (Mm)")
        self.ax.set_aspect("equal")
        self.ax.set_facecolor('black')

       
        # Create planet patches for animation
        self.planet_patches = [plt.Circle((self.pos[0,i,0], self.pos[0,i,1]), 0.1, 
                                        color=planet.plColour(), animated=True) 
                             for i, planet in enumerate(planets)]
        for patch in self.planet_patches:
            self.ax.add_patch(patch)
   
    def animate(self, i):
        """Update planet positions for each frame"""
        for n in range(len(self.planets)):

            self.planet_patches[n].center = (self.pos[i,n,0], self.pos[i,n,1])

        return self.planet_patches

    def run(self):
        """Run the animation"""
        self.anim = FuncAnimation(self.fig, self.animate, 
                                frames=self.timestamps, 
                                interval=20, blit=True)
        plt.show()
                    
    def calculate_periods(self):  
        """ Calculate orbital periods for each planet."""
        orbital_periods = []
        for i in range(len(self.planets)):
    
             if np.array_equal(self.pos[0, i], np.array([0.0, 0.0])):
                orbitedPos = self.pos[:,i] 
                
             
        for i in range(len(self.planets)):
             t = 0
             if not(np.array_equal(self.pos[:,i] , orbitedPos)):

                while (self.pos[t,i,1] - orbitedPos[t,1] < 0 and self.pos[t+1,i,1] - orbitedPos[t+1,1] > 0) == False:

                        t+=1

             orbital_periods.append(t*self.dt)
             print(self.planets[i].name, "has period ", orbital_periods[i])        


    def total_energy(self,pos,vel):
        total_energy = []
        


        num_planets = len(self.planets)  

        for i in range(self.timestamps-1):
            sum_pot = 0
            sum_kin = 0

            for j in range(num_planets):

                sum_kin = sum_kin + 1/2 * self.planets[j].mass * np.linalg.norm(vel[i,j]) ** 2

                for k in range(num_planets):

                    if k != j:
                        r = pos[i+1, j] - pos[i+1, k]
                        sum_pot = sum_pot + (G* self.planets[k].mass * self.planets[j].mass/ np.linalg.norm(r))
                        

            total_energy += [-1/2 * sum_pot + sum_kin]

        return total_energy
    
    def graph_energy(self):
        b_pos, b_vel = beeman_positions(self.planets,self.timestamps,self.dt)
        ec_pos, ec_vel = eulercromer_positions(self.planets,self.timestamps,self.dt)
        e_pos, e_vel = euler_positions(self.planets,self.timestamps,self.dt)

        b_tot = self.total_energy(b_pos,b_vel)
        ec_tot = self.total_energy(ec_pos,ec_vel)
        e_tot = self.total_energy(e_pos,e_vel)


        stamps = np.linspace(0,self.timestamps-1,self.timestamps-1)

        plt.plot(stamps, b_tot, label = 'Beeman Method', color = 'blue', linestyle='-')
        plt.plot(stamps, ec_tot, label = 'Euler-Cromer Method', color = 'red', linestyle='--')
        plt.plot(stamps, e_tot, label = 'Euler Method', color = 'green', linestyle=':')
        plt.title('Total Energy for Different Methods')
        plt.xlabel('Timestamp')
        plt.ylabel('Total Energy')
        plt.legend()
        plt.show()

            


def main():
    M = 332946.0

    with open('parameters_solar.json') as f:
        parameters_solar = json.load(f)



        #name,colour,mass,iniVel,iniPos
        planets  = [Planet(body['name'],body['colour'],body['mass'],[0,np.sqrt(G*M/body['orbital_radius'])] if body['orbital_radius'] != 0 else [0, 0],[body['orbital_radius'], 0]) for body in parameters_solar['bodies']]
        

    
    
    
    #self,method,planets,timestamps,dt
    a = Animation(euler_positions,planets,parameters_solar['num_iterations'],parameters_solar['timestep'])
    a.run()
    a.graph_energy()
    
    
main()