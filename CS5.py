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

def positions(planets,timestamps,dt):
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

        return pos

           
        
    

class Animation:

    def __init__(self,planets,timestamps,dt):

        self.planets = planets
        self.timestamps = timestamps
        self.dt = dt
        self.fig, self.ax = plt.subplots()


        self.pos = positions(planets,timestamps,dt)
        
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
                    
    '''def calculate_periods(self):  
        """ Calculate orbital periods for each planet."""
        orbital_periods = []
        for i in range(len(self.planets)):
             t = 0
             # since we can't calucalte the orbit of the sun around the sun
             if (self.pos[0,i,0]) == float(0) and (self.pos[0,i,1] == float(0)):
                 t = 0
             else:
                 while (self.pos[t,i,1] < 0 and self.pos[t+1,i,1] > 0) == False:
                     t+=1
             orbital_periods.append(t*self.dt)
             print(planet.name[i], "has period ", orbital_periods[i])     '''      


def main():
    M = 332946.0

    with open('parameters_solar.json') as f:
        parameters_solar = json.load(f)



        #name,colour,mass,iniVel,iniPos
        planets  = [Planet(body['name'],body['colour'],body['mass'],[0,np.sqrt(G*M/body['orbital_radius'])] if body['orbital_radius'] != 0 else [0, 0],[body['orbital_radius'], 0]) for body in parameters_solar['bodies']]
        

    
    
    
    #self,planets,timestamps,dt
    a = Animation(planets,parameters_solar['num_iterations'],parameters_solar['timestep'])
    a.run()
    a.calculate_periods()
    
    
main()