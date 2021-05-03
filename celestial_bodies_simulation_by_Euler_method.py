import os
import sys
import re
import math

class CelestialBodiesSimulation():

    def __init__(self, input_file):

        self.G = 0.100

        self.input_file = input_file
        self.output_file = os.path.splitext(input_file)[0]+'.out'
        self.f = open(self.output_file, 'w')

    def read_input_file(self, input_file):

        f = open(input_file, 'r')
        init_conditions = f.readlines()

        body_types = []
        ms = []
        xs = []
        vs = []

        for init_i in init_conditions:
            
            init_i = re.split(' [ ]*', init_i)

            body_types.append(init_i[0])

            ms.append(float(init_i[1]))

            xs_i = []
            for x in init_i[2:5]:
                xs_i.append(float(x))
            xs.append(xs_i)

            vs_i = []
            for v in init_i[5:8]:
                vs_i.append(float(v))
            vs.append(vs_i)

        f.close()

        print(body_types)
        print(ms)
        print(xs)
        print(vs)
        
        self.nb_bodies = len(ms)
        
        return body_types, ms, xs, vs

    def calc_fs(self, ms, xs):

        fs = [[0.0,0.0,0.0] for i in range(self.nb_bodies)]

        for i in range(self.nb_bodies):
            for j in range(i+1, self.nb_bodies):

                r3 = math.sqrt((xs[i][0]-xs[j][0])**2.0 + \
                               (xs[i][1]-xs[j][1])**2.0 + \
                               (xs[i][2]-xs[j][2])**2.0)

                r3 = r3 ** 3.0

                force_ijx = self.G * (ms[i]*ms[j])*(xs[i][0]-xs[j][0]) / r3
                force_ijy = self.G * (ms[i]*ms[j])*(xs[i][1]-xs[j][1]) / r3
                force_ijz = self.G * (ms[i]*ms[j])*(xs[i][2]-xs[j][2]) / r3

                fs[i][0] -= force_ijx
                fs[i][1] -= force_ijy
                fs[i][2] -= force_ijz

                fs[j][0] += force_ijx
                fs[j][1] += force_ijy
                fs[j][2] += force_ijz
                
        return fs

    def calc_vs(self, ms, fs, vs, deltaT):

        for i in range(self.nb_bodies):

            a_ix = fs[i][0] / ms[i]
            a_iy = fs[i][1] / ms[i]
            a_iz = fs[i][2] / ms[i]

            vs[i][0] += a_ix * deltaT
            vs[i][1] += a_iy * deltaT
            vs[i][2] += a_iz * deltaT

        return vs

    def calc_xs(self, vs, xs, deltaT):

        for i in range(self.nb_bodies):

            xs[i][0] += vs[i][0] * deltaT
            xs[i][1] += vs[i][1] * deltaT
            xs[i][2] += vs[i][2] * deltaT

        return xs

    def time_evolution(self, total_time, deltaT, output_interval):
    
        body_types, ms, xs, vs = self.read_input_file(self.input_file)

        nb_max_iters = math.ceil(total_time / deltaT)
        output_interval = int(output_interval / deltaT)

        for iters in range(nb_max_iters):
    
            fs = self.calc_fs(ms, xs)
            xs = self.calc_xs(vs, xs, deltaT)
            vs = self.calc_vs(ms, fs, vs, deltaT)
    
            if iters % output_interval == 0:
                self.write_trajectories_to_output_file(xs, body_types)
    
    
        self.f.close()
        
    def write_trajectories_to_output_file(self, xs, body_types):

       self.f.write(str(self.nb_bodies) + '\n')
       self.f.write('\n')
    
       for i in range(self.nb_bodies):
           traj = body_types[i] + ' ' + str(round(xs[i][0],4)) \
                                + ' ' + str(round(xs[i][1],4)) \
                                + ' ' + str(round(xs[i][2],4)) + '\n' 
           self.f.write(traj)

if __name__ == '__main__':

    input_file = sys.argv[1]
    cbs = CelestialBodiesSimulation(input_file)
    cbs.time_evolution(total_time=100, deltaT=1.00e-4, output_interval=1.0)

