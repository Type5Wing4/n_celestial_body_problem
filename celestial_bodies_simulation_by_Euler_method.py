import os
import sys
import re
import math

def time_evolution():

    G = 0.10

    input_file = sys.argv[1]

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
    
    nb_bodies = len(ms)

    deltaT = 1.0e-4
    total_time = 100.0
    nb_max_iters = math.ceil(total_time / deltaT)

    output_interval = int(nb_max_iters / 4000)

    output_file = os.path.splitext(input_file)[0]+'.out'
    f = open(output_file, 'w')


    for iters in range(nb_max_iters):

        fs = [[0.0,0.0,0.0] for i in range(nb_bodies)]

        for i in range(nb_bodies):
            for j in range(i+1, nb_bodies):

                r3 = math.sqrt((xs[i][0]-xs[j][0])**2.0 + \
                               (xs[i][1]-xs[j][1])**2.0 + \
                               (xs[i][2]-xs[j][2])**2.0)

                r3 = r3 ** 3.0

                force_ijx = G * (ms[i]*ms[j])*(xs[i][0]-xs[j][0]) / r3
                force_ijy = G * (ms[i]*ms[j])*(xs[i][1]-xs[j][1]) / r3
                force_ijz = G * (ms[i]*ms[j])*(xs[i][2]-xs[j][2]) / r3

                fs[i][0] -= force_ijx
                fs[i][1] -= force_ijy
                fs[i][2] -= force_ijz

                fs[j][0] += force_ijx
                fs[j][1] += force_ijy
                fs[j][2] += force_ijz

        for i in range(nb_bodies):

            a_ix = fs[i][0] / ms[i]
            a_iy = fs[i][1] / ms[i]
            a_iz = fs[i][2] / ms[i]

            vs[i][0] += a_ix * deltaT
            vs[i][1] += a_iy * deltaT
            vs[i][2] += a_iz * deltaT

            xs[i][0] += vs[i][0] * deltaT
            xs[i][1] += vs[i][1] * deltaT
            xs[i][2] += vs[i][2] * deltaT

        if iters % output_interval == 0:

            f.write(str(nb_bodies) + '\n')
            f.write('\n')

            for i in range(nb_bodies):
                traj = body_types[i] + ' ' + str(round(xs[i][0],4)) \
                                     + ' ' + str(round(xs[i][1],4)) \
                                     + ' ' + str(round(xs[i][2],4)) + '\n' 
                f.write(traj)

    f.close()

if __name__ == '__main__':

    time_evolution()
