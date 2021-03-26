# Garrett Giddings
# University of Alabama AAML
# 5-11-2020
# Antenna design using AIS BPSO algorithm


import math
import random
import ScriptEnv
import csv
import re
import os
import time

#Naming folder for this particular simulation
folder = '2021-03-26 Vivaldi Opt'
#Location where all results will be stored
output_loc = "C:/Users/jalnas.STUDENT/Documents/BPSO/Vivalid_outputs/"+folder
#Project name
project_name="2_6 GHz_Vivaldi_single"
#Design name
design="HFSSDesign1"
#Prefix of parasitic component models in HFSS
parasitic_prefix = "Rectangle"
#integer offset for parasitic names
parasitic_offset = 126
#label for HFSS boundary condition of parasitic elements
parasitic_boundary_name = "Finite_para"

#HFSS scripting object initialization
ScriptEnv.Initialize("Ansoft.ElectronicsDesktop") #repeated initializations are unnecessary/possibly harmful
oProject = oDesktop.SetActiveProject(project_name)
oDesign = oProject.SetActiveDesign(design)
oEditor = oDesign.SetActiveEditor("3D Modeler")
oModule = oDesign.GetModule("BoundarySetup")

class Particle:
    def __init__(self, position, velocity, cost, best_position, best_cost):
        self.position = position
        self.velocity = velocity
        self.cost = cost
        self.best_position = best_position
        self.best_cost = best_cost
        self.local_best_position = position
        self.local_best_cost = float('inf')

def main():
    # Creating output folder
    if not os.path.exists(output_loc):
        os.makedirs(output_loc)

    # Parameters
    max_it = 50  # Maximum num of iterations
    pop = 10  # Population size
    n_vars = 30  # Number of patches
    Vmax = 10  # AIS velocity maximum
    lbest = -1  # 0 -> Global best ||| 1 -> Local Best ||| (0,1) -> Hybrid ||| -1 -> Dynamic Hybrid
    last_i = 0  # The last iteration performed by the previous optimization

    #Storing Parameters
    with open(output_loc+"/Optimization_Params.txt", "w") as f:
        str_out = optim_param_header(max_it, pop, n_vars, Vmax, lbest)
        f.write(str_out)

    # Calling PSO
    PSO(max_it ,pop ,n_vars ,Vmax, lbest, last_i)

    oProject.Save()
    return

def particle_to_str(p):
    out = "\nPosition: " + str(p.position) + "\nVelocity: " + str(p.velocity) + "\nCost:" + str(
        p.cost) + "\nBest Position:" + str(p.best_position) + "\nBest Cost:" + str(p.best_cost) + "\n"
    return out

def PSO(max_it, pop, n_vars, Vmax, lbest_c, last_i):
    global times_list
    Positions_dic = {} #dictionary has cost associated with each distinct particle position /JA
    max_pos = 2 ** n_vars - 1

    if last_i:
        global folder
        with open(output_loc + "/swarm_iteration_" + str(last_i) + ".txt", "r") as f:
            report = f.read()
            #TO-DO: Insert routine to read n_vars, pop, from file /JA
            load_costs = re.findall(r"\nCost:\s*(-?\d+.?\d*)", report)
            load_velocities = re.findall(r"\nVelocity:\s*(\d+)", report)
            load_positions = re.findall(r"\nPosition:\s*(\d+)", report)
            load_best_positions = re.findall(r"\nBest Position:\s*(\d+)", report)
            load_best_costs = re.findall(r"\nBest Cost:\s*(-?\d+.?\d*)", report)
            load_global_cost = re.findall(r"\nGlobal Best Cost:\s*(-?\d+.?\d*)", report)
            load_global_position = re.findall(r"\nGlobal Best Position:\s*(\d+)", report)
        global_best_cost = float(load_global_cost[0])
        global_best_position = int(load_global_position[0])

        # Initialize population members
        swarm = []
        for i in range(pop):
            start_position = int(load_positions[i])
            start_velocity = int(load_velocities[i])
            start_cost = float(load_costs[i])
            start_bp = int(load_best_positions[i])
            start_bc = float(load_best_costs[i])
            swarm.append(Particle(start_position, start_velocity, start_cost, start_bp, start_bc))
    else:
        # Intitialization
        # Initialize to worst possible solution (infinity for minimization problems)
        global_best_cost = float('inf')
        global_best_position = -1

        # Initialize population members
        swarm = []
        hfss_message("Starting BPSO iteration 0",0)
        for i in range(pop):
            # Generate random starting locs and velocities
            start_velocity_b = []
            start_velocity_s = ''
            ones = [] #contains indices of 1 bits

            for n in range(n_vars):
                start_velocity_b.append(random.randint(0, 1))
                if start_velocity_b[n] == 1:
                    ones.append(n)
            while len(ones) > Vmax: #restrict 1 bits in velocity to "Vmax" 1 bits /JA
                n = random.randint(0, len(ones) - 1)    #randomly select a 1 bit /JA
                start_velocity_b[ones[n]] = 0           #remove from bit array /JA
                ones.pop(n)                             #remove index entry /JA

            for bit in start_velocity_b: #convert bit array to string
                start_velocity_s += str(bit)

            start_position = random.randint(0, max_pos)
            start_velocity = int(start_velocity_s, 2)
            swarm.append(Particle(start_position, start_velocity, float('inf'), 0, float('inf')))

            # Evaluation
            swarm[i].cost = simulate(swarm[i].position, n_vars)
            Positions_dic[swarm[i].position] = swarm[i].cost

            # Update personal best
            if swarm[i].cost < swarm[i].best_cost:
                swarm[i].best_position = swarm[i].position
                swarm[i].best_cost = swarm[i].cost
            # Update global best
            if swarm[i].best_cost < global_best_cost:
                global_best_cost = swarm[i].best_cost
                global_best_position = swarm[i].best_position
    # Update local best
    for i in range(pop):
        for n in [-1,0,1]:
            ind = (i+n) % pop #note that (-1)%x = x-1 /JA
            if swarm[ind].best_cost < swarm[i].local_best_cost:
                swarm[i].local_best_cost = swarm[ind].best_cost
                swarm[i].local_best_position = swarm[ind].best_position

    #Store swarm Data
    if last_i==0:
        with open(output_loc+"/swarm_iteration_" + str(0) + ".txt", "w") as swarm_data:
            str_out = optim_param_header(max_it, pop, n_vars, Vmax, lbest_c)
            str_out += "\n\n"
            str_out += swarm_data_string(swarm, 0, lbest_c, global_best_cost, global_best_position, 0)
            swarm_data.write(str_out)

    # Main loop of pso
    temp = lbest_c #stores lbest_c user parameter in temporary variable /JA
    for j in range(last_i, max_it):
        #Notify new iteration /JA
        message = "Starting BPSO iteration " + str(j+1)
        hfss_message(message,0)

        iter_start = time.time()  #iteration start time /JA

        convergence_test = 1
        best_costs = 0
        for i in range(pop):
            # Update Velocity
            if temp == -1:
                lbest_c = 1 - j / max_it #lbest_c inversely proportional to iteration number /JA

            # Combine gbest and lbest values
            mask_b = []
            for n in range(n_vars):
                if random.random() > lbest_c: #lbest_c is a probability threshold for a 1 bit /JA
                    mask_b.append(1)
                else:
                    mask_b.append(0)
            m_str = ''
            for bit in mask_b:
                m_str += str(bit)
            mask = int(m_str, 2)
            gbest = global_best_position & mask #apply mask to global position /JA
            mask = mask ^ (2 ** n_vars - 1) #invert mask /JA
            lbest = swarm[i].local_best_position & mask #apply inverted mask to local position /JA
            swarm_vector = gbest | lbest

            c1 = random.randint(0, max_pos) #random weighting /JA
            c2 = random.randint(0, max_pos) #random weighting /JA
            d1 = swarm[i].position ^ swarm[i].best_position
            d2 = swarm[i].position ^ swarm_vector
            velocity = (d1 & c1) | (d2 & c2)

            # Apply velocity limits
            v = [int(x) for x in bin(velocity)[2:]]
            ones = []
            for n, bit in enumerate(v):
                if bit == 1:
                    ones.append(n)

            while len(ones) > Vmax: #limit number of 1 bits to Vmax /JA
                n = random.randint(0, len(ones) - 1)
                v[ones[n]] = 0
                ones.pop(n)

            v_str = ''
            for bit in v:
                v_str += str(bit)
            v_new = int(v_str, 2)

            if v_new == 0: #if the new velocity is 0, set random bit to 1 /JA
                v_new = 2 ** (random.randint(0, n_vars - 1))

            swarm[i].velocity = v_new

            # Update Position
            swarm[i].position = swarm[i].position ^ swarm[i].velocity

            # Update Cost
            if swarm[i].position in Positions_dic: #if position has already been simulated, get results /JA
                swarm[i].cost = Positions_dic[swarm[i].position]

            # if results exists, read and analyze /JA
            elif os.path.exists(output_loc+"/Parasitic_S11_" + str(swarm[i].position) + ".csv"):
                swarm[i].cost = analyze(swarm[i].position)
                Positions_dic[swarm[i].position] = swarm[i].cost

            else:
                swarm[i].cost = simulate(swarm[i].position, n_vars)
                Positions_dic[swarm[i].position] = swarm[i].cost

            # Update Personal Best
            if swarm[i].cost < swarm[i].best_cost:
                swarm[i].best_position = swarm[i].position
                swarm[i].best_cost = swarm[i].cost
                # Update global best
                if swarm[i].best_cost < global_best_cost:
                    global_best_cost = swarm[i].best_cost
                    global_best_position = swarm[i].best_position

            # Convergence check
            if i == 0:
                best_costs = swarm[i].best_cost
            if swarm[i].best_cost != best_costs: #All particles must have encountered the best cost for convergence /JA
                convergence_test = 0

        # Update local best
        for i in range(pop):
            for n in [-1,0,1]:
                ind = (i+n) % pop
                if swarm[ind].best_cost < swarm[i].local_best_cost:
                    swarm[i].local_best_cost = swarm[ind].best_cost
                    swarm[i].local_best_position = swarm[ind].best_position

        iter_end = time.time() #iteration end time

        oProject.Save() #Save Project; hypothesis that saving clears RAM and speeds up execution speed /JA

        #Calculate and display iteration duration
        #TO-DO: log iteration duration to swarm data
        iter_duration = iter_end - iter_start #iteration duration in seconds /JA
        msg = "Iteration " + str(j+1) + " duration: " + str(iter_duration) + " seconds"
        hfss_message(msg, 0)

        #Store swarm Data
        with open(output_loc+"/swarm_iteration_" + str(j+1) + ".txt", "w") as swarm_data:
            str_out = optim_param_header(max_it, pop, n_vars, Vmax, temp)
            str_out += "\n\n"
            str_out += swarm_data_string(swarm, j+1, temp, global_best_cost, global_best_position, convergence_test)
            swarm_data.write(str_out)

        # Exit if converged
        convergence_test = 0 # execute all iterations
        if convergence_test == 1:
            hfss_message("Optimization converged", 0)
            break

    # After optimization load best cost design and simulate
    config_parasitics(particle=global_best_position, n_vars=n_vars)
    oDesign.Analyze("Setup1")
    return
def config_3D_parasitics(particle, n_vars):
    config_start = time.time()

    # pad argument to bit length of n_vars
    model_on = [int(x) for x in bin(particle)[2:]]
    while len(model_on) < n_vars:
        model_on.insert(0, 0)

    for j, bit in enumerate(model_on):
        oEditor.ChangeProperty(
            [
                "NAME:AllTabs",
                [
                    "NAME:Geometry3DAttributeTab",
                    [
                        "NAME:PropServers",
                        parasitic_prefix + str(j)
                    ],
                    [
                        "NAME:ChangedProps",
                        [
                            "NAME:Transparent",
                            "Value:="	, 1-bit
                        ],
                        [
                            "NAME:Model",
                            "Value:="		, bool(bit)
                        ]
                    ]
                ]
            ])
    return

def config_parasitics(particle, n_vars):
    # assumes parasitics are sheets assigned a copper boundary condition
    config_start = time.time()

    # pad argument to bit length of n_vars
    model_on = [int(x) for x in bin(particle)[2:]]
    while len(model_on) < n_vars:
        model_on.insert(0, 0)

    # check if desired boundary exists
    existing_boundaries = oModule.GetBoundaries()
    if parasitic_boundary_name in existing_boundaries:
        oModule.DeleteBoundaries([parasitic_boundary_name])

    objects = []
    for j, bit in enumerate(model_on):
        oEditor.ChangeProperty(
            [
                "NAME:AllTabs",
                [
                    "NAME:Geometry3DAttributeTab",
                    [
                        "NAME:PropServers",
                        parasitic_prefix + str(j + 1 + parasitic_offset)
                    ],
                    [
                        "NAME:ChangedProps",
                        [
                            "NAME:Model",
                            "Value:=", bool(bit)
                        ],
                        [
                            "NAME:Transparent",
                            "Value:=", 1 - bit
                        ]
                    ]
                ]
            ]
        )
        if bit == 1:
            objects.append(parasitic_prefix + str(j + 1 + parasitic_offset))

    #objects.append("IFA_arm")
    oModule.AssignFiniteCond(
        [
            "NAME:" + parasitic_boundary_name,
            "Objects:=", objects,
            "UseMaterial:="	, True,
            "Material:="		, "copper",
            "UseThickness:="	, False,
            "Roughness:="		, "0um",
            "InfGroundPlane:="	, False,
            "IsTwoSided:="		, False,
            "IsInternal:="		, True
        ])

    config_end = time.time()
    config_time = config_end - config_start
    hfss_message(str(config_time) + " seconds to configure",0)

def simulate(particle, n_vars):
    config_parasitics(particle, n_vars)

    sim_start = time.time()
    oDesign.Analyze("Setup1")
    sim_end = time.time()
    sim_duration  = sim_end - sim_start
    hfss_message("Particle Configuration "+str(particle)+" simulation duration: "+str(sim_duration), 0)

    oReportModule = oDesign.GetModule("ReportSetup")
    oReportModule.ExportToFile('Output Variables Table 1', output_loc +"/Parasitic_S11_" + str(particle) + ".csv")
    cost = analyze(particle)
    return cost

def analyze(particle):
    # Read simulation results file
    #global output_loc

    c1 = 1 #weight factor; minimize 2GHz reflection
    c2 = 1 #weight factor; minimize 3GHz reflection
    c3 = 1 #weight factor; minimize 4GHz reflection
    c4 = 1 #weight factor; minimize 6GHz reflection

    #Reads in the first data row which has the form:
    #[Frequency (GHz)],[S11 @ 2GHz],[S11 @ 3GHz],[S11 @ 4GHz],[S11 @ 6GHz]
    with open(output_loc + "/Parasitic_S11_" + str(particle) + ".csv") as csvFile:
        reader = csv.reader(csvFile)
        next(reader) # skip header
        S11_dB = next(reader) #save first row

    S11_2G_dB = float(S11_dB[1])
    S11_3G_dB = float(S11_dB[2])
    S11_4G_dB = float(S11_dB[3])
    S11_6G_dB = float(S11_dB[4])

    hfss_message("{}, {},{},{}".format(S11_2G_dB, S11_3G_dB, S11_4G_dB, S11_6G_dB), 0)
    cost = c1 * S11_2G_dB + c2 * S11_3G_dB + c3 * S11_4G_dB + c4 * S11_6G_dB
    '''
    bandwidth = float(bw[1])
    if bw[0] == 'bandwidth [GHz]':
        bandwidth *= 1000
    '''
    return cost

def optim_param_header(max_it, pop, n_vars, Vmax, lbest):
    str_out = "Optimization Parameters"
    str_out += "\n-------------------------"
    str_out += "\nmax_it: " + str(max_it)
    str_out += "\npop: " + str(pop)
    str_out += "\nn_vars: " + str(n_vars)
    str_out += "\nVmax: " + str(Vmax)
    str_out += "\nlbest: " + str(lbest)

    return str_out

def swarm_data_string(swarm, iteration, lbest, global_best_cost, global_best_position, convergence):
    str_out = "Swarm Data"
    str_out += "\nIteration: " + str(iteration)
    str_out += "\n------------"
    for n,p in enumerate(swarm):
        str_out += "\nParticle "+str(n)+":"
        str_out += particle_to_str(p)
        if lbest:
            str_out+= "Local Best Cost: " + str(p.local_best_cost)
        str_out += "\n"
    str_out += "\n\nGlobal Best Cost: " + str(global_best_cost)
    str_out += "\nGlobal Best Position: " + str(global_best_position)
    str_out += "\nConvergence: "
    if convergence:
        str_out += "True"
    else:
        str_out += "False"

    return str_out

def hfss_message(msg, priority):
    """
    Displays a string in HFSS; relies on project name and design name in global variables
    :param msg: string to display in HFSS Message Manager
    :param priority: integer in [0,3]; 0 is lowest priority
    :return: none
    """
    #displays a message in the HFSS message manager
    #priority = [0,3] where 0 is lowest priority
    oDesktop.AddMessage(project_name,design,priority,msg)
    return

if __name__ == "__main__": #allows functions to be defined at end-of-file
    main()
