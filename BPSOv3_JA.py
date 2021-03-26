#Garrett Giddings
#University of Alabama AAML
#5-11-2020
#Antenna design using AIS BPSO algorithm

import math
import random
import ScriptEnv
import csv
import re
import os

#Naming folder for this particular simulation
folder = 'IFA_5-12'
#Location where all results will be stored
output_loc = "C:/Users/grgiddings/Documents/Ansoft/ALL_Outputs/"+folder
#Project name
project_name="IFA_Garrett"
#Design name
design="IFA_para1"

class Particle:
    def __init__(self, position, velocity, cost, best_position, best_cost):
        self.position = position
        self.velocity = velocity
        self.cost = cost
        self.best_position = best_position
        self.best_cost = best_cost
        self.local_best_position = position
        self.local_best_cost = float('inf')

def particle_to_str(p):
	out = "\nPosition: " + str(p.position) + "\nVelocity: " + str(p.velocity) + "\nCost:" + str(p.cost) + "\nBest Position:" + str(p.best_position) + "\nBest Cost:" + str(p.best_cost) + "\n"
	return out

def PSO(max_it,pop,n_vars,Vmax,lbest_c, last_i):
    Positions_dic = {}
    max_pos = 2**n_vars - 1

    if last_i:
        global folder
        with open("C:/Users/grgiddings/Documents/Ansoft/ALL_OUTPUTS/"+folder+"/swarm_iteration_"+str(last_i)+".txt", "r") as f:
            report = f.read()
            #Insert routine to read n_vars, pop, from file
            load_costs = re.findall(r"\nCost:(-\d+.\d+)",report)
            load_velocities = re.findall(r"\nVelocity: (\d+)",report)
            load_positions = re.findall(r"\nPosition: (\d+)",report)
            load_best_positions = re.findall(r"\nBest Position:(\d+)",report)
            load_best_costs = re.findall(r"\nBest Cost:(-\d+.\d+)",report)
            load_global_cost = re.findall(r"\nGlobal Best Cost: (-\d+.\d+)",report)
            load_global_position = re.findall(r"\nGlobal Best Position: (\d+)",report)
        global_best_cost = float(load_global_cost[0])
        global_best_position = int(load_global_position[0])
        
        #Initialize population members
        swarm = []
        for i in range(pop):
            start_position = int(load_positions[i])
            start_velocity = int(load_velocities[i])
            start_cost = float(load_costs[i])
            start_bp = int(load_best_positions[i])
            start_bc = float(load_best_costs[i])
            swarm.append(Particle(start_position,start_velocity,start_cost,start_bp,start_bc))
    else:
        #Intitialization
        #Initialize to worst possible solution (infinity for minimization problems)
        global_best_cost = float('inf')
        global_best_position = -1

        #Initialize population members
        swarm = []
        for i in range(pop):
            #Generate random starting locs and velocities
            start_velocity_b = []
            start_velocity_s = ''
            ones = []
            for n in range(n_vars):
                start_velocity_b.append(random.randint(0,1))
                if start_velocity_b[n] == 1:
                    ones.append(n)
            while len(ones) > Vmax:
                n = random.randint(0,len(ones)-1)
                start_velocity_b[ones[n]] = 0
                ones.pop(n)
            for bit in start_velocity_b:
                start_velocity_s += str(bit)
            start_position = random.randint(0,max_pos)
            start_velocity = int(start_velocity_s,2)
            swarm.append(Particle(start_position,start_velocity, float('inf'),0,float('inf')))

            #Evaluation
            swarm[i].cost = simulate(swarm[i].position, n_vars)
            Positions_dic[swarm[i].position] = swarm[i].cost
            #Update personal best
            if swarm[i].cost < swarm[i].best_cost:
                swarm[i].best_position = swarm[i].position
                swarm[i].best_cost = swarm[i].cost
            #Update global best
            if swarm[i].best_cost < global_best_cost:
                global_best_cost = swarm[i].best_cost
                global_best_position = swarm[i].best_position

    #Update local best
    for i in range(pop):
        for n in [-1,0,1]:
            ind = (i+n) % pop
            if swarm[ind].best_cost < swarm[i].local_best_cost:
                swarm[i].local_best_cost = swarm[ind].best_cost
                swarm[i].local_best_position = swarm[ind].best_position

    #Store swarm Data
    if last_i==0:
        with open(output_loc+"/swarm_iteration_" + str(0) + ".txt", "w") as swarm_data:
            str_out = optParamHeader(max_it, pop, n_vars, Vmax, lbest)
            str_out += "\n"
            str_out += swarmDataString(swarm, global_best_cost, global_best_position)
            swarm_data.write(str_out)

    #Main loop of pso
    temp=lbest_c
    for j in range(last_i,max_it):
        convergence_test = 1
        best_costs = 0
        for i in range(pop):
            #Update Velocity
            if temp == -1:
                lbest_c = 1-j/max_it
            #Combine gbest and lbest values
            mask_b = []
            for n in range(n_vars):
                if random.random() > lbest_c:
                    mask_b.append(1)
                else:
                    mask_b.append(0)
            m_str = ''
            for bit in mask_b:
                m_str += str(bit)
            mask = int(m_str,2)
            gbest = global_best_position & mask
            mask = mask ^ (2**n_vars-1)
            lbest = swarm[i].local_best_position & mask
            swarm_vector = gbest | lbest

            c1 = random.randint(0,max_pos)
            c2 = random.randint(0,max_pos)
            d1 = swarm[i].position ^ swarm[i].best_position
            d2 = swarm[i].position ^ swarm_vector
            velocity = (d1 & c1) | (d2 & c2)

            #Apply velocity limits
            v = [int(x) for x in bin(velocity)[2:]]
            ones = []
            for n,bit in enumerate(v):
                if bit == 1:
                    ones.append(n)
            while len(ones) > Vmax:
                n = random.randint(0,len(ones)-1)
                v[ones[n]] = 0
                ones.pop(n)
            v_str = ''
            for bit in v:
                v_str += str(bit)
            v_new = int(v_str,2)
            if v_new == 0:
                v_new = 2**(random.randint(0,n_vars-1))
            swarm[i].velocity = v_new

            #Update Position
            swarm[i].position = swarm[i].position ^ swarm[i].velocity
                                    
            #Update Cost
            if swarm[i].position in Positions_dic:
                swarm[i].cost = Positions_dic[swarm[i].position]
            else:
                if os.path.exists(output_loc+"/Parasitic_S11_" + str(swarm[i].position) + ".csv"):
                    swarm[i].cost = analyze(swarm[i].position)
                    Positions_dic[swarm[i].position] = swarm[i].cost
                else:
                    swarm[i].cost = simulate(swarm[i].position, n_vars)
                    Positions_dic[swarm[i].position] = swarm[i].cost
            
            #Update Personal Best
            if swarm[i].cost < swarm[i].best_cost:
                swarm[i].best_position = swarm[i].position
                swarm[i].best_cost = swarm[i].cost 
                #Update global best
                if swarm[i].best_cost < global_best_cost:
                    global_best_cost = swarm[i].best_cost
                    global_best_position = swarm[i].best_position
            #Convergence check
            if i == 0:
                best_costs = swarm[i].best_cost
            if swarm[i].best_cost != best_costs:
                convergence_test = 0
        #Update local best
        for i in range(pop):
            for n in [-1,0,1]:
                ind = (i+n) % pop
                if swarm[ind].best_cost < swarm[i].local_best_cost:
                    swarm[i].local_best_cost = swarm[ind].best_cost
                    swarm[i].local_best_position = swarm[ind].best_position

        #Store swarm Data 
        with open(output_loc+"/swarm_iteration_" + str(j+1) + ".txt", "w") as swarm_data:
            str_out = optParamHeader(max_it, pop, n_vars, Vmax, lbest)
            str_out += "\n"
            str_out += swarmDataString(swarm, global_best_cost, global_best_position)
            swarm_data.write(str_out)

        #Exit if converged
        if convergence_test == 1:
            return
    return

def simulate(particle, n_vars):	
    global project_name
    global design
    ScriptEnv.Initialize("Ansoft.ElectronicsDesktop")
    oDesktop.RestoreWindow()
    oProject = oDesktop.SetActiveProject(project_name)
    oDesign = oProject.SetActiveDesign(design)
    oEditor = oDesign.SetActiveEditor("3D Modeler")
    oModule = oDesign.GetModule("BoundarySetup")

    model_on = [int(x) for x in bin(particle)[2:]]
    while len(model_on) < n_vars:
        model_on.insert(0,0)
    for i in range(len(model_on)):
        oEditor.ChangeProperty(
                    [
                        "NAME:AllTabs",
                        [
                            "NAME:Geometry3DAttributeTab",
                            [
                                "NAME:PropServers", 
                                "parasitic_Detach"+str(i+1)
                            ],
                            [
                                "NAME:ChangedProps",
                                [
                                    "NAME:Model",
                                    "Value:="		, False
                                ],
                                [
                                    "NAME:Transparent",
                                    "Value:="		, 1
                                ]
                            ]
                        ]
                    ]
                )

    for j,bit in enumerate(model_on):
        if bit == 1:
            oEditor.ChangeProperty(
                [
                    "NAME:AllTabs",
                    [
                        "NAME:Geometry3DAttributeTab",
                        [
                            "NAME:PropServers", 
                            "parasitic_Detach"+str(j+1)
                        ],
                        [
                            "NAME:ChangedProps",
                            [
                                "NAME:Model",
                                "Value:="		, True
                                ],
                            [
                                "NAME:Transparent",
                                "Value:="		, 0
                            ]
                        ]
                    ]
                ]
            )
            oModule.AssignFiniteCond(
            [
                "NAME:FiniteCond_para"+str(j+1),
                "Objects:="		, ["parasitic_Detach"+str(j+1)],
                "UseMaterial:="		, True,
                "Material:="		, "copper",
                "UseThickness:="	, False,
                "Roughness:="		, "0um",
                "InfGroundPlane:="	, False,
                "IsTwoSided:="		, False,
                "IsInternal:="		, True
            ])
    oDesign.AnalyzeAll()
    oModule = oDesign.GetModule("ReportSetup")
    global output_loc
    oModule.ExportToFile('S Parameter Plot 1', output_loc+"/Parasitic_S11_" + str(particle) + ".csv")
    cost = analyze(particle)
    return cost 

def analyze(particle):
    #Read simulation results file
	global output_loc
	frequencies = []
	db = []
	with open(output_loc+"/Parasitic_S11_" + str(particle) + ".csv") as csvFile:
		reader = csv.reader(csvFile)
		for row in reader:
			frequencies.append(row[0])
			db.append(row[1])
	csvFile.close()

	# Removing first row of data (Column Labels)
	frequencies.pop(0)
	db.pop(0)

	# Converting data from strings to floats
	for n in range(len(frequencies)):
		frequencies[n] = float(frequencies[n])
		db[n] = float(db[n])
    
    #Determining the maximum S value
	max_s = -1000
	for S in db:
		if S > max_s:
			max_s = S
	cost = max_s
	return cost

def optParamHeader(max_it, pop, n_vars, Vmax, lbest):
    str_out = "Optimization Parameters"
    str_out += "\n-------------------------"
    str_out += "\nmax_it: " + str(max_it)
    str_out += "\npop: " + str(pop)
    str_out += "\nn_vars: " + str(n_vars)
    str_out += "\nVmax: " + str(Vmax)
    str_out += "\nlbest: " + str(lbest)

    return str_out

def swarmDataString(swarm, global_best_cost, global_best_position):
    str_out = "Swarm Data"
    str_out += "\n------------\n"
    for n,p in enumerate(swarm):
                str_out += "Particle "+str(n)+":"
                str_out += particle_to_str(p)
                if lbest:
                    str_out+="Local Best Cost: " + str(p.local_best_cost)
            str_out += "\n\nGlobal Best Cost: " + str(global_best_cost)
            str_out += "\nGlobal Best Position: " + str(global_best_position)
            str_out += "\n\n"
    
    return str_out

def main():
    #Creating output folder
    global output_loc
    if not os.path.exists(output_loc):
        os.makedirs(output_loc)

    #Parameters
    max_it = 40         #Maximum num of iterations
    pop = 10            #Population size
    n_vars = 24         #Number of patches
    Vmax = 15           #AIS velocity maximum
    lbest = -1          #0 -> Global best ||| 1 -> Local Best ||| (0,1) -> Hybrid ||| -1 -> Dynamic Hybrid
    last_i = 0          #The last iteration performed by the previous optimization

    #Storing Parameters
    with open(output_loc+"/Optimization_Params.txt", "w") as f:
        str_out = optParamHeader(max_it, pop, n_vars, Vmax, lbest)
        f.write(str_out)

    #Calling PSO
    PSO(max_it,pop,n_vars,Vmax, lbest, last_i)
    return



if __name__ == "__main__":
    main()
