__author__ = "A. Dey, MSc"
__license__ = "GNU General Public License v3.0"
__version__ = "1.0"
__maintainer__ = "A. Dey"
__status__ = "Pre-Production"

import NPROC 
from sklearn.model_selection import train_test_split


#compute the distance of each of the solutions from the origin
def L2_Dist(points):
    l2 = []
    for i in range(len(points)):
        l2.append((points[i][0]**2+points[i][1]**2+points[i][2]**2+points[i][3]**2)**0.5)
    #print(np.where(l2 == min(l2))[0][0])
    index = np.where(l2 == min(l2))[0][0]
    return index

#compute the distance of each of the solutions from the g-best
def GlobalL2_Dist(points,gbest):
    return ((gbest[0]-points[0])**2+(gbest[1]-points[1])**2+(gbest[2]-points[2])**2+(gbest[3]-points[3])**2)**0.5


# compute the fitness function

def fitness(Train_data_zScore,y,POS_INDEX,kernel_summer,df1_pval,df1_foldchng):
        #feature selected training data.
        model = NPROC.npc() 
        X_ = Train_data_zScore[:,POS_INDEX]
        
        #calculate the p-value for particular patricle
        
        pav1 = df1_pval[['p Value_With_Lesion_Vs_Normal']].T[POS_INDEX]
        pav2 = df1_pval[['p Value_Without_Lesion_Vs_Normal']].T[POS_INDEX]
        pval = (pav1.T[['p Value_With_Lesion_Vs_Normal']].mean().values[0]+
              pav2.T[['p Value_Without_Lesion_Vs_Normal']].mean().values[0])/2
        
        #calculate the fold-change for particular patricle
        
        fc_Nor = (df1_foldchng[['FC (abs)_With_Lesion_Vs_Normal']].T[POS_INDEX].mean().values[0]+
              df1_foldchng['FC (abs)_Without_Lesion_Vs_Normal'].T[POS_INDEX].mean())/2
        
        
        #calculate the accuracy using MKL and count the no. of genes
        
        X_train , X_test , y_train, y_test = train_test_split(X_, y, random_state=0)
        clfr = model.npc(X_train, y_train, kernel_summer, method = "mkl")
        fitted_score = model.predict(clfr,X_test)
        Accuracy = fitted_score[0]
        gene_no = len(POS_INDEX) #No. of genes

        return [Accuracy, gene_no, fc_Nor, pval]

def Euclidian(points):
    return (points[0]**2+points[1]**2+points[2]**2+points[3]**2)**0.5

def particle_distance(x,y):
    return ((x[0]-y[0])**2+(x[1]-y[1])**2+(x[2]-y[2])**2+(x[3]-y[3])**2)**0.5

def comp_rad():
    return 4


# A function that returns a vector of dimention of 4X50 which represent the lbest of a particle
def lbest_particles(Cost,radius):
    pdm = []
    pd = []
    lbest = []
    neighbour = []
    #radius = 5

    for i in range(len(cost)):
        for j in range(len(cost)):
            if i != j:
                pd.append(particle_distance(cost[i],cost[j]))
        pdm.append(pd)
        pd = []


    Swarm = [] # this contains the nearest particles for each of the particle in the swarm 
    for j in range(NOP):
        for i in np.sort(pdm[j])[:radius]:
            neighbour.append(np.where(pdm[j]==i)[0][0])
        Swarm.append(neighbour)
        neighbour = []


    for i in Swarm:
        E = []
        for j in i:
            E.append(Euclidian(cost[j]))
        lbest.append(cost[np.where(E==min(E))[0][0]])
        #lbest.append(pareto.eps_sort(E)[L2_Dist(pareto.eps_sort(E))])

    return np.array(lbest)
def cost_gen(Pos,Train_data_zScore,position):
    #get the position of 1's in each of the gene exp. of a particle
    POS_INDEX = np.where(Pos==1)
    POS_INDEX = POS_INDEX[0]
    POS_INDEX = np.ndarray.tolist(POS_INDEX)
    #....... FITNESS FUNCTION ............
    return fitness(Train_Data,POS_INDEX,[position[i,19950], position[i,19951], position[i,19952]])