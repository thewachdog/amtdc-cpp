# for calculation
import numpy as np 
from scipy import interpolate

# for plotting
from bokeh.plotting import figure, show
from bokeh.embed import components

# for pantograph
import turtle
from PIL import Image

# for generating output
import json


valueDict = {
    'messenger_crossSec': '120', 
    'messenger_youngs': '97', 
    'messenger_linearDensity': '1.080', 
    'messenger_inertia': '1.1459', 
    'messenger_tension': '16000', 
    'contact_crossSec': '150', 
    'contact_youngs': '100', 
    'contact_inertia': '1.7905', 
    'contact_linearDensity': '1.35', 
    'contact_tension': '22000', 
    'contact_preSag': '0.055', 
    'dropper_crossSec': '100', 
    'dropper_linearDensity': '0.117', 
    'dropper_youngsModulus': '2', 
    'dropper_inertia': '0.23', 
    'dropper_count': '9', 
    'dropper_location0': '4.5',
    'dropper_location1': '10.25',
    'dropper_location2': '16',
    'dropper_location3': '21.75',
    'dropper_location4': '27.5',
    'dropper_location5': '33.25',
    'dropper_location6': '39',
    'dropper_location7': '44.75',
    'dropper_location8': '50.5',
    'pantograph_count': '1', 
    'pantograph_distance': '1', 
    'pantograph_alpha': '1', 
    'pantograph_beta': '1', 
    'pantograph_model': 'lumpedMass', 
    'degOfFreedomCount': '2', 
    'pantograph_mass0': 1,
    'pantograph_stiffness0': 2,
    'pantograph_damping0': 3,
    'pantograph_mass1': 1,
    'pantograph_stiffness1': 2,
    'pantograph_damping1': 3,
    'catenary_type': 'simple', 
    'miscellaneous_spanLength': '55', 
    'miscellaneous_spanCount': '1', 
    'miscellaneous_massRegArm': '1', 
    'miscellaneous_stiffnessRegArm': '1', 
    'miscellaneous_youngsModulusRegArm': '4', 
    'miscellaneous_crossSecRegArm': '4', 
    'miscellaneous_inertiaRegArm': '0', 
    'miscellaneous_massSteadyArm': '1', 
    'miscellaneous_stiffnessSteadyArm': '1', 
    'miscellaneous_youngsModulusSteadyArm': '6', 
    'miscellaneous_crossSecSteadyArm': '6', 
    'miscellaneous_inertiaSteadyArm': '0', 
    'miscellaneous_encumbrance': '1.2', 
    'miscellaneous_clampMassContact': '0.165', 
    'miscellaneous_clampMassMessenger': '0.195', 
    'csrfmiddlewaretoken': ['XGIgMYqoPE38th4ZXUGDLE4vOWkdkv31m6UEh4kbSBNOqSfxVXw5xPWT5gbjbcjC']
}


# ------------------------------- For Static Deflection of OHE ------------------------------------- #
'''
Some Sample Data for future reference:

length of one span = 55m
length of one element = 0.25m
total number of elements = 55 * 0.25
number of droppers per span = 9
number of dropper points in 2d graph = 9 * 2 = 18
Reg arm = contact wire support
Steady arm = messenger wire support
'''

# %% System parameters
# % rho=linear density C: contact wire, M: Messenger wire, D: droppers.
# % A= cross-sectional area  C: contact wire, M: Messenger wire, D: droppers.
# % E=Young's modulus, C: contact wire, M: Messenger wire, D: droppers.
# % I=Area moment of Inertia, C: contact wire, M: Messenger wire, D: droppers.
# % T=Tension in the wire, C: contact wire, M: Messenger wire

# % Mesenger wire
rho_M = valueDict['messenger_linearDensity']
A_M = valueDict['messenger_crossSec'] * (10**-6)
E_M = valueDict['messenger_youngs'] * (10**9)
I_M = valueDict['messenger_inertia'] * (10**-9)
T_M = - valueDict['messenger_tension']
EA_M = round(E_M * A_M)
EI_M = E_M * I_M

# Contact Values

rho_C = valueDict['contact_linearDensity']
A_C = valueDict['contact_crossSec'] * (10**-6)
E_C = valueDict['contact_youngs'] * (10**9)
I_C = valueDict['contact_inertia'] * (10**-9)
T_C = - valueDict['contact_tension']
EA_C = round(E_C * A_C)
EI_C = E_C * I_C

# Dropper

rho_D = valueDict['dropper_linearDensity']
EA_D = valueDict['dropper_youngsModulus'] * valueDict['dropper_crossSec'] * (10**3)
EI_D = valueDict['dropper_youngsModulus'] * valueDict['dropper_inertia']

No_span = valueDict['miscellaneous_spanCount'] # Number of Spans
# % contact wire support
EA_S_CW = valueDict['miscellaneous_youngsModulusRegArm'] * (10 ** 9) * valueDict['miscellaneous_crossSecRegArm'] * (10 ** -6)
EI_S_CW = valueDict['miscellaneous_youngsModulusRegArm'] * (10 ** 9) * valueDict['miscellaneous_inertiaRegArm'] * (10 ** -6)

# % Messenger wire support
EA_S_MW = valueDict['miscellaneous_youngsModulusSteadyArm'] * (10 ** 9) * valueDict['miscellaneous_crossSecSteadyArm'] * (10 ** -6)
EI_S_MW = valueDict['miscellaneous_youngsModulusSteadyArm'] * (10 ** 9) * valueDict['miscellaneous_inertiaSteadyArm'] * (10 ** -6)

lengthSpan = valueDict['miscellaneous_spanLength'] # Length of one span
T_c = - T_C
T_m = - T_M
w_c = rho_C * 9.8 # % weight density of the contact wire
w_m = rho_M * 9.8 # % weight density of the messenger wire
w_dm = valueDict['miscellaneous_clampMassMessenger'] * 9.8 # %weight of dropper clamp on messenger wire
w_dc = valueDict['miscellaneous_clampMassContact'] * 9.8 # %weight of dropper clamp on contact wire
No_droppers = valueDict['dropper_count']

dropper_location = []
for val in valueDict:
    if val.startswith('dropper_location'):
        dropper_location += [valueDict[val]]
d_c = valueDict['contact_preSag'] # % Maximum presag on the contact wire
encumbrance = valueDict['miscellaneous_encumbrance']

# %% Dropper force and length
# % presag calculated at each droppper location
D_dropper = (4 * d_c / ((dropper_location[No_droppers - 1] - dropper_location[0]) ** 2)) * np.array([dropper_location[x] - dropper_location[0] for x in range(No_droppers)]) * np.array([(dropper_location[No_droppers - 1] - dropper_location[0]) - dropper_location[x] + dropper_location[0] for x in range(No_droppers)])

# % Forces at the droppers
F = [0 for _ in range(No_droppers)]
f1 = dropper_location[0]
f2 = (dropper_location[1] - dropper_location[0])/2
f3 = (T_c * (D_dropper[1] - D_dropper[0])) / (w_c * (dropper_location[1] - dropper_location[0]))
F_1_9 = (f1+f2+f3) * w_c + w_dc
F[0] = F_1_9
F[8] = F_1_9

for i in range(1, No_droppers - 1):
    f4 = (dropper_location[i + 1] - dropper_location[i-1]) / 2
    f5 = (T_c * (D_dropper[i] - D_dropper[i-1])) / (w_c * (dropper_location[i] - dropper_location[i-1]))
    f6 = (T_c * (D_dropper[i+1] - D_dropper[i])) / (w_c * (dropper_location[i+1] - dropper_location[i]))
    F[i] = (f4 - f5 + f6) * w_c + w_dc

r1 = []
# % Reaction force at the messenger wire support
for i in range(No_droppers):
    r1 += [(F[i] * (lengthSpan - dropper_location[i])) / lengthSpan]

r2 = sum(r1)
R_a = ((w_m * lengthSpan) / 2) + r2

F_droppers_m = [0 for _ in range(No_droppers)]
F_droppers_m1 = [0 for _ in range(No_droppers)]

for i in range(1, No_droppers):
    for j in range(i):
        F_droppers_m1[j] = F[j] * 5.75 * (i - j)
    F_droppers_m[i] = sum(F_droppers_m1)
    if i > (No_droppers // 2):
        F_droppers_m[i] = F_droppers_m[No_droppers - i - 1]

c_m = []
# # % Sag calculation at the messenger wire at the dropper locations
for i in range(No_droppers):
    c1 = R_a * dropper_location[i]
    c2 = (w_m * dropper_location[i] ** 2) / 2
    c3 = F_droppers_m[i]
    c_m += [(c1 - c2 - c3) / T_m]
    if i > (No_droppers // 2):
        c_m[i] = c_m[No_droppers - i - 1]

length_dropper = []
# % Calculation of dropper lengths 
for i in range(No_droppers):
    length_dropper += [encumbrance - c_m[i] + D_dropper[i]]

# %% Node coordinates

EL = 0.25 # % length of each element in the wire
numberNodes_w = round(lengthSpan * 4 + 1) # % number of elements for contact wire and messenger wire for a single span
numberElements = No_span * ((numberNodes_w  - 1) * 2 + No_droppers) + 2 * No_span + 2 # % total number of elements (220+220+9) contact wire, messenger wire, droppers
numberElements_W = int(No_span * lengthSpan // EL) # % Total Number of elements for contact wire and messneger wire
numberNodes_W = int(No_span * lengthSpan // EL) + 1 # % Total Number of nodes for contact wire and messneger wire

NC_CW = [[0, 0] for _ in range(numberNodes_W)] #  % nodal coordinates of the contact wire
NC_MW = [[0, 0] for _ in range(numberNodes_W)] # % nodal coordinates of the messenger wire
mw_x = list(np.arange(0, lengthSpan + EL, EL)) # % Discretization of the wire into finite elements for sigle span
Nodes_effective = (2 * numberNodes_W + 2 * No_span + 2) # % Total number of nodes of the OHE system

# %dropper position and support on the contact wire (x coordinates)
cW_x = [0 for _ in range(No_droppers+2)]
for i in range(1, No_droppers + 1):
    cW_x[i] = dropper_location[i - 1]
    cW_x[No_droppers + 1] = lengthSpan

# % dropper position and support on the contact wire (y coordinates)
cW_y = [0 for _ in range(No_droppers + 2)]
for i in range(1, No_droppers+1):
    cW_y[i] = D_dropper[i - 1]
interp1_cy = interpolate.interp1d(cW_x, cW_y)
CW_Y = interp1_cy(mw_x)

# %Interpolation of the dropper point to find out the configuration of the entire contact wire
# CW_Y = interp1(cW_x,cW_y,mw_x);
CW_Y = -CW_Y
MW_x = cW_x # % messenger wire x coordinates
MW_y = [0 for _ in range(No_droppers+2)] # % messenger wire y coordinates for dropper positions
for i in range(No_droppers + 2):
    if i == 0 or i == No_droppers + 1:
        MW_y[i] = 1.2
    else:
        MW_y[i] = length_dropper[i - 1]

# MW_Y = interp1(MW_x,MW_y,mw_x,'spline'); % Interpolation to find out the y coordinates of the messenger wire for finding its initial position
interp1_my = interpolate.InterpolatedUnivariateSpline(MW_x, MW_y)
MW_Y = interp1_my(mw_x)
for i in range(numberNodes_w): # % Assigning the nodal coordinates to the messenger wire
    NC_MW[i][0] = mw_x[i]

# % Node coordinates for contact wire
for j in range(No_span): 
    for i in range(numberNodes_w):
        if j == 0:
            NC_CW[i][0] = mw_x[i]
            NC_CW[i][1] = CW_Y[i]
        else:
            NC_CW[(i + j * numberNodes_w) - j][0] = mw_x[i] + j * lengthSpan
            NC_CW[(i + j * numberNodes_w) - j][1] = CW_Y[i]

# % Node coordinates for Messenger wire
for j in range(No_span):
    for i in range(numberNodes_w):
        if j == 0:
            NC_MW[i][0] = mw_x[i]
            NC_MW[i][1] = MW_Y[i]
        else:
            NC_MW[(i + j * numberNodes_w) - j][0] = mw_x[i] + j * lengthSpan
            NC_MW[(i + j * numberNodes_w) - j][1] = MW_Y[i]

# %% Droppers Node coordinates
NC_droppers_single = [[0, 0] for _ in range(2 * No_droppers)]
for i in range(No_droppers):
    NC_droppers_single[2 * i][0] = dropper_location[i]
    NC_droppers_single[2 * i][1] = - D_dropper[i]
    NC_droppers_single[2 * i + 1][0] = dropper_location[i]
    NC_droppers_single[2 * i + 1][1] = length_dropper[i]     

# %Dropper node coordinates for multiple span
NC_droppers = [[0, 0] for _ in range(No_span * No_droppers * 2)]

for j in range(No_span):
    for i in range(No_droppers * 2):
        if j == 0:
            NC_droppers[i][0] = NC_droppers_single[i][0]
            NC_droppers[i][1] = NC_droppers_single[i][1]
        else:
            NC_droppers[i + j * No_droppers * 2][0] = NC_droppers_single[i][0] + j * lengthSpan
            NC_droppers[i + j * No_droppers * 2][1] = NC_droppers_single[i][1]

# %% Supports at the wires

NC_supports_CW = [[0, 0] for _ in range(2*No_span + 2)]
NC_supports_MW = [[0, 0] for _ in range(2*No_span + 2)]

for i in range(0, len(NC_supports_CW), 2):
    NC_supports_CW[i][1] = 0
    NC_supports_CW[i][0] = (i/2)*lengthSpan
    NC_supports_MW[i][1] = 1.2
    NC_supports_MW[i][0] = (i/2)*lengthSpan

for i in range(1, len(NC_supports_CW), 2):
    NC_supports_CW[i][1] = -1
    NC_supports_CW[i][0] = ((i-1)/2)*lengthSpan
    NC_supports_MW[i][1] = 2.2
    NC_supports_MW[i][0] = ((i-1)/2)*lengthSpan

# %% Combining different matrices
nodeCoordinates = [[0, 0] for _ in range(Nodes_effective)]

# % node coordinates for Contact wire
for i in range(len(NC_CW)):
    for j in range(2):
        nodeCoordinates[i][j] = NC_CW[i][j]

# % node coordinates for Messenger wire
for i in range(len(NC_CW), 2*len(NC_MW)):
    for j in range(2):
        nodeCoordinates[i][j] = NC_MW[i-len(NC_CW)][j]

# % node coordinates for Contact wire Supports
for i in range((2*len(NC_MW)), (2*len(NC_MW)) + (No_span+1)):
    for j in range(2):
        nodeCoordinates[i][j] = NC_supports_CW[2 * (i - (2*len(NC_MW))) + 1][j]

for i in range((2*len(NC_MW)) + (No_span+1), (2*len(NC_MW)) + (No_span+1)+(No_span+1)):
    for j in range(2):
        nodeCoordinates[i][j] = NC_supports_MW[2*(i - ((2*len(NC_MW))+(No_span+1))) + 1][j]

# %% Connection between different elements
elementNodes = []
# % Connection between contact wire elements nodes
for i in range(numberElements_W):
    elementNodes += [[i, i+1]]

# % Connection between messenger wire elements nodes
for i in range(numberElements_W, 2*numberElements_W):
    elementNodes += [[i+1, i+2]]

# % Connection between droppers element nodes
for i in range(2*numberElements_W, 2*numberElements_W+No_span * No_droppers):
    elementNodes += [[int(NC_droppers[2*(i-(2*numberElements_W))][0] // EL), int(NC_droppers[2*(i-(2*numberElements_W))][0]//EL)+1 + int((No_span*lengthSpan)//EL)]]

# %Connections between contact wire support
for i in range((No_span * No_droppers + 2 * numberElements_W), (No_span * No_droppers + 2 * numberElements_W) + No_span + 1):
    elementNodes += [[ int(((i - ((No_span * No_droppers + 2*numberElements_W)))*lengthSpan)//EL), (i - (No_span * No_droppers - 2))]]

# % Connection between Messenger wire support element nodes 
for i in range((No_span * No_droppers + 2 * numberElements_W) + 1 + No_span + 1, (No_span * No_droppers + 2 * numberElements_W) + 3 + 2 * No_span):
    elementNodes += [[int(((i - ((No_span * No_droppers + 2 * numberElements_W) + 1 + No_span + 1)) * lengthSpan)//EL)  + numberNodes_W, i - (No_span * No_droppers - 1)]]
xx = [_[0] for _ in nodeCoordinates]
yy = [_[1] for _ in nodeCoordinates]

# %% Number of degrees of freedom and effective number of nodes

GDof = 3 * int(2*((No_span*lengthSpan/EL)+1)+(2*No_span+2))
Nodes_effective = (2*numberNodes_W+2*No_span+2)

# %% external force
# % Self weight
force = [0 for _ in range(GDof)]
Length_Element = 0.25

# % Self weight of the contact wire
for i in range(Nodes_effective, Nodes_effective + numberElements_W):
    force[i] = force[i] - (9.81 * rho_C * Length_Element / 4)

# %% Stiffness formation
stiffness = np.zeros((GDof, GDof)).tolist()

# % Contact wire
for e in range(numberElements_W):
    # % elementDof: element degrees of freedom (Dof)
    indice = elementNodes[e]
    elementDof = np.concatenate((indice, np.array(indice)+ Nodes_effective, np.array(indice)+ 2*Nodes_effective), axis=None)
    nn = len(indice)
    xa = xx[indice[1]] - xx[indice[0]]
    ya = yy[indice[1]] - yy[indice[0]]
    length_element = np.sqrt(xa*xa+ya*ya)
    cosa = round(xa/length_element)
    sena = ya/length_element
    ll = round(length_element, 2)
    L = np.zeros((6, 6)).tolist()
    oneu = [[1, -1], [-1, 1]]
    oneu2 = [[1, -1], [1, -1]]
    oneu3 = [[1, 1], [-1, -1]]
    oneu4 = [[4, 2], [2, 4]]
    A=[[round(4*EI_C/ll, 4) + round(4*ll*T_C/30, 4), round(2*EI_C/ll, 4) - round(ll*T_C/30, 4)], [round(2*EI_C/ll, 4) - round(ll*T_C/30, 4), round(4*EI_C/ll, 4) + round(4*ll*T_C/30, 4)]]
    k1 = np.zeros((6, 6)).tolist()
    for i in range(6):
        for j in range(6):
            if i > 3:
                if i == j:
                    L[i][j] = 1
            else:
                if j < 4:
                    identity = [[1, 0], [0, 1]]
                    t = identity[i%2][j%2]
                    if (i < 2 and j < 2) or (i > 1 and j > 1):
                        L[i][j] = t * cosa
                    elif j < 2:
                        L[i][j] = - t * sena
                    else:
                        L[i][j] = t * sena
            if i < 2 and j < 2:
                    k1[i][j] = round(EA_C/ll*oneu[i][j])
            else:
                if i >= 2 and j >= 2:
                    if i in range(2, 4) and j in range(2, 4):
                        k1[i][j] = round(((36*T_C/30*ll)+(12*EI_C/ll**3))*oneu[i - 2][j - 2])
                    elif i >= 4 and j >= 4:
                        k1[i][j] = A[i-4][j-4]
                    elif i > j:
                        k1[i][j] = round(((3*T_C/30)+(6*EI_C/ll**2))*oneu2[i - 4][j - 2])
                    elif j > i:
                        k1[i][j] = round(((3*T_C/30)+(6*EI_C/ll**2))*oneu3[i - 2][j - 4])
    L = np.matrix(L)
    temp = (L.transpose()* k1 * L).tolist()
    for i in range(len(elementDof)):
        for j in range(len(elementDof)):
            stiffness[elementDof[i]][elementDof[j]] += temp[i][j]

# % Messenger wire
for e in range(numberElements_W, 2*numberElements_W):
    # % elementDof: element degrees of freedom (Dof)
    indice = elementNodes[e]
    elementDof1 = np.concatenate((indice, np.array(indice)+ Nodes_effective, np.array(indice)+ 2*Nodes_effective), axis=None)
    nn = len(indice)
    xa = xx[indice[1]] - xx[indice[0]]
    ya = yy[indice[1]] - yy[indice[0]]
    length_element = np.sqrt(xa*xa+ya*ya)
    cosa = round(xa/length_element)
    sena = ya/length_element
    ll = length_element
    L = np.zeros((6, 6)).tolist()
    oneu = [[1, -1], [-1, 1]]
    oneu2 = [[1, -1], [1, -1]]
    oneu3 = [[1, 1], [-1, -1]]
    oneu4 = [[4, 2], [2, 4]]
    A=[[round(4*EI_M/ll, 4) + round(4*ll*T_M/30, 4), round(2*EI_M/ll, 4) - round(ll*T_M/30, 4)], [round(2*EI_M/ll, 4) - round(ll*T_M/30, 4), round(4*EI_M/ll, 4) + round(4*ll*T_M/30, 4)]]
    k1 = np.zeros((6, 6)).tolist()
    for i in range(6):
        for j in range(6):
            if i > 3:
                if i == j:
                    L[i][j] = 1
            else:
                if j < 4:
                    identity = [[1, 0], [0, 1]]
                    temp = identity[i%2][j%2]
                    if (i < 2 and j < 2) or (i > 1 and j > 1):
                        L[i][j] = temp * cosa
                    elif j < 2:
                        L[i][j] = - temp * sena
                    else:
                        L[i][j] = temp * sena
            if i < 2 and j < 2:
                    k1[i][j] = EA_M/ll*oneu[i][j]
            else:
                if i >= 2 and j >= 2:
                    if i in range(2, 4) and j in range(2, 4):
                        k1[i][j] = ((36*T_M/30*ll)+(12*EI_M/ll**3))*oneu[i - 2][j - 2]
                    elif i >= 4 and j >= 4:
                        k1[i][j] = A[i - 4][j - 4]
                    elif i > j:
                        k1[i][j] = ((3*T_M/30)+(6*EI_M/ll**2))*oneu2[i - 4][j - 2]
                    elif j > i:
                        k1[i][j] = ((3*T_M/30)+(6*EI_M/ll**2))*oneu3[i - 2][j - 4]
    L = np.matrix(L)
    temp = (L.transpose() * k1 * L).tolist()
    for i in range(len(elementDof1)):
        for j in range(len(elementDof1)):
            stiffness[elementDof1[i]][elementDof1[j]] += temp[i][j]

# % Dropper
for e in range(2*numberElements_W, 2 * numberElements_W + No_droppers * No_span):
    # % elementDof: element degrees of freedom (Dof)
    indice = elementNodes[e]
    elementDof2 = np.concatenate((indice, np.array(indice)+ Nodes_effective, np.array(indice)+ 2*Nodes_effective), axis=None)
    nn = len(indice)
    xa = xx[indice[1]] - xx[indice[0]]
    ya = yy[indice[1]] - yy[indice[0]]
    length_element = np.sqrt(xa*xa+ya*ya)
    cosa = round(xa/length_element)
    sena = ya/length_element
    ll = length_element
    L = np.zeros((6, 6)).tolist()
    oneu = [[1, -1], [-1, 1]]
    oneu2 = [[1, -1], [1, -1]]
    oneu3 = [[1, 1], [-1, -1]]
    oneu4 = [[4, 2], [2, 4]]
    T = 0
    A=[[(4*EI_D/ll)+(4*ll*T/30), (2*EI_D/ll)-(ll*T/30)], [(2*EI_D/ll)-(ll*T/30), (4*EI_D/ll)+(4*ll*T/30)]]
    k1 = np.zeros((6, 6)).tolist()
    for i in range(6):
        for j in range(6):
            if i > 3:
                if i == j:
                    L[i][j] = 1
            else:
                if j < 4:
                    identity = [[1, 0], [0, 1]]
                    temp = identity[i%2][j%2]
                    if (i < 2 and j < 2) or (i > 1 and j > 1):
                        L[i][j] = temp * cosa
                    elif j < 2:
                        L[i][j] = - temp * sena
                    else:
                        L[i][j] = temp * sena
            if i < 2 and j < 2:
                    k1[i][j] = EA_D/ll*oneu[i][j]
            else:
                if i >= 2 and j >= 2:
                    if i in range(2, 4) and j in range(2, 4):
                        k1[i][j] = ((36*T/30*ll)+(12*EI_D/ll**3))*oneu[i - 2][j - 2]
                    elif i >= 4 and j >= 4:
                        k1[i][j] = A[i-4][j-4]
                    elif i > j:
                        k1[i][j] = ((3*T/30)+(6*EI_D/ll**2))*oneu2[i - 4][j - 2]
                    elif j > i:
                        k1[i][j] = ((3*T/30)+(6*EI_D/ll**2))*oneu3[i - 2][j - 4]
    L = np.matrix(L)
    temp = (L*k1*L.transpose()).tolist()
    for i in range(len(elementDof2)):
        for j in range(len(elementDof2)):
            stiffness[elementDof2[i]][elementDof2[j]] += temp[i][j]

# % contact wire support
for e in range(2 * numberElements_W + No_droppers * No_span, (No_span * No_droppers + 2 * numberElements_W) + No_span + 1):
    # % elementDof: element degrees of freedom (Dof)
    indice = elementNodes[e]
    elementDof3 = np.concatenate((indice, np.array(indice)+ Nodes_effective, np.array(indice)+ 2*Nodes_effective), axis=None)
    nn = len(indice)
    xa = xx[indice[1]] - xx[indice[0]]
    ya = yy[indice[1]] - yy[indice[0]]
    length_element = np.sqrt(xa*xa+ya*ya)
    cosa = round(xa/length_element)
    sena = ya/length_element
    ll = length_element
    L = np.zeros((6, 6)).tolist()
    oneu = [[1, -1], [-1, 1]]
    oneu2 = [[1, -1], [1, -1]]
    oneu3 = [[1, 1], [-1, -1]]
    oneu4 = [[4, 2], [2, 4]]
    T = 0
    A=[[(4*EI_S_CW/ll)+(4*ll*T/30), (2*EI_S_CW/ll)-(ll*T/30)], [(2*EI_S_CW/ll)-(ll*T/30), (4*EI_S_CW/ll)+(4*ll*T/30)]]
    k1 = np.zeros((6, 6)).tolist()
    for i in range(6):
        for j in range(6):
            if i > 3:
                if i == j:
                    L[i][j] = 1
            else:
                if j < 4:
                    identity = [[1, 0], [0, 1]]
                    temp = identity[i%2][j%2]
                    if (i < 2 and j < 2) or (i > 1 and j > 1):
                        L[i][j] = temp * cosa
                    elif j < 2:
                        L[i][j] = - temp * sena
                    else:
                        L[i][j] = temp * sena
            if i < 2 and j < 2:
                    k1[i][j] = EA_S_CW/ll*oneu[i][j]
            else:
                if i >= 2 and j >= 2:
                    if i in range(2, 4) and j in range(2, 4):
                        k1[i][j] = ((36*T/30*ll)+(12*EI_S_CW/ll**3))*oneu[i - 2][j - 2]
                    elif i >= 4 and j >= 4:
                        k1[i][j] = A[i-4][j-4]
                    elif i > j:
                        k1[i][j] = ((3*T/30)+(6*EI_S_CW/ll**2))*oneu2[i - 4][j - 2]
                    elif j > i:
                        k1[i][j] = ((3*T/30)+(6*EI_S_CW/ll**2))*oneu3[i - 2][j - 4]

    L = np.matrix(L)
    temp = (L*k1*L.transpose()).tolist()
    for i in range(len(elementDof3)):
        for j in range(len(elementDof3)):
            stiffness[elementDof3[i]][elementDof3[j]] += temp[i][j]

#  % Messenger wire support

for e in range((No_span * No_droppers + 2 * numberElements_W) + 1 + No_span, (No_span * No_droppers + 2 * numberElements_W) + 2 + 2 * No_span):
    # % elementDof: element degrees of freedom (Dof)
    indice = elementNodes[e]
    elementDof4 = np.concatenate((indice, np.array(indice)+ Nodes_effective, np.array(indice)+ 2*Nodes_effective), axis=None)
    nn = len(indice)
    xa = xx[indice[1]] - xx[indice[0]]
    ya = yy[indice[1]] - yy[indice[0]]
    length_element = np.sqrt(xa*xa+ya*ya)
    cosa = round(xa/length_element)
    sena = ya/length_element
    ll = length_element
    L = np.zeros((6, 6)).tolist()
    oneu = [[1, -1], [-1, 1]]
    oneu2 = [[1, -1], [1, -1]]
    oneu3 = [[1, 1], [-1, -1]]
    oneu4 = [[4, 2], [2, 4]]
    T = 0
    A=[[(4*EI_S_MW/ll)+(4*ll*T/30), (2*EI_S_MW/ll)-(ll*T/30)], [(2*EI_S_MW/ll)-(ll*T/30), (4*EI_S_MW/ll)+(4*ll*T/30)]]
    k1 = np.zeros((6, 6)).tolist()
    for i in range(6):
        for j in range(6):
            if i > 3:
                if i == j:
                    L[i][j] = 1
            else:
                if j < 4:
                    identity = [[1, 0], [0, 1]]
                    temp = identity[i%2][j%2]
                    if (i < 2 and j < 2) or (i > 1 and j > 1):
                        L[i][j] = temp * cosa
                    elif j < 2:
                        L[i][j] = - temp * sena
                    else:
                        L[i][j] = temp * sena
            if i < 2 and j < 2:
                    k1[i][j] = round(EA_S_MW/ll*oneu[i][j])
            else:
                if i >= 2 and j >= 2:
                    if i in range(2, 4) and j in range(2, 4):
                        k1[i][j] = ((36*T/30*ll)+(12*EI_S_MW/ll**3))*oneu[i - 2][j - 2]
                    elif i >= 4 and j >= 4:
                        k1[i][j] = A[i-4][j-4]
                    elif i > j:
                        k1[i][j] = ((3*T/30)+(6*EI_S_MW/ll**2))*oneu2[i - 4][j - 2]
                    elif j > i:
                        k1[i][j] = ((3*T/30)+(6*EI_S_MW/ll**2))*oneu3[i - 2][j - 4]
    L = np.matrix(L)
    temp = (L*k1*L.transpose()).tolist()
    for i in range(len(elementDof4)):
        for j in range(len(elementDof4)):
            stiffness[elementDof4[i]][elementDof4[j]] += temp[i][j]

# %% Mass Matrix formulation

mass = [[0 for _ in range(GDof)] for _ in range(GDof)]

# % Contact wire
for e in range(numberElements_W):
    # % elementDof: element degrees of freedom (Dof)
    indice = elementNodes[e]
    elementDof = np.concatenate((indice, np.array(indice)+ Nodes_effective, np.array(indice)+ 2*Nodes_effective), axis=None)
    xa = xx[indice[1]] - xx[indice[0]]
    ya = yy[indice[1]] - yy[indice[0]]
    length_element = np.sqrt(xa*xa+ya*ya)
    cosa = int(xa//length_element)
    sena = int(ya/length_element)
    ll = length_element
    L = [[0 for _ in range(len(elementDof))] for i in range(len(elementDof))]

    oneu = [[2/6, 1/6], [1/6, 2/6]]
    oneu2 = [[26/70, 9/70], [9/70, 26/70]]
    oneu3 = [[22 * ll/420, -13 * ll/420],[13 * ll/420, -22 * ll/420]]
    oneu4 = [[4 * ll*ll/420, -3 * ll*ll/420], [-3 * ll*ll/420, 4 * ll*ll/420]]
    oneu3_t = [[oneu3[j][i] for j in range(len(oneu3))] for i in range(len(oneu3[0]))]

    m1 = [[0 for _ in range(len(elementDof))] for i in range(len(elementDof))]
    for i in range(len(elementDof)):
        for j in range(len(elementDof)):
            if i > 3 and i == j:
                    L[i][j] = 1
            else:
                if j < 4:
                    identity = [[1, 0], [0, 1]]
                    temp = identity[i%2][j%2]
                    if (i < 2 and j < 2) or (i > 1 and j > 1):
                        L[i][j] = temp * cosa
                    elif j < 2:
                        L[i][j] = - temp * sena
                    else:
                        L[i][j] = temp * sena

            if i < 2 and j < 2:
                    m1[i][j] = oneu[i][j] * rho_C * ll
            else:
                if i >= 2 and j >= 2:
                    if i in range(2, 4) and j in range(2, 4):
                        m1[i][j] = oneu2[i - 2][j - 2] * rho_C * ll
                    elif i >= 4 and j >= 4:
                        m1[i][j] = oneu4[i - 4][i - 4] * rho_C * ll
                    elif i > j:
                        m1[i][j] = oneu3_t[i - 4][j - 2] * rho_C * ll
                    elif j > i:
                        m1[i][j] = oneu3[i - 2][j - 4] * rho_C * ll
    L = np.matrix(L)
    temp = (L*m1*L.transpose()).tolist()
    for i in range(len(elementDof)):
        for j in range(len(elementDof)):
            mass[elementDof[i]][elementDof[j]] += temp[i][j]

# % Messenger wire
for e in range(numberElements_W, 2*numberElements_W):
    # % elementDof: element degrees of freedom (Dof)
    indice = elementNodes[e]
    elementDof = np.concatenate((indice, np.array(indice)+ Nodes_effective, np.array(indice)+ 2*Nodes_effective), axis=None)
    xa = xx[indice[1]] - xx[indice[0]]
    ya = yy[indice[1]] - yy[indice[0]]
    length_element = np.sqrt(xa*xa+ya*ya)
    cosa = int(xa/length_element)
    sena = int(ya/length_element)
    ll = length_element
    L = [[0 for _ in range(len(elementDof))] for i in range(len(elementDof))]

    oneu = [[2/6, 1/6], [1/6, 2/6]]
    oneu2 = [[26/70, 9/70], [9/70, 26/70]]
    oneu3 = [[22 * ll/420, -13 * ll/420],[13 * ll/420, -22 * ll/420]]
    oneu4 = [[4 * ll*ll/420, -3 * ll*ll/420], [-3 * ll*ll/420, 4 * ll*ll/420]]
    oneu3_t = [[oneu3[j][i] for j in range(len(oneu3))] for i in range(len(oneu3[0]))]

    m1 = [[0 for _ in range(len(elementDof))] for i in range(len(elementDof))]
    for i in range(len(elementDof)):
        for j in range(len(elementDof)):
            if i > 3 and i == j:
                    L[i][j] = 1
            else:
                if j < 4:
                    identity = [[1, 0], [0, 1]]
                    temp = identity[i%2][j%2]
                    if (i < 2 and j < 2) or (i > 1 and j > 1):
                        L[i][j] = temp * cosa
                    elif j < 2:
                        L[i][j] = - temp * sena
                    else:
                        L[i][j] = temp * sena

            if i < 2 and j < 2:
                    m1[i][j] = oneu[i][j] * rho_M * ll
            else:
                if i >= 2 and j >= 2:
                    if i in range(2, 4) and j in range(2, 4):
                        m1[i][j] = oneu2[i - 2][j - 2] * rho_M * ll
                    elif i >= 4 and j >= 4:
                        m1[i][j] = oneu4[i - 4][i - 4] * rho_M * ll
                    elif i > j:
                        m1[i][j] = oneu3_t[i - 4][j - 2] * rho_M * ll
                    elif j > i:
                        m1[i][j] = oneu3[i - 2][j - 4] * rho_M * ll
    L = np.matrix(L)
    temp = (L*m1*L.transpose()).tolist()
    for i in range(len(elementDof)):
        for j in range(len(elementDof)):
            mass[elementDof[i]][elementDof[j]] += temp[i][j]

# % Droppers
for e in range(2 * numberElements_W, 2 * numberElements_W + No_droppers * No_span):
    # % elementDof: element degrees of freedom (Dof)
    indice = elementNodes[e]
    elementDof = np.concatenate((indice, np.array(indice)+ Nodes_effective, np.array(indice)+ 2*Nodes_effective), axis=None)
    xa = xx[indice[1]] - xx[indice[0]]
    ya = yy[indice[1]] - yy[indice[0]]
    length_element = np.sqrt(xa*xa+ya*ya)
    cosa = int(xa/length_element)
    sena = int(ya/length_element)
    ll = length_element
    L = [[0 for _ in range(len(elementDof))] for i in range(len(elementDof))]

    oneu = [[2/6, 1/6], [1/6, 2/6]]
    oneu2 = [[26/70, 9/70], [9/70, 26/70]]
    oneu3 = [[22 * ll/420, -13 * ll/420],[13 * ll/420, -22 * ll/420]]
    oneu4 = [[4 * ll*ll/420, -3 * ll*ll/420], [-3 * ll*ll/420, 4 * ll*ll/420]]
    oneu3_t = [[oneu3[j][i] for j in range(len(oneu3))] for i in range(len(oneu3[0]))]

    m1 = [[0 for _ in range(len(elementDof))] for i in range(len(elementDof))]
    for i in range(len(elementDof)):
        for j in range(len(elementDof)):
            if i > 3 and i == j:
                    L[i][j] = 1
            else:
                if j < 4:
                    identity = [[1, 0], [0, 1]]
                    temp = identity[i%2][j%2]
                    if (i < 2 and j < 2) or (i > 1 and j > 1):
                        L[i][j] = temp * cosa
                    elif j < 2:
                        L[i][j] = - temp * sena
                    else:
                        L[i][j] = temp * sena

            if i < 2 and j < 2:
                    m1[i][j] = oneu[i][j] * rho_D * ll
            else:
                if i >= 2 and j >= 2:
                    if i in range(2, 4) and j in range(2, 4):
                        m1[i][j] = oneu2[i - 2][j - 2] * rho_D * ll
                    elif i >= 4 and j >= 4:
                        m1[i][j] = oneu4[i - 4][i - 4] * rho_D * ll
                    elif i > j:
                        m1[i][j] = oneu3_t[i - 4][j - 2] * rho_D * ll
                    elif j > i:
                        m1[i][j] = oneu3[i - 2][j - 4] * rho_D * ll
    L = np.matrix(L)
    temp = (L*m1*L.transpose()).tolist()
    for i in range(len(elementDof)):
        for j in range(len(elementDof)):
            mass[elementDof[i]][elementDof[j]] += temp[i][j]

# %% Prescribed Degree of freedom
# % Prescribed Degree of freedom for supports
prescribedDof = np.zeros(6 * (2 * No_span + 2)).tolist()
i = 0
for e in range((2 * numberElements_W + 9 * No_span), (No_span * 9 + 2 * numberElements_W) + 2 + 2 * No_span):
    indice = elementNodes[e]
    prescribedDof1 = [*indice, *[Nodes_effective + _ for _ in indice], *[2*Nodes_effective + _ for _ in indice]]
    for j in range(6):
        prescribedDof[j+i] = prescribedDof1[j]
    i += 6

prescribedDof2 = [0 for _ in range(No_droppers * No_span)]
# % Prescribed Degree of freedom for droppers
for e in range(2 * numberElements_W, 2 * numberElements_W + No_droppers * No_span):
    indice = elementNodes[e]
    prescribedDof2[e - 2 * numberElements_W] = indice[0] + Nodes_effective
del prescribedDof[2::6]

# % Prescribed Degree of freedom for droppers  
q = len(prescribedDof)
for i in range(len(prescribedDof), len(prescribedDof) + No_droppers * No_span):
    prescribedDof += [prescribedDof2[i - q]]

# %% Solution
# % solution
activeDof = []
for val in range(GDof):
    if val not in prescribedDof:
        activeDof += [val]

s_activ = np.zeros((len(activeDof), len(activeDof)))
f_activ = []
for i in range(len(activeDof)):
    f_activ += [force[activeDof[i]]]
    for j in range(len(activeDof)):
        s_activ[i][j] = stiffness[activeDof[i]][activeDof[j]]
s_activ = np.linalg.pinv(s_activ)
f_activ = np.matrix(f_activ)

U = np.matmul(s_activ, f_activ.transpose()).tolist()
displacements = np.zeros(GDof)
for val in range(len(activeDof)):
    displacements[activeDof[val]] = U[val][0]

# % reactions
F1 = np.matmul(np.matrix(stiffness), np.matrix(displacements).transpose())
reactions = []
for val in prescribedDof:
    reactions += [F1[val]]
q = []
for val in range(Nodes_effective+numberNodes_W, ((Nodes_effective+numberNodes_W)+1)+numberElements_W):
    q += [displacements[val]]
MW_selfweight = []
for i in range(numberNodes_W):
    MW_selfweight += [NC_MW[i][1] + q[i]]

#%% Plotting the results for the static condition after applying self weight

#% plot for contact wire

#% adding the initial configuration of contact wire to the results obtained
for j in range(No_span):
    for i in range(numberNodes_w):
        if j == 0:
            displacements[Nodes_effective + i] += CW_Y[i]
        else:
            displacements[Nodes_effective + i + (j - 1) * numberNodes_w - j - 1] += CW_Y[i]
