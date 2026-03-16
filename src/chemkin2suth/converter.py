import numpy as np 
from scipy.optimize import curve_fit

def omega(idx, T_star, delta_star):
    #calculate collision integrals
    #idx: collision integral index
    #T_star: normalized temperature
    #delta_star: reduced dipole moument
    #omega_: collision integral
    if idx == 1:
        #- Table 12.1 "Chemically reacting flow" (2003)
        a = [1.0548, 0.15504, 0.55909, 2.1705, 0.093193, 1.5]
        #- Equation (12.19) "Chemically reacting flow" (2003)
        omega_ = a[0]*(T_star**(-a[1]))+(T_star+a[2])**(-a[3]) #Equation (12.6)
        f = 1.0 + (np.exp(a[4]/T_star)-np.exp(-a[5]/T_star))*(delta_star)**2/(2.0+2.5*delta_star)

        omega_ *= f
    elif idx ==2:
        #- Table 12.1 "Chemically reacting flow" (2003) 
        b = [1.0413, 0.11930, 0.43628, 1.6041, 0.095661, 2.0]
        #- Equation (12.20) "Chemically reacting flow" (2003)
        omega_ = b[0]*(T_star**(-b[1]))+(T_star+b[2])**(-b[3])#Equation (12.7)
        f = 1.0 + (np.exp(b[4]/T_star)-np.exp(-b[5]/T_star))*(delta_star)**2/(2.0+2.5*delta_star)

        omega_ *= f
    else:
        raise ValueError('invalid input for collision integrals')
    
    return omega_

def molecularWeight(elem_list):
    #calculate molecular weight of species
    #elem_list: list of elements of species elem_list = [(element, number),....]
    # molecularWeight_: molecular wieght of species (kg/kmol)
    element_dic_ori = {
    "E":    0,
    "e":    5.45e-4,
    "H":    1.00797,
    "D":    2.01410,
    "T":    3.01604,
    "He":   4.00260,
    "Li":   6.93900,
    "Be":   9.01220,
    "B":   10.81100,
    "C":   12.01115,
    "N":   14.00670,
    "O":   15.99940,
    "F":   18.99840,
    "Ne":  20.18300,
    "Na":  22.98980,
    "Mg":  24.31200,
    "Al":  26.98150,
    "Si":  28.08600,
    "P":   30.97380,
    "S":   32.06400,
    "Cl":  35.45300,
    "Ar":  39.94800,
    "K":   39.10200,
    "Ca":  40.08000,
    "Sc":  44.95600,
    "Ti":  47.90000,
    "V":   50.94200,
    "Cr":  51.99600,
    "Mn":  54.93800,
    "Fe":  55.84700,
    "Co":  58.93320,
    "Ni":  58.71000,
    "Cu":  63.54000,
    "Zn":  65.37000,
    "Ga":  69.72000,
    "Ge":  72.59000,
    "As":  74.92160,
    "Se":  78.96000,
    "Br":  79.90090,
    "Kr":  83.80000,
    "Rb":  85.47000,
    "Sr":  87.62000,
    "Y":   88.90500,
    "Zr":  91.22000,
    "Nb":  92.90600,
    "Mo":  95.94000,
    "Tc":  99.00000,
    "Ru": 101.07000,
    "Rh": 102.90500,
    "Pd": 106.40000,
    "Ag": 107.87000,
    "Cd": 112.40000,
    "In": 114.82000,
    "Sn": 118.69000,
    "Sb": 121.75000,
    "Te": 127.60000,
    "I":  126.90440,
    "Xe": 131.30000,
    "Cs": 132.90500,
    "Ba": 137.34000,
    "La": 138.91000,
    "Ce": 140.12000,
    "Pr": 140.90700,
    "Nd": 144.24000,
    "Pm": 145.00000,
    "Sm": 150.35000,
    "Eu": 151.96000,
    "Gd": 157.25000,
    "Tb": 158.92400,
    "Dy": 162.50000,
    "Ho": 164.93000,
    "Er": 167.26000,
    "Tm": 168.93400,
    "Yb": 173.04000,
    "Lu": 174.99700,
    "Hf": 178.49000,
    "Ta": 180.94800,
    "W":  183.85000,
    "Re": 186.20000,
    "Os": 190.20000,
    "Ir": 192.20000,
    "Pt": 195.09000,
    "Au": 196.96700,
    "Hg": 200.59000,
    "Tl": 204.37000,
    "Pb": 207.19000,
    "Bi": 208.98000,
    "Po": 210.00000,
    "At": 210.00000,
    "Rn": 222.00000,
    "Fr": 223.00000,
    "Ra": 226.00000,
    "Ac": 227.00000,
    "Th": 232.03800,
    "Pa": 231.00000,
    "U":  238.03000,
    "Np": 237.00000,
    "Pu": 242.00000,
    "Am": 243.00000,
    "Cm": 247.00000,
    "Bk": 249.00000,
    "Cf": 251.00000,
    "Es": 254.00000,
    "Fm": 253.00000
    }
    element_dic = {k.casefold(): v for k, v in element_dic_ori.items()}
    molecularWeight_ = 0
    for element, cnt in elem_list:
        try:
            molecularWeight_ += element_dic[element.casefold()]*cnt
        except:
            raise ValueError('unknown element {:s}'.format(element))
    return molecularWeight_

def read_thermo(f):
    #read thermo data file in chemkin format (composition data only)
    #species_dic: composition dic of species {specie name: [(element name, element count),...]}
    species_dic = {}
    lines = f.readlines()
    for line in lines:
        line = line.strip()
        if len(line) != 80:
            continue
        if line[0] == '!':
            continue
        if line[-1] != '1':
            continue 
        species_name = line[0:16].strip()
        atomics = ''.join(line[24:44].split())
        element_list = []
        buff = ''
        number = ''
        for ch in atomics:
            if ch.isalpha():
                if number:
                    element_list.append((buff,int(number)))
                    number = ''
                    buff = ''
                    buff += ch
                else:
                    buff += ch
            elif ch.isdigit() and buff:
                number += ch
            else:
                raise ValueError('incorrect thermo format in line: {:s}'.format(line))
        element_list.append((buff,int(number)))
            
        species_dic[species_name] = element_list
    return species_dic



def read_trans(f):
    #read tranportation data file in chemkin format
    #f: file object
    #trans_dic: tranportation data dictionray. trans = {species: [geometry, welldepth, diameter, dipole, polar, relaxation],...}
    trans_dic = {}
    lines = f.readlines()
    for line in lines:
        line = line.split('!')[0].strip()
        if len(line.split()) == 7:
            line = line.split()
            trans_dic[line[0]] = [float(val) for val in line[1:]]

    return trans_dic

def mu(trans, W, T):
    #calculate dynamic viscosity
    #trans: transportation data trans = [speices name, geometry, welldepth, diameter, dipole, polar, relaxation]
    #W: molecular weight (kg/kmol)
    #T: temperature (K)
    #mu_: dynamic viscosity (kg/m/s)

    kb = 1.380649e-23#boltzmann constant (m^2kg/K/s^2)
    Na = 6.02214076e23#Avogadro constant (mol^-1)
    Pi = 3.1415926

    geometry, welldepth, diameter, dipole, polar, relaxation= trans
    sigma = 1e-10 * diameter#collision diameter (m)
    #- Equation (12.8) "Chemically reacting flow" (2003)
    T_star = T / welldepth 
    epsilon = kb * welldepth#(J)
    miu = 1e-21*(10**(-3.5))*dipole#(m^(3/2)J^(1/2))
    delta_star  = 1/2*miu**2/(epsilon*sigma**3)
    mu_ = 5.0/16.0*(W/1000.0/Na*kb*T/Pi)**0.5/(sigma**2*omega(2,T_star,delta_star))
    return mu_

def sutherland(T, As, Ts):
    #sutherland expression for curve fitting
    #\mu = A_s \frac{\sqrt{T}}{1 + T_s / T}
    #As, Ts: model parameters
    #mu_: dynamic viscosity (kg/m/s)
    mu_ = As * ((T)**0.5) / (1+Ts/T)
    return mu_

def convert_transport_file(trans_input_file, thermo_input_file, output_file):

    with open(trans_input_file,'r') as f:
        trans_dic = read_trans(f)
    with open(thermo_input_file,'r') as f:
        thermo_dic = read_thermo(f)

    sutherland_trans_dic = {}
    for species, trans in trans_dic.items():
        T_arr = np.arange(300,5000,100)
        try:
            elements = thermo_dic[species]
        except:
            raise ValueError('species {:s} undefined in thermo file'.format(species))

        W = molecularWeight(elements)#kg/kmol
        mu_arr = np.array([mu(trans,W,T) for T in T_arr])

        params ,_ = curve_fit(sutherland,T_arr, mu_arr)
        mu_arr_fit = sutherland(T_arr, *params)
        error = np.mean((mu_arr_fit-mu_arr)**2)**0.5
        print('species name:{:>10s}    fitting error:    {:.4E}'.format(species, error))

        sutherland_trans_dic[species] = params

    with open(output_file,'w') as f:
        for species, params in sutherland_trans_dic.items():
            f.write(species+'\n')
            f.write('{\n')
            f.write('    transport\n')
            f.write('    {\n')
            f.write('        '+'As'.ljust(16)+str(params[0])+';\n')
            f.write('        '+'Ts'.ljust(16)+str(params[1])+';\n')
            f.write('    }\n')
            f.write('}\n\n')




