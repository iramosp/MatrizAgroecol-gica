import copy, random, math
import numpy as np
from scipy.integrate import odeint


def matriz_agroecologica(poblacion_0, tipo, t_total, n_especies, f, r_alea, a_alea, m_milpa, m_intensivo, D, vecinos):
    """matriz_agroecologica iteraciones migración y muerte y luego lotka:  
    poblacion = [tiempo] [x][y] [especieA][especieB][...]

    """
    #corre simulacion en 2D
    poblacion = [poblacion_0] #inicializa array poblacion de 3 dimensiones con forma (x,y,n_especies)
    for t in range(t_total):
        temp = np.zeros_like(poblacion[-1])        
        for i in range(iter_difymuerte):
            for i in range(x_celdas): #para todo x y
                for j in range(y_celdas):
                    if tipo[i][j] == 'm': #milpa
                        temp[i][j] = muerte(poblacion[-1][i][j], m_milpa)
                    elif tipo[i][j] == 'i': #intensivo
                        temp[i][j] = muerte(poblacion[-1][i][j], m_intensivo)
                    elif tipo[i][j] == 'b':
                        temp[i][j] = poblacion[-1][i][j]
            #migracion
            #for i in range(n_especies):
                #si quieres variar la taza de migracion por especie aqui es donde debes de variarla            
            temp = migracion(temp, tipo, Disp)
            poblacion.append(temp)
        #interacciones ecologicas y muerte
        for i in range(x_celdas): #para todo x y
            for j in range(y_celdas):
                if tipo[i][j] == 'b': #interacciones ecologicas
                    temp[i][j] = odeint(f, poblacion[-1][i][j], [0,1], args=(r_alea,a_alea))[-1]
        poblacion.append(temp)
    return np.array(poblacion)  


def migracion(X, tp, L):
    """
    Funcion que asigna cuánta poblacion migra
    dependiendo del parche en el que esté.
    Los valores de migración están guardados en el 
    diccionario L
    L['b'] = valor bosque
    L['m'] = valor milpa 
    L['i'] = valor intensivo
    X es la matriz de dimensión 3 X(x,y,i) donde 
    la especie i-esima tiene su representación para 
    todos los parches
    t es la distribución de tipos de parche en todo el 
    espacio
    """
    s = X.shape
    t = np.array(tp)
    P = np.zeros(s, dtype=float)
    G = np.zeros(s, dtype=float)
    R = np.zeros(s, dtype=float) 

    for idx in range(s[2]):
        esp = X[:,:,idx]
        dm = esp.shape
        xesp, yesp = dm

        loss_e = np.zeros( dm, dtype=float ) 
        gain_e = np.zeros ( dm, dtype=float )       
        for x in range(xesp):
            for y in range(yesp):
                loss_e[x,y] = esp[x,y] * L[t[x,y]]
        for x in range(-1,dm[0]):
            for y in range(-1,dm[1]):
                gain_e[x,y] = (loss_e[(x-1)%dm[0],(y-1)%dm[1]]+loss_e[(x-1)%dm[0],y%dm[1]]+loss_e[(x-1)%dm[0],(y+1)%dm[1]]+loss_e[x%dm[0],(y-1)%dm[1]]+loss_e[x%dm[0],(y+1)%dm[1]]+loss_e[(x+1)%dm[0],(y-1)%dm[1]]+loss_e[(x+1)%dm[0],y%dm[1]]+loss_e[(x+1)%dm[0],(y+1)%dm[1]])/8
        P[:,:,idx] = loss_e
        G[:,:,idx] = gain_e
        R[:,:,idx] = esp +( gain_e - loss_e )
      
    return R   


def lotka(x,t,r_alea, a_alea):
    """ecuacion de lotka volterra generalizada
    """
    dx = x*(r_alea+np.dot(a_alea,x))
    return dx 


def genera_poblacion_inicial(tipo_matriz_agroecologica, n_especies, p0_bosque=0, p0_milpa=0, p0_intensivo=0): #Poblacion inicial
    if type(p0_bosque)==float: #all same value
        p0_bosque = [p0_bosque for i in range(n_especies)]
    elif type(p0_bosque)==int: #users
        p0_bosque = [p0_bosque for i in range(n_especies)]
    elif p0_bosque=="eq_caos":p0_bosque = [ 0.3013,  0.4586,  0.1307,  0.3557]
    
    if type(p0_milpa)==float: #all same value
        p0_milpa = [p0_milpa for i in range(n_especies)]
    elif type(p0_milpa)==int: #users
        p0_milpa = [p0_milpa for i in range(n_especies)]
        
    if type(p0_intensivo)==float: #all same value
        p0_intensivo = [p0_intensivo for i in range(n_especies)]
    elif type(p0_intensivo)==int: #users
        p0_intensivo = [p0_intensivo for i in range(n_especies)]
    
    poblacion_0 = copy.deepcopy(tipo)
    for x in range(len(tipo)): #inicializar poblaciones
        for y in range(len(tipo[0])):
            if tipo[x][y] == 'b': 
                if p0_bosque=="random": poblacion_0[x][y] = [1/random.random() for i in range(n_especies)]
                else: poblacion_0[x][y] = p0_bosque
            if tipo[x][y] == 'm': 
                if p0_milpa=="random": poblacion_0[x][y] = [random.random() for i in range(n_especies)]
                else: poblacion_0[x][y] = p0_milpa
            if tipo[x][y] == 'i': 
                if p0_intensivo=="random": poblacion_0[x][y] = [random.random() for i in range(n_especies)]              
                else: poblacion_0[x][y] = p0_intensivo
                
    #poblacion_0[1][1] = [0,0,0,0,0,0,0,0,0,0] #para vaciar los 4 bosques de las esquinas 
    #poblacion_0[1][8] = [0,0,0,0,0,0,0,0,0,0]
    #poblacion_0[8][1] = [0,0,0,0,0,0,0,0,0,0]
    #poblacion_0[8][8] = [0,0,0,0,0,0,0,0,0,0]    
                
    return np.array(poblacion_0) #array de 3 dimensiones con forma (x,y,n_especies)


def muerte(x, m):
   """recibe x = poblacion   
        m = taza muerte cte o np.array
   regresa x = poblacion superviviente  np.array
   """
   x = x - x*m
   return x