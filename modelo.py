import copy, random, math
import numpy as np
from scipy.integrate import odeint


def matriz_agroecologica(paisaje, matriz_interacciones, tasas_reproduccion, condiciones_iniciales, t_total, Dispersion, Mortalidad, pasos_mm=2):
    """Modelando la matriz agroecológica. Los pasos que sigue el modelo son: dinámica de poblaciones Lotka-Volterra -> (migración -> muerte) x pasos_mm. La población se guarda en cada iteración de LV y en cada iteración de migración+muerte. Regresa un arreglo de numpy con la poblacion en cada paso, con forma:
    poblacion = [(pasos_mm +1 ) * t_total + 1, x_celdas, y_celdas, n_especies]
    """

    x_celdas = len(paisaje)
    y_celdas = len(paisaje[1])
    n_especies = len(matriz_interacciones)

    # Inicializar un arreglo de numpy de la forma: poblacion = [tiempo] [x][y] [especieA][especieB][...]
    poblacion = np.zeros(((pasos_mm +1 ) * t_total + 1, x_celdas, y_celdas, n_especies))

    # Población inicial con las condiciones iniciales en las celdas de bosque   
    poblacion[0,:,:,:] = genera_poblacion_inicial(paisaje, n_especies, p0_bosque = condiciones_iniciales)

    for t in range(1, t_total+1):
        T = (pasos_mm + 1) * t - pasos_mm

        poblacion[T, :, :, :] = copy.deepcopy(poblacion[T-1, :, :, :])

        # Dinámica de poblaciones en bosque
        for i in range(x_celdas): #para todo x y
            for j in range(y_celdas):
                if paisaje[i][j] == "b":
                    poblacion[T, i, j, :] = odeint(lotka, poblacion[T, i, j, :], [0, 0.5], 
                                            args=(tasas_reproduccion, matriz_interacciones))[-1]
    
        for k in range(1, pasos_mm + 1):
            # Migración
            poblacion[T + k, :, :, :] = migracion(poblacion[T + k -1, :, :, :] , paisaje, Dispersion)

            # Muerte
            for i in range(x_celdas): #para todo x y
                for j in range(y_celdas):
                    tasa_muerte = Mortalidad[paisaje[i][j]]
                    poblacion[T + k, i, j, :] = (1-tasa_muerte) * poblacion[T + k, i, j, :]

    return poblacion

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


def lotka(x, t, r, a):
    """Ecuacion de lotka volterra generalizada
    """
    dx = x * (r + np.dot(a, x))
    return dx 


def genera_poblacion_inicial(paisaje, n_especies, p0_bosque, p0_milpa=0, p0_intensivo=0):
    """
    Genera la poblacion inicial. Se especifican condiciones iniciales para cada especie, para las celdas de "bosque". Las celdas "milpa" e "intensivo" inician vacías.
    Recibe: 
        paisaje: lista de x por y, especifica el tipo de celda ("b" = bosque, "m" = milpa, "i" = intensivo)
        n_especies: int 
        p0_bosque: condiciones iniciales para las celdas "bosque". Es numpy array de 1D con forma (n_especies, ).

    Regresa:
        poblacion_0: numpy array de 3D con forma (x, y, n_especies).
    """

    p0_milpa = [p0_milpa for i in range(n_especies)]
    p0_intensivo = [p0_intensivo for i in range(n_especies)]
    
    poblacion_0 = copy.deepcopy(paisaje)
    for x in range(len(paisaje)): #inicializar poblaciones
        for y in range(len(paisaje[0])):
            if paisaje[x][y] == 'b': 
                poblacion_0[x][y] = p0_bosque

            if paisaje[x][y] == 'm': 
                poblacion_0[x][y] = p0_milpa
           
            if paisaje[x][y] == 'i': 
                poblacion_0[x][y] = p0_intensivo
                
    return np.array(poblacion_0)