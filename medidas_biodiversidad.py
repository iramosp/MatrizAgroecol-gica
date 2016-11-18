import numpy as np

def riqueza_paisaje(poblacion, t=-1, biomasa_min = 5.):
    """ Entrada: un arreglo poblacion = [tiempo] [x][y] [especieA][especieB][...]
        Salida: la biomasa y la riqueza de especies en UN tiempo, en todo el paisaje;
        si no se especifica el tiempo se toma la última iteración.
    """
    riqueza = []
    for idx in range(poblacion.shape[3]):
        if np.sum(poblacion[t,:,:,idx]) > biomasa_min: #suma de cada especie en la iteración t
            riqueza.append(1) 

    biomasa = np.sum(poblacion[t,:,:,:]) #suma total de individuos al final de n iteraciones
    return biomasa, np.sum(riqueza)  #riqueza de especies al final de n iteraciones


def medidas_por_celda(poblacion, t=-1, biomasa_min = 5.):
    """
    Entra un arreglo 4D con la forma: 
    poblacion = [tiempo] [x] [y] [especies]
    Produce arreglos 2D que indican la biomasa y la riqueza 
    de cada celda de poblacion en UN tiempo; 
    si no se indica el tiempo, se toma la última iteración
    """
    poblacion = poblacion[t]
    biomasa = np.zeros((poblacion.shape[0], poblacion.shape[1]))
    riqueza = np.zeros((poblacion.shape[0], poblacion.shape[1]))

    for x in range(poblacion.shape[0]):
        for y in range(poblacion.shape[1]):
            biomasa[x] [y] = np.sum(poblacion[x, y, :])
            riqueza[x] [y] = np.sum(poblacion[x, y, :] > biomasa_min)
              
    return biomasa, riqueza


def riqueza_agricola(poblacion, paisaje, t=-1, biomasa_min = 5.):
    """ Entrada: un arreglo poblacion = [tiempo] [x][y] [especieA][especieB][...],
        el paisaje y el tiempo (si no se indica el tiempo se toma la última iteración).
        Salida: la biomasa y la riqueza de especies en UN tiempo, en las celdas que no son bosque;
        si no se especifica el tiempo se toma la última iteración.
        Adaptada del programa original.
    """
    x_celdas = len(paisaje)
    y_celdas = len(paisaje[1])
    
    riqueza = np.zeros(poblacion.shape[3])
    for idx in range(poblacion.shape[3]):
        
        for i in range(x_celdas): #para todo x y
            for j in range(y_celdas):
                if paisaje[i][j] != "b":
                    riqueza[idx] += poblacion[t,i,j,idx]  
    
    biomasa = np.sum(riqueza) #suma total de individuos al final de n iteraciones
    return biomasa, len(riqueza[riqueza > biomasa_min])  #riqueza de especies al final de n iteraciones



def shannon_wiener(poblacion, paisaje, t=-1, biomasa_min = 5.):
    """ Entrada: un arreglo poblacion = [tiempo] [x][y] [especieA][especieB][...],
        el paisaje y el tiempo (si no se indica el tiempo se toma la última iteración).
        Salida: la biomasa y la riqueza de especies en UN tiempo, en las celdas que no son bosque;
        si no se especifica el tiempo se toma la última iteración.
        Adaptada del programa original.
    """
    x_celdas = len(paisaje)
    y_celdas = len(paisaje[1])
    
    riqueza = np.zeros(poblacion.shape[3])
    for idx in range(poblacion.shape[3]):
        
        for i in range(x_celdas): #para todo x y
            for j in range(y_celdas):
                if paisaje[i][j] != "b":
                    riqueza[idx] += poblacion[t,i,j,idx]  
    
    biomasa = np.sum(riqueza) #suma total de individuos al final de n iteraciones - biomasa
    riqueza = riqueza[riqueza > biomasa_min]
    p = riqueza / biomasa
    sw = p * np.log(p)
    sw = -1 * np.sum(sw)

    return sw

def medida_area(poblacion, paisaje, t=-1, biomasa_min = 5.):
    """ Entrada: un arreglo poblacion = [tiempo] [x][y] [especieA][especieB][...],
        el paisaje y el tiempo (si no se indica el tiempo se toma la última iteración).
        Salida: la biomasa y la riqueza de especies en UN tiempo, en las celdas que no son bosque;
        si no se especifica el tiempo se toma la última iteración.
        Adaptada del programa original.
        Vivos si están vivos con más de 0.05 en al menos 10% del paisaje + biomasa min
    """
    x_celdas = len(paisaje)
    y_celdas = len(paisaje[1])
    
    riqueza = np.zeros(poblacion.shape[3])
    area = np.zeros(poblacion.shape[3])
    vivos = np.zeros(poblacion.shape[3])


    for idx in range(poblacion.shape[3]):
        
        for i in range(x_celdas): #para todo x y
            for j in range(y_celdas):
                if paisaje[i][j] != "b":
                    riqueza[idx] += poblacion[t,i,j,idx]  

                    if poblacion[t,i,j,idx] >= biomasa_min/50:
                        area[idx] += 1

        if riqueza[idx] > biomasa_min and area[idx]>=30:
            vivos[idx] = 1
    return len(vivos[vivos>0])
