#!/usr/bin/env python
# -*- coding: utf-8 -*- 

# 30 agosto 2015


# poblacion = [tiempo] [x][y] [especieA][especieB][...]



import sys, copy, random, os, glob, math
import numpy as np
import time
from itertools import combinations, product
from scipy.integrate import odeint
from scipy.signal import convolve2d
from matplotlib import pyplot as plt
from matplotlib import colorbar as cbar
from matplotlib import colors as colors


"""
CORRIDAS
"""

def guardarexperimento(momento, m1, v1, ci):
    """
    Este metodo guarda tres archivos con el nombre
    momento-{m, v, c}.txt
    usando el metodo savetxt de numpy. 
    
    Para cargarlos hay que usar el metodo, también de numpy, loadtxt
    M = np.loadtxt('ejemplo.txt')
    entonces, en la variable M contiene lo leído
    """
    pref = '-'.join(momento.split())
    print 'Abriendo archivo para guardar datos'
    narch = pref+'-m.txt'
    fh = open( narch, 'w' )
    fh.write( '#Matriz \n')
    np.savetxt(fh, m1)
    fh.close()
    narch = pref+'-v.txt'
    fh = open( narch, 'w' )
    fh.write( '#Tazas de reproduccion \n' )
    np.savetxt(fh, v1)
    fh.close()
    narch = pref+'-c.txt'
    fh = open( narch, 'w' )
    fh.write( '#Condiciones iniciales\n' )
    np.savetxt(fh, ci)
    fh.close()
    print 'Datos guardados en ', narch
    
def migracion_esp(X, tp, L):
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
    #print "Especies "
    for idx in range(s[2]):
        esp = X[:,:,idx]
        dm = esp.shape
        xesp, yesp = dm
        #print "Espacio de la especie ", idx, ": ", xesp,":",yesp
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
    
    
def correr_aleatorizado_sin_matriz( n_especies, t_total, f, r_alea, a_alea, plot=False): #integra sistema con la ecuación lotka-volterra
    global x_0
#    x_0 = 1/(np.random.random(n_especies)) #pob inicial aleatoria mayor a 1 EL NORMAL
    x_0 = Nombres[actual][0] #iniciando con condiciones guardadas
#    x_0 = puntofijo #inciando desde los puntos fijos
    t = np.linspace(0,t_total,t_total+1)
    x = odeint(f, x_0, t, args=(r_alea,a_alea))
    if plot:
        fig = plt.figure()
        fig.add_subplot(111)
        plt.plot(t, x)
        plt.show()  
    return x[-1]
    
#correr_2DMM iteraciones migración y muerte y luego lotka:  
def correr_2DMM(poblacion_0, tipo, t_total, n_especies, f, r_alea, a_alea, m_milpa, m_intensivo, D, vecinos):
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
            temp = migracion_esp(temp, tipo, Disp)
            poblacion.append(temp)
        #interacciones ecologicas y muerte
        for i in range(x_celdas): #para todo x y
            for j in range(y_celdas):
                if tipo[i][j] == 'b': #interacciones ecologicas
                    temp[i][j] = odeint(f, poblacion[-1][i][j], [0,1], args=(r_alea,a_alea))[-1]
        poblacion.append(temp)
    return np.array(poblacion)  

#correr_2DL lotka y luego iteraciones migración y muerte:
def correr_2DL(poblacion_0, tipo, t_total, n_especies, f, r_alea, a_alea, m_milpa, m_intensivo, D, vecinos):
    #corre simulacion en 2D
    poblacion = [poblacion_0] #inicializa array poblacion de 3 dimensiones con forma (x,y,n_especies)
    for t in range(t_total):
        temp = np.zeros_like(poblacion[-1])
        #interacciones ecologicas y muerte
        for i in range(x_celdas): #para todo x y
            for j in range(y_celdas):
                if tipo[i][j] == 'b': #interacciones ecologicas
                    temp[i][j] = odeint(f, poblacion[-1][i][j], [0,1], args=(r_alea,a_alea))[-1]
        for i in range(iter_difymuerte):
            for i in range(x_celdas): #para todo x y
                for j in range(y_celdas):
                    if tipo[i][j] == 'm': #milpa
                        temp[i][j] = muerte(poblacion[-1][i][j], m_milpa)
                    elif tipo[i][j] == 'i': #intensivo
                        temp[i][j] = muerte(poblacion[-1][i][j], m_intensivo)
            #migracion
            #for i in range(n_especies):
                #si quieres variar la taza de migracion por especie aqui es donde debes de variarla
        
            temp = migracion_esp(temp, tipo, Disp)
            poblacion.append(temp)
    return np.array(poblacion)

#correr_2D un paso de lotka por un paso de migración y muerte:
def correr_2D(poblacion_0, tipo, t_total, n_especies, f, r_alea, a_alea, m_milpa, m_intensivo, D, vecinos):
    #corre simulacion en 2D
    poblacion = [poblacion_0] #inicializa array poblacion de 3 dimensiones con forma (x,y,n_especies)
    for t in range(t_total):
        temp = np.zeros_like(poblacion[-1])
        #interacciones ecologicas y muerte
        for i in range(x_celdas): #para todo x y
            for j in range(y_celdas):
                if tipo[i][j] == 'b': #interacciones ecologicas
                    temp[i][j] = odeint(f, poblacion[-1][i][j], [0,1], args=(r_alea,a_alea))[-1]
                elif tipo[i][j] == 'm': #milpa
                    temp[i][j] = muerte(poblacion[-1][i][j], m_milpa)
                elif tipo[i][j] == 'i': #intensivo
                    temp[i][j] = muerte(poblacion[-1][i][j], m_intensivo)
        #migracion
         
        #for i in range(n_especies):
            #si quieres variar la taza de migracion por especie aqui es donde debes de variarla
         
        temp = migracion_esp(temp, tipo, Disp)
        poblacion.append(temp)
    return np.array(poblacion)       

"""
VALORES INICIALES
"""

def genera_tipo_matriz_agroecologica(x_celdas, y_celdas, n_bosque=0, posicion_bosque=[], n_milpa=0, posicion_milpa=[]): #Distribucion de matriz agroecologica
    #recibe una lista de tuples o un comando
    #ej: [(2,3),(4,1),(2,2)]
    tipo = [['i' for i in range(y_celdas)] for j in range(x_celdas)] #inicializa todo con intensivo
    
    if posicion_bosque == "extremos":
        posicion_bosque = [(0,0),(x_celdas-1, y_celdas-1)]
    if posicion_milpa == "stepping": #method for square matrix
        posicion_milpa = [(i,i) for i in range(n_milpa,x_celdas-1,n_milpa)]
    
    if posicion_bosque == "random" or posicion_milpa == "random":
        #genera todas las posibles coordenadas
        pairs = [(x,y) for x in [i for i in range(x_celdas)] for y in [j for j in range(y_celdas)]]
        random.shuffle(pairs) #randomisa
        #quita los bosque y milpas ya declarados
        if type(posicion_bosque) == list:
            for p in posicion_bosque:
                try: pairs.remove(p)
                except: pass
        if type(posicion_milpa) == list:
            for p in posicion_milpa:
                try: pairs.remove(p)
                except: pass
        #selecciona n_bosques y n_milpas
        if posicion_bosque == "random":
            posicion_bosque = random.sample(pairs,n_bosque)
            for p in posicion_bosque: #quita para no confundir a random milpa
                try: pairs.remove(p)
                except: pass
        if posicion_milpa == "random":
            posicion_milpa = random.sample(pairs,n_milpa)
    
    for p in posicion_bosque: #escribe bosques en matriz
        tipo[p[0]][p[1]] = 'b'
    for p in posicion_milpa: #escribe bosques en matriz
        tipo[p[0]][p[1]] = 'm'
    return tipo

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

def matriznicho(s_especies=10, c_conectancia= 0.2): #hace la matriz de interacciones con 0, 1 y -1s
    
    f_matrizbuena = False
    while not f_matrizbuena:        
        ni = [random.random() for _ in xrange(s_especies)]
        array_ni = np.array(ni)
#        print array_ni        
        ci = [np.random.random()*x for x in ni]
        array_ci = np.array(ci)
#        print array_ci        
        ripre = np.random.beta(1,(1.0/(2*c_conectancia))-1,s_especies) #c_conectancia <= .5 
        array_ripre = np.array(ripre)
#        print array_ripre        
        ri = array_ripre * array_ni
#        print ri        
        rimedios = ri/2
#        print rimedios        
        extizq = array_ci - rimedios
#        print extizq        
        extder = array_ci + rimedios
#       print extder        
        redtrof = np.zeros(shape=(s_especies,s_especies)) #inicializa matriz de especies vs especies en ceros        
        for i in range(len(array_ni)): 
            for j in range(len(extizq)):
                if array_ni[i] >= extizq[j] and array_ni[i] <= extder[j]: 
                    redtrof[i,j] = -1 
                    redtrof[j,i] = 1 #aij=-1 en la matriz trófica si j se come a i, aji=-1
                else: pass                
        np.fill_diagonal(redtrof, -1)        #pone diagonal negativa         
#        print redtrof #devuelve la matriz trófica         
        #comparar filas y columnas, devuelve los pares idénticos        
        f_identica = False        
        for i in range(len(redtrof)): #generate pairs
            for j in range(i+1,len(redtrof)): 
                if np.array_equal(redtrof[i],redtrof[j]): #compare rows
                    if np.array_equal(redtrof[:,i],redtrof[:,j]): #compare columns
#                        print (i, j)
                        f_identica = True
                else: pass                
        #devuelve especies disconexas        
        f_disconexa = False        
        for i in range(len(redtrof)): 
                if sum(np.absolute(redtrof[i])) <= 1:
                    if sum(np.absolute(redtrof[:,i])) <= 1:
#                        print i,
                        f_disconexa = True
#                else: print "cool"                
        #calcula conectancia        
        f_conectancia = False        
        links = (np.sum(np.absolute(redtrof))/2)-s_especies
#        print "links existentes", links
        conectanciae = links/s_especies**2
#        print "conectancia empirica", conectanciae        
        if np.absolute(conectanciae - c_conectancia) > .03:
            f_conectancia = True            
        if not f_identica and not f_disconexa and not f_conectancia:
#            print redtrof
            return redtrof
            f_matrizbuena=True

def a_aleatorizada(n_especies, a):   #toma la matriz trofica "a" creada con el modelo de nicho y aleatoriza sus magnitudes, conservando el signo
    global amodificada
    amodificada = np.zeros([n_especies, n_especies]) #para hacer interacciones negativas menos fuertes cambiando -1 por algo negativo mayor
    for i in range(len(a)):
        for j in range(len(a)):
            if a[i,j]==-1:
                amodificada[i,j]=-0.1   #valor por el cual se sustituyen los -1
            else:
                amodificada[i,j]=a[i,j]  #lo demás se queda igual
    M = np.zeros([n_especies, n_especies])
    for i in range(len(M)):
        for j in range (len(M)):
            M[i,j] = np.random.uniform(0,2) #aleatorizados de [0,2]
#            M[i,j] = np.absolute(np.random.normal())   #aleatorizados con el valor absoluto de un aleatorio con dist normal       
    matrizaleatorizada = M*amodificada #a_alea con interacc negativas menos fuertes
#    matrizaleatorizada = M*a   #a_alea con interacciones negativas y positivas aleatorias en el mismo intervalo
    return matrizaleatorizada
    
def nicho_taza_crecimiento(a): #establecer taza crecimiento dependiendo de si la especie es basal (+1) o no (-1)
    r = np.empty(len(a))
    matrizr = a*-1  #ahora nos interesan los 1s en las columnas
    r.fill(-1)    #inicializa array de rs con -1s
    matrizr[matrizr<=0] = 0  # hace 0 todo lo que no son unos
    #print matrizr
    for i in range (len(matrizr)):
    #print sum(matrizr[:,i])
        if sum(matrizr[:,i])==1:    #suma los unos en cada columna
            r[i]=1                  #si sólo hay un uno, la especie es basal, su r se hace 1
    return r        
    
def r_aleatorizado(r):   #toma el vector de r y aleatoriza sus magnitudes, conservando el signo
    raleatorio = np.empty(len(r))
    for i in range (len(r)):
        if r[i]==1:
            raleatorio[i]=r[i]*np.random.uniform(0,10) #r positivas ente 0 y 2.5
        else:
            raleatorio[i]=r[i]*np.random.random() #r negativas entre -1 y 0
#    raleatorio = [np.random.uniform(0,2.5)*x for x in r]
#    raleatorio = [np.random.normal(2.5)*x for x in r] #aleatorizado con dist normal con centro no cero
    return raleatorio
      

def vector_identidades(n_especies, a): #hace un vector que indica si cada especie es basal, intermedia o top
    identidades = np.zeros(n_especies) #inicia vector con ceros, terminará 0=intermedia, 1=basal, 2=top
    #Basales (sólo un -1 en la matriz de nicho original)    
    matrizbasales = a*-1  #ahora sólo no simportan los 1    
    matrizbasales[matrizbasales<=0] = 0  
    # print matrizr    
    for i in range (len(matrizbasales)):
    #    print sum(matrizr[:,i])
        if sum(matrizbasales[:,i])==1: #si son basales (solo tienen un uno)
            #print "la especie", i, "es basal"
            identidades[i]=1
        else: pass    
    #Top (al menos dos -1 y ningún 1 en la original)    
    matriztop = a*1
    matriztop[matriztop >= 0] = 0        
    for i in range (len(a)):
        if sum(a[:,i]) == sum(matriztop[:,i]):
            #print "la especie", i, "es top"
            identidades[i]=2
        else: pass    
    return identidades   


"""
FUNCIONES
"""
   
def d_lotkavolterra_alea(x,t,r_alea, a_alea): #ecuacion de lotka volterra generalizada
    dx = x*(r_alea+np.dot(a_alea,x))
    return dx 

def migracion(x, D, tipo='linea',limite='fill'):
    # recibe matriz de celdas con poblaciones y coeficiente difusion
    # regresa nueva matriz de celdas con poblaciones  np.array
    # D es una lista de valores tal que cada entrda indica lo siguiente
    # D[0] = el tipo de la celda actual, es decir, puede valer 'm','b','i'
    # D[1] = es la tasa de difusion para ese mismo tipo de celda
    # D[2] = es la tasa de difusion para cualquier otro (bosque)
    # tipo determina el numero de vecinos
    #     linea, vecinos4, vecinos8, y especial
    # limite es cerrado o circular
    #     fill o wrap
    
    if tipo == 'linea': # Migracion linea 1D
        if limite == 'fill': #cerrado
            #calcular perdida
            loss = np.array([1.] + [2. for i in range(len(x[0])-2)] + [1.])
            loss = (D/2.0) * loss * x[0]
            #calcular ganancia
            D = np.array([D/2.0,0,D/2.0])
            gain = np.convolve(x[0], D, 'same')
        if limite == 'wrap': #circular
            #calcular perdida
            loss = D * x[0]
            #calcular ganancia
            D = np.array([D/2.0,0,D/2.0])
            gain = np.convolve(x, D, 'full') #np.convolve(x, D, 'full')
            gain[1]+=gain[-1] #wrap
            gain[-2]+=gain[0] #wrap
            gain = gain[1:-1] #reduce
        return  [x[0] + gain - loss]
   
    # Migracion matriz 2D
    elif tipo == 'vecinos4':
        #calcular perdida
        if limite == 'fill': #cerrado
            loss = np.array([[2] + [3 for i in range(len(x[0])-2)] + [2]] + [[3] + [4 for i in range(len(x[0])-2)] + [3] for j in range(len(x)-2)] + [[2] + [3 for i in range(len(x[0])-2)] + [2]])
            loss = (D/4.0) * loss * x
        if limite == 'wrap': #circular
            loss = D * x
        print loss
        print sum(sum(loss))
        #calcular ganancia
        D = np.array([[0,D/4.0,0],[D/4.0,0,D/4.0],[0,D/4.0,0],]) #matriz dispercion
        gain = convolve2d(x, D, mode='same', boundary=limite)
        print gain
        print sum(sum(gain))
        #calcula total
        return  x + gain - loss
   
    elif tipo == 'vecinos8':
        ca = D[0]
        tda = D[1]
        tdb = D[2]
        #calcular perdida
        if limite == 'fill': #cerrado, este mundo no lo vamos a usar porque son demasiadas interpretaciones
            loss = np.array([[3] + [5 for i in range(len(x[0])-2)] + [3]] + [[5] + [8 for i in range(len(x[0])-2)] + [5] for j in range(len(x)-2)] + [[3] + [5 for i in range(len(x[0])-2)] + [3]])
            loss = (D/8.0) * loss * x
            #if ca=='i' or ca=='m':
            #    loss = (tda/8.0) * loss * x
            #elif ca=='b':
            #    loss = (tdb/8.0) * loss * x
        if limite == 'wrap': #circular
            loss = tda*x if ca=='m' or ca=='i' else tdb*x
        #calcular ganancia
        D = np.array([[D/8.0,D/8.0,D/8.0],[D/8.0,0,D/8.0],[D/8.0,D/8.0,D/8.0],]) #matriz dispercion
        gain = convolve2d(x, D, mode='same', boundary=limite)
        #calcula total
        return x + gain - loss
   
    # Si se pasa directamente un array de dispercion
    elif tipo == 'especial':
        #calcular perdida
        loss = sum(sum(D)) * x
        #calcular ganancia
        gain = convolve2d(x,D,mode='same', boundary='wrap')
        #calcula total
        return x + gain - loss

def muerte(x, m):
   # recibe x = poblacion   
   #     m = taza muerte cte o np.array
   # regresa x = poblacion superviviente  np.array
   x = x - x*m
   return x



"""
GRAFICAS
"""
def plot_xy(name, x, y, title='',labels=[], x_label='', y_label=''): #grafica x contra y
    #recibe x como numpy.array
    #Ejemplo imprimir poblacion vs tiempo
    #plot_xy(name, x, t, title,labels, x_label='Tiempo', y_label='Poblacion')
    plt.clf()
    try:
        for i in range(len(x[0])): #grafica cada x contra y
        #print i, labels[i]
            plt.plot(y, x[:,i], label=labels[i])
    except: plt.plot(y, x)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(title)
    plt.legend(loc=0)
    plt.plot()
    #save plot
    f_format = name.split('.')[-1]
    name = name.split('.')[0]
    plt.savefig(name+'.'+f_format, format=f_format, bbox_inches='tight')
    plt.show()

def plotHeatmap(name, data , x_label='' , y_label='', x_tick_labels=[], y_tick_labels=[]):
    #data is a 2x2 array normalized [0,1]
    plt.clf()
    fig, ax = plt.subplots()
    #delete top and right axis
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    ## put the major ticks at the middle of each cell
    ax.set_xticks(np.arange(data.shape[1])+0.5, minor=False)
    ax.set_yticks(np.arange(data.shape[0])+0.5, minor=False)
    ## want a more natural, table-like display
    ##ax.invert_yaxis()
    ##ax.xaxis.tick_top()
    ax.set_xticklabels(x_tick_labels, rotation=90, minor=False)
    ax.set_yticklabels(y_tick_labels, minor=False)
    #set colorbar
    cdict = {'red':   [(0.0,  1.0, 1.0),(0.01,  0.5, 0.5),(0.5,  0.0, 0.0),(1.0,  0.0, 0.0)],
        'green': [(0.0,  1.0, 1.0),(0.1, 1.0, 1.0),(1.0,  0.0, 0.0)],
        'blue':  [(0.0,  1.0, 1.0),(0.5,  1.0, 1.0),(1.0,  0.5, 0.5)]}
    my_cmap=colors.LinearSegmentedColormap('my_colormap',cdict,256)
    #heatmap = ax.pcolor(data, cmap=plt.cm.Blues)
    heatmap = ax.pcolor(data, cmap=my_cmap, vmin=0, vmax=1)
    cbar = plt.colorbar(heatmap)
    plt.title(name)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.plot()
    #save plot
    f_format = name.split('.')[-1]
    name = name.split('.')[0]
    plt.savefig(name+'.'+f_format, format=f_format, bbox_inches='tight')
    #plt.show()

def plot_lineas_especies(name, x, t_total, labels=[]):
    #graficar sum xy para todas las especies
    #print np.shape(labels)
    t = np.linspace(0,t_total+1,t_total+1)
    plot_xy(name, x, t, 'especies_vs_tiempo', labels)

def plot_diagramas_fase(name, x, t_total, labels=[]):
    t = np.linspace(0,t_total+1,t_total+1)
    for i,j in combinations ([n for n in range(len(x[0]))],2):
        plot_xy(name.split('.')[0] + '_'+labels[i]+labels[j]+'.' + name.split('.')[-1], x[:,i], x[:,j], 'fase_'+labels[i]+labels[j])
    
def heatmaps_especie_en_tiempo(name,poblacion, n=[], labels=[]):
    ter = '.'+name.split('.')[-1]
    name = name.split('.')[0]
    #graficar heatmap especie en un tiempo dado
    #n debe ser un array
    if type(n) == int: n = [n]
    for i in range(n_especies):
        plotHeatmap(name+label[i]+'_t0'+ter, poblacion[0,:,:,i]) #estado inicial
        plotHeatmap(name+label[i]+'_tf'+ter, poblacion[-1,:,:,i]) #estado final
        for j in n: #grafica vector de tiempos
            plotHeatmap(name+label+'_t'+str(j)+ter, poblacion[j,:,:,i])

def heatmap_tipo(name,tipo,m_milpa=.5,m_intensivo=1):
    for i in range(len(tipo)):
        for j in range(len(tipo[0])):
            if tipo[i][j] == 'b': tipo[i][j] = 1.
            if tipo[i][j] == 'm': tipo[i][j] = 1-m_milpa
            if tipo[i][j] == 'i': tipo[i][j] = 1-m_intensivo
    tipo = np.array(tipo)
    plotHeatmap(name, tipo) #estado inicial    
    
    

"""
MEDIDAS DE BIODIVERSIDAD
"""

def shannonwiener(x): #índice de shannon para el sistema sin espacio o para una sola celda
    #recibe poblacion x que es un vector con 10 entradas
    tot=np.sum(x)
    pi=[]
    for i in range(len(x)):
        p_i=x[i]/tot
        if p_i>0:
            pi.append(p_i)
    lnpi=np.array(np.log(pi))
    mult=pi*lnpi
    suma=np.sum(mult)
    return suma*-1

def nefectivo(sh): #num efectivo para el sistema sin espacio
    #recibe índice de shannon
    return math.exp(sh)

def alpha_shannon(pob): #índice de shannon para el sistema completo
    #recibe poblacion[-1] con forma [x=10][y=10][sp=10]
    #regresa float
    #H' = -\sum(p_i * ln(p_i))
    sp = np.sum(np.sum(pob, axis=0),axis=0) #sumar valores en xy. valor total de cada especie
    tot = np.sum(sp) #total sp
    pi=[]
    for i in range(len(sp)):
        p_i=sp[i]/tot
        if p_i>0:
            pi.append(p_i)       
    lnpi=np.array(np.log(pi))
    mult=pi*lnpi
    suma=np.sum(mult)
    return suma*-1
        




"""
M   M     A  III  N   N
MM MM    AA   I   NN  N
M M M   A A   I   N N N
M   M  AAAA   I   N  NN
M   M A   A  III  N   N
"""


"""
PARAMETROS
"""

n_especies = 10 # numero especies
x_celdas = 10 # numero celdas en x
    #para 1D x=1
y_celdas = 10 # numero celdas en y
              #        si 1 modelo en 1D

a = matriznicho(n_especies)
a_alea = a_aleatorizada(n_especies, a)
r = nicho_taza_crecimiento(a)
r_alea = r_aleatorizado(r)

numerodemilpas = 0

#iden = vector_identidades(n_especies, a)
#puntofijo = -(np.dot(np.linalg.inv(a_alea), r_alea)) #para calcular puntos fijos       

#D = 0.3 #coeficiente de difusión/migracion total cte o np.array para función migración
Disp = {'b':0.3, 'm':1.0, 'i':1.0} #tasas dipersión para nueva función migracion_esp

m_milpa = 0.3 #taza muerte negra cte o np.array
m_intensivo = 0.6 #taza muerte blanca cte o np.array

h = 0.001 #diferencial de cambio en t (euler y graficas)
t_total = 1000 #tiempo total de simulacion
iter_difymuerte = 1 #iteraciones de dif y muerte entre cada una de lotka volterra

now = time.strftime('%c')

"""CORRIDAS"""

"""para correr valores guardados"""
os.chdir('/Users/personal2/Desktop/equilibrio2')
CI = glob.glob('*-c.txt')
Nombres = dict()
idx=1
for nci in CI:
    pn = nci[:-5]
    ntc = pn+"v.txt"
    nma = pn+"m.txt"
    ci = np.loadtxt(nci)
    tc = np.loadtxt(ntc)
    ma = np.loadtxt(nma)
    Nombres[idx] = [ci,tc,ma]
    idx += 1   
actual= int(sys.argv[1])  #número de corrida a utilizar. va de 1 a n.
print CI[actual-1]

"""corre el sistema aleatorizado sin matriz espacial"""
 
x= correr_aleatorizado_sin_matriz(n_especies, t_total, d_lotkavolterra_alea, r_alea, a_alea, True) #corre la ecuación con valores recién generados
#x= correr_aleatorizado_sin_matriz(n_especies, t_total, d_lotkavolterra_alea, Nombres[actual][1], Nombres[actual][2], True) #la corre con valores guardados
#sh=shannonwiener(x)
#print ((sum([i>0.0001 for i in x]))/(n_especies*1.0)) #cuantas sobrevivieron con x mayor a 0.0001
#print nefectivo(sh) # num efectivo de especies
#print sh # índice de shannon-wiener
#for i in range(len(x)): print x[i]    #total de cada especie
#print sum(x)     #suma total

     
"""características de la matriz de interacciones"""
#Para ver CONDICIONANTE de la matriz:
#print " de a_alea", np.linalg.cond(a_alea)
#print " de amodificada", np.linalg.cond(amodificada)

#para ver valores sigulares:
#_,s1,_ = np.linalg.svd(amodificada)
#_,s2,_ = np.linalg.svd(a_alea)
#print "Valores singulares de amodificada ", s1,"\n"
#print "Valores singulares de a_alea ", s2, "\n"

#print"determinante", np.linalg.det(a_alea)
#print "eigenvalores", np.linalg.eigvals(a_alea) 

"""Para GUARDAR experimento si alguna condición se cumple"""
#if ((sum([i>0.0000 for i in x]))/(n_especies*1.0))*100 <100: #qué porcentaje de las especies sobrevivieron
#now = time.strftime('%c')
#guardarexperimento(now, a_alea, r_alea, x_0) 


"""correr muchos y GUARDAR SISTEMAS DONDE SOBREVIVE  el x% """
  
#for n in range(100):
#    n_especies = 10
#    a = matriznicho(n_especies)
#    r = nicho_taza_crecimiento(a)
#    a_alea = a_aleatorizada(n_especies, a)
#    r_alea = r_aleatorizado(r)
#    h = 0.001 
#    t_total = 500
#    x= correr_aleatorizado_sin_matriz(n_especies, t_total, d_lotkavolterra_alea, r_alea, a_alea, False)
#    if ((sum([i>0.0000 for i in x]))/(n_especies*1.0))*100 == 100:
#        print "se guarda"
#        now = time.strftime('%c')
#        guardarexperimento(now, a_alea, r_alea, x_0)
               
                                                  
"""Distribucion de matriz agroecologica"""

#tipo = genera_tipo_matriz_agroecologica(x_celdas, y_celdas, n_bosque=5, posicion_bosque=[(1,1),(1,8),(5,5),(8,1),(8,8)], n_milpa=numerodemilpas, posicion_milpa="random") #con 5 bosques
#tipo = genera_tipo_matriz_agroecologica(x_celdas, y_celdas, n_bosque=10, posicion_bosque=[(1,1),(1,8),(5,5),(8,1),(8,8),(2,1),(7,1),(2,8),(7,8),(5,4)], n_milpa=numerodemilpas, posicion_milpa="random") #con 10 bosques
#print tipo
                          

"""Poblacion inicial"""
#poblacion_0 = genera_poblacion_inicial(tipo, n_especies, p0_bosque = Nombres[actual][0], p0_milpa=0, p0_intensivo=0)
#poblacion_0 = genera_poblacion_inicial(tipo, n_especies, p0_bosque = 'random', p0_milpa=0, p0_intensivo=0)
#poblacion_0 = genera_poblacion_inicial(tipo, n_especies, p0_bosque = x_0, p0_milpa=0, p0_intensivo=0)
#print CI[actual-1] #escribe nombre dle archivo abierto
#print "numero de milpas", numerodemilpas
              
#poblacion = correr_2D(poblacion_0, tipo, t_total, n_especies, d_lotkavolterra_alea, r_alea, a_alea, m_milpa, m_intensivo, Disp, 'vecinos8')
#poblacion = correr_2DMM(poblacion_0, tipo, t_total, n_especies, d_lotkavolterra_alea, Nombres[actual][1], Nombres[actual][2], m_milpa, m_intensivo, Disp, 'vecinos8')
#hay correr_2D, correr_2DMM o correr_2DL
#for idx in range(n_especies):
#    print "Especie ",idx
#    print poblacion[-1,:,:,idx]

"""medidas para sistemas que llegan a punto fijo"""
#riqueza=[]
#for idx in range(n_especies):
#    print np.sum(poblacion[-1,:,:,idx]) #suma de cada especie en la última iteración
#    if np.sum(poblacion[-1,:,:,idx])>0.0001: riqueza.append(1)  #vivos arriba de 0.0001   
#print "sum(x)", np.sum(poblacion[-1,:,:,:]) #suma total en la última iteración
#sumbosques=np.sum(poblacion[-1,1,1,:])+np.sum(poblacion[-1,1,8,:])+np.sum(poblacion[-1,5,5,:])+np.sum(poblacion[-1,8,1,:])+np.sum(poblacion[-1,8,8,:]) #suma total en los 5 bosques en la última iteración
#print "sum bosques", sumbosques #sum(x) bosques
#print "sum matriz", np.sum(poblacion[-1,:,:,:])-sumbosques #sum(x) matriz
#print "riqueza", np.sum(riqueza) 
#shannon = alpha_shannon(poblacion[-1])
#print "shannon", shannon
#print "nefec", nefectivo(shannon)
#
#
##riqueza promedio en bosques y matriz y alfa y beta 
#riqb=[]
#riqm=[]
#for i in range(len(tipo)):
#    for j in range(len(tipo)):
#        if tipo[i][j] == 'b':
#            riquezas=(poblacion[-1,i,j,:]>0.0001).sum()
#            riqb.append(riquezas)
#        else:
#            riquezas2=(poblacion[-1,i,j,:]>0.0001).sum()
#            riqm.append(riquezas2)
#print "riq bosques", np.mean(riqb) #riqueza promedio en bosques
#print "riz matriz", np.mean(riqm) #riqueza promedio en matriz
#promtot=float(np.sum(riqb)+np.sum(riqm))/(x_celdas*y_celdas)
#print "riq promedio", promtot #riqueza promedio total
#print "riq beta", np.sum(riqueza)/promtot #riqueza beta
#
##shannon promedio en bosques y matriz y alfa y beta
#shanpromb=[]
#shanpromm=[]
#for i in range(len(tipo)):
#    for j in range(len(tipo)):
#        if tipo[i][j] == 'b':
#            shanpromb.append(shannonwiener(poblacion[-1,i,j,:]))
#        else:
#            shanpromm.append(shannonwiener(poblacion[-1,i,j,:]))            
#print "shannon bosques", np.mean(shanpromb) #shannon promedio en bosques
#print "shannon mat", np.mean(shanpromm) #shannon promedio en matriz
#print "shannon prom", (np.sum(shanpromb)+np.sum(shanpromm))/(x_celdas*y_celdas) #shannon promedio total
#print "shannon beta", shannon-((np.sum(shanpromb)+np.sum(shanpromm))/(x_celdas*y_celdas)) #shannon beta

"""medidas para sistemas oscilantes"""

#riqueza=[]
#promedios=[]
#for idx in range(n_especies):
#    promedio=np.mean(np.sum(np.sum(poblacion[-200:,:,:,idx],1),1)) #promedio de cada especie en las últimas 200 iteraciones
#    print promedio
#    promedios.append(promedio) #pone el promedio de cada especie en una lista
#    if np.mean(np.sum(np.sum(poblacion[-200:,:,:,idx],1),1))>0.0001: riqueza.append(1)  #si el promedio de cada sp está arriba de 0.0001, agrega 1 a la lista "riqueza"
#print "sum(x) promedio", np.sum(promedios) #suma promedio de las especies con las últimas 200 iteraciones
#sumbosques=np.mean(np.sum(poblacion[-200:,1,1,:],1))+np.mean(np.sum(poblacion[-200:,1,8,:],1))+np.mean(np.sum(poblacion[-200:,5,5,:],1))+np.mean(np.sum(poblacion[-200:,8,1,:],1))+np.mean(np.sum(poblacion[-200:,8,8,:],1)) #suma total en los 5 bosques en la última iteración
#print "sum(x) bosques promedio", sumbosques #sum(x) bosques
#print "sum(x) matriz promedio", np.sum(promedios)-sumbosques #sum(x) matriz
#print "riqueza promedio", np.sum(riqueza) 
#poblacionpromedio=np.mean(poblacion[-200:,:,:,:],0) #promedia sobre eje tiempo: arroja array de forma (10,10,10) con las poblaciones promedio de las últimas 200 iteraciones
#shannon = alpha_shannon(poblacionpromedio)
#print "shannon", shannon
#print "nefectivo", nefectivo(shannon)
#
##riqueza promedio en bosques y matriz y alfa y beta 
#riqb=[] #aquí irá la riqueza en cada celda bosque para depués promediar
#riqm=[] #aquí irá la riqueza en cada celda matriz para depués promediar
#for i in range(len(tipo)):
#    for j in range(len(tipo)):
#        if tipo[i][j] == 'b':
#            riquezas=(np.mean(poblacion[-200:,i,j,:],0)>0.001).sum() #usa array promediado en eje tiempo
#            riqb.append(riquezas)
#        else:
#            riquezas2=(np.mean(poblacion[-200:,i,j,:],0)>0.001).sum() #usa array promediado en eje tiempo
#            riqm.append(riquezas2)
#print np.mean(riqb) #riqueza promedio en bosques
#print np.mean(riqm) #riqueza promedio en matriz
#promtot=float(np.sum(riqb)+np.sum(riqm))/(x_celdas*y_celdas)
#print promtot #riqueza promedio total
#print np.sum(riqueza)/promtot #riqueza beta
#
##shannon promedio en bosques y matriz y alfa y beta
#shanpromb=[]
#shanpromm=[]
#for i in range(len(tipo)):
#    for j in range(len(tipo)):
#        if tipo[i][j] == 'b':
#            shanpromb.append(shannonwiener(np.mean(poblacion[-200:,i,j,:],0))) #usa array promediado en eje tiempo
#        else:
#            shanpromm.append(shannonwiener(np.mean(poblacion[-200:,i,j,:],0))) #usa array promediado en eje tiempo            
#print np.mean(shanpromb) #shannon promedio en bosques
#print np.mean(shanpromm) #shannon promedio en matriz
#print (np.sum(shanpromb)+np.sum(shanpromm))/(x_celdas*y_celdas) #shannon promedio total
#print shannon-((np.sum(shanpromb)+np.sum(shanpromm))/(x_celdas*y_celdas)) #shannon beta

"""Graficas varias"""

#heatmap_tipo('tipo.png',tipo)
#heatmaps_especie_en_tiempo('especie.png',poblacion, labels = map(str, range(n_especies)))

#sum_xy = np.sum(np.sum(poblacion, axis=1),axis=1) #calcula sum en xy para cada especie
#plot_lineas_especies('sum_xy.png', sum_xy, t_total, labels=map(str, range(n_especies)))
#print "sumxy:", sum_xy[-1]


"""para correr desde 95 hasta 0"""

#n_especies = 10 # numero especies
#x_celdas = 10 # numero celdas en x
#    #para 1D x=1
#y_celdas = 10 # numero celdas en y
#              #        si 1 modelo en 1D
#
#a = matriznicho(n_especies)
#a_alea = a_aleatorizada(n_especies, a)
#r = nicho_taza_crecimiento(a)
#r_alea = r_aleatorizado(r)
#
#numerodemilpas = 90
#
#Disp = {'b':0.3, 'm':1.0, 'i':1.0} #tasas dipersión para nueva función migracion_esp
#
#m_milpa = 0.3 #taza muerte cte o np.array
#m_intensivo = 0.6 #taza muerte cte o np.array
#
#h = 0.001 #diferencial de cambio en t (euler y graficas)
#t_total = 1000 #tiempo total de simulacion
#iter_difymuerte = sys.argv[1] #input("num iteraciones dif y muerte?") #iteraciones de dif y muerte entre cada una de lotka volterra
#print "iteraciones", iter_difymuerte
#
#now = time.strftime('%c')
#
##para correr valores guardados
#os.chdir('/Users/personal2/Desktop/ciclo')
#CI = glob.glob('*-c.txt')
#Nombres = dict()
#idx=1
#for nci in CI:
#    pn = nci[:-5]
#    ntc = pn+"v.txt"
#    nma = pn+"m.txt"
#    ci = np.loadtxt(nci)
#    tc = np.loadtxt(ntc)
#    ma = np.loadtxt(nma)
#    Nombres[idx] = [ci,tc,ma]
#    idx += 1   
#actual= sys.argv[2] #input ("actual?" ) #número de corrida a utilizar. va de 1 a n.
#print "actual", actual
##print CI[actual-1] #escribe nombre dle archivo abierto
#print "de 90 a 0"
#


#for n in range(19):
#
#    tipo = genera_tipo_matriz_agroecologica(x_celdas, y_celdas, n_bosque=10, posicion_bosque=[(1,1),(1,8),(4,5),(8,1),(8,8),(2,1),(7,1),(2,8),(7,8),(5,4)], n_milpa=numerodemilpas, posicion_milpa="random")
#    print numerodemilpas
    #poblacion_0 = genera_poblacion_inicial(tipo, n_especies, p0_bosque = Nombres[actual][0], p0_milpa=0, p0_intensivo=0)
                
    #poblacion = correr_2DMM(poblacion_0, tipo, t_total, n_especies, d_lotkavolterra_alea, Nombres[actual][1], Nombres[actual][2], m_milpa, m_intensivo, Disp, 'vecinos8')

    #aquí pegar medidas para equilibrio o ciclo, según el caso  
      
    #riqueza=[]
    #for idx in range(n_especies):
    #    #print np.sum(poblacion[-1,:,:,idx]) #suma de cada especie en la última iteración
    #    if np.sum(poblacion[-1,:,:,idx])>0.0001: riqueza.append(1)  #vivos arriba de 0.0001   
    #print  np.sum(poblacion[-1,:,:,:]) #suma total en la última iteración
    #sumbosques=np.sum(poblacion[-1,1,1,:])+np.sum(poblacion[-1,1,8,:])+np.sum(poblacion[-1,5,5,:])+np.sum(poblacion[-1,8,1,:])+np.sum(poblacion[-1,8,8,:]) #suma total en los 5 bosques en la última iteración
    #print sumbosques #sum(x) bosques
    #print np.sum(poblacion[-1,:,:,:])-sumbosques #sum(x) matriz
    #print np.sum(riqueza) 
    ##shannon = alpha_shannon(poblacion[-1])
    ##print "shannon", shannon
    ##print "nefec", nefectivo(shannon)
    #
    #
    ##riqueza promedio en bosques y matriz y alfa y beta 
    #riqb=[]
    #riqm=[]
    #for i in range(len(tipo)):
    #    for j in range(len(tipo)):
    #        if tipo[i][j] == 'b':
    #            riquezas=(poblacion[-1,i,j,:]>0.0001).sum()
    #            riqb.append(riquezas)
    #        else:
    #            riquezas2=(poblacion[-1,i,j,:]>0.0001).sum()
    #            riqm.append(riquezas2)
    #print np.mean(riqb) #riqueza promedio en bosques
    #print np.mean(riqm) #riqueza promedio en matriz
    #promtot=float(np.sum(riqb)+np.sum(riqm))/(x_celdas*y_celdas)
    #print promtot #riqueza promedio total
    ##print "riq beta", np.sum(riqueza)/promtot #riqueza beta
    
    #riqueza=[]
    #promedios=[]
    #for idx in range(n_especies):
    #    promedio=np.mean(np.sum(np.sum(poblacion[-200:,:,:,idx],1),1)) #promedio de cada especie en las últimas 200 iteraciones
    #    #print promedio
    #    promedios.append(promedio) #pone el promedio de cada especie en una lista
    #    if np.mean(np.sum(np.sum(poblacion[-200:,:,:,idx],1),1))>0.0001: riqueza.append(1)  #si el promedio de cada sp está arriba de 0.0001, agrega 1 a la lista "riqueza"
    #print  np.sum(promedios) #suma promedio de las especies con las últimas 200 iteraciones
    #sumbosques=np.mean(np.sum(poblacion[-200:,1,1,:],1))+np.mean(np.sum(poblacion[-200:,1,8,:],1))+np.mean(np.sum(poblacion[-200:,5,5,:],1))+np.mean(np.sum(poblacion[-200:,8,1,:],1))+np.mean(np.sum(poblacion[-200:,8,8,:],1)) #suma total en los 5 bosques en la última iteración
    #print sumbosques #sum(x) bosques
    #print np.sum(promedios)-sumbosques #sum(x) matriz
    #print np.sum(riqueza) 
    ##poblacionpromedio=np.mean(poblacion[-200:,:,:,:],0) #promedia sobre eje tiempo: arroja array de forma (10,10,10) con las poblaciones promedio de las últimas 200 iteraciones
    ##shannon = alpha_shannon(poblacionpromedio)
    ##print "shannon", shannon
    ##print "nefectivo", nefectivo(shannon)
    #
    ##riqueza promedio en bosques y matriz y alfa y beta 
    #riqb=[] #aquí irá la riqueza en cada celda bosque para depués promediar
    #riqm=[] #aquí irá la riqueza en cada celda matriz para depués promediar
    #for i in range(len(tipo)):
    #    for j in range(len(tipo)):
    #        if tipo[i][j] == 'b':
    #            riquezas=(np.mean(poblacion[-200:,i,j,:],0)>0.001).sum() #usa array promediado en eje tiempo
    #            riqb.append(riquezas)
    #        else:
    #            riquezas2=(np.mean(poblacion[-200:,i,j,:],0)>0.001).sum() #usa array promediado en eje tiempo
    #            riqm.append(riquezas2)
    #print np.mean(riqb) #riqueza promedio en bosques
    #print np.mean(riqm) #riqueza promedio en matriz
    #promtot=float(np.sum(riqb)+np.sum(riqm))/(x_celdas*y_celdas)
    #print promtot #riqueza promedio total
    ##print np.sum(riqueza)/promtot #riqueza beta
    #
    #
    
#    numerodemilpas-=5



