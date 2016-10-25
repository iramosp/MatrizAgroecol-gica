import random
def genera_tipo_matriz_agroecologica(x_celdas, y_celdas, n_bosque=0, posicion_bosque=None, n_milpa=0, posicion_milpa=None): 
    """
    Nota: intercambié posicion_bosque=[] por posicion_bosque
    y posicion_milpa=[] por posicion_milpa
    (mutable data types should not be used as default vaules)

    Genera el paisaje
    Distribucion de matriz agroecologica
    Recibe una lista de tuples o un comando
    ej: [(2,3),(4,1),(2,2)]

    Recibe: x_celdas: número de celdas en x
            y_celdas: número de celdas en y
            n_bosque: cantidad de celdas con vegetación primaria (bosque)
            posicion_bosque: una lista de tuples o un comando, ej: [(2,3),(4,1),(2,2)], que corresponde a las coordenadas de las celdas de vegetación primaria (bosque)
            n_milpa: cantidad de celdas que representan milpa (manejo agroecológico)
            posicion_milpa: cómo ubicar las milpas
    """
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