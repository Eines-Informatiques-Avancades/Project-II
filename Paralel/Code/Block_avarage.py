import numpy as np
import matplotlib.pyplot as plt
import sys

# parameters file
params_file=sys.argv[1]

# read the name of the data file
def read_file(filename):
    with open(filename, "r") as f:
        lines = f.readlines()

        i = 0
        for line in lines:

            if "Time Evolution:" in line:
                file1 = lines[i + 1]
                file1_split = file1.split('"')[1]
                break
            i = i + 1

    return file1_split



input_file=read_file(params_file)
results_file = open(f"{input_file}_mean_values.dat", "w")
#print('A información do método atópase en:  https://doi.org/10.1063/1.457480')
def blockAverage(datos,nome,tamaño_bloque_maximo=False,grafica=True):
    "Método para calcular incertidumes con datos correlacionados temporalmente, block-avarage. Sería ideal engadir o método de bootstraping para comparar"
    Nobs = len(datos) # número total de datos, dar un array ou lista
    tamaño_bloque_minimo=1 # mínimo tamaño de bloque, corresponde a non facer block_avarage
    if tamaño_bloque_maximo==False:
        tamaño_bloque_maximo=int(Nobs/4)#Criterio por se non se selecciona un dato óptimo de bloque
    Numero_de_bloques = tamaño_bloque_maximo-tamaño_bloque_minimo # total number of block sizes
    Media_bloque = np.zeros(Numero_de_bloques) 
    co = np.zeros(Numero_de_bloques) # definese no seguinte papper, pero é a varianza
    bloqueCtr=0
    incertidume=np.zeros(Numero_de_bloques)
    barra_erros=np.zeros(Numero_de_bloques)
    for datos_bloque in range(tamaño_bloque_minimo,tamaño_bloque_maximo):
        Nbloque=int(Nobs/datos_bloque) # En que bloque estou según o numero de observacions
        almacenamento_temporal=np.zeros(Nbloque) # isto é para almacenar a media según en que bloque esté, é temporal
        # Loop to chop datastream into blocks
        # and take average
        for i in range(1,Nbloque+1):
            #Esta é a parte máis complicada do código pois (fago de 1 a N+1, para empregar logo as definicións matemáticas)
            #E como vou calculando as medias según o bloque no que me atope, e dicir o número de datos por bloque
            comeza = (i-1)*datos_bloque
            remata = comeza+datos_bloque
            almacenamento_temporal[i-1] = np.mean(datos[comeza:remata])
        Media_bloque[ bloqueCtr] = np.mean(almacenamento_temporal)
        co[bloqueCtr]  = np.var(almacenamento_temporal)
        incertidume[ bloqueCtr]=np.sqrt(np.var(almacenamento_temporal)/((Nbloque-1)))
        barra_erros[ bloqueCtr]=np.sqrt(np.var(almacenamento_temporal)/(Nbloque-1))*1/(np.sqrt(2*(Nbloque-1)))
        bloqueCtr=bloqueCtr+1
    tamaño_bloque= np.arange(tamaño_bloque_minimo,tamaño_bloque_maximo)
    if grafica:
        plt.title(nome)
        plt.errorbar(tamaño_bloque, incertidume,barra_erros, marker='o',ls='None',markersize=8, capsize=6)
        plt.xlabel('Block length')#(Número de datos por bloque, que deben ter todos o mesmo tamaño)
        plt.ylabel('Statistical error')
        plt.xticks()
        plt.yticks()
        plt.savefig(f"{input_file}_{nome}_error.png",dpi=300,bbox_inches = 'tight')
        #plt.show()
        plt.clf()
        plt.errorbar(tamaño_bloque, Media_bloque, incertidume,marker='o',ls='None',markersize=8, capsize=6,color='Orange')
        plt.ylabel('Mean value')
        plt.title(nome)
        plt.xticks()
        plt.yticks()
        plt.xlabel('Block length')
        plt.savefig(f"{input_file}_{nome}_mean_value.png",dpi=300,bbox_inches = 'tight')
        plt.clf()
        #plt.show()
        #print(f"Mean value {nome} = {Media_bloque[-1]:.3f} +/- {incertidume[-1]:.3f}")
        results_file.write(f"Mean value {nome} = {Media_bloque[-1]:.3f} +/- {incertidume[-1]:.3f}" + "\n")
    return tamaño_bloque, Media_bloque, incertidume
    
##########################################
#Añadir a posteriori
print('The statistical analysis started')
print('Do you want to select the number of blocks? [yes/no]')
condicion=input()
if str(condicion)=='yes' or str(condicion)=='si' or str(condicion)=='Yes'or str(condicion)=='YES' or str(condicion)=='SI':
    print("How many blocks do you want to select? (an integer, please)")
    bloque=int(input())
else:
    bloque=False


###########################################
dato = np.loadtxt(input_file,skiprows=2)
#print(dato)
Time=dato[:,0] # Primera columna)
Kinetic=dato[:,1]
Potential=dato[:,2]
Total_energy=dato[:,3]
Temperature=dato[:,4]
Pressure=dato[:,5]

blockAverage(Kinetic,'Kinetic Energy',bloque,grafica=True)
blockAverage(Potential,'Potential Energy',bloque,grafica=True)
blockAverage(Total_energy,'Total Energy',bloque,grafica=True)
blockAverage(Temperature,'Temperature',bloque,grafica=True)
blockAverage(Pressure,'Pressure',bloque,grafica=True)
results_file.close()

print("The statistical analysis ended.")
