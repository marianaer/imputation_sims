#!/usr/bin/env python


import string
from sys import argv
import numpy as np

filein = argv[1]
filenameout= argv[2]

fileout = open(filenameout+'.ped','w')
file = open(filein,'r')

nb_ind= argv[3]
nb_ind=int(nb_ind)
x = file.readline()

allele_dic = {0:'A', 1:'T', 2:'C', 3:'G'}


loci=[]
ind_counter=0
ind_list=[]


while x!='': #mientras exista la linea (o sea, hasta que se acabe el archivo)
    x=x.replace('\r','') #quitar saltos de lnea
    if '[Profile]' in x[:-1]: #Si es la linea que dice profile
        x = file.readline() #  saltarse dos lineas
        x = file.readline() # arriba ^
        x=x.replace('\r','') #quitar saltos de line
        nb_samps=int(string.split(x[:-1],'=')[1]) #Leer renglon de # de samples

        sample_nb=0 #contador de numero de muestras que hay

################################################################################
    #if 'Num linked loci ' in x[:-1]: # en fsc2 debe ser Total number of polymorphic sites
    if 'Total number of polymorphic sites' in x[:-1]:
        x=x.replace(' ','')
        n_snps=int(string.split(x[:-1],':')[1])
        loci.append(n_snps)
        geno_array=np.zeros((n_snps,nb_ind),dtype='int32')
        
    if 'Reporting status of a maximum of' in x[:-1]:
        loci=[]
        n_snps=int(string.split(x[:-1],)[6])#se queda solo con el # y lo guarda
        loci.append(n_snps)
        geno_array=np.zeros((n_snps,nb_ind),dtype='int32')
        
    
################################################################################

    if 'SampleName=' in x[:-1]: #Si llega al samplename
        x = x.replace('\r', '') #quita salto
        x = x.replace(' ','')
        t22 = x.rsplit('=', -1)
        t23 = str(t22[1:])
        t23 = t23.replace('\\n', '')
        #t23.translate(str.maketrans('', '', string.punctuation))
        t23 = t23.translate(None, string.punctuation) # nombre de la sample




################################################################################

    if 'SampleData= {' in x[:-1]: #Si encuentra esta parte
        sample_nb+=1 # Agrega al contador de # de muestras que hay
        x=file.readline() # se salta a la sig linea
        x=x.replace('\r','') #quita saltos de linea
        #### Aqui tendria que ir lo de poner espacios




        while x!='\n': # Mientras no sea salto de linea: mientras siga la linea

        ## FSpecific for Fastsimcoal2 output format


            x=x.split()
            x[2]=list(x[2])
            ini=x[0]+'\t'+x[1]+'\t '
            x=ini+str(list(x[2]))
            x=x.replace(",","") # Limpiando caracteres de lista
            x=x.replace("[","")
            x=x.replace("]","")
            x=x.replace("\'","")
            x=x.replace('"', '')
            x=x.replace('\r','') #quita saltos de linea

            y1=string.split(x) # splitea la line
            x=file.readline() # se salta a la sig linea??
            x=x.split()
            x[2]=list(x[2])
            ini=x[0]+'\t'+x[1]+'\t '
            x=ini+str(list(x[2]))
            x=x.replace(",","") # Limpiando caracteres de lista
            x=x.replace("[","")
            x=x.replace("]","")
            x=x.replace("\'","")
            x=x.replace('"', '')
            x=x.replace('\r','') #quita saltos de linea
            y2=string.split(x) # guarda linea again (que deberia ser del mismo tipo??)
            ttt25 = y1[0] #Guarda ID 1_1 etc

 
            

            if len(y1)==len(y2): #Si las primeras 2 lines de genotype daa miden igual
                genos1=y1[2:] #Guarda solo datos de geno (los 0 y 1)
                genos2=y2[2:] #Guarda solo datos de geno (los 0 y 1)

            elif len(y1)==len(y2)+2: #Si es same leng +2( esto es si el ouput de fsc es diploide)
                genos1=y1[2:] #Guarda solo datos de geno (los 0 y 1)
                genos2=y2[:] #Guarda toda la linea

            geno_array[:,ind_counter]=np.array(genos1) #crea arr con genotipos
            ind_counter+=1 
            geno_array[:,ind_counter]=np.array(genos2) #crea arr con genotipos
            ind_counter+=1


            if x=='\n': # Si por fin llega a una linea en blanco
                break # dejar de iterar sobre las lineas


            x=file.readline()
            x=x.replace('\r','')

    x=file.readline()
#fileout.close() # cerrar archivo de salida
mask=np.zeros(n_snps,dtype='int32')


for g in range(n_snps):
    if len(np.unique(geno_array[g]))>2: # si tiene mas de 2 snp
        mask[g]=1

usable_snps=np.where(mask==0)
geno_usable=geno_array[usable_snps]

g1=geno_usable[:,::2] ## Odd indexes (1,3,5...)

g2=geno_usable[:,1::2] ## Even indexes (2,4,6...)

conc=np.vstack((g1,g2))


### PED AND MAP FILE CREATION WORK FINE. DONT CHANGE



fileout = open(filenameout+'.ped','w')

ids=[]
for i in range(1, nb_samps+1):
    for j in range( 0, nb_ind/2):
        if (j) % 2 == 0:
            ids.append((str(i)+'_'+str(j+1)))

##for i in range(1,nb_samps+1):
for j in range(0, nb_ind/2):
       # if (j) % 2 == 0:

    nt=np.vectorize(allele_dic.get)(conc[:,j]) 
    nt=np.array2string(nt, threshold=np.inf)
        #nt= nt.replace('\n', '')
    nt= nt.replace('[', '')
    nt= nt.replace("'", '')
    nt= nt.replace(']', '')
    #nt= nt.replace(' ', '') #esto es p contar cuantos snps, si quitas dice el # de snps
    nt=nt.replace('\n', '')
    nt=str(nt)
    ind=ids[j]

    out=ind+' ' + ind+ ' 0 0 0 '+'Sample_'+str(ind[0])+' '+nt+'\n'
        #print(out)
    fileout.write(out)
fileout.close()




### PED AND MAP FILE CREATION WORK FINE. DONT CHANGE



# fileout = open(filenameout+'.ped','w')

# for i in range(1,nb_samps+1):
#     for j in range(0, nb_ind/2):
#         print(j)
#         #if (j) % 2 != 0:

#         nt=np.vectorize(allele_dic.get)(geno_usable[:,j]) 
#         nt=np.array2string(nt, threshold=np.inf)
#         #nt= nt.replace('\n', '')
#         nt= nt.replace('[', '')
#         nt= nt.replace("'", '')
#         nt= nt.replace(']', '')
#         #nt= nt.replace(' ', '') #esto es p contar cuantos snps, si quitas dice el # de snps
#         nt=nt.replace('\n', '')
#         nt=str(nt)+' '+str(nt)+'\n'
#         ind=str(i)+'_'+str(j+1)
#         out=ind+' ' + ind+ ' 0 0 0 '+'Sample'+str(i)+' '+nt
#         #print(out)
#         fileout.write(out)
# fileout.close()



### MAP FILE WORKS FINE NO MOVERLE A ESTA MADRE

salida=open(filenameout+'.map','w')

#Al chile esta parte Sya no c que brga

# either generate snps on chromosome with x number of linkage blocks or
# generate them on x number of chromosomes with equal number of snps on each
cat=1
i=0
for g in usable_snps[0]:
    i=i+1
    if cat ==23: 
        cat=1
    pos=str((g*10000))
    if g ==0:
        pos='1'
    #if g in usable_snps[0]:
    out = str(cat)+' '+'snp_'+str(g)+' 0 '+pos+'\n'
    cat =  cat +1
    salida.write(out)

fileout.close()
fileout2 = open(filenameout+'.par','w')

genotypename = 'genotypename:\t' + filenameout+ '.ped'
snpname = 'snpname:\t' + filenameout + '.map'
indivname = 'indivname:\t' + filenameout + '.ped'
evecoutname = 'evecoutname:\t' + filenameout + '.evec'
evaloutname = 'evaloutname:\t' + filenameout + '.eval'

fileout2.write(genotypename+'\n')
fileout2.write(snpname+'\n')
fileout2.write(indivname+'\n')
fileout2.write(evecoutname+'\n')
fileout2.write(evaloutname+'\n')

fileout2.close()
