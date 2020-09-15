
#!/usr/bin/env python


import string
from sys import argv

filein = argv[1]
filenameout= argv[2]

fileout = open(filenameout+'.ped','w')
file = open(filein,'r')


x = file.readline()

allele_dic = {}
allele_dic['0']='A'
allele_dic['1']='T'
allele_dic['2']='C'
allele_dic['3']='G'



loci=[] #counts number of loci/snps, as a list of 1's

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
      x=x.replace(' ','') #
      n_snps=int(string.split(x[:-1],':')[1])#se queda solo con el # y lo guarda
      loci.append(n_snps) #agrega a la lista el # de snps

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
        x=x.split()
        ## FSpecific for Fastsimcoal2 output format
        x[2]=list(x[2])
        ini=x[0]+'\t'+x[1]+'\t '
        x=ini+str(list(x[2]))
        x=x.replace(",","") # Limpiando caracteres de lista
        x=x.replace("[","")
        x=x.replace("]","")
        x=x.replace("\'","")
        x=x.replace('"', '')
        x=x.replace('\r','') #quita saltos de linea



        while x!='\n': # Mientras no sea salto de linea: mientras siga la linea

#            y1=string.split(x[:-1],'\t')
            y1=string.split(x[:-1]) # splitea la line
            x=file.readline() # se salta a la sig linea??
            #x=x.replace('\r','') #Quita salto de linea
            x=x.split()
        ## FSpecific for Fastsimcoal2 output format
            x[2]=list(x[2])
            ini=x[0]+'\t'+x[1]+'\t '
            x=ini+str(list(x[2]))
            x=x.replace(",","") # Limpiando caracteres de lista
            x=x.replace("[","")
            x=x.replace("]","")
            x=x.replace("\'","")
            x=x.replace('"', '')
            x=x.replace('\r','') #quita saltos de linea
            #x=x.replace(' ','')
#            y2=string.split(x[:-1],'\t')
            y2=string.split(x[:-1]) # guarda linea again (que deberia ser del mismo tipo??)
            ttt25 = y1[0] #Guarda ID 1_1 etc

            
           # with open('ys', 'w') as f:
            #	print >> f, 'Filename', y1, y2
            
            

            if len(y1)==len(y2): #Si las primeras 2 lines de genotype daa miden igual
                genos1=y1[2:] #Guarda solo datos de geno (los 0 y 1)
                genos2=y2[2:] #Guarda solo datos de geno (los 0 y 1)

            elif len(y1)==len(y2)+2: #Si es same leng +2(pasa a partir de 1_10)
                genos1=y1[2:] #Guarda solo datos de geno (los 0 y 1)
                genos2=y2[:] #Guarda toda la linea


            # Construyendo map (o ped?) file
            out=y1[0]+' '+y1[0]+' 0 0 0 ' + t23

            if x=='\n': # Si por fin llega a una linea en blanco
                break # dejar de iterar sobre las lineas


            for g in range(len(genos1)): # itera sobre geno data
                #print(g)
                # imprime out (arriba) y traduce 0 y 1 a A y T
                out=out+' '+allele_dic[genos1[g]]+' '+allele_dic[genos2[g]]
                #print(out)
            fileout.write(out+'\n') # escribir en archivo de salida de PED file

            x=file.readline()
            x=x.replace('\r','')

    x=file.readline()
fileout.close() # cerrar archivo de salida

### MAP FILE WORKS FINE NO MOVERLE A ESTA MADRE

salida=open(filenameout+'.map','w')

#Al chile esta parte ya no c que brga

# either generate snps on chromosome with x number of linkage blocks or
# generate them on x number of chromosomes with equal number of snps on each
cat =1
snp_nb=0 # # of snps is 0
for g in range(len(loci)): # itera lista con # de snps
    print('g:', g)
    for gg in range(loci[g]): # por cada # de snps en la lista de # de snps
    	#print('gg:', gg)
        if cat == 23:
            cat = 1 #Reiniciar cat al llegar a 23
        out = str(cat)+' '+str(snp_nb+1)+' 0 '+str((gg+1)*1000000)+'\n'       
        #print(out) 
        cat = cat + 1
        snp_nb+=1
        salida.write(out)

#fileout.close()
