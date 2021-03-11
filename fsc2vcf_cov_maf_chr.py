#!/usr/bin/env python

import datetime
import string
from sys import argv
import numpy as np
import  re
import math


#from datetime import datetime
startTime = datetime.datetime.now()

np.set_printoptions(threshold=np.inf)


argv=['1','/home/mescobar/Escritorio/fsc_sims/scripts/awa.arp', 'test3', 'vcf','cov=1.0','maf=0.05']



filein = argv[1]
filenameout= argv[2]

#fileout = open(filenameout+'.ped','w')
file = open(filein,'r')

###nb_ind= argv[3]
###nb_ind=int(nb_ind)
x = file.readline()

allele_dic = {0:'A', 1:'T', 2:'C', 3:'G'}


loci=[]
ind_counter=0
ind_list=[]


while x!='': # while theres a line
    x=x.replace('\r','') # remove
    if '[Profile]' in x[:-1]: # if its the profile line
        x = file.readline() #  read next two lines
        x = file.readline() # 
        x=x.replace('\r','') 
        nb_samps=int(string.split(x[:-1],'=')[1]) # read how many samps

        sample_nb=0 # counts # of samples (populations)

################################################################################

    if 'Total number of polymorphic sites' in x[:-1]:
        x=x.replace(' ','')
        n_snps=int(string.split(x[:-1],':')[1])
        loci.append(n_snps)
       # geno_array=np.zeros((n_snps,nb_ind),dtype='int32')
        
    if 'Reporting status of a maximum of' in x[:-1]:
        loci=[]
        n_snps=int(string.split(x[:-1],)[6])#se queda solo con el # y lo guarda
        loci.append(n_snps)
        #geno_array=np.zeros((n_snps,nb_ind),dtype='int32')

    if 'polymorphic positions on chromosome' in x[:-1]:
        x=file.readline()
        x=x.replace('#', '')#se queda solo con el # y lo guarda
        x=x.replace(' ','')
        x=x.replace("'","")
        x=x.split(",")
        snp_pos=[int (i) for i in x]        
    
################################################################################

    if 'SampleName=' in x[:-1]: 
        x = x.replace('\r', '') 
        x = x.replace(' ','')
        t22 = x.rsplit('=', -1)
        t23 = str(t22[1:])
        t23 = t23.replace('\\n', '')
        #t23.translate(str.maketrans('', '', string.punctuation))
        t23 = t23.translate(None, string.punctuation) # name of the sample
	
# Keep # of inds
    if 'SampleSize=' in x[:-1]:
    	x = x.replace('\r', '')
    	inds = (string.split(x[:],"="))
    	inds = int(inds[1])
    	nb_ind=inds*nb_samps
        if ind_counter==0:
            geno_array=np.zeros((n_snps,nb_ind),dtype='int32')       
################################################################################


    if 'SampleData= {' in x[:-1]: 
        sample_nb+=1 # Adds to sample counter
        x=file.readline() # reads next line
        x=x.replace('\r','') 


        while x!='\n': # as long as its reading the line

        ## Parse genotype line


            x=x.split()
            x[2]=list(x[2])
            ini=x[0]+'\t'+x[1]+'\t '
            x=ini+str(list(x[2]))
            x=x.replace(",","") 
            x=x.replace("[","")
            x=x.replace("]","")
            x=x.replace("\'","")
            x=x.replace('"', '')
            x=x.replace('\r','')

            y1=string.split(x)
            x=file.readline() 
            x=x.split()
            x[2]=list(x[2])
            ini=x[0]+'\t'+x[1]+'\t '
            x=ini+str(list(x[2]))
            x=x.replace(",","") 
            x=x.replace("[","")
            x=x.replace("]","")
            x=x.replace("\'","")
            x=x.replace('"', '')
            x=x.replace('\r','') 
            y2=string.split(x) 
            ttt25 = y1[0] #Guarda ID 1_1 etc


            if len(y1)==len(y2): 
                genos1=y1[2:] # stores genotype data (los 0 y 1)
		genos1=genos1[0:n_snps]
                genos2=y2[2:] # stores genotype data (los 0 y 1)
		genos2=genos2[0:n_snps]		

            elif len(y1)==len(y2)+2: 
                genos1=y1[2:] #Guarda solo datos de geno (los 0 y 1)
                genos2=y2[:] #Guarda toda la linea


                    
            geno_array[:,ind_counter]=np.array(genos1) # stores genotype data
            ind_counter+=1 
         
            geno_array[:,ind_counter]=np.array(genos2)
            ind_counter+=1
            #print(ind_counter)




            if x=='\n': # If its not sample data anymore
                break # stop looping


            x=file.readline()
            x=x.replace('\r','')

    x=file.readline()
#fileout.close() 


     


mask=np.zeros(n_snps,dtype='int32') # create mask array for triallelic snps


for g in range(n_snps):
    if len(np.unique(geno_array[g]))!=2: # if there are more than 2 alleles at a given snp
        mask[g]=1

usable_snps=np.where(mask==0) # Keep only biallelic snps
geno_usable=geno_array[usable_snps] 



#####################
#### MAF filter #####
#####################

maf_arg = re.compile("maf*") 
#look for "missing" in argv array
maf = filter(maf_arg.match, argv)
maf=str(maf)
maf=maf.split("=")
maf=str(maf[1])
maf=maf.replace("'", '')
maf=maf.replace("]", '')
maf=maf #Stores missing %

#nb_ind = # of alleles
counts_list=[]

for i in range(geno_usable.shape[0]):

    c=((np.unique(geno_usable[i,:], return_counts=True))[1])

    if min(c) > (float(maf)*nb_ind):
        counts_list.append(i)

geno_usable=geno_usable[counts_list]

g1=np.vectorize(allele_dic.get)(geno_usable[:,::2]) ## Odd indexes (1,3,5...) genotypes form ind 1_1, 1_3...

g2=np.vectorize(allele_dic.get)(geno_usable[:,1::2]) ## Even indexes (2,4,6...) genotypes from 1_2, 1_4...

conc=np.vstack((g1,g2)) #concatenate 1_1 and 1_2 (assumes thyre both copies of the same individual)


### PED AND MAP FILE CREATION WORK FINE

if 'plink' in argv:
    fileout = open(filenameout+'.ped','w')

    ids=[]
    for i in range(1, nb_samps+1):
        for j in range( 0, nb_ind/nb_samps):
            if (j) % 2 == 0:
                ids.append((str(i)+'_'+str(j+1)))

    ##for i in range(1,nb_samps+1):
    for j in (range(nb_ind/2)):
           # if (j) % 2 == 0:

        #nt=np.vectorize(allele_dic.get)(conc[:,j]) # translate numbers to nucleotides 
        nt=conc[:,j]
        nt=np.array2string(nt)
            #nt= nt.replace('\n', '')
        nt= nt.replace('[', '')
        nt= nt.replace("'", '')
        nt= nt.replace(']', '')
        nt=nt.replace('\n', '')
        nt=str(nt)
        ind=ids[j]

        out=ind+' ' + ind+ ' 0 0 0 '+'Sample'+str(ind[0])+' '+nt+'\n'
            #print(out)
        fileout.write(out)
    fileout.close()




 ### PED AND MAP FILE CREATION 

#us=np.array(usable_snps[0])
ps=np.array(snp_pos)
usable_pos=ps[counts_list] # Usable positions filtered by maf

if 'plink' in argv:
    salida=open(filenameout+'.map','w')

    cat=1
    i=0
    for g in usable_pos:
        i=i+1
        if cat ==23: 
            cat=1
        pos=str(g)
        if g ==0:
            pos='1'
        #if g in usable_snps[0]:
        out = str(cat)+' '+'snp_'+str(i)+' 0 '+pos+'\n'
        cat =  cat +1
        salida.write(out)

    fileout.close()



 #####################
 ### VCF #############
 #####################


length=np.zeros(len(usable_pos),  dtype='int32') # stores the length
c=1
for i in (range(len(usable_pos))): # Only keep the last 22 snps of all
 	if i >= len(usable_pos)-22:
   		length[i]=c # appendes the position (length-1) of the last 22 snps (one belongs to each chr)
    			#length.append(c)
   	c=c+1
   	if c==23:
        		#print(c)
       		c=1
chrm=length[-22:]
        		
        		
contig_length=usable_pos[-22:]+1 # the length is the largest position of a snp in each chr
	
contig_length=contig_length[chrm.argsort()] #sorts the position of each chr from chr 1 to 22
	
inds=[]
for i in range(1, nb_samps): #We only care about pop 1 and 2
	for j in range( 0, nb_ind/nb_samps):
   		if (j) % 2 == 0:
   			inds.append((str(i)+'_'+str(j+1)+'_'+str(i)+'_'+str(j+1)))

inds_pop1=inds[:len(inds)/2]
inds_pop2=inds[len(inds)/2:]

phase_arr=np.zeros((len(geno_usable), nb_ind/3), dtype='object') # Divide by 6 bc we only want to output pop1 and pop2 asdiploid organisms (we have 3 pops)

for i in range(nb_ind/3):
    for g in range(len(geno_usable)): # Divide by 3 bc we only want to output pop1 and pop2

        g1_allele= g1[g,i]        
        g2_allele= g2[g,i]
        ref_alt=np.unique(np.concatenate((g1[g,:], g2[g,:])))
        if len(ref_alt)<2: #if the site is monomorphic. creo que no pasa, ya quitamos los monomorficos antes # print('0 | 0', np.unique(g1[g,:]))
            phase_arr[g,i]='0|0'

        else:        
            if g1_allele==ref_alt[0] and g2_allele==ref_alt[0]: # if they're the same BUT the site is  not monomorphic= hom ref
                    #    print('0 | 0', g1_allele, g2_allele)
                phase_arr[g,i]='0|0'


            elif g1_allele==ref_alt[1] and g2_allele==ref_alt[1]: # if they're the same BUT the site is  not monomorphic= hom alt
                   #   print('1 | 1', g1_allele, g2_allele)
                phase_arr[g,i]='1|1'

            elif g1_allele==ref_alt[0] and g2_allele==ref_alt[1]: #het
                      #  print('0 | 1', g1_allele, g2_allele)
                phase_arr[g,i]='0|1'

            elif g1_allele==ref_alt[1] and g2_allele==ref_alt[0]: #het
                      #  print('1 | 0', g1_allele, g2_allele)
                phase_arr[g,i]='1|0'

			#elif g1_allele==np.unique(np.concatenate((g1[g,:], g2[g,:])))[1] and g2_allele==np.unique(np.concatenate((g1[g,:], g2[g,:])))[1]:
                       # print('1 | 0', g1_allele, g2_allele)
			#	phase_arr[g,i]='1|0'

chr_ord=np.resize(range(1,23), len(usable_pos))# assigns a chr to each site




if 'vcf' in argv:

   # inds_pop2=inds[len(inds)/2:]
	d = datetime.datetime.now()
	date=(d.strftime("%Y"),d.strftime("%x").split('/')[0],d.strftime("%x").split('/')[1])
	date=str(date)
	date=date.replace("'", '')
	date=date.replace(",", '')
	date=date.replace("(", '')
	date=date.replace(")", '')
	date=date.replace(" ", '')
	date=('##fileDate='+date)

for i in range(1,3):
    for j in range(1,23):
        chr_file=open(filenameout+'_'+maf.split(".")[1]+'maf'+'_Pop'+str(i)+'_chr'+str(j)+'_phased'+'.vcf', 'w')
        chr_file.write(date+'\n')
        chr_file.write('##contig=<ID='+str(j)+',length='+str(contig_length[j-1])+'>\n')
        chr_file.write('##INFO=<ID=PR,Number=0,Type=Flag,Description="Provisional reference allele, may not be based on real reference genome">\n')
        chr_file.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        header=('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t')

        inds=str(inds)
        inds=inds.replace("'", '')
        inds=inds.replace("[", '')
        inds=inds.replace("]", '')
        inds=inds.replace(",", '')
        inds=inds.replace(" ", '\t')

        if i == 1: # Pop1
            chr_file.write(header+inds[:len(inds)/2]+'\n')# only outputs individuals from pop1
            chr_sites=np.where(chr_ord==j) # Keeps the sites in this chr
            geno_pos=usable_pos[chr_sites]
            phase_chr=phase_arr[(chr_sites)]
            chr_genos=g1[chr_sites]

            for l in range(len(geno_pos)): # For each snp in this chr
                pos=str(geno_pos[l]) #keep theposition
                ref=str(np.unique(chr_genos[l,:])[0]) #keep the ref allele

                if len(np.unique(chr_genos[l,:]))==2: 
                    alt=str(np.unique(chr_genos[l,:])[1]) #Keep the alt allele

                    out= str(j)+ '\t'+pos+'\t'+'snp_'+pos+'\t'+ref+'\t'+alt+'\t.\t.\tPR\t'+'GT'+'\t'+str(str(phase_chr[l, :len(inds_pop1)]).split())+'\n'
                    out=out.replace('[','')
                    out=out.replace(']','')
                    out=out.replace('"','')
                    out=out.replace("'",'')
                    out=out.replace(",",'')
                    out=out.replace(" ",'\t')
                    chr_file.write(out)

        if i == 2: # Pop2
            chr_file.write(header+inds[len(inds)/2:]+'\n')# only outputs individuals from pop2
            chr_sites=np.where(chr_ord==j) # Goes chr by chr
            geno_pos=usable_pos[chr_sites]
            phase_chr=phase_arr[(chr_sites)]
            chr_genos=g1[chr_sites]

            for l in range(len(geno_pos)): # For each snp
                pos=str(geno_pos[l]) #keep theposition
                ref=str(np.unique(chr_genos[l,:])[0]) #keep the ref allele

                if len(np.unique(chr_genos[l,:]))==2: 
                    alt=str(np.unique(chr_genos[l,:])[1]) #Keep the alt allele

                    out= str(j)+ '\t'+pos+'\t'+'snp_'+pos+'\t'+ref+'\t'+alt+'\t.\t.\tPR\t'+'GT'+'\t'+str(str(phase_chr[l, len(inds_pop2):]).split())+'\n'
                    out=out.replace('[','')
                    out=out.replace(']','')
                    out=out.replace('"','')
                    out=out.replace("'",'')
                    out=out.replace(",",'')
                    out=out.replace(" ",'\t')
                    chr_file.write(out)

#############################
### SIMULATE MISSING RATE ###
#############################

inds=[]
for i in range(1, nb_samps): #We only care about pop 1 and 2
    for j in range( 0, nb_ind/nb_samps):
        if (j) % 2 == 0:
            inds.append((str(i)+'_'+str(j+1)+'_'+str(i)+'_'+str(j+1)))  
            
#inds

def phred2prob(x):
    return 10.0**(-x/10.0) #error

def prob2phred(x):
    return -10*math.log10(x)

all_dic={}
all_dic['A']=0
all_dic['C']=1
all_dic['G']=2
all_dic['T']=3

def geno_caller_3GT(X,ref,alt,all_dic):
    #diploid caller assuming that only assesses likelihood for three possible genotypes (ref/ref,ref/alt,alt/alt)           
    GL=[0.0,0.0,0.0]

    count=0
    for g in range(len(X)):
        if all_dic.has_key(X[g][0])==False:
            continue
        err=phred2prob(X[g][1]) #el input es el error en phred
        tru=1-phred2prob(X[g][1])
      
        if X[g][0]==ref:
            GL[0]=GL[0]+math.log10(tru)
            GL[1]=GL[1]+math.log10((tru+err)/2)
            GL[2]=GL[2]+math.log10(err)           
        elif X[g][0]==alt:
            GL[0]=GL[0]+math.log10(err)
            GL[1]=GL[1]+math.log10((tru+err)/2)
            GL[2]=GL[2]+math.log10(tru)  
        else:
            GL[0]=GL[0]+math.log10(err) 
            GL[1]=GL[1]+math.log10(err)
            GL[2]=GL[2]+math.log10(err)
        count+=1

    if count==0:
        GL=[-9.0,-9.0,-9.0]
    return GL 



### If user specified to simulate a misingness rate:
r = re.compile("cov*") #look for "missing" in argv array
m = filter(r.match, argv)
if m: #If user did specify a coverage:
    m=str(m)
    m=m.split("=")
    m=str(m[1])
    m=m.replace("'", '')
    m=m.replace("]", '')
    m=m #Stores coverage %

    d = datetime.datetime.now()
    date=(d.strftime("%Y"),d.strftime("%x").split('/')[0],d.strftime("%x").split('/')[1])
    date=str(date)
    date=date.replace("'", '')
    date=date.replace(",", '')
    date=date.replace("(", '')
    date=date.replace(")", '')
    date=date.replace(" ", '')
    date=('##fileDate='+date)
    # VCF genotypes
    #m_rt=np.random.uniform(0, 1, (len(geno_usable), nb_ind/3)) #how much data we got
    #unknown_genos = np.where(m_rt>1-float(m))

    #phase_arr[unknown_genos]="./." # Add missing genotypes

    
    for indx, kk in enumerate(inds):
        if indx >= len(inds)/2: ### We only simulate misisng rate in pop2
            for j in range(1,23):
                chr_file=open(filenameout+'_'+maf.split(".")[1]+'maf_'+str(m)+'cov_'+'ind_'+kk+'_chr'+str(j)+'_unphased'+'.vcf', 'w')

                chr_file.write(date+'\n')
                chr_file.write('##contig=<ID='+str(j)+',length='+str(contig_length[j-1])+'>\n')
                chr_file.write('##INFO=<ID=PR,Number=0,Type=Flag,Description="Provisional reference allele, may not be based on real reference genome">\n')
                chr_file.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
                chr_file.write('##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">\n')
                chr_file.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read Depth">\n')
                chr_file.write('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">\n')
                header=('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t')
                chr_file.write(header+kk+'\n')


                chr_sites=np.where(chr_ord==j)
                geno_pos=usable_pos[chr_sites] # keep genomic positions of sites in chr1 
                phase_chr=phase_arr[(chr_sites)] #keep phased genotype info of sites in chr1
                chr_genos=g1[chr_sites]


                cov_arr=np.zeros((len(phase_chr), (nb_ind/3)), dtype='object')
                a1_dp=np.zeros((len(phase_chr), (nb_ind/3)), dtype='object')
                a2_dp=np.zeros((len(phase_chr), (nb_ind/3)), dtype='object')
        


                for inx, g in enumerate(phase_chr[:,indx]):
                #print(g)
                    cov=np.random.poisson(float(m)) #coverage for this site
    
                    if g=='0|0':
                        a1=cov     
                        a2=0
                        flip1_p=np.random.uniform(0, 1, (a1))
                        flip_a1=len(np.where(flip1_p>1-0.010)[0]) #1% chance of flip error. How many reads will flip
                        flip_a2=0
                        a1_err=a1-flip_a1
                        a2_err=flip_a2+flip_a1
        
                        cov_arr[inx, indx]=cov
                        a1_dp[inx, indx]=a1_err
                        a2_dp[inx, indx]=a2_err
        
        
                    elif g=='1|1':
                        a1=0
                        a2=cov
                        flip2_p=np.random.uniform(0, 1,(a2))
                        flip_a2=len(np.where(flip2_p>1-0.010)[0])
                        flip_a1=0
                        a1_err=flip_a1+flip_a2
                        a2_err=a2-flip_a2
        
                        cov_arr[inx, indx]=cov
                        a1_dp[inx,indx]=a1_err
                        a2_dp[inx,indx]=a2_err
            
        
                    elif g=='1|0' or g=='0|1':
                        a1=np.random.binomial(cov, 0.5)
                        a2=cov-a1
                        flip1_p=np.random.uniform(0, 1, (a1))
                        flip2_p=np.random.uniform(0, 1, (a2))
                        flip_a1=len(np.where(flip1_p>1-0.010)[0])
                        flip_a2=len(np.where(flip2_p>1-0.010)[0])
                        a1_err=a1-flip_a1+flip_a2
                        a2_err=a2-flip_a2+flip_a1    
    
                        cov_arr[inx,indx]=cov
                        a1_dp[inx,indx]=a1_err
                        a2_dp[inx,indx]=a2_err        
        
        
        
                    #for l in range(len(geno_pos)):

                    pos=str(geno_pos[inx])
                    ref=str(np.unique(chr_genos[inx,:])[0])
                
                    if len(np.unique(chr_genos[inx,:]))==2:
                        alt=str(np.unique(chr_genos[inx,:])[1]) #Keep the alt allele
                
                # For GL estimation
                    if (phase_chr[inx,indx])[0]=='0':
                        a1=ref
                    elif (phase_chr[inx,indx])[0]=='1':
                        a1=alt
                
                    if (phase_chr[inx,indx])[2]=='0':
                        a2=ref

                    elif (phase_chr[inx,indx])[2]=='1':
                        a2=alt
                    
                    
                
                    if a2_dp[inx,indx]==0 and a1_dp[inx,indx]!=0: # if all the cov goes to one allele
                        gg1=[a1,40]
                        d1=a1_dp[inx,indx]
                        d2=a2_dp[inx,indx]
                        log_GL=geno_caller_3GT([gg1 for i in range(d1)], ref, alt, all_dic)

                        #print("alavergaaaa")

                        if a1==ref and a2==alt:
                        	dp_ref=d1
                        	dp_alt=0
                        elif a1==alt and a2==ref: #1/0 toda la cvg va pal alt
                            dp_ref=0
                            dp_alt=d1
                        elif a1==ref and a2==ref:
                        	dp_ref=d1
                        	dp_alt=0
                        elif a1==alt and a2==alt:
                        	dp_ref=0
                        	dp_alt=d1

                    
                    elif a1_dp[inx,indx]==0 and a2_dp[inx,indx]!=0: #if all the cov goes to the other allele

                        gg2=[a2,40]
                        d1=a2_dp[inx,indx]
                        d2=a2_dp[inx,indx]
                        log_GL=geno_caller_3GT([gg2 for i in range(d2)], ref, alt, all_dic)

                        if a1==ref and a2==alt: #0/1
                            dp_ref=0
                            dp_alt=d2
                        elif a1==alt and a2==ref: #1/0
                        	dp_ref=d2
                        	dp_alt=0
                        elif a1==ref and a2==ref: #0/0
                        	dp_ref=d2
                        	dp_alt=0
                        elif a1==alt and a2==alt: #1/1
                        	dp_ref=0
                        	dp_alt=d2

                    elif a1_dp[inx,indx]!=0 and a2_dp[inx,indx]!=0:
                    	gg1=[a1,40]
                    	d1=a1_dp[inx,indx]
                    	gg2=[a2,40]
                    	d2=a2_dp[inx,indx]
                    	lista1=[gg1 for i in range(d1)]
                    	lista2=[gg2 for i in range(d2)]

                    	if a1==alt and a2==ref:
                    		log_GL=geno_caller_3GT(lista1+lista2, a2, a1, all_dic)
                    		dp_ref=d2
                    		dp_alt=d1
                    	elif a1==ref and a2==alt:
                    		log_GL=geno_caller_3GT(lista1+lista2, a1, a2, all_dic)
                    		dp_ref=d1
                    		dp_alt=d2
                    	elif a1==ref and a2==ref:
                    		a2=alt
                    		gg1=[a1,40]
                    		gg2=[a2,40]
                    		d1=a1_dp[inx,indx]
                    		d2=a2_dp[inx,indx]
                    		dp_ref=d1
                    		dp_alt=d2
                    		lista1=[gg1 for i in range(d1)]
                    		lista2=[gg2 for i in range(d2)]
                    		log_GL=geno_caller_3GT(lista1+lista2,a1,a2,all_dic)
                    	elif a1==alt and a2==alt: #1/1 ambos con cov (bc of seq error)
                    		a1=ref 
                    		gg1=[a1,40]
                    		gg2=[a2,40]
                    		d1=a1_dp[inx,indx]
                    		d2=a2_dp[inx,indx]
                    		dp_ref=d1
                    		dp_alt=d2
                    		lista1=[gg1 for i in range(d1)]
                    		lista2=[gg2 for i in range(d2)]
                    		log_GL=geno_caller_3GT(lista1+lista2,a1,a2,all_dic)
                   
                    
                    if cov_arr[inx,indx]!=0:

                        raw_PL=[math.log10(10**i)*-10 for i in log_GL]
                        norm_PL=[(round(i-min(raw_PL))) for i in raw_PL]
                        GQ=int(sorted(norm_PL)[1])
                        norm_PL=str(norm_PL).replace('.0', '')
                        norm_PL=norm_PL.replace(' ', '')

                    elif cov_arr[inx,indx]==0:
                        phase_chr[inx,indx]='./.'
                        GQ='.'
                        cov_arr[inx,indx]='0'
                        a1_dp[inx,indx]='0'
                        a2_dp[inx,indx]='0'
                        dp_ref='0'
                        dp_alt='0'
                        norm_PL='0,0,0' 

                    if len(np.unique(chr_genos[inx,:]))==2:
                        alt=str(np.unique(chr_genos[inx,:])[1])
                        out= str(j)+ '\t'+pos+'\t'+'snp_'+pos+'\t'+ref+'\t'+alt+'\t.\t.\tPR\t'+'GT:GQ:DP:AD:PL'+'\t'+str(str(phase_chr[inx,indx]).split())+':'+str(GQ)+':'+str(cov_arr[inx,indx])+':'+str(dp_ref)+','+str(dp_alt)+':'+norm_PL+'\n'
                        out=out.replace('[','')
                        out=out.replace(']','')
                        out=out.replace('"','')
                        out=out.replace("'",'')
                        out=out.replace("|",'/')
                        out=out.replace(" ",'\t')
                        out=out.replace(".0", '')
                        chr_file.write(out)


print datetime.datetime.now() - startTime 
