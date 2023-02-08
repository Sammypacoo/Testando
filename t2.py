

# -*- coding: utf-8 -*-
'''
# Modules to import
from collections import defaultdict
import gzip
import pandas as pd
import re
import numpy as np
import time
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

# Files used

# This is the file that has Gene, CCDSCode, Strand, location in CCDS,CCDSid
arq1=pd.read_csv("./CCDS.current.txt", sep='\t')

# This is the file that has CCDSid, NP, NM
arq2=arquivo=pd.read_csv("./CCDS2Sequence.current.txt", sep='\t')

# This is the file was Produced from the file
# CCDS_protein_exons.current.faa, was produced using
#the exon code it has the ptns Sequences and the CCDSid
arq3=pd.read_csv("./exon_seq1",  sep='\t')

# Table with the transcripts present in NCBI but not CCDDS, created by the script "RemoveCCDS"/ Translated
arq4=pd.read_csv("./tabela1.csv",header= None, sep=',')

# Table with the transcripts present in NCBI with CCDS/Translated
arq5=pd.read_csv("./GCF_000001405.40_GRCh38.p14_translated_cds.faa",header= None, sep='\t')

#The file with the data about NM id
arq6=pd.read_csv("./GCF_000001405.40_GRCh38.p14_rna_from_genomic.fna", header= None, sep='\t')

# This file is the table NM and NP
arq7=pd.read_csv("./LRG_RefSeqGene.txt",sep='\t')

# This file is the PTMs database that was downloaded
ptm=pd.read_csv("./phosphositeplus_allptms_REFSEQ1", sep='\t')

arq5=arq5.dropna()
index=[]
for i in range(len(arq5)):
  index.append(i)
index
arq5=arq5.set_axis(index, axis=0)
arq5.tail()

def join_arq1_arq2_ncbi_id(gene):
   # Rename the columns
    arq1.rename(columns={"ccds_id":"#ccds"}, inplace = True)

    #Select the rows related to the gene that you are working on
    app=arq1[arq1["gene"]==gene]

    # Performs a Merge between files arq 1 and arq2 using the CCDSid column 
    #o join the two tables (arq 1 and arq2), so now we have both access to the information of arq1 and arq2  
    result = pd.merge(app, arq2, on="#ccds",how="left")

    #Remove duplicates rows, just to guarantee that is not repeating
    result = result.drop_duplicates()

    #Selects only the rows with the NCBI data, removing others like those from the essemble
    result=result[result["source"]=="NCBI"]

    #Keeps only the accepted version of the transcript, filtering the status_in_CCDS column to only what is accepted
    result=result[result["status_in_CCDS"]=="Accepted"]
    return result

#teste1=join_arq1_arq2_ncbi_id("APP")
#print(teste1)

# Find NP
def NP(s1):
  texto1="NP_"
  texto2="XP_"
  v=s1.split(' ')
  g=v[0]
  p=re.findall(r"NP_",g)
  if len(p)!=0:
    g1=g.split(texto1)
    g1=g1[1].split("_")
    g2=texto1+g1[0]
    
  else:
    p=re.findall(r"XP_",g)
    if len(p)!=0:
        g1=g.split(texto2)
        g1=g1[1].split("_")
        g2=texto2+g1[0]
    else:
        g2="erro"
  return g2

#Exons Limits
import re
def exon1 (s):
  l=[]
  texto="join"
  texto2="complement"
  texto3="location="
  try:
      l=exon(texto,s)
      return l
  except:
    try:
      l=exon(texto2,s)
      return l
    except:
       try:
          l=exon(texto3,s)
          return l
       except:
          #print(s)
          #print ("erro")
          l='erro'
          return l

def exon(texto,s):
  l=[]
  p=s.split(texto)
  p=p[1]
  p=p.replace("<","")
  p=p.replace(">","")
  p1=p.split(']')
  p1=p1[0]
  p1=p1.replace("(","")
  p1=p1.replace(")","")
  p2=p1.split(",")
  p2
  for i in range (len(p2)):
    c=p2[i]
    c=c.split("..")
    c1=int(c[0])
    c2=int(c[1])
    c=[]
    c.append(c1)
    c.append(c2)
    l+=[c]
  return l

  
def exon3(s): 
  try:    
    u=s.split("join")
    u1=u[1]
  except:
    try:
      u=s.split("complement")
      u1=u[1]
    except:
      u=s.split("location")
      u1=u[1]

  u2=u1.split(")")
  u2=u2[0]
  u2=u2.replace("(","")
  u2=u2.replace("[","")
  u2=u2.replace("=","")
  u2=u2.replace(":","")
  u2=u2.replace(">","")
  u2=u2.replace("<","")
  p2=u2.split(",")
  p2
  l=[]
  #print(p2)
  i=0
  while i <len(p2):
    #print("entrei")
    c=p2[i]
    #print(c)
    if c.find('..')!=-1:    
      c=c.split("..")
      c1=int(c[0])
      c2=int(c[1])
      c=[]
      c.append(c1)
      c.append(c2)
      l+=[c]
    i+=1    
  return l  

def exon4(s): 
  try:    
    u=s.split("join")
    u1=u[1]
  except:
    try:
      u=s.split("complement")
      u1=u[1]
    except:
      u=s.split("location")
      u1=u[1]
  u2=u1.split("]")
  u2=u2[0]
  u2=u2.replace("(","")
  u2=u2.replace("[","")
  u2=u2.replace("=","")
  u2=u2.replace(":","")
  u2=u2.replace(">","")
  u2=u2.replace("<","")
  p2=u2.split(",")
  p2
  l=[]
  i=0
  while i <len(p2):
    c=p2[i]
    if c.find('..')!=-1:    
      c=c.split("..")
      c1=int(c[0])
      c2=int(c[1])
      c=[]
      c.append(c1)
      c.append(c2)
      l+=[c]
    else:
      continue
    i+=1
  return l


#Preparar a tabela
def tabela_arq5(arq5,item):
  table=pd.DataFrame(columns=["Gene","NP","Exon","Protein"])
  t2=""
  lista=[]
  bol=False
  veri=0
  i=item
  while i<len(arq5):
    s=arq5.iloc[i,0]
    if s[0]==">":
        bol=True
        if len(lista)!=0 :
          lista.append(t2)   
          #print(lista)
          table.loc[len(table)] = lista 
        lista=[]
        #Gene
        u=s.split('gene=')
        u1=u[1]
        u2=u1.split(']')
        u3=u2[0]
        # NP
        v2=NP(s)
        #print(v2)
        if v2=="erro":
          bol=False
        #Exons Limits
        l=exon1(s)
        if l=="erro":
            try: 
              l=exon3(s)
              #print(l)
            except:
              bol=False
              print ("ERRO")
        if bol==True :
          lista.append(u3)
          lista.append(v2)
          lista.append(l)
        veri+=1  
    else:
      if bol==True:
            t2=t2+s
    if veri ==2:
      break
    i+=1
  return table

# Itens da linah importante para adicionar no NCBI
def NP_NM_items(table,tabela,arq7,i1,gene,l):

  item4=f'new_transcript_NCBI_{gene}_{i1}' # colocar f' para add nome diferente
  item7=table.iloc[l,2]
  item7=item7[0][0]
  item8=table.iloc[l,2]
  item8=item8[-1][1]
  item9=table.iloc[l,2]
  item9
  t=''
  
  for i in range(len(item9)):
    x=item9[i]
    a=int(x[0])
    b=int(x[1])
    u=f'{a}-{b}'
    if i==0:
      t=f'{u}'
    else:
      t=f'{t}, {u}'
  t=f'[{t}]'
  item9=t

  np=table.iloc[l,1]
  linha_np=arq7[arq7['Protein']==np]
  #print(np)
  #print(linha_np)
  nm=linha_np.iloc[0,5]
  #print(nm)
  item14=nm
  item15=np
  item14
  linha_ant=list(tabela.iloc[0,:])
  linha_new=[]
  date=[4,7,8,9,14,15]
  dic={4:item4,7:item7,8:item8,9:item9,14:item14,15:item15}
  for i3 in range(len(linha_ant)):
    #print(i3)
    if i3 in date:
      #print("entrei")
      item=dic[i3]
      linha_new+=[item]
    else:
      #print("entrei1")
      item=linha_ant[i3]
      linha_new+=[item]
  return linha_new

def join_table(gene,arq5):
  # Juntando os dados arq5 com a tabela produto das tabelas arq1 e arq2
  #nm_np=NM_NP_NCBI(arq5,arq6,gene)
  tabela=join_arq1_arq2_ncbi_id(gene)
  index=[]
  for i1 in range(len(tabela)):
    index.append(i1)
  tabela=tabela.set_axis(index, axis=0)
  gene1=f'gene={gene}]'
  contain_values = arq5[arq5[0].str.contains(gene1)].index.to_numpy()  
  contain_values=list(contain_values)
  x1=contain_values
  #print(len(x1))
  #print(x1)
  if len(x1) > 0:
    table1=pd.DataFrame()
    j=0
    for i2 in range(len(x1)):
      #print(i1)
      item=x1[i2]
      #print(item)
      table=tabela_arq5(arq5,item)     
      table1=pd.concat([table1,table])
  index=[]
  for i in range(len(table1)):
    index.append(i)
  table1=table1.set_axis(index, axis=0)
  for j1 in range(len(table1)):
    np=table1.iloc[j1,1]
    #print(j1)
    #print(np)
    testando= tabela[tabela["protein_ID"].str.contains(np)]
    if  (len(testando)) ==0:   
      #print("Nao tem no ccds")
      linha_new=NP_NM_items(table1,tabela,arq7,j,gene,j1)
      tabela.loc[len(tabela)] = linha_new
      #print(tabela)
      j+=1
  return tabela,table1

#teste2,teste3=join_table("APOE",arq5)
#print(teste2)
#print(teste3)

#Check the anotation

## CCDS Strand-
def aa_position_CCDS_strand_neg(ccdsid,arq3,tabela):
  table_genomic_aa=pd.DataFrame(columns=["condon","aminoácido"])
  #CCDS STRAND -
  contain_values1 = arq3[arq3["CCDSid"].str.contains(ccdsid)]
  seq=''
  for lines in range(len(contain_values1)):
    seq1=contain_values1.iloc[lines,4]
    seq=seq+seq1
  seq
  exons_lista=[]
  contain_values = tabela[tabela["#ccds"].str.contains(ccdsid)]
  contain_values
  exons=contain_values.iloc[0,9]
  exons=exons.replace("[","")
  exons=exons.replace("]","")
  exons=exons.split(",")
  bo=False
  exons.reverse()
  exons
  for ex in range(len(exons)):
    ran=exons[ex]
    #print(ran)
    ran=ran.split("-") 
    inicio=ran[1]
    inicio=int(inicio)
    inicio+=1          
    fim=ran[0]
    fim=int(fim)
    fim+=1
    if ex==len(exons)-1:
      fim=fim+3
    if bo==True:
      #print("Entrei_condon_truncado")
      if len(extra)==2:
        extra.append(inicio)
        inicio-=1
        #print(extra)
      if len(extra)==1:
        extra.append(inicio)
        extra.append(inicio-1)
        inicio-=2
        #print(extra)
      aa=seq[0]
      codon=extra
      new_row = {'condon':codon, 'aminoácido':aa}
      seq=seq[1:]
      table_genomic_aa = table_genomic_aa.append(new_row, ignore_index=True)
    bo=False
    while inicio-2>=fim:
      #print("entrei_condon_normal")
      codon=[inicio,inicio-1,inicio-2]
      #print(codon)
      inicio-=3
      aa=seq[0]
      new_row = {'condon':codon, 'aminoácido':aa}
      seq=seq[1:]
      table_genomic_aa = table_genomic_aa.append(new_row, ignore_index=True)
    if  fim-inicio == 0 or fim-inicio==-1 :
        bo=True
        if fim-inicio ==-1:
          extra=[inicio,inicio-1]
        else:
          extra=[inicio]
          #print(extra)
  return  table_genomic_aa

def aa_position_CCDS_strand_posi(ccdsid,arq3,tabela):

  ## CCDS Strand+
  table_genomic_aa=pd.DataFrame(columns=["condon","aminoácido"])
  #CCDS STRAND -
  contain_values1 = arq3[arq3["CCDSid"].str.contains(ccdsid)]
  seq=''
  for lines in range(len(contain_values1)):
    seq1=contain_values1.iloc[lines,4]
    seq=seq+seq1
  seq
  exons_lista=[]
  #print(tabela)
  contain_values = tabela[tabela["#ccds"].str.contains(ccdsid)]
  #print(contain_values)
  exons=contain_values.iloc[0,9]
  exons=exons.replace("[","")
  exons=exons.replace("]","")
  exons=exons.split(",")
  bo=False
  for ex in range(len(exons)):
    #print("oie")
    ran=exons[ex]
    #print(ran)
    ran=ran.split("-") 
    inicio=ran[0]
    inicio=int(inicio)
    inicio+=1          
    fim=ran[1]
    fim=int(fim)
    fim+=1
    #print(inicio)
    if ex==len(exons)-1:
      fim=fim-3
    #print(fim)
    if bo==True:
      #print("entre2")
      if len(extra)==2:
        extra.append(inicio)
        inicio+=1
        #print(extra)
      if len(extra)==1:
        extra.append(inicio)
        extra.append(inicio+1)
        inicio+=2
        #print(extra)
      aa=seq[0]
      codon=extra
      new_row = {'condon':codon, 'aminoácido':aa}
      seq=seq[1:]
      table_genomic_aa = table_genomic_aa.append(new_row, ignore_index=True)
    bo=False
    while inicio+2<=fim:
      codon=[inicio,inicio+1,inicio+2]
      inicio+=3
      aa=seq[0]
      new_row = {'condon':codon, 'aminoácido':aa}
      seq=seq[1:]
      table_genomic_aa = table_genomic_aa.append(new_row, ignore_index=True)
    if  fim-inicio < 2 and fim-inicio>-1 :
      if inicio-1 !=fim:
        bo=True
        if fim-inicio ==1:
          extra=[inicio,inicio+1]
        else:
          #print("entrei1")
          extra=[inicio]
          #print(extra)

  return  table_genomic_aa

def aa_position_NCBI_strand_posi(ccdsid,arq3,table1,tabela):
  ## NCBI Strand+
  table_genomic_aa=pd.DataFrame(columns=["condon","aminoácido"])
  np=tabela[tabela["#ccds"].str.contains(ccdsid)]
  np=np.iloc[0,15]
  table1_a=table1[table1["NP"].str.contains(np)]
  seq=table1_a.iloc[0,3]
  seq
  exons_lista=[]
  contain_values = tabela[tabela["#ccds"].str.contains(ccdsid)]
  contain_values
  exons=contain_values.iloc[0,9]
  exons=exons.replace("[","")
  exons=exons.replace("]","")
  exons=exons.split(",")
  bo=False
  for ex in range(len(exons)):
    #print("oie")
    ran=exons[ex]
    #print(ran)
    ran=ran.split("-") 
    inicio=ran[0]
    inicio=int(inicio)    
    fim=ran[1]
    fim=int(fim)
    #print(inicio)
    if ex==len(exons)-1:
      fim=fim-3
    #print(fim)
    if bo==True:
      #print("entre2")
      if len(extra)==2:
        extra.append(inicio)
        inicio+=1
        #print(extra)
      if len(extra)==1:
        extra.append(inicio)
        extra.append(inicio+1)
        inicio+=2
        #print(extra)
      aa=seq[0]
      codon=extra
      new_row = {'condon':codon, 'aminoácido':aa}
      seq=seq[1:]
      table_genomic_aa = table_genomic_aa.append(new_row, ignore_index=True)
    bo=False
    while inicio+2<=fim:
      codon=[inicio,inicio+1,inicio+2]
      inicio+=3
      aa=seq[0]
      new_row = {'condon':codon, 'aminoácido':aa}
      seq=seq[1:]
      table_genomic_aa = table_genomic_aa.append(new_row, ignore_index=True)
    if  fim-inicio < 2 and fim-inicio>-1 :
      if inicio-1 !=fim:
        bo=True
        if fim-inicio ==1:
          extra=[inicio,inicio+1]
        else:
          #print("entrei1")
          extra=[inicio]
          #print(extra)

  return  table_genomic_aa

## NCBI Strand-
def aa_position_NCBI_strand_neg(ccdsid,arq3,table1,tabela):
  table_genomic_aa=pd.DataFrame(columns=["condon","aminoácido"])
  np=tabela[tabela["#ccds"].str.contains(ccdsid)]
  np=np.iloc[0,15]
  table1_a=table1[table1["NP"].str.contains(np)]
  seq=table1_a.iloc[0,3]
  seq
  #lista1
  exons_lista=[]
  contain_values = tabela[tabela["#ccds"].str.contains(ccdsid)]
  contain_values
  exons=contain_values.iloc[0,9]
  exons=exons.replace("[","")
  exons=exons.replace("]","")
  exons=exons.split(",")
  bo=False
  exons.reverse()
  exons
  for ex in range(len(exons)):
    ran=exons[ex]
    #print(ran)
    ran=ran.split("-") 
    inicio=ran[1]
    inicio=int(inicio)       
    fim=ran[0]
    fim=int(fim)
    if ex==len(exons)-1:
      fim=fim+3
    if bo==True:
      #print("Entrei_condon_truncado")
      if len(extra)==2:
        extra.append(inicio)
        inicio-=1
        #print(extra)
      if len(extra)==1:
        extra.append(inicio)
        extra.append(inicio-1)
        inicio-=2
        #print(extra)
      aa=seq[0]
      codon=extra
      new_row = {'condon':codon, 'aminoácido':aa}
      seq=seq[1:]
      table_genomic_aa = table_genomic_aa.append(new_row, ignore_index=True)
    bo=False
    while inicio-2>=fim:
      #print("entrei_condon_normal")
      codon=[inicio,inicio-1,inicio-2]
      #print(codon)
      inicio-=3
      aa=seq[0]
      new_row = {'condon':codon, 'aminoácido':aa}
      seq=seq[1:]
      table_genomic_aa = table_genomic_aa.append(new_row, ignore_index=True)
    if  fim-inicio == 0 or fim-inicio==-1 :
        bo=True
        if fim-inicio ==-1:
          extra=[inicio,inicio-1]
        else:
          extra=[inicio]
          #print(extra)
  return  table_genomic_aa


def convertion_nucle_aa(tabela,table1,id_transcript):
  if tabela.iloc[0,6]=="+":
    bol=True
  else:
    bol=False
  x=id_transcript.find("CCDS")
  if x!=-1:
    if bol==True:
      #print("Strand+")
      aa_rela=aa_position_CCDS_strand_posi(id_transcript,arq3,tabela)
      #print(aa_rela)
    else:
      #print("Strand-")
      aa_rela=aa_position_CCDS_strand_neg(id_transcript,arq3,tabela)
      #print(aa_rela)
  else:
    if bol==True:
      #print("Strand+")
      aa_rela=aa_position_NCBI_strand_posi(id_transcript,arq3,table1,tabela)
      #print(aa_rela)
    else:
      #print("Strand-")
      aa_rela=aa_position_NCBI_strand_neg(id_transcript,arq3,table1,tabela)
      #print(aa_rela)
  aa_rela
  pos=[]
  for i in range (1,len(aa_rela)+1):
    pos+=[i]
  aa_rela["Positon"]=pos
  return aa_rela

#teste2,teste3=join_table("APOE",arq5)
#teste4=convertion_nucle_aa(teste2,teste3,"new_transcript_NCBI_APOE_0")
#print(teste4)
#print(teste4.iloc[0,0])


#1. Identificar qual gene ele ese refere
# Irei utilizar a tabela CCDS
#21 25891740
def gene_ccds(chrom, exon1, arq1,exon2=None):
  if exon2==None:
    exon2=exon1
  exon1=int(exon1)
  exon2=int(exon2)
  chrom=str(chrom)
  arq1a=arq1[arq1['#chromosome']==chrom]
  arq1a=arq1a.dropna()
  arq1a = arq1a[arq1a["cds_from"].str.contains("-") == False]
  arq1a['cds_from'] = arq1a['cds_from'].astype(float)
  arq1b=arq1a[arq1a["cds_from"]<=exon1]
  arq1b = arq1b[arq1b["cds_to"].str.contains("-") == False]
  arq1b['cds_to'] = arq1b['cds_to'].astype(float)
  arq1c=arq1b[arq1b["cds_to"]>=exon2]
  lista=arq1c["gene"]
  set_res = set(lista) 
  lista= (list(set_res))
  lista=list(lista)
  return lista
x=gene_ccds("21",25891740,arq1)
x

def gene_ncbi(chrom, exon1, arq5,exon2=None):
  arq5=arq5.dropna()
  lista1=[]
  #print(exon1)
  if exon2==None:
      exon2=exon1
  chrom=int(chrom)
  exon1=int(exon1)
  exon2=int(exon2)
  #print("Exon2")
  #print(exon2)
  if chrom <10:
    chrom=f'0{chrom}'
  j = arq5[arq5[0].str.contains(">lcl")]
  chrom=f'NC_0000{chrom}'
  j = j[j[1].str.contains(chrom)]
  #print(j)
  for i in range (len(j)):
    #print("Entrei")
    s=j.iloc[i,1]
    #print(s)
    try:
      #print("Entrei1")
      exon=exon3(s)
    except:
        exon=exon4(s)
    if len(exon)!=0:
      inicio=exon[0][0]
      inicio=int(inicio)
      fim=exon[-1][1]
      fim=int(fim)
      #print(inicio,fim)
      if inicio <= exon1 and fim >= exon2:
        #gene
        u=s.split('gene=')
        u1=u[1]
        u2=u1.split(']')
        u3=u2[0]
        lista1+=[u3]
  set_res = set(lista1) 
  lista1= (list(set_res))
  lista1=list(lista1)
  return lista1

### Pesquisa por mutação gênica
# Primeira parte do projeto 
# Parte 1
def gene_finder(lista):
  lista_erros=[]
  chrom=lista[0]
  chrom=int(chrom)
  exon1=lista[1]
  exon1=int(exon1)
  exon2=lista[2]
  exon2=int(exon2)
  resul=gene_ccds(chrom,exon1,arq1,exon2)
  if len(resul)==0:
    #print("Entrei")
    resul=gene_ncbi(chrom,exon1,arq5,exon2)
  if len(resul)>1 or len(resul)==0 :
    resul="It not a part of a gene"
  return resul

def find_genomic_range1(tabela,lista):

    chrom=lista[1]
    tabela1=tabela.iloc[:,[4,9]]
    lista1=[]
    chrom1=lista[1]
    chrom2=lista[2]
    for c in range(chrom1,chrom2+1,1):
      chrom=c
      for i2 in range(len(tabela1)):
        #print(i2)
        xa=tabela1.iloc[i2,0]
        x=xa.find("CCDS")
        x1=tabela1.iloc[i2,1]
        x1=x1.replace("[","")
        x1=x1.replace("]","")
        x1=x1.split(",")
        #print(x1)
        if x!=-1:
          #print("CCDS")
          for ccds in range(len(x1)):
            exon=x1[ccds]
            exon=exon.split("-")
            inicio=exon[0]
            inicio=int(inicio)
            inicio+=1          
            fim=exon[1]
            fim=int(fim)
            fim+=2
            r=range(inicio,fim)
            #print(r)
            if chrom in r:
              lista1+=[xa]
        else:
          #print("NCBI")
          for ncbi in range(len(x1)):
            exon=x1[ncbi]
            exon=exon.split("-")
            inicio=exon[0]                    
            inicio=int(inicio)
            fim=exon[1]
            fim=int(fim)
            fim+=1
            r=range(inicio,fim)
            if chrom in r:
              lista1+=[xa]
    return lista1

def LAMP_anotador(lista):
  variant=lista[1]
  x=gene_finder(lista)
  #print(x)
  if x[0]!="Erro":
    gene=x[0]
    tabela,table1=join_table(gene,arq5)
  #print(tabela)
  lista1=find_genomic_range1(tabela,lista)
  tabela_posi=pd.DataFrame(columns=["NM","Aminoacid","Site"])
  string=""
  for j in range(len(lista1)):
    ccdsid=lista1[j]
    aa_rela=convertion_nucle_aa(tabela,table1,ccdsid)
    #print(aa_rela)
    find_nm=tabela[tabela["#ccds"]==ccdsid]
    if len(find_nm)==1:
      row=[]
      nm=find_nm.iloc[0,14]
      nm=nm.split(".")
      nm=nm[0]
      variant=lista[1]
      variant=[variant]
      xee=aa_rela[pd.DataFrame(aa_rela.condon.tolist()).isin(variant).any(1).values]
      aa=xee.iloc[0,1]
      posi=xee.iloc[0,2]
      row.append(nm)
      row.append(aa)
      row.append(posi)
      tabela_posi.loc[len(tabela_posi)] = row
      if j ==len(lista1)-1:
        string=f'{string}{nm}:p.{aa}{posi}'
      else:
        string=f'{string}{nm}:p.{aa}{posi};'
      
    else:
      for j1 in range(len(find_nm)):
        row=[]
        nm=find_nm.iloc[j1,14]
        nm=nm.split(".")
        nm=nm[0]
        variant=lista[1]
        variant=[variant]
        xee=aa_rela[pd.DataFrame(aa_rela.condon.tolist()).isin(variant).any(1).values]
        aa=xee.iloc[0,1]
        posi=xee.iloc[0,2]
        row.append(nm)
        row.append(aa)
        row.append(posi)
        tabela_posi.loc[len(tabela_posi)] = row
        if j1 ==len(find_nm)-1:
            string=f'{string}{nm}:p.{aa}{posi}'
        else:
            string=f'{string}{nm}:p.{aa}{posi};'
  return tabela_posi,string

def table_annovar(variante,annovar_resul):
  tabela_posi=pd.DataFrame(columns=["NM","Aminoacid","Site"])
  posi=annovar_resul.iloc[variante,9]  
  posi=posi.split(",")
  string=""
  for i in range(len(posi)):
    x=posi[i]
    row=[]
    if "p." in x:
      x=x.split(":")
      transcrito=x[1]
      aa_posi=x[4]
      aa_posi=aa_posi.split("p.")
      aa_posi=aa_posi[1]
      aa=aa_posi[0]
      aa_posi=aa_posi[1:]
      posicao=aa_posi[:-1]
      posicao=int(posicao)
      row.append(transcrito)
      row.append(aa)
      row.append(posicao)
      tabela_posi.loc[len(tabela_posi)] = row      
    if i ==len(posi)-1:
        string=f'{string}{transcrito}:p.{aa}{posicao}'
    else:
        string=f'{string}{transcrito}:p.{aa}{posicao};'
     
  return tabela_posi,string

def annovar_variant(variant):
  test=annovar_resul.iloc[variant,:3]
  test=test.tolist()
  test
  resul=[]
  for p in range(len(test)):
    valor=test[p]
    valor=int(valor)
    resul.append(valor)
  return resul
annovar_resul=pd.read_csv("./teste.hg38_multianno.csv", sep=',')

## Validador utilizando o ANNOVAR revelando que eu consigo informar o quais transcritos e a posição correta do aa
annovar_resul=annovar_resul.drop([0])
#teste5,s=table_annovar(0,annovar_resul)
#teste6=annovar_variant(0)
#teste7,s2=LAMP_anotador(teste6)
#print(teste5)
#print(teste6)
#print(teste7)

def verify(x,lamp):
  for j in range(len(x)):
    nm=x.iloc[j,0]
    #print(nm)
    x1=x[x["NM"]==nm]
    x1=x1.set_axis([0], axis=0)
    lamp1=lamp[lamp["NM"]==nm]
    lamp1=lamp1.set_axis([0], axis=0)
    bol=x1.equals(lamp1)
    if bol==False:
      final=False
      break
    if j==len(x)-1:
      if bol==True:
        final=True
  return final
#teste8=verify(teste5,teste7)
#print(teste8)

#Verificar
annovar_resul.head()
t=annovar_resul.iloc[:5,[0,1,2,3,4,6]]
lista_resul=[]
lista_structure_annovar=[]
lista_structure_lamp=[]
print(len(annovar_resul))
for i in range(5):
  print(i)
  x,string_anovar=table_annovar(i,annovar_resul)
  #print(x)
  lista=annovar_variant(i)
  #print(lista)
  lamp,string_lamp=LAMP_anotador(lista)
  #print(lamp)
  resul=verify(x,lamp)
  lista_resul.append(resul)
  lista_structure_annovar.append(string_anovar)
  lista_structure_lamp.append(string_lamp)

teste9=lista_resul
print(teste9)
print(len(teste9))
t["Annovar Structure"]=lista_structure_annovar
t["Lamp Structure"]=lista_structure_lamp
t["Check"]=teste9
teste10=t
print(t)
'''



import pandas as pd
import gzip

def get_vcf_names(vcf_path):
    with gzip.open(vcf_path, "rt") as ifile:
          for line in ifile:
            if line.startswith("#CHROM"):
                  vcf_names = [x for x in line.split('\t')]
                  break
    ifile.close()
    return vcf_names


names = get_vcf_names("./h1.vcf.gz")
vcf = pd.read_csv("./h1.vcf.gz", compression='gzip', comment='#', chunksize=10000, delim_whitespace=True, header=None, names=names)
print(vcf[20])