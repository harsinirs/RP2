import pandas as pd
import numpy as np
import re as re
import time
from tabulate import tabulate
from pyteomics import parser
from tabulate import tabulate
from pyteomics import parser
from collections import Counter
from itertools import permutations
import time

def readData(filename):
    from pyteomics import parser

    """
    Opens and reads in a csv formatted file exported from Mascot 
    This function takes a single argument, the name of an external LC-MS/MS search result file
    and reads it into a dataframe
    It expects proteins mapped to their peptides, masses, sequences and scores in columns and no additional search information
    It extracts distinct protein names and sequences into a list in order.
    """
    protname = []
    protseq=[]
    peptides=[]
    ppeptides=[]
    protdictname={}
    protdictseq={}
    msdata = pd.read_csv(filename)
    
    #saves unique entries of protein from search results
    protname=msdata.prot_desc.unique()
    protseq=msdata.prot_seq.unique()
    
    #cleaves proteins into peptides following tryptic cleavage rules
    for i in protseq:
        pep=parser.cleave(str(i), parser.expasy_rules['trypsin'], 2)
        peptides.append(list(pep))
    
    #filters peptides to exclude small hydrophillic and long hydrophobic peptides
    pp= [[ele for ele in sub if len(ele) in range(5,25)]for sub in peptides]
    
    #links proteins and peptides into a dictionary
    for i in range (0,len(protname)):
        protdictname[protname[i]]=pp[i]
        protdictseq[protname[i]]=protseq[i]
        
    return(pp,protdictname,protdictseq)

def peptidemass(pp):
     """
Calculates mass of peptides using their monoisotopic values as a linear sum
This function takes a single argument, the list of peptides
and using a mass dictionary, it calculates the mass of the peptides
it also records the possible variants of the peptide due to modifications by identifying the residues susceptible to modification 
The maximum number of modifications allowed per peptide is limited to 5 and all the possible permutations and combinations of modified sites is recorded as a new entry.
It creates a new dataframe of peptides and thier mass along with the number, type and site of modification.
    """
    #dictionary of amino acids mass
    mono = {'A' :71.0371, 'B': 114.5349,'C' :160.0306, 'D' :115.0269, 'E' :129.0426, 'F' :147.0684, 'G' :57.02146, 'H' :137.0589, 'I' :113.0841, 'K' :128.0950, 'L' :113.0841, 'M' :131.0405, 'N' :114.0429, 'P' :97.0528, 'Q' :128.0586, 'R' :156.1011, 'S' :87.0320, 'T' :101.0477, 'U': 150.9536, 'V' :99.0684, 'W' :186.0793, 'X':111, 'Y' :163.0633, '*' :0.0}
    masspep=[]
    
    #List of residues that can be modified by oxidation or deamination
    mod=['P','K','M']
    modnq=['N','Q']
    
    n=0
    
    for i in pp:
        for j in i:
            c=0
            m=0
            cn=0
            cc=0
            for k in range(0,len(j)):
                m=m+mono[j[k]]
                #counts number of modifiable residues
                if j[k] in mod:
                    c=c+1
                elif j[k] in modnq:
                    cn=cn+1
                elif j[k]=='C':
                    cc=c+1
                #list of position of modified residues
                pos = [k for k, ltr in enumerate(j) if ltr in mod]
                posnq = [k for k, ltr in enumerate(j) if ltr in modnq]
            #encoding modification presence, number and site
            for (a,b) in [(a,b) for a in range(c+1) for b in range(cn+1)]:
                if a+b<=5:
                    n=n+1
                    perm = list(permutations(pos, a))
                    permnq=list(permutations(posnq,b))
                    for (x,y) in [(set(x),set(y)) for x in perm for y in permnq]:
                        if x == set() and y==set():
                            masspep.append([n,m+18.0104+a*(16)+b*(0.984),j,a,b,' '.join(str(p) for p in x),' '.join(str(q) for q in y)])
                        else:
                            masspep.append([n,m+18.0104+cc*1+a*(16)+b*(0.984),j,a,b,' '.join(str(p) for p in x),' '.join(str(q) for q in y)])
    peptidemass = pd.DataFrame(masspep, columns=['pep no','mass','peptide','oxi','deam','oxipos','deampos'])
    peptidemass=peptidemass.drop_duplicates()
    peptidemass = peptidemass.reset_index(drop=True)#fix index after deletion of duplicates
    return(peptidemass)

def fragmentation(peptidemass):
    import time

    """
    This function fragments the peptides into b and y fragment ions
    using the peptide information, it create 2 dataframes for b and y ions 
    and encodes the modification information to each fragment
    """
    import time 
    dfb = pd.DataFrame(columns=['pep no.','peptide','pepmass', 'oxi','deam','b ion','bn','oxipos','deampos'])
    fragb=[]
    dfy = pd.DataFrame(columns=['pep no.','peptide','pepmass', 'oxi','deam','y ion','yn','oxipos','deampos'])
    fragy=[]
    timer = time.time() # to monitor execution time
    st=time.time()
    #fragment each peptide to b and y ions
    for i in range(0,len(peptidemass['peptide'])):
        pep=peptidemass['peptide'][i]
        pepy=pep[::-1]#reverse string to fragment from C terminal
        pmass=peptidemass['mass'][i]
        pno=peptidemass['pep no'][i]
        modox=peptidemass['oxi'][i]
        modam=peptidemass['deam'][i]
        pos1=peptidemass['oxipos'][i]
        pos2=peptidemass['deampos'][i]

        for k in range(1,len(pep)):
            bion=pep[:k]
            yion=pepy[:k]
            fragb=[pno,pep,pmass,modox,modam,bion,k,pos1,pos2]
            fragy=[pno,pep,pmass,modox,modam,yion,k,pos1,pos2]
            dfy.loc[len(dfy)] = fragy
            dfb.loc[len(dfb)] = fragb
        # when programs runs for long report the progress 90 sec
        if time.time() - timer > 90:
            print(f"fragmenting peptide {pep} at {pno}")
            timer = time.time()
    et = time.time()
    print("elapsed time: ", (et-st))
    return(dfb,dfy)

def fragmass(df,iontyp):
    """
    Fragment mass is calcuated using the monoisotopic mass and presence of modification on residue.
    it creates 2 other dataframe or searching with peptide, fragment, mass and modification information
    """
    mass=[]
    #massy=[]
    m=0
    n=1
    mono = {'A' :71.0371, 'B': 114.5349,'C' :161.0, 'D' :115.0269, 'E' :129.0426, 'F' :147.0684, 'G' :57.02146, 'H' :137.0589, 'I' :113.0841, 'K' :128.0950, 'L' :113.0841, 'M' :131.0405, 'N' :114.0429, 'P' :97.0528, 'Q' :128.0586, 'R' :156.1011, 'S' :87.0320, 'T' :101.0477, 'U': 150.9536, 'V' :99.0684, 'W' :186.0793, 'X':111, 'Y' :163.0633,'Z':128.55059, '*' :0.0}
    
    for i in range(0,len(df['pep no.'])):
        pep=df['peptide'][i]
        pno=df['pep no.'][i]
        if iontyp=='b':
            ion=df['b ion'][i]
        elif iontyp=='y':
            ion=df['y ion'][i]
        pepmass=df['pepmass'][i]
        modd1=df['oxi'][i]
        modd2=df['deam'][i]
        pos1=df['oxipos'][i]
        pos2=df['deampos'][i]
        c=0
        cn=0
        cc=0
        for j in range(0,len(ion)):
            m = m+mono[ion[j]]
            if pos1!=pos2:
                if ion[j]=='C':
                    cc=cc+1
            if pos1!='':
                pos= pos1.split(' ')
                for x in tuple(pos):
                    if j ==int(x):
                        m=m+16
            if pos2!='':
                pos= pos2.split(' ')
                for y in pos:
                    if j== int(y):
                        m=m+0.984
        mass.append([pno,pep,pepmass,ion,(m+1)+cc*(1),j+1,modd1,modd2,pos1,pos2])
    #massy.append([pno,pep,pepmass,iony,(m+19)+cc*(1),j+1,modd1,modd2,pos1,pos2])    
        m=0
        #print(mass)
    #link all information about precurssor and fragment ion
    #linkedy=pd.DataFrame(massy,columns=['pep no.','peptide', 'pepmass','y ion', 'y m/z','yn','oxi','deam','oxipos','deampos'])
    if iontyp=='b':
        linked=pd.DataFrame(mass,columns=['pep no.','peptide', 'pepmass','b ion', 'b m/z','bn','oxi','deam','oxipos','deampos'])
    elif iontyp=='y':
        linked=pd.DataFrame(mass,columns=['pep no.','peptide', 'pepmass','y ion', 'y m/z','bn','oxi','deam','oxipos','deampos'])
    return(linked)

def searchquery(queryn):
    querydict={}
    for i in range(0,int(queryn)):
        fname=input("Enter the file name (eg.: 1214.973peaks.txt): ")
        qpep =fname.split('p')[0]#assumes name entered as in example
        querydf = pd.read_csv(fname, sep='\t', header=None, names=['mass', 'intensity']) 
        querydict[qpep]=querydf['mass'].to_list()#dictionary of precurssor mass keyed to list of fragment mass values
    return(querydict)

def searchdb(linkedb,linkedy,query):
    """
    this function searches database created for query peptide and its b and y fragments are matched with peak list data
    saves all the matched peptides and fragment data in a dataframe
    """
    result=pd.DataFrame(columns=['queryfragm','ionmass','Delmass','iontype','ionseq','ion_n','pepseq','pep no.','pepmass', 'oxipos','deampos'])
    for key, value in query.items():
        for j in range(0,len(linkedy['y m/z'])):
            if float(key)-0.5<linkedy['pepmass'][j]<float(key)+0.5:
                for k in range(0,len(value)):
                    if linkedy['y m/z'][j]-0.5< value[k] < linkedy['y m/z'][j]+0.5:#falls within range
                        delta=value[j]-linkedy['y m/z'][j]#difference between obv and cal
                        pepdelta=float(key)-float(linkedy['pepmass'][j])
                        result.loc[len(result)] = [round(value[j],4),round(linkedy['y m/z'][j],4),delta,"y ion", linkedy['y ion'][j],linkedy['yn'][j],linkedy['peptide'][j],linkedy['pep no.'][j],round(linkedy['pepmass'][j],4), linkedy['oxipos'][j],linkedy['deampos'][j]]
                    if linkedb['b m/z'][j]-0.5< i < linkedb['b m/z'][j]+0.5:#falls within range
                        delta=value[j]-linkedb['b m/z'][j]
                        pepdelta=float(key)-linkedb['pepmass'][j]
                        result.loc[len(result)] = [round(value[j],4),round(linkedb['b m/z'][j],4),delta,"b ion", linkedb['b ion'][j],linkedb['bn'][j],linkedb['peptide'][j],linkedb['pep no.'][j],round(linkedb['pepmass'][j],4), linkedb['oxipos'][j],linkedb['deampos'][j]]
                
    return(result)

def printresults(result,protdictname,protdictseq):
    """
    print output of protein hit, protein sequence covered, peptide matched, modification present,
    and matched fragment data
    it also computes a score based on number of fragments matched for identified peptide and consecutive atches identified
    normalised to the length of peptide
    it prints on terminal
    """
    #dictionary of number of peptide matches obtained
    a=dict(result['pep no.'].value_counts())
    
    for i in a:
        n_consec=0 #number of consecutive matches counter
        res=result[result["pep no."]==i]
        res = res.reset_index(drop=True)#fix index after deletion of duplicates
        res['consec'] = res.ion_n.diff()
        maxmatch=2*(len(res['pepseq'][0])-1)#maximum number of fragment matches for peptide possible
        
        for b in res.consec[1:]:
            if (b == 1):
                n_consec += 1
                
        for key, value in protdictname.items():
            for j in range(0,len(value)):
                if result['pepseq'][0] in value[j]:
                    print(tabulate([['protein_name','protein _seq','protein_coverage'],[key,protdictseq[key],parser.coverage(protdictseq[key],result['pepseq'][0])]],tablefmt="grid"))
                    print(tabulate([['pep no.','pepseq','pepmass','oxipos','deampos','frag match','match','consec match'],[i,res['pepseq'][0],res['pepmass'][0], res['oxipos'][0],res['deampos'][0],a[i],round((a[i]-1)/maxmatch,3), n_consec]],tablefmt="grid")) #score calcuated based on length and maximum matches possible
                    print(tabulate((res.iloc[:,:6]),headers=['index','queryfragm','ionmass','fragmassdel','ion type','ionseq','ion_n'],tablefmt="grid"))
    
    
if __name__ == "__main__":

    import argparse
    
    parser = argparse.ArgumentParser(description="Create target dataabse and search against MALDI MS2 data for confident peptide and protein identification")
    parser.add_argument("input_file", help="Input file containing LC-MS/MS data exported from Mascot search")
    parser.add_argument("queryn", help="Input number of query files to search")

    # assign required input arguments to variables
    input_file = args.input_file
    queryn = args.queryn
    
    # read the input file into  dataframe
    peplist,protname,protseq=readData(input_file)
    pepmass=peptidemass(peplist)
    datab,datay=fragmentation(pepmass)
    databaseb=fragmass(datab,'b')
    databasey=fragmass(datay,'y')
    query=searchquery(queryn)
    search_op=searchdb(databaseb,databasey,query)
    printresults(search_op,protname,protseq)
    
    