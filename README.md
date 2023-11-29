# RP2
Workflow that takes LC-MS/MS search results from MASCOT to create a targeted reduced database and
searches MALDI MSI peaklists for peptide assignment and isobaric peptide differenciation with a rudimentary score assigned to matches 

#To run program:
#step 1: ensure script and files are in same directory
step 2: open terminal and open directory 
          >> cd "file path"
step 3: run program
          >> python RP2.py <.csv file name> <n. of peaklists to search>
          eg: python RP2.py test1.csv 2
step 4: enter peaklist name when prompted eg:783.647peaks.txt 





