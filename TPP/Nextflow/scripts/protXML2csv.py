# converter of protxml files to csv INCLUDING details from StPeter
from pyteomics import protxml,pepxml
from math import nan
import pandas as pd
pp = protxml.read("StPeterOut.prot.xml")
pepp = pepxml.read("Sample.pep.xml")
pdata = protxml.DataFrame(pp)
pepdata = pepxml.DataFrame(pepp)
pepdata["Quantification_SI"] = nan

def extract_stpeter(analysis):
    if isinstance(analysis, list):    
        a = analysis[0]
        if (len(a) == 2):
            b = a.get("StPeterQuant")
            SI = b.get("SI")
            SIn = b.get("SIn")
            peps = b.get("StPeterQuant_peptide")
            peptable = pd.DataFrame(peps)
            pepseqs = peptable.get("sequence").to_string()
            # adding SI values to peptide table
            for index, r in peptable.iterrows():
                pepdata.at[pepdata["modified_peptide"] == r.get("sequence"), "Quantification_SI"] = r.get("SI") 
                
            pepSI = peptable.get("SI").to_string()
            return({"analysis":a.get("analysis"),"SI":SI, "SIn":SIn, "peptides": pepseqs, "pepSIs": pepSI})
    return({"analysis":nan,"SI":nan, "SIn":nan, "peptides": nan, "pepSIs": nan})
    
StPeter = pdata.apply(lambda row: extract_stpeter(row.get("analysis_result")), axis=1)

StPeter = pd.DataFrame(StPeter.tolist())

All = pd.concat([pdata, StPeter], axis=1)

All.to_csv("StPeterProts.csv")
pepdata.to_csv("StPeterPeps.csv")
