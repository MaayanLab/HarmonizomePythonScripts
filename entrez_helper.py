from Bio import Entrez
import sys
import time

class entrez:
  def __init__(self, email, api_key):
    Entrez.email = email
    Entrez.api_key = api_key


  def get_pmid(self, term, retmax=2000):
    handle = Entrez.esearch(db="pubmed", retmax=retmax, term=term)
    record = Entrez.read(handle)
    handle.close()
    if "IdList" in record:
      return record["IdList"]
    else:
      return None

  def get_all_pmids_of_list(self, term_list, term_dict=None, timeout1=0.3, timeout2=5):
    terms_pmids = {}
    for term in term_list:
        # Continue where it failed
        # sys.stdout.write("Processing disease %s\r" % (disease))
        if term not in terms_pmids:
            sys.stdout.write("Processing term %s\r" % (term))
            if term_dict and term in term_dict:
              pmids = term_dict[term] # Check if the term was previously discovered on a different library
            else:
              time.sleep(timeout1)
              for i in range(5):
                  try:
                      pmids = self.get_pmid(term)
                  except Exception as e:
                      if i == 4:
                          raise e
                      else:
                          time.sleep(timeout2)
                          continue
            terms_pmids[term] = pmids
    return terms_pmids