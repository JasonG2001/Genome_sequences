from Bio import SeqIO
from collections import defaultdict
from typing import Dict
import pandas as pd


def read_csv(name: str) -> pd.DataFrame:

    """Reads in the csv file, extracting only the column under the name 'KH ID' and returns this as a pandas dataframe

    Parameters
    ----------
    name: str
        name of the csv file if in the same directory, otherwise, enter the relative path

    Returns
    -------
    pd.DataFrame
        the names under 'KH ID' are returned as a pandas dataframe
    """

    df: pd.DataFrame = pd.read_csv(name)
    df: pd.DataFrame = df.loc[:, "KH ID"]
    return df


def make_dictionary(fasta_name: str) -> Dict[str,str]:

    """Converts a fasta file into the form of a dictionary where the title is the key and sequence is the value

    Parameters
    ----------
    fasta_name: str
        name of the fasta file if in the same directory, otherwise, enter the relative path

    Returns
    -------
    Dict[str,str]
        converts the fasta file into a python dictionary format with title as the key and sequence as the value
    """
    
    sequences: Dict[str,str] = defaultdict(str)
    
    fasta_sequences = SeqIO.parse(open(fasta_name),'fasta')
    with open(fasta_name) as _:
        for fasta in fasta_sequences:
            name, sequence = fasta.id, str(fasta.seq)
            
            sequences[name] = sequence
            
    return sequences


def get_result(df: pd.DataFrame, sequences: Dict[str,str], total_sequence: str) -> Dict[str,str]:
    
    """Returns the resultant dictionary ready for final conversion to a fasta file.

    Paramters
    ---------
    df: pd.Dataframe
        dataframe consisting of only the names of interest
    sequences: Dict[str,str]
        dictionary of all the names with their corresponding sequences obtained from a general fasta file
    total_sequence: str
        the overall genome sequence which we plan to search for the specific gene

    Returns
    -------
    Dict[str,str]
        a dictioanry with our wanted names and their corresponding sequences in the genome
    """

    res: Dict[str,str] = defaultdict(str)

    wanted_name: str
    for wanted_name in df:
        wanted_sequence: str = get_wanted_sequence(wanted_name, sequences, total_sequence)
        
        res[wanted_name]: str = wanted_sequence
        
    return res
    

def get_wanted_sequence(sequence_name: str, sequences: Dict[str,str], total_sequence: str) -> str: 

    """Searches a dataset of gene names and their corresponding sequence on record and returns the wanted section of the genome

    Parameters
    ----------
    sequence_name: str
        name of the sequence that we want to search for
    sequences: str
        the dataset on record containing all models and names of a gene with their corresponding sequence
    total_sequence: str
        the entire genome sequence which we plan on checking our wanted sequence against

    Returns
    -------
    str
        if either the normal sequence or the reverse sequence is found in the genome, promoters of length 
        2000 are returned. If the normal sequence is found, the 2000 bases are taken before the match. 
        Else if the inverse sequence is found, the 2000 bases after the match is returned.

        if the normal and reverse sequence are not found in the genome, then 'not available' is returned
    """

    name: str
    for name in sequences:
        if sequence_name in name:
            normal_sequence_to_check: str = sequences[name]
            reverse_sequence_to_check: str = reverse_sequence(normal_sequence_to_check)
            if check_presence(normal_sequence_to_check, total_sequence):
                promoter: str = get_2000(normal_sequence_to_check, total_sequence, before=True)
                
                return promoter           
              
            if check_presence(reverse_sequence_to_check, total_sequence): 
                reverse_promoter: str = get_2000(reverse_sequence_to_check, total_sequence, before=False)
                reversed_sequence: str = reverse_sequence(reverse_promoter)   
                
                return reversed_sequence
            
    return "not_available"


def get_2000(sequence: str, total_sequence: str, before: bool) -> str:

    """Obtains the 2000 base sequence wanted for each gene name

    Parameters
    ----------
    sequence: str
        the sequence being located somewhere in the genome, obtained from the overall dataset of names to 
        sequence
    total_sequence: str
        the total sequence of the genome
    before: bool
        True or False values will determine if the 2000 base sequence is taken before or after the match

    Returns
    -------
    str
        returns the 2000 base sequence before or after. Except block catches instance if there are no 2000
        base sequence before of after. This case returns "out_of_genome"
    """

    try:
        start_index: int = total_sequence.find(sequence)
        if before:
            return total_sequence[start_index - 2000:start_index]
        else:
            return total_sequence[start_index + len(sequence):start_index + len(sequence) + 2000]
    except:
        return "out_of_genome"


def reverse_sequence(sequence: str) -> str:
    
    """Function 'reverses' the sequence through changing the bases to their corresponding base and
    reversing the entire order

    Parameters
    ----------
    sequence: str
        gene sequence that we wish to obtain the corresponding sequence for

    Returns
    -------
    str
        the corresponding sequence
    """

    corresponding_sequence: str = ""
    
    base_pairs: Dict[str,str] = {
        "A": "T",
        "C": "G",
        "T": "A",
        "G": "C",
        "N": "N"
    }
    
    base: str
    for base in sequence:
        corresponding_sequence += base_pairs[base]
        
    corresponding_sequence: str = corresponding_sequence[::-1]
    
    return corresponding_sequence
                    
            
def convert_to_one_sequence(check_sequences: Dict[str,str]) -> str:
    
    """The entire genome fasta file contains names which break up each section of the genome. To simplfy
    the process of finding a specific gene inside the genome, the segmented sections of the genome are 
    combined together to give a long sequence

    Parameters
    ---------
    check_sequences: Dict[str,str]
        the whole genome fasta file converted to the form of a dictionary

    Returns
    -------
    str
        the whole sequence chain for the genome
    """

    total_sequence: str = ""
    
    sequence: str
    for sequence in check_sequences.values():
        total_sequence += sequence
        
    return total_sequence
            
            
def check_presence(sequences_to_check: str, total_sequence: str) -> bool:
    
    """Quick check to determine if the sequence is within the genome

    Parameters
    ----------
    sequence_to_check: str
        the gene sequence to look for
    total_sequence: str
        the entire genome sequence

    Returns
    -------
    bool
        True is gene is found, else false if gene is not found in the genome
    """

    if sequences_to_check in total_sequence:
        return True
    else:
        return False
    
    
def convert_to_fasta(df: pd.DataFrame, sequences: Dict[str,str], total_sequence: str, name: str) -> None:
    
    """Converts the final dictionary being processed into a .txt file in the format of a fasta file

    Parameters
    ----------
    df: pd.DataFrame
        dataframe of just the wanted gene names
    sequences: Dict[str,str]
        sequence of the names of all genes and their models to their gene sequence
    total_sequence: str
        total genome sequence
    name: str
        name of the that's wanted to be generated
    """

    final_dictionary: Dict[str,str] = get_result(df, sequences, total_sequence)
    
    text_file = open(name, "w")
    
    name: str
    for name in final_dictionary:
        text_file.write(">" + name + "\n" + final_dictionary[name] + "\n")
        
    text_file.close()
