from Bio import SeqIO

def pseudo_fasta_coordinates(path_to_fasta, restriction=None):
    '''
    Takes in FASTA with separate chr entries and returns a table with the pseudo entry borders of each entry after the collapsing of the entries. Currently filters atypical chromosomes.
    '''
    with open(path_to_fasta, 'r') as fasta_file:
        entries={}
        seqlist=[]
        updated_length=0
        for record in SeqIO.parse(fasta_file, 'fasta'):
            if restriction == "unrestricted":
                #this part makes sure now downstream errors occur while doing splits by _
                record.id = record.id.replace("_", "-")
                entries[record.id] = [updated_length, updated_length + len(record.seq)]
                updated_length = updated_length + len(record.seq)
                seqlist.append(str(record.seq))
            elif '_' not in record.id and 'M' not in record.id:
                entries[record.id] = [updated_length, updated_length + len(record.seq)]
                updated_length = updated_length + len(record.seq)
                seqlist.append(str(record.seq))
    return len(''.join(seqlist)), entries
