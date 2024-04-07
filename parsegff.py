# python3 -m pip install numpy
import numpy as np
from Bio import SeqIO

# Open the GFF3 file
gff_file = "example.gff3"

# Open the genome FASTA file
fasta_file = "genome.fa"

# Parse the genome FASTA file
genome_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

# Parse the GFF3 file
for record in SeqIO.parse(gff_file, "gff"):
    print(f"Record ID: {record.id}")
    
    # Get the corresponding genome sequence
    genome_seq = genome_dict[record.id].seq
    
    # Iterate over the features in the record
    for feature in record.features:
        print(f"Feature Type: {feature.type}")
        print(f"Location: {feature.location}")
        
        # Extract the feature sequence from the genome
        feature_seq = feature.extract(genome_seq)
        print(f"Sequence: {feature_seq}")
        
        print(f"Qualifiers:")
        for key, value in feature.qualifiers.items():
            print(f"  {key}: {value}")
        print("---")