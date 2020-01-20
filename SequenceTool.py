from Bio import SeqIO
from Bio import pairwise2
from Bio.SubsMat.MatrixInfo import blosum62
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio import AlignIO
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceCalculator
from Bio.SeqUtils import GC
from Bio.SeqUtils import six_frame_translations

def nj_tree_constructor(x):
    constructor = DistanceTreeConstructor()
    calculator = DistanceCalculator("identity")
    dm = calculator.get_distance(x)
    njtree = constructor.nj(dm)
    print(njtree)
    Phylo.draw_ascii(njtree)

def upgma_tree_constructor(x):
    constructor = DistanceTreeConstructor()
    calculator = DistanceCalculator("identity")
    dm = calculator.get_distance(x)
    upgmatree = constructor.upgma(dm)
    print(upgmatree)
    Phylo.draw_ascii(upgmatree)

num_of_sequences = int(input("Enter the number of sequences to analyse... "))

if num_of_sequences == 0:
    print("What are you even doing?")

if num_of_sequences == 1:
    file = SeqIO.read(input("Enter the path to the FASTA file "), "fasta")
    complete_sequence = input("Is your sequence complete? yes or no ")
    if complete_sequence == str("yes"):
        print("GC %")
        gc_values = GC(file.seq)
        print(gc_values)
        translation = file.seq.translate()
        print(translation)
    if complete_sequence == str("no"):
        print("Six Frame Translation:")
        six_frame = six_frame_translations(file.seq)
        print(six_frame)

if num_of_sequences == 2:
    file1 = SeqIO.read(input("Enter the path to the first FASTA file "), "fasta")
    file2 = SeqIO.read(input("Enter the path to the second FASTA file "), "fasta")
    alignment = pairwise2.align.globalds(file1.seq, file2.seq, blosum62, -10, -0.5)
    print("Pairwise Alignment")
    print(pairwise2.format_alignment(*alignment[0]))

if num_of_sequences > 2:
    print("Get ready for some MSA")
    in_file = input("Enter the path to the FASTA file containing multiple sequences ")
    out_file = "aligned.fasta"
    clustalomega_cline = ClustalOmegaCommandline(infile=in_file, outfile=out_file, verbose=True, auto=True)
    clustalomega_cline()
    align = AlignIO.read(out_file, "fasta")
    print("Multiple Sequence Alignment")
    print(align)
    next_step = input("do you want to view the phylogenetic tree? yes or no ")
    if next_step == "yes":
        choice = input("UPGMA, Maximum Likelihood or Neighbor Joining? ")
        if choice == str("UPGMA"):
            upgma_tree_constructor(align)
        if choice == str("Neighbour Joining"):
            nj_tree_constructor(align)
        else:
            print("Error, input unknown")
    if next_step == "no":
        print("Thank you for using The Kailignment Tool")
    if next_step != "yes" and next_step != "no":
        print("Error, input unknown.")

