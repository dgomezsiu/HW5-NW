# Import NeedlemanWunsch class and read_fasta function
from align import read_fasta, NeedlemanWunsch

def main():
    """
    This function should
    (1) Align all species to humans and print species in order of most similar to human BRD
    (2) Print all alignment scores between each species BRD2 and human BRD2
    """
    hs_seq, hs_header = read_fasta("./data/Homo_sapiens_BRD2.fa")
    gg_seq, gg_header = read_fasta("./data/Gallus_gallus_BRD2.fa")
    mm_seq, mm_header = read_fasta("./data/Mus_musculus_BRD2.fa")
    br_seq, br_header = read_fasta("./data/Balaeniceps_rex_BRD2.fa")
    tt_seq, tt_header = read_fasta("./data/tursiops_truncatus_BRD2.fa")

    # TODO Align all species to humans and print species in order of most similar to human BRD
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    nw = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", -10.0, 01.0)

    # gg
    gg_score, gg_hs_align, hs_gs_align = nw.align(hs_seq, gg_seq)

    # mm
    mm_score, mm_hs_align, hs_mm_align = nw.align(hs_seq, mm_seq)

    # br
    br_score, br_hs_align, hs_br_align = nw.align(hs_seq, br_seq)

    # tt
    tt_score, tt_hs_align, hs_tt_align = nw.align(hs_seq, tt_seq)

    # sort species in order of most similar to least similar compared to hs
    species = {
        'gg': (gg_score, gg_hs_align),
        'mm': (mm_score, mm_hs_align),
        'br': (br_score, br_hs_align),
        'tt': (tt_score, tt_hs_align)
    }

    sorted_species = {i: species[i] for i in sorted(species, key = species.get, reverse = True)}

    # TODO print all of the alignment score between each species BRD2 and human BRD2
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    print('Gallus gallus: ', gg_hs_align, " ", gg_score)
    print('Mus musculus: ', mm_hs_align, " ", mm_score))
    print('Balaeniceps rex: ', br_hs_align, " ", br_score))
    print('Tursiops truncatus: ', tt_hs_align, " ", tt_score))
    

if __name__ == "__main__":
    main()
